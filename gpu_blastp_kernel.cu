#include <algo/blast/core/blast_aalookup.h>

#include <algo/blast/core/gpu_cpu_common.h>
#include "gpu_blastp_kernel.h"


//shared and constant memory
extern __shared__ int array[];
__constant__ char d_Query_const[16384];
__constant__ int d_x_dropoff_const[4096];
__constant__ int d_cutoff_score_const[4096];
__constant__ int d_query_offset_const[4096];


/**
 * Each thread scans on of the subjet sequences and
 * it returns the score and coordinats of each successful alignment
 * @param d_Database2Dpadded                   the formatted sequence database[in]
 * @param d_RepeatedSubstitutionMatrix         multiple copies of the substitution matrix (one for each block for parallel loading from global to shared) [in]
 * @param SubstitutionMatrix_length            the substitution matrix [in]
 * @param d_RepeatedPV                         multiple copies of the presence vector (one for each block for parallel loading from global to shared) [in]
 * @param d_ThickBackbone                      ThickBackbone as taken from NCBI_BLAST [in]
 * @param d_overflow                           overflow array as taken from NCIB_BLAST [in]
 * @param *d_RepeatedSequenceMaxlength_vector  multiple copies of the vector that has the maximum length for each group [in]
 * @param SequenceMaxlength_vector_stride      the stride of the array that holds the maximum lengths of each group [in]
 * @param NumSequences                         number of sequences the GPU works on [in]
 * @param pv_length                            length of the presence vector [in]
 * @param Query_length                         length of the query [in]
 * @param diag_mask                            mask for the diagonal array [in]
 * @param diag_array_length                    length of the diagonal array used by each thread [in]
 * @param num_quries                           number of queries [in]
 * @param word_length                          word length [in]
 * @param charsize                             char size of the lookup table [in]
 * @param index_mask                           mask for the iniatialization of the index [in]
 * @param window                               window for the two-hit alignment [in]
 * @param diag_offset                          diagonal offset [in] 
 * @param d_diag_array                         allocated global memory used for storing the diagonal arrays [in]
 * @param d_Hits returns                       returns hits found for each sequence [out]
 * @param d_GPUBlastInitHitList                returns the statistics of the successful extensions [out]
*/

__global__ void
GPU_BLASTP_kernelTwoHit( const int* d_Database2Dpadded,
                         const int* d_RepeatedSubstitutionMatrix,
                         const int SubstitutionMatrix_length,
                         const PV_ARRAY_TYPE* d_RepeatedPV,
                         const AaLookupSmallboneCell* d_ThickBackbone,                         
                         const unsigned short* d_overflow,                         
                         const int* d_RepeatedSequence_length_vector,
                         const int Sequence_length_vector_stride,
                         const int NumSequences,
                         const int pv_length,
                         const int Query_length,
                         const int diag_mask,
                         const int diag_array_length,
                         const int num_queries,
                         const int word_length,
                         const int charsize,
                         const int index_mask,
                         const int window,
                         int diag_offset,
                         int* d_diag_array,
                         int* d_Hits,
                         GPUBlastInitHitList* d_GPUBlastInitHitList)
{
    int bx = blockIdx.x;
    int tx = threadIdx.x;
    int blockDimx = blockDim.x;

    //define the limints of each array in the shared memory 

    int* SubstitutionMatrix_shared_int = (int*)array;
    int* Sequence_length_vector_shared = (int*)&SubstitutionMatrix_shared_int[SubstitutionMatrix_length>>2];//The substitution matrix is 640 bytes or 160 integers long
    char* Sequence_AA = (char*) &Sequence_length_vector_shared[Sequence_length_vector_stride];
    PV_ARRAY_TYPE* pv_shared = (PV_ARRAY_TYPE*) &Sequence_AA[2* blockDimx * 4];

    //pv array from global to shared memory
    int pv_Index_global = bx*pv_length + tx;
    int pv_Index_shared = tx;
    while( pv_Index_shared < pv_length )//SubstitutionMatrix_length is a multiple of 64 which includes some padded zeros so there is some waste of memory space. Can it be recovered?
    {
        pv_shared[ pv_Index_shared ] = d_RepeatedPV[ pv_Index_global ];

        pv_Index_global += blockDimx;
        pv_Index_shared += blockDimx;
    }

    __syncthreads();

    //Substitution Matrix from global to shared memory
    int SubstitutionMatrix_Index_global = bx*(SubstitutionMatrix_length>>2) + tx;
    int SubstitutionMatrix_Index_shared = tx;
    while( SubstitutionMatrix_Index_shared < (SubstitutionMatrix_length>>2) )//SubstitutionMatrix_length is a multiple of 64 which includes some padded zeros so there is some waste of memory space. Can it be recovered?
    {
        SubstitutionMatrix_shared_int[ SubstitutionMatrix_Index_shared ] = d_RepeatedSubstitutionMatrix[ SubstitutionMatrix_Index_global ];

        SubstitutionMatrix_Index_global += blockDimx;
        SubstitutionMatrix_Index_shared += blockDimx;
    }

    __syncthreads();

    //SequenceMaxlength_vector from global to shared memory
    int Sequence_length_vector_Index_global = bx*Sequence_length_vector_stride + tx;
    int Sequence_length_vector_Index_shared = tx;
    while( Sequence_length_vector_Index_shared < Sequence_length_vector_stride )
    {
        Sequence_length_vector_shared[ Sequence_length_vector_Index_shared ] = d_RepeatedSequence_length_vector[ Sequence_length_vector_Index_global ];

        Sequence_length_vector_Index_global += blockDimx;
        Sequence_length_vector_Index_shared += blockDimx;
    }

    __syncthreads();


    char* SubstitutionMatrix_shared_char = (char*)SubstitutionMatrix_shared_int;

    int Sequence = (bx * blockDimx + tx);
    int Sequence_Index_base = 0;
    int Sequence_Index_base2 = 0;
    int Sequence_Index = 0;

    int Sequence_AAI, Sequence_AAII;
    int counter = 0;
    GPUBlastInitHitList data;

    int Group_size = gridDim.x*blockDimx;
    int diag_array_base = ((bx*blockDimx+tx)*diag_array_length);
    //tx>>4 = divide by 16 (each half warp has 16 threads) to find in which half warp does this thread belong in
    //<<7 = multiply by 128 because each half warp occupies 128 bytes (16 (threads) * (4 bytes (for Sequence_AAI)+ 4 bytes(for Sequence_AAII))
    //(tx&0xF) = Module 16 to find the exact thread of the half warp that the current (tx) thread belongs in
    //<<2 = Multiply by 4 to  move to the bank that we store Sequence_AAI. For Sequence_AAII add 64 (see below)
    int Shared_AA_Index = ((tx>>4)<<7) + ((tx&0xF)<<2);

    while( Sequence < NumSequences )
    {
        int Hits = 0;
        int score = 0;
        int s_last_off = 0;
        int hits_extended = 0;
        int successful_extensions = 0;
        int Sequence_length = Sequence_length_vector_shared[ counter ];

        Sequence_Index = Sequence_Index_base + bx*blockDimx + tx;

        //When the lower triangle of the matris is scanned the Sequence is always scanned from the first element thus initially always the first byte is loaded1~
        Sequence_AAI = d_Database2Dpadded[ Sequence_Index + 0 ];
        Sequence_AA[ Shared_AA_Index + 0 ] = Sequence_AAI >>  0;
        Sequence_AA[ Shared_AA_Index + 1 ] = Sequence_AAI >>  8;
        Sequence_AA[ Shared_AA_Index + 2 ] = Sequence_AAI >> 16;
        Sequence_AA[ Shared_AA_Index + 3 ] = Sequence_AAI >> 24;

        int index = (Sequence_AA[ Shared_AA_Index + 0 ] << charsize) |  Sequence_AA[ Shared_AA_Index + 1 ];

        for(int s_off = 0; (s_off <= Sequence_length - word_length); ++s_off)
        {
            //In every iteration eight amino acids are stored in shared memory and they trasnferred from global only when necessary
            //The last condition in the "if" statement is to exclude the extra loading after the end of the subject
            //E.g. if the length of the subject is 64 characters, then when i == 56 the last loading should happen by the "else" statemenet below
            //and when i == 60 ( (60>>2) == (64>>2)-1 ) there should be not loading

            if( ((s_off&0x3) == 0) && (((s_off >> 2) & 0x1) == 1) && ((s_off>>2) < ((Sequence_length>>2) - 1)) )
            {
                Sequence_AAI = d_Database2Dpadded[ Sequence_Index + (((s_off>>2) + 1)*Group_size) ];//Division by 4 (>>2) because the sequence is accessed as int
                Sequence_AA[ Shared_AA_Index + 0 ] = Sequence_AAI >>  0;
                Sequence_AA[ Shared_AA_Index + 1 ] = Sequence_AAI >>  8;
                Sequence_AA[ Shared_AA_Index + 2 ] = Sequence_AAI >> 16;
                Sequence_AA[ Shared_AA_Index + 3 ] = Sequence_AAI >> 24;

            }else if( ((s_off&0x3) == 0) && ( ((s_off >> 2) & 0x1) == 0) ) {

                Sequence_AAII = d_Database2Dpadded[ Sequence_Index + (((s_off>>2) + 1)*Group_size) ];
                Sequence_AA[ Shared_AA_Index + 64 ] = Sequence_AAII >>  0;
                Sequence_AA[ Shared_AA_Index + 65 ] = Sequence_AAII >>  8;
                Sequence_AA[ Shared_AA_Index + 66 ] = Sequence_AAII >> 16;
                Sequence_AA[ Shared_AA_Index + 67 ] = Sequence_AAII >> 24;

            }
            //&0x7 = Modulo 8
            //>>2 = divide by 4 to find which of the two sections (s_off+2) belongs in (lower 64 or the upper 64)
            //<<6 = multiply by 64 to get into the right section
            //(s_off+2)&0x3 = Modulo 4 to go to the exact position
            index = ((index<<charsize) | Sequence_AA[ Shared_AA_Index + ((((s_off + 2)&0x7)>>2)<<6) + (((s_off+2)&0x3)) ]) & index_mask;

            int numhits = 0;
            if( PV_TEST(pv_shared, index, PV_ARRAY_BTS) ) //check via the precence vector (PV) if this database subject word exists in the query
            {
                numhits = d_ThickBackbone[index].num_used;
//                Hits += numhits; TO DO:
                Hits += d_ThickBackbone[index].num_used;
                
                const unsigned short *src;
                if( numhits <= AA_HITS_PER_CELL ) //check there are less than AA_HITS_PER_CELL occurrences of this word in the query
                    src = d_ThickBackbone[index].payload.entries;
                else
                    src = &(d_overflow[d_ThickBackbone[index].payload.overflow_cursor]); //if not fetch the occurrences from the overflow array
                
                for(int i = 0; i < numhits; ++i) //for each occurrence expand and calculate score
                {

                    ExtendTwoHitReturn ReturnValue;

                    int query_offset = src[i];
                    int subject_offset = s_off;

                    int diag_coord = (query_offset - subject_offset) & diag_mask;

                    int diagonal = d_diag_array[diag_array_base + diag_coord];
                    int last_hit = diagonal & 0x7FFFFFFF;
                    int flag = diagonal & 0x80000000;


                    /* If the reset bit is set, an extension just happened. */
                    if (flag)
                    {
                        /* If we've already extended past this hit, skip it. */
                        if ( (subject_offset + diag_offset) < last_hit)
                            continue;
                        else /* Otherwise, start a new hit. */
                            d_diag_array[diag_array_base + diag_coord] = (subject_offset + diag_offset) & 0x7FFFFFFF;//set the MSB to zero and the rest to (subject_offset + diag_offset)
                        
                    } else { /* If the reset bit is cleared, try to start an extension. */
                        
                        /* find the distance to the last hit on this diagonal */
                        last_hit = last_hit - diag_offset;
                        int diff = subject_offset - last_hit;
                        
                        if (diff >= window)
                        {
                            /* We are beyond the window for this diagonal; start a new hit */
                            /* If the difference is less than the wordsize (i.e. last hit and this hit overlap), give up */
                            d_diag_array[diag_array_base + diag_coord] = (subject_offset + diag_offset) & 0x7FFFFFFF; //set only the 31 LSB's (the flag is zero since we are here)
                            continue;
                        }

                        /* If the difference is less than the wordsize (i.e. last hit and this hit overlap), give up */
                        if (diff < word_length)
                            continue;

                        /* Extend this pair of hits. The extension to the left must reach the end of the first word in order for extension to the right to proceed.
                           To use the cutoff and X-drop values appropriate for this extension, the query context must first be found */
                        int k = 0;
                        for(k = num_queries-1; k >= 0; --k)
                            if( query_offset >= d_query_offset_const[k] )
                                break;
                            
                        
                        ReturnValue = s_BlastAaExtendTwoHit(SubstitutionMatrix_shared_char,
                                                            d_Database2Dpadded + Sequence_Index,
                                                            Sequence_length,
                                                            d_Query_const,
                                                            Query_length,
                                                            last_hit + word_length,
                                                            subject_offset, query_offset,
                                                            d_x_dropoff_const[k],
                                                            word_length);

                        

                        score = ReturnValue.MAX_score;
                        s_last_off = ReturnValue.s_last_off;

                        ++hits_extended;

                        /* if the hsp meets the score threshold, report it */
                        if ( score >= d_cutoff_score_const[k] ){

                            data.hsp_q = ReturnValue.hsp_q;
                            data.hsp_s = ReturnValue.hsp_s;
                            data.query_offset = query_offset;
                            data.subject_offset = subject_offset;
                            data.hsp_len = ReturnValue.hsp_len;
                            data.score = score;

                            //Check to see if the successful extended seed can fit in the dedicated memory for this thread
                            //The successfully extended seeds are stored only up to the number of 32, after that they are overwritten and in that case the
                            //"succssful_extensions" variable will have a value larger than 32 and this will signal the CPU to carry out that alignment
                            d_GPUBlastInitHitList[ Sequence*MAX_SUCCESSFUL_HITS_PER_SEQUENCE + (successful_extensions&MOD_MAX_SUCCESSFUL_HITS_PER_SEQUENCE) ] = data;

                            successful_extensions++;

//                            if( successful_extensions > MAX_SUCCESSFUL_HITS_PER_SEQUENCE ) {
//
//                                hits_extended = 1000000;
//                                goto end_of_batch;
                                
//                            }
                                
                            
                        }

                        /* If an extension to the right happened, reset the last hit so that future hits to this diagonal must start over. */

                        if (ReturnValue.right_extend) {

                            d_diag_array[diag_array_base + diag_coord] = (s_last_off - (word_length - 1) + diag_offset) & 0x7FFFFFFF;//set the 31 LSB's (the flag is zero size we are here)
                            d_diag_array[diag_array_base + diag_coord] |= 0x80000000; //set the MSB which holds the flag (this can be merged in the previous assignment)

                        } else  /* Otherwise, make the present hit into the previous hit for this diagonal */

                            d_diag_array[diag_array_base + diag_coord] = (subject_offset + diag_offset) & 0x7FFFFFFF;//set only the 31 LSB's (the flag is zero sinze we are here)
                        
                    } /* 
                         if (flag) */
                } /* end for - done with this batch of hits*/
            } //if (PV_TEST...)
            
        }//for Sequence_length loop

//      end_of_batch:
        
        diag_offset += Sequence_length + window;

        d_Hits[ Sequence*3 + 0 ] = Hits;
        d_Hits[ Sequence*3 + 1 ] = successful_extensions;
        d_Hits[ Sequence*3 + 2 ] = hits_extended;

        Sequence += Group_size;
        Sequence_Index_base  += Group_size * (Sequence_length>>2);
        Sequence_Index_base2 += Group_size * (Sequence_length>>2);

        counter++;
    } //while loop
}



/** Perform a two-hit extension. Given two hits L and R, begin
 * at R and extend to the left. If we do not reach L, abort the extension.
 * Otherwise, begin at R and extend to the right.
 *
 * @param matrix          substitution matrix [in]
 * @param subject         subject sequence [in]
 * @param subject_length  subject length [in]
 * @param query           query sequence [in]
 * @param query           query length [in]
 * @param s_left_off      left subject offset [in]
 * @param s_right_off     right subject offset [in]
 * @param q_right_off     right query offset [in]
 * @param dropoff         X dropoff parameter [in]
 * @param word_size       number of letters in one word [in]
 * @return the scores and coordinates of the alignment of the hsp.
 */

static ExtendTwoHitReturn
__device__ s_BlastAaExtendTwoHit(const char* matrix,
                                 const int* subject,
                                 const int subject_length,
                                 const char* query,
                                 const int query_length,
                                 const int s_left_off,
                                 int s_right_off,
                                 int q_right_off,
                                 const int dropoff,
                                 const int word_size)
{
    int left_d = 0, right_d = 0;       /* left and right displacements */
    int left_score = 0, right_score = 0;       /* left and right scores */
    int i, score = 0;
    ExtendTwoHitReturn ReturnValue;

    /* find one beyond the position (up to word_size-1 letters to the right)
       that gives the largest starting score */
    /* Use "one beyond" to make the numbering consistent with how it's done
       for BlastAaExtendOneHit and the "Extend" functions called here. */
    for (i = 0; i < word_size; ++i)
    {
        int index = ((((s_right_off+i)>>2)*gridDim.x*blockDim.x)<<2) + ((s_right_off + i)&0x3);

        score += matrix[ query[q_right_off + i] + ((char*)subject)[index]*ALPHABET_SIZE ];

        if (score > left_score) {
            left_score = score;
            right_d = i + 1;    /* Position is one beyond the end of the word. */
        }
    }

    q_right_off += right_d;
    s_right_off += right_d;

    right_d = 0;
    int right_extend = FALSE;

    ReturnValue.s_last_off = s_right_off;

    /* first, try to extend left, from the second hit to the first hit. */
    ExtendLeftReturn ReturnLeftValue = s_BlastAaExtendLeft(matrix, subject, query, s_right_off-1, q_right_off-1, dropoff, 0);
    left_score = ReturnLeftValue.maxscore;
    left_d = ReturnLeftValue.length;

    /* Extend to the right only if left extension reached the first hit. */
    ExtendRightReturn ReturnRightValue;
    if (left_d >= (s_right_off - s_left_off))
    {
        right_extend = TRUE;

        ReturnRightValue = s_BlastAaExtendRight(matrix, subject, query, s_right_off, q_right_off, dropoff, left_score, subject_length, query_length);
        right_score = ReturnRightValue.maxscore;
        right_d = ReturnRightValue.length;

    }



    ReturnValue.s_last_off = ReturnRightValue.s_last_off;
    ReturnValue.hsp_q = q_right_off - left_d;
    ReturnValue.hsp_s = s_right_off - left_d;;
    ReturnValue.hsp_len = left_d + right_d;
    ReturnValue.MAX_score = MAX(left_score, right_score);
    ReturnValue.right_extend = right_extend;

    ReturnValue.left_score = left_score;
    ReturnValue.left_disp = left_d;
    ReturnValue.right_disp = right_d;

    return ReturnValue;


}


/*
 * Beginning at s_off and q_off in the subject and query, respectively,
 * extend to the left until the cumulative score drops by at least
 * 'dropoff', or the end of at least one sequence is reached.
 *
 * @param matrix     the substitution matrix [in]
 * @param subject    subject sequence [in]
 * @param query      query sequence [in]
 * @param s_off      subject offset [in]
 * @param q_off      query offset [in]
 * @param dropoff    the X dropoff parameter [in]
 * @param length     the length of the computed extension [out]
 * @param maxscore   the best running score from a previous extension;
 *                   the running score of the current extension must exceed
 *                   this value if the extension is to be of nonzero size [in]
 * @return The score and the length of the extension
 */
__device__ ExtendLeftReturn s_BlastAaExtendLeft(const char* matrix,
                                                const int* subject,
                                                const char* query,
                                                const int s_off,
                                                const int q_off,
                                                const int dropoff,
                                                int maxscore)
{
    int n, best_i;
    int score = maxscore;

    n = MIN(s_off, q_off);
    best_i = n + 1;

    char query_temp = 0, sequence_temp = 0;

    for (int i = n; i >= 0; i--)
    {
        int index = ((((s_off - n + i)>>2)*gridDim.x*blockDim.x)<<2) + ((s_off - n + i)&0x3);

        sequence_temp = ((char*)subject)[index];

        query_temp = query[ q_off - n + i ];

        score += matrix[ query_temp*ALPHABET_SIZE + sequence_temp ];

        if (score > maxscore)
        {
            maxscore = score;
            best_i = i;
        }

//         if( (blockIdx.x == 93) && (threadIdx.x == 21) )
//             printf("\nLeft4  bx = %4d, tx = %4d, i = %4d, q_off = %3d, s_off = %3d: S(%c,%c)=%4d, score = %4d, maxscore = %4d, dropoff = %4d, best_i = %4d, ",blockIdx.x, threadIdx.x,i, q_off-n+i, s_off - n + i,
//                    inversemap(query_temp), inversemap(sequence_temp), matrix[ query_temp*ALPHABET_SIZE + sequence_temp], score, maxscore, maxscore-score,  best_i);


        if ((maxscore - score) >= dropoff)
            break;
    }

    ExtendLeftReturn ReturnValues;
    ReturnValues.maxscore = maxscore;
    ReturnValues.length = n - best_i + 1;
    return ReturnValues;

}



/*
 * Beginning at s_off and q_off in the subject and query, respectively,
 * extend to the right until the cumulative score drops by at least
 * 'dropoff', or the end of at least one sequence is reached.
 *
 * @param matrix                      the substitution matrix [in]
 * @param subject                     subject sequence [in]
 * @param query                       query sequence [in]
 * @param s_off                       subject offset [in]
 * @param q_off                       query offset [in]
 * @param dropoff                     the X dropoff parameter [in]
 * @param maxscore                    the best running score from a previous extension;
 *                                    the running score of the current extension must exceed
 *                                    this value if the extension is to be of nonzero size [in]
 * @param Sequence_length             the length of the subject sequence [in]
 * @param Query_length                the length of the query sequence [in]
 * @return The score, lenght, and coordinates of the extension
 */
__device__ ExtendRightReturn s_BlastAaExtendRight(const char* matrix,
                                                  const int* subject,
                                                  const char* query,
                                                  const int s_off,
                                                  const int q_off,
                                                  const int dropoff,
                                                  int maxscore,
                                                  const int Sequence_length,
                                                  const int Query_length)
{
    int i, n, best_i = -1;
    int score = maxscore;

    n = MIN(Sequence_length - s_off, Query_length - q_off);

    char query_temp = 0, sequence_temp = 0;

    //These two initial loadings of Sequence_AAI and Query_AAI can be omitted because i starts from 0 and
    //at the beginning of the for loop they will be loaded whatsoever.

    for (i = 0; i < n; i++)
    {
        int index = ((((s_off + i)>>2)*gridDim.x*blockDim.x)<<2) + ((s_off + i)&0x3);

        sequence_temp = ((char*)subject)[index];

        query_temp = query[q_off+i];

        score += matrix[ query_temp*ALPHABET_SIZE + sequence_temp ];

        if (score > maxscore)
        {
            maxscore = score;
            best_i = i;
        }

        if (score <= 0 || (maxscore - score) >= dropoff)
            break;
    }

    ExtendRightReturn ReturnValues;
    ReturnValues.length = best_i + 1;
    ReturnValues.s_last_off = s_off + i;
    ReturnValues.maxscore = maxscore;

    return ReturnValues;
}


