#include "gpu_blastp.h"
#include "gpu_blastp_kernel.cu"
#include "gpu_blastp_kernel.h"

#include <algo/blast/core/gpu_cpu_common.h>
#include <algo/blast/core/blast_aalookup.h>

#include <cuda.h>
#include <cuda_runtime.h>

#include <stdio.h>
#include <string.h>
#include <math.h>


void  GPU_BLASTP(const BLAST_SequenceBlk* query, const BlastQueryInfo* query_info,
                 const LookupTableWrap* lookup_wrap,
                 const Blast_ExtendWord* ewp, const BlastInitialWordParameters* word_params,
                 const BlastGPUOptions* gpu_options, const int* Sequence_length_vector,
                 const int h_Database2Dpadded_bytes, const char* h_Database2Dpadded,
                 Int4** d_Hits, const Int4 h_Hits_bytes, GPUBlastInitHitList** d_GPUBlastInitHitList, const Int4 h_GPUBlastInitHitList_bytes,
                 Int4** d_Database2Dpadded, Int4** d_RepeatedSubstitutionMatrix, Int4** d_RepeatedSequence_length_vector,
                 PV_ARRAY_TYPE** d_RepeatedPV, AaLookupSmallboneCell** d_ThickBackbone, Uint2** d_overflow, Int4** d_RepeatedDiag_array
                 )
{

    const Int4 num_threadsy = 1, num_blocksy = 1;
    const Int4 num_blocksx = gpu_options->num_blocksx, num_threadsx = gpu_options->num_threadsx,  num_sequences = gpu_options->num_sequences_to;

    if( GPU_VERBOSE )
        fprintf(stderr,"threadId = 0: Number of processed sequences by the GPU = %5d\n", num_sequences);

    //Group_number is the number of groups that the sequences are divided into, and is determined by the number of threads, since each thread processes one sequence at a time.
    //However, the number of threads iterate to cover all the sequences in the case that the number of threads is smaller than the number of sequences.
    int Group_number = (num_sequences / (num_threadsx*num_blocksx) ) + (( (num_sequences % (num_threadsx*num_blocksx) ) == 0)?0:1);

    int Sequence_length_vector_stride = ( Group_number / PITCH + ((Group_number % PITCH)?1:0) ) * PITCH;
    int h_RepeatedSequence_length_vector_integers = Sequence_length_vector_stride * num_blocksx;
    unsigned long int h_RepeatedSequence_length_vector_bytes = Sequence_length_vector_stride * num_blocksx * sizeof(int);
    int* h_RepeatedSequence_length_vector = (int*) calloc( h_RepeatedSequence_length_vector_integers, sizeof(int)  );
    MultipleCopyMaxlength(h_RepeatedSequence_length_vector, Sequence_length_vector, Group_number, Sequence_length_vector_stride, num_blocksx);

    BlastAaLookupTable *lookup = (BlastAaLookupTable*) lookup_wrap->lut;

    PV_ARRAY_TYPE *h_pv = h_pv = lookup->pv;
    int pv_length = (lookup->backbone_size >> PV_ARRAY_BTS) + 1;

    Round2Multiple(&pv_length);

    unsigned long int h_RepeatedPV_bytes = num_blocksx * pv_length * sizeof(PV_ARRAY_TYPE);
    PV_ARRAY_TYPE* h_RepeatedPV = (PV_ARRAY_TYPE*) calloc( h_RepeatedPV_bytes, sizeof(char)); 
    MultipleCopyPV(h_RepeatedPV, h_pv, pv_length, num_blocksx);

    AaLookupSmallboneCell* h_ThickBackbone = (AaLookupSmallboneCell*) lookup->thick_backbone;

    Uint2* h_overflow = (Uint2*)lookup->overflow;

    unsigned long int h_ThickBackbone_bytes = lookup->backbone_size * sizeof(AaLookupSmallboneCell);
    unsigned long int h_overflow_bytes = lookup->overflow_size * sizeof(Uint2);

    unsigned int diag_array_length = ewp->diag_table->diag_array_length;
    int* h_diag_array = (int*) calloc(diag_array_length, sizeof(int));
    unsigned long long int h_RepeatedDiag_array_bytes = (unsigned long long int)diag_array_length * num_blocksx*num_threadsx*sizeof(int);
    int diag_mask = ewp->diag_table->diag_mask;
    int window = ewp->diag_table->window;
    int diag_offset = ewp->diag_table->offset;



    int Query_length = query->length;
    unsigned char* h_Query = query->sequence;

    int* h_x_dropoff = (int*)malloc(query_info->num_queries * sizeof(int) );
    int* h_cutoff_score = (int*)malloc(query_info->num_queries * sizeof(int) );
    int* h_query_offset = (int*)malloc(query_info->num_queries * sizeof(int) );

    for(int i = 0; i < query_info->num_queries; ++i)
    {
        h_x_dropoff[i] = word_params->cutoffs[i].x_dropoff;
        h_cutoff_score[i] = word_params->cutoffs[i].cutoff_score;
        h_query_offset[i] = query_info->contexts[i].query_offset;
    }

    //Read Substitution Matrix
    char* h_SubstitutionMatrix = (char*) calloc( SUBSTITUTION_MATRIX_LENGTH * sizeof(char), sizeof(char) );
    unsigned int h_RepeatedSubstitutionMatrix_bytes = SUBSTITUTION_MATRIX_LENGTH * num_blocksx * sizeof(char);
    char* h_RepeatedSubstitutionMatrix = (char*) calloc ( h_RepeatedSubstitutionMatrix_bytes, sizeof(char) );
    ReadSubstitutionMatrix( h_SubstitutionMatrix, h_RepeatedSubstitutionMatrix, SUBSTITUTION_MATRIX_LENGTH, num_blocksx );

    if( GPU_VERBOSE )
        fprintf(stderr,"threadId = 0: Global memory allocated memory size = %lf Mb\n",(double)(h_Database2Dpadded_bytes +
                                                                                               h_RepeatedSubstitutionMatrix_bytes +
                                                                                               h_RepeatedSequence_length_vector_bytes +
                                                                                               h_Hits_bytes +
                                                                                               h_RepeatedPV_bytes+
                                                                                               h_ThickBackbone_bytes+
                                                                                               h_overflow_bytes+
                                                                                               h_RepeatedDiag_array_bytes +
                                                                                               h_GPUBlastInitHitList_bytes)/1e6);


    unsigned long long int BytesTransferredToGlobalMemory =
        h_Database2Dpadded_bytes +
        h_RepeatedSubstitutionMatrix_bytes +
        h_RepeatedSequence_length_vector_bytes +
        h_RepeatedPV_bytes+
        h_ThickBackbone_bytes+
        h_overflow_bytes;

    if( GPU_VERBOSE )
    {
        fprintf(stderr,"\nthreadId = 0: Bytes transferred to global memory\n");
        fprintf(stderr,"\n\nthreadId = 0: Bytes allocated on host\n");
        fprintf(stderr,"threadId = 0: Database           = %d bytes\n",h_Database2Dpadded_bytes );
        fprintf(stderr,"threadId = 0: SubstitutionMatrix = %d bytes\n",h_RepeatedSubstitutionMatrix_bytes );
        fprintf(stderr,"threadId = 0: SequenceMaxlength  = %d bytes\n",h_RepeatedSequence_length_vector_bytes );
        fprintf(stderr,"threadId = 0: Hits               = %d bytes\n",h_Hits_bytes );
        fprintf(stderr,"threadId = 0: PV                 = %d bytes\n",h_RepeatedPV_bytes );
        fprintf(stderr,"threadId = 0: ThickBackbone      = %d bytes\n",h_ThickBackbone_bytes );
        fprintf(stderr,"threadId = 0: Overflow           = %d bytes\n",h_overflow_bytes );
        fprintf(stderr,"threadId = 0: RepeatedDiagArray  = %d bytes\n",h_RepeatedDiag_array_bytes) ;
        fprintf(stderr,"threadId = 0: InitHitList        = %d bytes\n",h_GPUBlastInitHitList_bytes ) ;

        fprintf(stderr,"\nthreadId = 0: Bytes transferred to global memory\n");
        fprintf(stderr,"threadId = 0: Database           = %10.1f Mb\n",(float)h_Database2Dpadded_bytes / 1e6);
        fprintf(stderr,"threadId = 0: SubstitutionMatrix = %10.1f Mb\n",(float)h_RepeatedSubstitutionMatrix_bytes / 1e6);
        fprintf(stderr,"threadId = 0: SequenceMaxlength  = %10.1f Mb\n",(float)h_RepeatedSequence_length_vector_bytes / 1e6);
        fprintf(stderr,"threadId = 0: PV                 = %10.1f Mb\n",(float)h_RepeatedPV_bytes / 1e6);
        fprintf(stderr,"threadId = 0: ThickBackbone      = %10.1f Mb\n",(float)h_ThickBackbone_bytes / 1e6);
        fprintf(stderr,"threadId = 0: Overflow           = %10.1f Mb\n",(float)h_overflow_bytes / 1e6);
        fprintf(stderr,"threadId = 0: Total              = %10.1f Mb\n",(double)BytesTransferredToGlobalMemory/1e6);
    }
    
    if( GPU_VERBOSE )
        fprintf(stderr,"threadId = 0: Allocating GPU memory...\n");

    cudaMalloc( (void**) &(*d_Database2Dpadded), h_Database2Dpadded_bytes ) ;
    cudaMalloc( (void**) &(*d_RepeatedSubstitutionMatrix), h_RepeatedSubstitutionMatrix_bytes ) ;
    cudaMalloc( (void**) &(*d_RepeatedSequence_length_vector), h_RepeatedSequence_length_vector_bytes ) ;
    cudaMalloc( (void**) &(*d_Hits), h_Hits_bytes ) ;
    cudaMalloc( (void**) &(*d_RepeatedPV), h_RepeatedPV_bytes ) ;
    cudaMalloc( (void**) &(*d_ThickBackbone), h_ThickBackbone_bytes ) ;
    cudaMalloc( (void**) &(*d_overflow), h_overflow_bytes ) ;
    cudaMalloc( (void**) &(*d_RepeatedDiag_array), h_RepeatedDiag_array_bytes ) ;
    cudaMemset( *d_RepeatedDiag_array, 0, h_RepeatedDiag_array_bytes) ;
    cudaMalloc( (void**) &(*d_GPUBlastInitHitList), h_GPUBlastInitHitList_bytes ) ;
    cudaMemset( *d_GPUBlastInitHitList,0, h_GPUBlastInitHitList_bytes ) ;
    
    if( GPU_VERBOSE )
        fprintf(stderr,"threadId = 0: Done with allocating GPU memory\n");

    // copy host memory to device
    if( GPU_VERBOSE )
        fprintf(stderr,"threadId = 0: Copying to GPU memory...\n");

    cudaMemcpy( *d_Database2Dpadded, h_Database2Dpadded, h_Database2Dpadded_bytes, cudaMemcpyHostToDevice) ;
    cudaMemcpy( *d_RepeatedSubstitutionMatrix, h_RepeatedSubstitutionMatrix, h_RepeatedSubstitutionMatrix_bytes, cudaMemcpyHostToDevice) ;
    cudaMemcpy( *d_RepeatedSequence_length_vector, h_RepeatedSequence_length_vector, h_RepeatedSequence_length_vector_bytes, cudaMemcpyHostToDevice) ;
    cudaMemcpy( *d_RepeatedPV, h_RepeatedPV, h_RepeatedPV_bytes, cudaMemcpyHostToDevice) ;
    cudaMemcpy( *d_ThickBackbone, h_ThickBackbone, h_ThickBackbone_bytes, cudaMemcpyHostToDevice) ;
    cudaMemcpy( *d_overflow, h_overflow, h_overflow_bytes, cudaMemcpyHostToDevice) ;

    if( GPU_VERBOSE )
        fprintf(stderr,"threadId = 0: Done with copying GPU memory\n");

    // setup execution parameters
    dim3  grid( num_blocksx, num_blocksy, 1);
    dim3  threads( num_threadsx, num_threadsy, 1);
    int mem_size = 0;

    mem_size = pv_length*sizeof(PV_ARRAY_TYPE) + SUBSTITUTION_MATRIX_LENGTH*sizeof(char) + Sequence_length_vector_stride*sizeof(int) + 2 * num_threadsx * 4 * sizeof(char);
    cudaMemcpyToSymbol( d_Query_const, h_Query, Query_length ) ;
    cudaMemcpyToSymbol( d_x_dropoff_const, h_x_dropoff, query_info->num_queries*sizeof(int) );
    cudaMemcpyToSymbol( d_cutoff_score_const, h_cutoff_score, query_info->num_queries*sizeof(int) );
    cudaMemcpyToSymbol( d_query_offset_const, h_query_offset, query_info->num_queries*sizeof(int) );

    cudaDeviceProp deviceProp;
    cudaGetDeviceProperties(&deviceProp, 0);
    if(mem_size > deviceProp.sharedMemPerBlock)
        fprintf(stderr,"ERROR: Shared memory requirements (%d bytes per block) bigger than available memory (16384 bytes)!!\n",mem_size);

    // execute the kernel
    if( GPU_VERBOSE )
        fprintf(stderr,"threadId = 0: Kernel execution...\n");

    GPU_BLASTP_kernelTwoHit<<< grid, threads, mem_size >>>( *d_Database2Dpadded,
                                                            *d_RepeatedSubstitutionMatrix,
                                                            SUBSTITUTION_MATRIX_LENGTH,
                                                            *d_RepeatedPV,
                                                            *d_ThickBackbone, *d_overflow,
                                                            *d_RepeatedSequence_length_vector,
                                                            Sequence_length_vector_stride,
                                                            num_sequences, pv_length,
                                                            Query_length, diag_mask, diag_array_length,
                                                            query_info->num_queries,
                                                            lookup->word_length, lookup->charsize,
                                                            lookup->mask, window, diag_offset,
                                                            *d_RepeatedDiag_array, *d_Hits, *d_GPUBlastInitHitList );
    if( GPU_VERBOSE )
        fprintf(stderr,"threadId = 0: Done with kernel execution\n");

    free(h_RepeatedSequence_length_vector);
    free(h_RepeatedPV);
    free(h_diag_array);
    free(h_x_dropoff);
    free(h_cutoff_score);
    free(h_query_offset);
    free(h_SubstitutionMatrix);
    free(h_RepeatedSubstitutionMatrix);

}

void GPU_BLASTP_get_data(const Int4* h_Hits,
                         Int4* d_Hits,
                         const Int4 h_Hits_bytes,
                         GPUBlastInitHitList* h_GPUBlastInitHitList,
                         GPUBlastInitHitList* d_GPUBlastInitHitList,
                         const Int4 h_GPUBlastInitHitList_bytes,
                         const Int4 num_sequences)
{


    // copy result from device to host
    if( GPU_VERBOSE )
        fprintf(stderr,"threadId = 0: Transferring data from GPU...\n");

    cudaMemcpy( (void*) h_Hits, d_Hits, h_Hits_bytes, cudaMemcpyDeviceToHost ) ;
    cudaMemcpy( h_GPUBlastInitHitList, d_GPUBlastInitHitList, h_GPUBlastInitHitList_bytes, cudaMemcpyDeviceToHost ) ;

    if( GPU_VERBOSE )
        fprintf(stderr,"threadId = 0: Done transferring data from GPU\n");

    if( GPU_VERBOSE )
    {
        long unsigned int Total_Hits = 0, Total_Ext = 0, Total_Succ_Ext = 0;
        for( int Sequence = 0; Sequence < num_sequences; ++Sequence )
        {
            Total_Hits += h_Hits[ Sequence*3 + 0 ];
            Total_Succ_Ext += h_Hits[ Sequence*3 + 1 ];
            Total_Ext += h_Hits[ Sequence*3 + 2 ];

//            fprintf(stderr,"GPU:Sequence = %6d: hits = %4d, Ext = %5d, Succ Ext = %5d\n",
//		Sequence, h_Hits[Sequence*3 ], h_Hits[Sequence*3 + 1], h_Hits[Sequence*3 + 2]);
//		      printf("GPU:subject->oid = %7d, hits = %3d, hits_extended = %3d, hits_succ_extended = %3d\n",
//			Sequence, h_Hits[Sequence*3 ], h_Hits[Sequence*3 + 1], h_Hits[Sequence*3 + 2]);


        }

        fprintf(stderr,"\nTotal_Hits = %lu, Total_Ext = %lu, Total_Succ_Ext = %lu\n",Total_Hits,Total_Ext, Total_Succ_Ext);

    }

}


void GPU_BLASTP_free_memory(Int4** d_Hits, GPUBlastInitHitList** d_GPUBlastInitHitList,
                            Int4** d_Database2Dpadded, Int4** d_RepeatedSubstitutionMatrix, Int4** d_RepeatedSequence_length_vector,
                            PV_ARRAY_TYPE** d_RepeatedPV, AaLookupSmallboneCell** d_ThickBackbone, Uint2** d_overflow, Int4** d_RepeatedDiag_array)
{

    cudaFree( *d_Hits ) ;
    cudaFree( *d_GPUBlastInitHitList ) ;
    cudaFree( *d_Database2Dpadded );
    cudaFree( *d_RepeatedSubstitutionMatrix );
    cudaFree( *d_RepeatedSequence_length_vector );
    cudaFree( *d_RepeatedPV );
    cudaFree( *d_ThickBackbone );
    cudaFree( *d_overflow );
    cudaFree( *d_RepeatedDiag_array );

    cudaThreadExit();

}

Boolean GPU_BLASTP_check_memory(const LookupTableWrap* lookup_wrap,
                                const Blast_ExtendWord* ewp,
                                const BlastGPUOptions* gpu_options,
                                const int h_Database2Dpadded_bytes,
                                const Int4 h_Hits_bytes, const Int4 h_GPUBlastInitHitList_bytes,
                                const Int4 Group_number, const Int4 num_queries,
                                const Int4 query_length
                                )
{

    if( num_queries > NUM_QUERIES_MAX )
        return FALSE;
    if( query_length > QUERY_LENGTH_MAX )
        return FALSE;

    int h_RepeatedSubstitutionMatrix_bytes = SUBSTITUTION_MATRIX_LENGTH * (gpu_options->num_blocksx) * sizeof(char);

    int Sequence_length_vector_stride = ( Group_number / PITCH + ((Group_number % PITCH)?1:0) ) * PITCH;
    long int h_RepeatedSequence_length_vector_bytes = Sequence_length_vector_stride * (gpu_options->num_blocksx) * sizeof(int);

    BlastAaLookupTable *lookup = (BlastAaLookupTable*) lookup_wrap->lut;

    int pv_length = (lookup->backbone_size >> PV_ARRAY_BTS) + 1;
    Round2Multiple(&pv_length);
    int h_RepeatedPV_bytes = (gpu_options->num_blocksx) * pv_length * sizeof(PV_ARRAY_TYPE);

    int h_ThickBackbone_bytes = lookup->backbone_size * sizeof(AaLookupSmallboneCell);

    int h_overflow_bytes = lookup->overflow_size * sizeof(Uint2);

    int diag_array_length = ewp->diag_table->diag_array_length;
    unsigned long long int h_RepeatedDiag_array_bytes =
        (unsigned long long int)diag_array_length * (gpu_options->num_blocksx)*(gpu_options->num_threadsx)*sizeof(int);


    unsigned long long int Total_global_memory_bytes  =
        h_Database2Dpadded_bytes +
        h_RepeatedSubstitutionMatrix_bytes +
        h_RepeatedSequence_length_vector_bytes +
        h_Hits_bytes +
        h_RepeatedPV_bytes +
        h_ThickBackbone_bytes +
        h_overflow_bytes +
        h_RepeatedDiag_array_bytes +
        h_GPUBlastInitHitList_bytes;

    cudaDeviceProp deviceProp;
    cudaGetDeviceProperties(&deviceProp, 0);
    unsigned long long int GPU_global_memory = (unsigned long long int) deviceProp.totalGlobalMem;
    if( Total_global_memory_bytes > GPU_global_memory) {
        return FALSE;
    }

    return TRUE;
}

Boolean GPU_BLAST_check_availability()
{
	int deviceCount = 0;
	if (cudaGetDeviceCount(&deviceCount) != cudaSuccess)
		return FALSE;
	if (deviceCount == 0)
		return FALSE;
	return TRUE;
}
