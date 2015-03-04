#ifndef GPU_BLASTP_KERNEL_H
#define GPU_BLASTP_KERNEL_H

typedef struct ExtendLeftReturn {

    int maxscore;
    int length;
    
} ExtendLeftReturn;


typedef struct ExtendRightReturn {

    int s_last_off;
    int maxscore;
    int length;
    
} ExtendRightReturn;


typedef struct ExtendTwoHitReturn {

    int s_last_off;
    int hsp_q;
    int hsp_s;
    int hsp_len;
    int MAX_score;
    int right_extend;
    
    int left_score;
    int left_disp;
    int right_disp;
    
} ExtendTwoHitReturn;


texture<int, 1, cudaReadModeElementType> texRef;
#define ALPHABET_SIZE  28
#define AA_HITS_PER_CELL 3 /**< maximum number of hits in one lookup table cell */

__global__ void
GPU_BLASTP_kernel_TwoHit( const int* d_Database2Dpadded,
                          const int *d_RepeatedSubstitutionMatrix,
                          const int SubstitutionMatrix_length,
                          const PV_ARRAY_TYPE* d_RepeatedPV,
                          const AaLookupSmallboneCell* d_ThickBackbone,
                          const unsigned short* d_overflow,
                          const int* d_RepeatedSequenceMaxlength_vector,
                          const int SequenceMaxlength_vector_stride,
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
                          GPUBlastInitHitList* d_GPUBlastInitHitList);


static ExtendTwoHitReturn
__device__ s_BlastAaExtendTwoHit(const char* matrix,
                                 const int* subject,
                                 const int subject_length,
                                 const char* query,
                                 const int query_length,
                                 const int s_left_off,
                                 const int s_right_off,
                                 const int q_right_off,
                                 const int dropoff,
                                 const int word_size);


__device__ ExtendLeftReturn s_BlastAaExtendLeft(const char* matrix,
                                                const int* subject,
                                                const char * query,
                                                const int s_off,
                                                const int q_off,
                                                const int dropoff,
                                                int maxscore);

__device__ ExtendRightReturn s_BlastAaExtendRight(const char* matrix,
                                                  const int* subject,
                                                  const char* query,
                                                  const int s_off,
                                                  const int q_off,
                                                  const int dropoff,
                                                  int maxscore,
                                                  const int Sequence_actual_length,
                                                  const int Query_length);

#endif /* GPU_BLASTP_KERNEL_H */
