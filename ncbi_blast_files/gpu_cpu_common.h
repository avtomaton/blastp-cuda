
#ifndef GPU_CPU_COMMON_H
#define GPU_CPU_COMMON_H

//#include <algo/blast/api/blast_api.h>
#include <algo/blast/core/blast_aalookup.h>
#include <algo/blast/core/blast_nalookup.h>
#include <algo/blast/core/lookup_wrap.h>
#include <algo/blast/core/blast_extend.h>

#define MAX_SUCCESSFUL_HITS_PER_SEQUENCE 2
#define MOD_MAX_SUCCESSFUL_HITS_PER_SEQUENCE (MAX_SUCCESSFUL_HITS_PER_SEQUENCE - 1)
#define NUM_SEQUENCES                   95
#define PERCENTAGE                      95
#define NUM_THREADSX                    64
#define NUM_THREADSX_MAX               512
#define NUM_BLOCKSX                    512
#define NUM_BLOCKSX_MAX              65535
#define GPU                              1
#define GPU_VERBOSE                      0
#define METHOD                           1
#define PITCH                           16
#define CPU_THREADS                      1

#define QUERY_LENGTH_MAX             16384
#define NUM_QUERIES_MAX               4096

#define SUBSTITUTION_MATRIX_LENGTH     640
typedef struct GPUBlastInitHitList {

  int hsp_q;
  int hsp_s;
  int query_offset;
  int subject_offset;
  int hsp_len;
  int score;

} GPUBlastInitHitList;

#ifdef __cplusplus
extern "C" {
#endif

void GPU_BLASTP(const BLAST_SequenceBlk* query,
		const BlastQueryInfo* query_info,
		const LookupTableWrap* lookup_wrap,
		const Blast_ExtendWord* ewp,
		const BlastInitialWordParameters* word_params,
		const BlastGPUOptions* gpu_options,
		const int* Sequence_length_vector,
		const int h_Database2Dpadded_bytes,
		const char* h_Database2Dpadded,
		Int4** d_Hits,
		const Int4 h_Hits_bytes,
		GPUBlastInitHitList** d_GPUBlastInitHitList,
		const Int4 h_GPUBlastInitHitList_bytes,
		Int4** d_Database2Dpadded,
		Int4** d_RepeatedSubstitutionMatrix,
		Int4** d_RepeatedSequence_length_vector,
		PV_ARRAY_TYPE** d_RepeatedPV,
		AaLookupSmallboneCell** d_ThickBackbone,
		Uint2** d_overflow,
		Int4** d_RepeatedDiag_array);

void GPU_BLASTP_get_data(const Int4* h_Hits,
			 Int4* d_Hits,
			 const Int4 h_Hits_bytes,
			 GPUBlastInitHitList* h_GPUBlastInitHitList,
			 GPUBlastInitHitList* d_GPUBlastInitHitList,
			 const Int4 h_GPUBlastInitHitList_bytes,
			 Int4 NumSequences);

void GPU_BLASTP_free_memory(Int4** d_Hits,
			    GPUBlastInitHitList** d_GPUBlastInitHitList,
			    Int4** d_Database2Dpadded,
			    Int4** d_RepeatedSubstitutionMatrix,
			    Int4** d_RepeatedSequence_length_vector,
			    PV_ARRAY_TYPE** d_RepeatedPV,
			    AaLookupSmallboneCell** d_ThickBackbone,
			    Uint2** d_overflow,
			    Int4** d_RepeatedDiag_array);

Boolean GPU_BLASTP_check_memory(const LookupTableWrap* lookup_wrap,
				const Blast_ExtendWord* ewp,
				const BlastGPUOptions* gpu_options,
				const int h_Database2Dpadded_bytes,
				const Int4 h_Hits_bytes,
				const Int4 h_GPUBlastInitHitList_bytes,
				const Int4 Group_number,
				const Int4 num_queries,
				const Int4 query_length);

Boolean GPU_BLAST_check_availability();

void Round2Multiple(Int4* MaxLength);


#ifdef __cplusplus
}
#endif

#endif /* GPU_CPU_COMMON_H */
