#ifndef GPU_BLASTP_H
#define GPU_BLASTP_H

#define WIDTH_MULTIPLE 64

/* --------------- protein blast defines -----------------------------*/

typedef uint PV_ARRAY_TYPE;


void ReadSubstitutionMatrix(char *SubstitutionMatrix, char *h_RepeatedSubstitutionMatrix, int SubstitutionMatrix_length, int num_blocksx);
void MultipleCopyPV(PV_ARRAY_TYPE* h_RepeatedPV, PV_ARRAY_TYPE* h_pv, int pv_length, int num_blocksx);
void MultipleCopyMaxlength(int* h_Repeated_SequenceMaxlength_vector, const int* Sequence_Maxlength_vector, int Group_number, int stride, int NumSequences);

#endif /*GPU_BLASTP_H */


