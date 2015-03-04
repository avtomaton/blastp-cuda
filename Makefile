.SUFFIXES : .cu 

CC	:= gcc
CXX 	:= g++
NVCC	:= $(NVCC_PATH)

COMMONFLAGS := -c
# the "--no-align-double" is necessary. Without it the nvcc compiler aligns structure sizes every 8 bytes which is 
#incompatible with some NCBI_BLAST data structures, i.e. the "BlastContextInfo"
NVCCFLAGS :=  --compiler-options -fno-strict-aliasing --compiler-options -fno-inline -DUNIX -O3 --no-align-double 

all: libgpublast.a

libgpublast.a : gpu_blastp.o gpu_blastp.cu.o gpu_blastp_kernel.h 
	ar cru libgpublast.a gpu_blastp.o gpu_blastp.cu.o
	ranlib libgpublast.a	
	cp libgpublast.a                          $(PROJ_DIR)/lib/
	cp $(CUDA_LIB)/libcudart.so               $(PROJ_DIR)/lib/

gpu_blastp.o : gpu_blastp.c gpu_blastp.h gpu_blastp_kernel.cu 
	$(CXX) $(COMMONFLAGS) gpu_blastp.c -o gpu_blastp.o -I$(CXX_DIR)/include/ -I$(PROJ_DIR)/inc/

gpu_blastp.cu.o : gpu_blastp.cu gpu_blastp.h gpu_blastp_kernel.cu gpu_blastp_kernel.h 
	$(NVCC) $(COMMONFLAGS) $(NVCCFLAGS) gpu_blastp.cu -o gpu_blastp.cu.o -I$(CXX_DIR)/include/ -I$(PROJ_DIR)/inc/


