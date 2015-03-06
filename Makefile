.SUFFIXES : .cu 

CC	:= gcc
CXX 	:= g++
NVCC	:= nvcc

COMMONFLAGS := -c
# the "--no-align-double" is necessary. Without it the nvcc compiler aligns structure sizes every 8 bytes which is 
#incompatible with some NCBI_BLAST data structures, i.e. the "BlastContextInfo"
NVCCFLAGS :=  --compiler-options -fno-strict-aliasing --compiler-options -fno-inline -DUNIX -O3 --no-align-double 

INCLUDES= -I$(BLAST_DIR)/include/ -I$(BUILD_DIR)/inc/ -I$(PATCH_DIR)/ncbi_blast_files

all: $(BUILD_DIR)/lib/libgpublast.a

$(BUILD_DIR)/lib/libgpublast.a : $(BUILD_DIR)/gpu_blastp.o $(BUILD_DIR)/gpu_blastp.cu.o gpu_blastp_kernel.h 
	ar cru $(BUILD_DIR)/lib/libgpublast.a $(BUILD_DIR)/gpu_blastp.o $(BUILD_DIR)/gpu_blastp.cu.o
	ranlib $(BUILD_DIR)/lib/libgpublast.a

$(BUILD_DIR)/gpu_blastp.o : gpu_blastp.c gpu_blastp.h gpu_blastp_kernel.cu 
	$(CXX) $(COMMONFLAGS) gpu_blastp.c -o $(BUILD_DIR)/gpu_blastp.o $(INCLUDES)

$(BUILD_DIR)/gpu_blastp.cu.o : gpu_blastp.cu gpu_blastp.h gpu_blastp_kernel.cu gpu_blastp_kernel.h
	$(NVCC) $(COMMONFLAGS) $(NVCCFLAGS) gpu_blastp.cu -o $(BUILD_DIR)/gpu_blastp.cu.o $(INCLUDES)


