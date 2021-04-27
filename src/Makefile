## Targets defined are:
#
# make			compile all source code files available
# make clean   	remove tmp files and executables

# Specify extensions of files to delete when cleaning
CLEANEXTS   = o a
LIBDIR 	   := ../libs
LIBDIRS    := $(LIBDIR)/libmmio
EXE_CUDA   := $(basename $(wildcard ./*.cu))
EXE_CPP    := $(basename $(wildcard ./*.cpp))
EXE		   := $(EXE_CUDA) $(EXE_CPP)
NVCC       =  nvcc
CPPC	   =  g++
NVCFLAGS   += -Wno-deprecated-gpu-targets -arch=sm_61 -DCUDA_DEBUG --include-path ../libs/libmmio #-O3
CPPFLAGS   += -march=native -Wall -g -Wpedantic -I../libs/libmmio #-O3
NVLDLIBS   += $(LIBDIRS)/libmmio.a
CPPLIBS    += $(LIBDIRS)/libmmio.a

.PHONY: all

ALL: $(EXE)

% : %.cpp
	$(CPPC) $(CPPFLAGS) $< -o $@ $(CPPLIBS)

% : %.cu
	$(NVCC) $(NVCFLAGS) $< -o $@ $(NVLDLIBS)

cuda: $(EXE_CUDA)

serial: $(EXE_CPP)

clean:
	for file in $(CLEANEXTS); do rm -f *.$$file; done

veryclean: clean
	rm -f $(EXE_CUDA) $(EXE_CPP)
