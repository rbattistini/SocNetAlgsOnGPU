## Targets defined are:
#
# make all		compile all source code files available
# make libs     compile static libraries in lib folder
# make clean   	remove tmp files and executables

PROJECT_NAME := sna_bc
SRC_DIR      := .
INC_DIR      := ../include
LIB_DIR      := ../lib

CPPCC        := g++
NVCC         := nvcc

NV_SRC       := $(wildcard $(SRC_DIR)/*.cu)
CPP_SRC      := $(wildcard $(SRC_DIR)/*.cpp)
OBJ_CPP      := $(patsubst %.cpp,%.o,$(CPP_SRC))
OBJ_CUDA     := $(patsubst %.cu,%.o,$(NV_SRC))

LDLIBS       := -lmmio -lzf_log
LDFLAGS      := -L$(LIB_DIR)/mmio -L$(LIB_DIR)/zf_log

CFLAGS       := -g

NVCFLAGS     := -arch=sm_61 -O3
CPPFLAGS     := -std=c++11
CPPFLAGS     += $(foreach includedir, $(INC_DIR), -I$(includedir))
CPPFLAGS     += -I../lib/zf_log -I../lib/mmio
ALLCFLAGS    :=  $(NVCFLAGS) $(CPPFLAGS) $(MPICFLAGS) $(CFLAGS)

all: $(PROJECT_NAME)

$(PROJECT_NAME): $(OBJ_CPP) $(OBJ_CUDA) libmmio.a libzf_log.a
	$(NVCC) $(ALLCFLAGS) $(LDFLAGS) $(LDLIBS) $(NV_SRC) $(OBJ_CPP) -o $(PROJECT_NAME)

$(OBJ_CUDA): $(NV_SRC)
	 $(NVCC) $(NVCFLAGS) -I../include -I../lib/zf_log -g -c $< -o $@

libs: libmmio.a libzf_log.a

libmmio.a:
	cd $(LIB_DIR)/mmio && $(MAKE)

libzf_log.a:
	cd $(LIB_DIR)/zf_log && $(MAKE)

.PHONY: all libs clean

clean:
	\rm -f $(EXE) *.o *~
