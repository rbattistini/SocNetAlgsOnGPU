## I principali target definiti da questo makefile sono:
##
# tag::maketargets[]
# make			compila tutti i sorgenti disponibili
# make clean   	cancella i file temporanei e gli eseguibili
# make graphs	compila i grafici dai .csv con gnuplot
# make clgraphs	rimuove i file .png generati
# make context	compila un pdf con context
# make docs		compila un pdf con asciidoc
# make cldocs	rimuove il pdf generato
# end::maketargets[]

SCRIPT_DIR := ./scripts
DATA_DIR   := ./datasets
SRC_DIR    := ./src
DOCS_DIR   := ./docs
IMG_DIR    := $(DOCS_DIR)/images
EXE_CUDA   := $(basename $(wildcard $(SRC_DIR)/*.cu))
EXE_CPP    := $(basename $(wildcard $(SRC_DIR)/*.cpp))
EXE		   := $(EXE_CUDA) $(EXE_CPP)
NVCC       =  nvcc
CPPC	   =  g++
NVCFLAGS   += -arch=sm_61 -DCUDA_DEBUG -DDEBUG #-O3
CFLAGS     += -std=c99 -march=native -Wall -Wpedantic #-O3
NVLDLIBS   += -lm
LDLIBS     += -lm

ALL: $(EXE)

% : %.cpp
	$(CPPC) $(CPPFLAGS) $< $(CPPLIBS) -o $@

% : %.cu
	$(NVCC) $(NVCFLAGS) $< $(NVLDLIBS) -o $@

.PHONY: clean clgraphs graphs cldocs vcldocs docs

context: $(DOCS_DIR)/*.tex
	context $(DOCS_DIR)/*.tex

docs: $(DOCS_DIR)/*.adoc clgraphs
	asciidoctor-pdf -r asciidoctor-diagram -r asciidoctor-bibtex $(DOCS_DIR)/*.adoc

graphs: $(SCRIPT_DIR)/*.plot $(DATA_DIR)/*.dat
	gnuplot $(SCRIPT_DIR)/*.plot

cldocs: clgraphs
	rm -f $(DOCS_DIR)/*.log $(DOCS_DIR)/*.tuc

vcldocs: cldocs
	rm - f $(DOCS_DIR)/*.pdf $(DOCS_DIR)/.asciidoc/*

clgraphs:
	\rm -f $(IMG_DIR)/*.png $(IMG_DIR)/*.svg

clean:
	\rm -f $(EXE) *.o *~
