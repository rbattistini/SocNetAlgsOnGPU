## Taken from: https://www.oreilly.com/library/view/c-cookbook/0596007612/ch01s17.html
# Specify extensions of files to delete when cleaning
CLEANEXTS   = o a

# Specify the target file and the install directory
OUTPUTFILE  = libmmio.a
INSTALLDIR  = .

.PHONY: all

all: $(OUTPUTFILE)

# Build mmio.a from mmio.o
$(OUTPUTFILE): mmio.o
	ar r $@ $^
	ranlib $@
	rm -f *.o

# No rule to build mmio.o from .c files is required;
# this is handled by make's database of implicit rules

clean:
	for file in $(CLEANEXTS); do rm -f *.$$file; done

# Indicate dependencies of .c files on .h files
mmio.o: mmio.h
