## Taken from: https://www.oreilly.com/library/view/c-cookbook/0596007612/ch01s17.html
# Specify extensions of files to delete when cleaning
CLEANEXTS   = o a

# Specify the target file and the install directory
OUTPUTFILE  = libzf_log.a
INSTALLDIR  = .

.PHONY: all

all: $(OUTPUTFILE)

# Build mmio.a from zf_log.o
$(OUTPUTFILE): zf_log.o
	ar r $@ $^
	ranlib $@
	rm -f *.o

# No rule to build zf_log.o from .c files is required;
# this is handled by make's database of implicit rules

clean:
	for file in $(CLEANEXTS); do rm -f *.$$file; done

# Indicate dependencies of .c files on .h files
zf_log.o: zf_log.h
