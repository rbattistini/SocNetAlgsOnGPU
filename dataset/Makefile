## Makefile to fetch and install graph data
#
# Based on the makefile used to build the dataset in Gunrock
#
# Targets defined are:
#
# make
# make fetch
# make clean
# make veryclean

# Each graph lives in its own directory
SUBDIRS := ca-AstroPh ca-GrQc ca-HepPh ca-HepTh Erdos02 as-22july06 USpowerGrid delaunay_n14

.PHONY: $(GRAPHS)

# Fetches all graphs, extracts and sets up files for test
all: recurse

# Only download the graphs, but do not proceed further
fetch: recurse

# Clean everything except the downloaded graphs
clean: recurse

veryclean: recurse

# Recurse into each subdirectory and try to build the provided targets
recurse:
	for subdir in $(SUBDIRS); do $(MAKE) -C $$subdir $(MAKECMDGOALS); done

$(GRAPHS):
	$(MAKE) -C $@
