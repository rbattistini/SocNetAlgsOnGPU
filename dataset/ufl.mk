## Taken from: https://github.com/gunrock/gunrock/blob/dynamic-graph/dataset/large/ufl.mk
# common make file fragment for ufl graph dataset
# just define GRAPH_NAME prior to including this fragment

GRAPH_TAR  = $(GRAPH_NAME).tar.gz

setup: $(GRAPH_NAME).mtx

$(GRAPH_NAME).mtx: $(GRAPH_TAR)
	tar xvfz $(GRAPH_TAR)
	cp $(GRAPH_NAME)/$(GRAPH_NAME).mtx $(GRAPH_NAME).mtx
	rm -rf $(GRAPH_NAME)

clean:
	rm $(GRAPH_NAME).mtx

veryclean: clean
	rm $(GRAPH_TAR)
