#The following variables must be defined prior to including this
#makefile fragment
#
#GRAPH_URL:  the url path to the file

GRAPH_FILE := $(notdir $(GRAPH_URL))

all: setup

fetch: $(GRAPH_FILE)

$(GRAPH_FILE):
	wget -N $(GRAPH_URL)
