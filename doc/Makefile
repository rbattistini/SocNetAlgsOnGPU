## Targets defined are:
#
# make			compile all source code files available
# make clean   	remove tmp files and executables
# make graphs	compila i grafici dai .csv con gnuplot
# make clgraphs	rimuove i file .png generati
# make doc		compila un pdf con asciidoc
# make cldocs	rimuove il pdf generato

SCRIPT_DIR := ../scripts
DATA_DIR   := ../datasets
IMG_DIR    := ./images
DOCS       := clgraphs graphs cldocs docs

.PHONY: $(DOCS)

docs: ./*.adoc clgraphs
	asciidoctor-pdf -r asciidoctor-diagram -r asciidoctor-bibtex \
 -r asciidoctor-mathematical ./*.adoc \
 -a pdf-themesdir=. -a pdf-theme=custom \
 -a pdf-fontsdir="fonts;GEM_FONTS_DIR"

graphs: $(SCRIPT_DIR)/*.plot $(DATA_DIR)/*.dat
	gnuplot $(SCRIPT_DIR)/*.plot

cldocs: clgraphs
	rm - f ./*.pdf ./.asciidoc/*

clgraphs:
	\rm -f $(IMG_DIR)/*.png $(IMG_DIR)/*.svg
