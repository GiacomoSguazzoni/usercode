# makefile for pdflatex
MAIN=TSnote
TEXSOURCES= $(wildcard *.tex)
PLOTFILES= $(wildcard fig/*.png fig/*.pdf)
STYLEFILES= $(wildcard *.sty *.cls )

all: $(MAIN).pdf
$(MAIN).pdf: $(TEXSOURCES) $(STYLEFILES) $(PLOTFILES)
	pdflatex $(MAIN)
	pdflatex $(MAIN)
force:
	pdflatex $(MAIN)
	pdflatex $(MAIN)

clean:
	rm -f $(MAIN).log $(MAIN).pdf $(MAIN).aux
