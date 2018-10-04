IMG_PATH=img
BIB_FILE=ref.bib
TEX_FILE=paper.tex
SI_FILE=SI.tex
REMOTE = ~/Dropbox/Ion\ rectification/TeX/
GS_TAG= -sDEVICE=pdfwrite -dQUIET -sBATCH -dNOPAUSE




all: main SI compress

main: $(TEX_FILE) $(BIB_FILE)
	latexmk -f -pdf -quiet -view=none -pdflatex='pdflatex -interaction=nonstopmode' $(TEX_FILE)

SI: $(SI_FILE) $(BIB_FILE) 
	latexmk -f -pdf -quiet -view=none -pdflatex='pdflatex -interaction=nonstopmode' $(SI_FILE)

compress:
	gs $(GS_TAG) -sOutputFile="paper-compak.pdf" paper.pdf
	gs $(GS_TAG) -sOutputFile="SI-compak.pdf" SI.pdf

transfer:
	cp $(TEX_FILE) $(SI_FILE) $(BIB_FILE) $(REMOTE) 
	cp paper-compak.pdf $(REMOTE)/paper.pdf
	cp SI-compak.pdf $(REMOTE)/SI.pdf
	mkdir -p $(REMOTE)/img
	cp img/fig* img/SI* $(REMOTE)/img

Clean:
	latexmk -c


clean-all:
	clean
	latexmk -C
	rm -f $(IMG_PATH)/*.pdf	#clean up pdf
