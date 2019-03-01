DOCX_REF=my_reference.docx
DOCX_NAME=paper.docx
CSL_FILE=aps.csl
IMG_PATH=img
BIB_FILE=ref.bib
TEX_FILE=paper.tex
SI_FILE=SI.tex
REMOTE = ~/Dropbox/Ion\ rectification/TeX/
GS_TAG= -sDEVICE=pdfwrite -dQUIET -sBATCH -dNOPAUSE
PDFLATEX_TAG=-interaction=nonstopmode -draftmode
LATEXMK_TAG= -f -pdf -quiet -view=none -pdflatex='pdflatex -interaction=nonstopmode'
DIFF_TAG= --exclude-textcmd="section,subsection" --config="PICTUREENV=(?:figure|section|DIFnomarkup)[*]*" --graphics-markup=0 --disable-citation-markup



all: main SI

convert:
	python convert_svg.py
main:
	pdflatex $(PDFLATEX_TAG) $(SI_FILE)
	latexmk $(LATEXMK_TAG) $(TEX_FILE)

SI:
	pdflatex $(PDFLATEX_TAG) $(TEX_FILE)
	latexmk $(LATEXMK_TAG) $(SI_FILE)

compress:
	gs $(GS_TAG) -sOutputFile="paper-compak.pdf" paper.pdf
	gs $(GS_TAG) -sOutputFile="SI-compak.pdf" SI.pdf

doc:
	pandoc --from=latex --to=docx $(TEX_FILE) -o $(DOCX_NAME) \
		--reference-doc=$(DOCX_REF) \
		--bibliography=$(BIB_FILE) \
		--csl=$(CSL_FILE) \
		$(VERBOSE)


diff:
	cp ./paper.tex ./paper.tex.bak
	git checkout 6f1f33 ./paper.tex #old version
	cp paper.tex paper_old.tex
	cp ./SI.tex ./SI.tex.bak	
	git checkout 6f1f33 ./SI.tex
	cp SI.tex SI_old.tex
	git checkout HEAD ./paper.tex
	git checkout HEAD ./SI.tex
	cp ./paper.tex.bak paper.tex #restore to current mode
	cp ./SI.tex.bak SI.tex
	latexdiff $(DIFF_TAG) paper_old.tex paper.tex > paper_change.tex
	latexdiff $(DIFF_TAG) SI_old.tex SI.tex > SI_change.tex
	latexmk $(LATEXMK_TAG) paper_change.tex
	latexmk $(LATEXMK_TAG) SI_change.tex

transfer:
	cp $(TEX_FILE) $(SI_FILE) $(BIB_FILE) $(REMOTE) 
	cp paper.pdf paper_change.pdf $(REMOTE)
	cp SI.pdf SI_change.pdf $(REMOTE)
	mkdir -p $(REMOTE)/img
	cp img/fig* img/SI* $(REMOTE)/img

clean:
	latexmk -c


clean-all:
	clean
	latexmk -C
	rm -f $(IMG_PATH)/*.pdf	#clean up pdf
