all: paper.pdf

paper.dvi: paper.tex bibliography.bib
	latex -shell-escape paper
	#bibtex paper
	latex -shell-escape paper
	latex -shell-escape paper

paper.pdf: paper.dvi
	dvipdf paper.dvi

clean:
	rm -f paper.aux paper.bbl paper.blg paper.dvi paper.log paper.pdf
