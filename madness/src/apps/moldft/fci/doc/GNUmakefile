
SOURCE   = matrixelements.tex \
           mspingraph.eps  orbitalgraph.eps spingraph.eps \
           orbitalarcs.eps spinarcs.eps

matrixelements.ps:	$(SOURCE)
	latex matrixelements.tex
	latex matrixelements.tex
	latex matrixelements.tex
	dvips -D 300 -N0 -o matrixelements.ps matrixelements.dvi	

clean:
	/bin/rm -f *toc *aux *log *dvi matrixelements.ps *~ core
