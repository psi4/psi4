all: post174-daily.pdf post174-daily.ps

post174-daily.pdf: post174-daily.ps
	ps2pdf post174-daily.ps

# the following dvips invocation was recommended on 
# http://www.acm.org/sigs/publications/sigfaq#a14
post174-daily.ps: post174-daily.dvi
	dvips -P pdf -t letter -o post174-daily.ps post174-daily.dvi

post174-daily.dvi: post174-daily.tex post174-daily.bib laplace.eps python.eps
	latex  post174-daily
	bibtex post174-daily
	latex  post174-daily
	latex  post174-daily
	latex  post174-daily

laplace.eps: laplace.png
	convert laplace.png laplace.eps

laplace.png: laplace.gnuplot laplace.dat
	gnuplot laplace.gnuplot

python.eps: python.png
	convert python.png python.eps

python.png: python.gnuplot python.dat
	gnuplot python.gnuplot

clean:
	rm -f post174-daily.blg
	rm -f post174-daily.bbl
	rm -f post174-daily.dvi
	rm -f post174-daily.log
	rm -f post174-daily.aux
	rm -f post174-daily.ps
	rm -f post174-daily.pdf
	rm -f python.png
	rm -f python.eps
	rm -f laplace.png
	rm -f laplace.eps
