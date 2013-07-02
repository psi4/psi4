   WARN = -Waddress \
          -Warray-bounds \
          -Wc++0x-compat \
          -Wchar-subscripts \
          -Wformat \
          -Wmain \
          -Wmissing-braces \
          -Wreorder \
          -Wreturn-type \
          -Wsequence-point \
          -Wsign-compare \
          -Wstrict-aliasing \
          -Wstrict-overflow=1 \
          -Wswitch \
          -Wtrigraphs \
          -Wuninitialized \
          -Wunknown-pragmas \
          -Wunused-function \
          -Wunused-label \
          -Wunused-value \
          -Wunused-variable \
          -Wvolatile-register-var


mra:	mra.yy mra.ll mra-driver.hh mra-driver.cc mra.cc 
	bison -d -ra -o mra-parser.cc mra.yy
	flex -o mra-scanner.cc mra.ll
	g++ $(WARN) -o $@ -I. mra.cc mra-driver.cc mra-scanner.cc mra-parser.cc

run:	mra
	./mra input
	pdflatex eqn.tex
	evince eqn.pdf

clean:
	/bin/rm -f mra *.o position.hh stack.hh location.hh mra-parser.* mra-scanner.* *~ *aux *dvi *pdf *log stuff prog.tex prog.cc prog.tree

# keep
# Makefile   mra.cc   mra-driver.cc   mra-driver.hh   mra.ll   mra.yy
