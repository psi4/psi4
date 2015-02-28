# Settle docs-specific dependencies
find_package(Perl QUIET)
if(NOT PERL_FOUND)
   message(STATUS "No Perl, no docs. Pre-built documentation at http://sirius.chem.vt.edu/psi4manual/latest/index.html")
endif()

find_package(Sphinx QUIET)
if(NOT SPHINX_FOUND)
   message(STATUS "No Sphinx, no docs. Pre-built documentation at http://sirius.chem.vt.edu/psi4manual/latest/index.html")
endif()

find_package(LATEX QUIET)
if(NOT (LATEX_COMPILER AND PDFLATEX_COMPILER))
   message(STATUS "No LaTeX (incl. pdflatex), no PDF docs. Pre-built documentation at http://sirius.chem.vt.edu/psi4manual/latest/index.html")
endif()

find_package(Doxygen QUIET)
if(NOT (LATEX_COMPILER AND PDFLATEX_COMPILER))
   message(STATUS "No Doxygen, no docs.")
endif()

