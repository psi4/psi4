# Makefile for Sphinx documentation
#

# You can set these variables from the command line.
SPHINXOPTS    = -w warnings.txt
SPHINXBUILD   = sphinx-build
PAPER         =
BUILDDIR      = build
SOURCEDIR     = source

# Internal variables.
PAPEROPT_a4     = -D latex_paper_size=a4
PAPEROPT_letter = -D latex_paper_size=letter
ALLSPHINXOPTS   = -d $(BUILDDIR)/doctrees $(PAPEROPT_$(PAPER)) $(SPHINXOPTS) $(SOURCEDIR)
# the i18n builder cannot share the environment and doctrees with the others
I18NSPHINXOPTS  = $(PAPEROPT_$(PAPER)) $(SPHINXOPTS) $(SOURCEDIR)

.PHONY: help clean html dirhtml singlehtml pickle json htmlhelp qthelp devhelp epub latex latexpdf text man changes linkcheck doctest gettext latexpdf-user latexpdf-prog index.rst distclean psi4docs html-user html-prog

DISCARD = expert_keywords.tex keywords.tex physconst.tex psi4bases.tex psi4bases.html tests_descriptions.tex variables_list.tex warnings.txt $(SOURCEDIR)/index.rst $(SOURCEDIR)/PsiMod.pyc $(SOURCEDIR)/appendices.rst

help:
	@echo "For PSI4, use one of"
	@echo "  html                to make html complete manual"
	@echo "  html-user           to make html user manual"
	@echo "  html-prog           to make html programmer manual"
	@echo "  latexpdf-user       to make LaTeX pdf user manual"
	@echo "  latexpdf-prog       to make LaTeX pdf programmer manual"
	@echo "  clean               to clean intermediate and autodoc files"
	@echo "  distclean           to clean additionally the build directory"
	@echo "  psi4docs            to make all three manual targets (NYI)"

clean:
	-rm -rf $(SOURCEDIR)/autodoc_* $(SOURCEDIR)/autodir_*/*.rst $(DISCARD)

distclean: clean
	-rm -rf $(BUILDDIR)/*
	touch $(BUILDDIR)/ignore

index.rst:
	./document_python_modules.py
	./document_bases.pl
	./document_options_tests_physconst_psivariables.pl
	./separate_userprogman.py

html: index.rst
	cp -p $(SOURCEDIR)/autodoc_index_html.rst $(SOURCEDIR)/index.rst
	cp -p $(SOURCEDIR)/autodoc_appendices_html.rst $(SOURCEDIR)/appendices.rst
	$(SPHINXBUILD) -b html $(ALLSPHINXOPTS) $(BUILDDIR)/html
	@echo
	@echo "Build finished. The HTML pages are in $(BUILDDIR)/html."

html-user: index.rst
	cp -p $(SOURCEDIR)/autodoc_index_htmluser.rst $(SOURCEDIR)/index.rst
	cp -p $(SOURCEDIR)/autodoc_appendices_htmluser.rst $(SOURCEDIR)/appendices.rst
	$(SPHINXBUILD) -b html $(ALLSPHINXOPTS) $(BUILDDIR)/html-user
	@echo
	@echo "Build finished. The HTML pages are in $(BUILDDIR)/html-user."

html-prog: index.rst
	cp -p $(SOURCEDIR)/autodoc_index_htmlprog.rst $(SOURCEDIR)/index.rst
	cp -p $(SOURCEDIR)/autodoc_appendices_htmlprog.rst $(SOURCEDIR)/appendices.rst
	$(SPHINXBUILD) -b html $(ALLSPHINXOPTS) $(BUILDDIR)/html-prog
	@echo
	@echo "Build finished. The HTML pages are in $(BUILDDIR)/html-prog."

dirhtml:
	$(SPHINXBUILD) -b dirhtml $(ALLSPHINXOPTS) $(BUILDDIR)/dirhtml
	@echo
	@echo "Build finished. The HTML pages are in $(BUILDDIR)/dirhtml."

singlehtml:
	$(SPHINXBUILD) -b singlehtml $(ALLSPHINXOPTS) $(BUILDDIR)/singlehtml
	@echo
	@echo "Build finished. The HTML page is in $(BUILDDIR)/singlehtml."

pickle:
	$(SPHINXBUILD) -b pickle $(ALLSPHINXOPTS) $(BUILDDIR)/pickle
	@echo
	@echo "Build finished; now you can process the pickle files."

json:
	$(SPHINXBUILD) -b json $(ALLSPHINXOPTS) $(BUILDDIR)/json
	@echo
	@echo "Build finished; now you can process the JSON files."

htmlhelp:
	$(SPHINXBUILD) -b htmlhelp $(ALLSPHINXOPTS) $(BUILDDIR)/htmlhelp
	@echo
	@echo "Build finished; now you can run HTML Help Workshop with the" \
	      ".hhp project file in $(BUILDDIR)/htmlhelp."

qthelp:
	$(SPHINXBUILD) -b qthelp $(ALLSPHINXOPTS) $(BUILDDIR)/qthelp
	@echo
	@echo "Build finished; now you can run "qcollectiongenerator" with the" \
	      ".qhcp project file in $(BUILDDIR)/qthelp, like this:"
	@echo "# qcollectiongenerator $(BUILDDIR)/qthelp/Psithon.qhcp"
	@echo "To view the help file:"
	@echo "# assistant -collectionFile $(BUILDDIR)/qthelp/Psithon.qhc"

devhelp:
	$(SPHINXBUILD) -b devhelp $(ALLSPHINXOPTS) $(BUILDDIR)/devhelp
	@echo
	@echo "Build finished."
	@echo "To view the help file:"
	@echo "# mkdir -p $$HOME/.local/share/devhelp/Psithon"
	@echo "# ln -s $(BUILDDIR)/devhelp $$HOME/.local/share/devhelp/Psithon"
	@echo "# devhelp"

epub:
	$(SPHINXBUILD) -b epub $(ALLSPHINXOPTS) $(BUILDDIR)/epub
	@echo
	@echo "Build finished. The epub file is in $(BUILDDIR)/epub."

latex:
	$(SPHINXBUILD) -b latex $(ALLSPHINXOPTS) $(BUILDDIR)/latex
	@echo
	@echo "Build finished; the LaTeX files are in $(BUILDDIR)/latex."
	@echo "Run \`make' in that directory to run these through (pdf)latex" \
	      "(use \`make latexpdf' here to do that automatically)."

latexpdf:
	$(SPHINXBUILD) -b latex $(ALLSPHINXOPTS) $(BUILDDIR)/latex
	@echo "Running LaTeX files through pdflatex..."
	$(MAKE) -C $(BUILDDIR)/latex all-pdf
	@echo "pdflatex finished; the PDF files are in $(BUILDDIR)/latex."

latexpdf-user: index.rst
	cp -p $(SOURCEDIR)/autodoc_index_latexuser.rst $(SOURCEDIR)/index.rst
	cp -p $(SOURCEDIR)/autodoc_appendices_latexuser.rst $(SOURCEDIR)/appendices.rst
	$(SPHINXBUILD) -b latex $(ALLSPHINXOPTS) $(BUILDDIR)/latex-user
	@echo "USERMAN Running LaTeX files through pdflatex..."
	$(MAKE) -C $(BUILDDIR)/latex-user all-pdf
	@echo "USERMAN pdflatex finished; the PDF files are in $(BUILDDIR)/latex-user."
	mv $(BUILDDIR)/latex-user/Psi4.pdf $(BUILDDIR)/latex-user/psi4_userman.pdf

latexpdf-prog: index.rst
	cp -p $(SOURCEDIR)/autodoc_index_latexprog.rst $(SOURCEDIR)/index.rst
	cp -p $(SOURCEDIR)/autodoc_appendices_latexprog.rst $(SOURCEDIR)/appendices.rst
	$(SPHINXBUILD) -b latex $(ALLSPHINXOPTS) $(BUILDDIR)/latex-prog
	@echo "PROGMAN Running LaTeX files through pdflatex..."
	$(MAKE) -C $(BUILDDIR)/latex-prog all-pdf
	@echo "PROGMAN pdflatex finished; the PDF files are in $(BUILDDIR)/latex-prog."
	mv $(BUILDDIR)/latex-prog/Psi4.pdf $(BUILDDIR)/latex-prog/psi4_progman.pdf

text:
	$(SPHINXBUILD) -b text $(ALLSPHINXOPTS) $(BUILDDIR)/text
	@echo
	@echo "Build finished. The text files are in $(BUILDDIR)/text."

man:
	$(SPHINXBUILD) -b man $(ALLSPHINXOPTS) $(BUILDDIR)/man
	@echo
	@echo "Build finished. The manual pages are in $(BUILDDIR)/man."

texinfo:
	$(SPHINXBUILD) -b texinfo $(ALLSPHINXOPTS) $(BUILDDIR)/texinfo
	@echo
	@echo "Build finished. The Texinfo files are in $(BUILDDIR)/texinfo."
	@echo "Run \`make' in that directory to run these through makeinfo" \
	      "(use \`make info' here to do that automatically)."

info:
	$(SPHINXBUILD) -b texinfo $(ALLSPHINXOPTS) $(BUILDDIR)/texinfo
	@echo "Running Texinfo files through makeinfo..."
	make -C $(BUILDDIR)/texinfo info
	@echo "makeinfo finished; the Info files are in $(BUILDDIR)/texinfo."

gettext:
	$(SPHINXBUILD) -b gettext $(I18NSPHINXOPTS) $(BUILDDIR)/locale
	@echo
	@echo "Build finished. The message catalogs are in $(BUILDDIR)/locale."

changes:
	$(SPHINXBUILD) -b changes $(ALLSPHINXOPTS) $(BUILDDIR)/changes
	@echo
	@echo "The overview file is in $(BUILDDIR)/changes."

linkcheck:
	$(SPHINXBUILD) -b linkcheck $(ALLSPHINXOPTS) $(BUILDDIR)/linkcheck
	@echo
	@echo "Link check complete; look for any errors in the above output " \
	      "or in $(BUILDDIR)/linkcheck/output.txt."

doctest:
	$(SPHINXBUILD) -b doctest $(ALLSPHINXOPTS) $(BUILDDIR)/doctest
	@echo "Testing of doctests in the sources finished, look at the " \
	      "results in $(BUILDDIR)/doctest/output.txt."

psi4docs: html latexpdf-user latexpdf-prog

