
srcdir = /Users/masaaki/psi4/doc/progman

MAN = progman

include ../MakeVars

TEX = progman.tex

SUBTEX = \
svn/svn.tex \
introduction/introduction.tex \
introduction/history.tex \
introduction/license.tex \
python/python.tex \
testing/testing.tex \
sample-codes/sample-codes.tex \
sample-codes/integrals/integrals.tex \
sample-codes/df-mp2/df-mp2.tex

DVI = $(TEX:%.tex=%.dvi)
PS = $(DVI:%.dvi=%.ps)
HTML = $(TEX:%.tex=%.html)

include ../MakeRules
