include $(top_srcdir)/config/MakeGlobal.am

bin_PROGRAMS = testpg moldft mcpfit

moldft_SOURCES = moldft.cc lda.cc  molecule.cc \
                 molecule.h molecularbasis.h molecularbasis.cc \
                 atomutil.h atomutil.cc \
                 corepotential.h corepotential.cc

testpg_SOURCES = testpg.cc pointgroup.h

mcpfit_SOURCES = mcpfit.cc atomutil.h atomutil.cc corepotential.h corepotential.cc

moldft_LDADD = $(MRALIBS)

testpg_LDADD = $(MRALIBS)

mcpfit_LDADD = $(MRALIBS)

