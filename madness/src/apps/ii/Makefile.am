include $(top_srcdir)/config/MakeGlobal.am

#SUBDIRS = xc


bin_PROGRAMS = testsystolic eigen_solver localize_boys tensor_training

#moldft_SOURCES = moldft.cc lda.cc  molecule.cc \
#                 \
#                molecule.h molecularbasis.h

#moldft_LDADD = xc/libMADxc.a $(MRALIBS)

testsystolic_SOURCES =  testsystolic.cc
testsystolic_LDADD = $(MRALIBS)

eigen_solver_SOURCES =  eigen_solver.cc
eigen_solver_LDADD = $(MRALIBS)

localize_boys_SOURCES = localize_boys.cc
localize_boys_LDADD = $(MRALIBS)

tensor_training_SOURCES = tensor_training.cc
tensor_training_LDADD = $(MRALIBS)

#moldft_SOURCES = moldft.cc lda.cc  molecule.cc
