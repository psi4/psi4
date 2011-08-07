srcdir = /Users/masaaki/psi4/src/bin/transqt
VPATH = /Users/masaaki/psi4/src/bin/transqt

include ../MakeVars

LDLIBS += $(LAPACK) $(BLAS)
CDEF += -DUSE_BLAS

PSITARGET=$(top_objdir)/lib/libPSI_transqt.a
#PSILIBS = -lPSI_iwl -lPSI_qt -lPSI_chkpt -lPSI_psio -lPSI_ciomr -lPSI_ipv1

CXXSRC = transqt.cc \
semicanonical.cc \
transform_one.cc \
transform_two.cc \
transform_two_mp2.cc \
transform_two_mp2r12a_gr.cc \
transform_two_mp2r12a_t.cc \
transform_two_uhf.cc \
transform_two_backtr_uhf.cc \
yoshimine.cc \
frozen_core.cc \
backsort.cc \
cleanup.cc

BINOBJ = $(CXXSRC:%.cc=%.o)

include ../MakeRules

ifneq ($(DODEPEND),no)
$(BINOBJ:%.o=%.d): $(DEPENDINCLUDE)
include $(BINOBJ:%.o=%.d)
endif

install_man:: transqt.1
	$(MKDIRS) $(DESTDIR)$(mandir)/man1
	$(INSTALL_INCLUDE) $^ $(DESTDIR)$(mandir)/man1

