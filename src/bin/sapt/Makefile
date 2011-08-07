srcdir = /Users/masaaki/psi4/src/bin/sapt
VPATH = /Users/masaaki/psi4/src/bin/sapt

include ../MakeVars

LDLIBS += $(LAPACK) $(BLAS)

PSITARGET = $(top_objdir)/lib/libPSI_sapt.a

CXXSRC = wrapper.cc

BINOBJ = $(CXXSRC:%.cc=%.o)

include ../MakeRules

ifneq ($(DODEPEND),no)
$(BINOBJ:%.o=%.d): $(DEPENDINCLUDE)
include $(BINOBJ:%.o=%.d)
endif
