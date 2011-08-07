srcdir = /Users/masaaki/psi4/src/bin/df-cc
VPATH = /Users/masaaki/psi4/src/bin/df-cc

include ../MakeVars

LDLIBS += $(LAPACK) $(BLAS)

PSITARGET = $(top_objdir)/lib/libPSI_dfcc.a

CXXSRC = wrapper.cc

BINOBJ = $(CXXSRC:%.cc=%.o)

include ../MakeRules

ifneq ($(DODEPEND),no)
$(BINOBJ:%.o=%.d): $(DEPENDINCLUDE)
include $(BINOBJ:%.o=%.d)
endif
