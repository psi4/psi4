srcdir = /Users/masaaki/psi4/src/bin/df-mp2
VPATH = /Users/masaaki/psi4/src/bin/df-mp2

include ../MakeVars

LDLIBS += $(LAPACK) $(BLAS)

PSITARGET = $(top_objdir)/lib/libPSI_dfmp2.a

CXXSRC = dfmp2.cc dfmp2_energy.cc wrapper.cc

BINOBJ = $(CXXSRC:%.cc=%.o)

include ../MakeRules

ifneq ($(DODEPEND),no)
$(BINOBJ:%.o=%.d): $(DEPENDINCLUDE)
include $(BINOBJ:%.o=%.d)
endif
