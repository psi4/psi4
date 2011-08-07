srcdir = /Users/masaaki/psi4/src/bin/mints
VPATH = /Users/masaaki/psi4/src/bin/mints

include ../MakeVars

PSITARGET=$(top_objdir)/lib/libPSI_mints_wrapper.a

CXXSRC = $(notdir $(wildcard $(srcdir)/*.cc))
DEPENDINCLUDE = $(notdir $(wildcard $(srcdir)/*.h))

BINOBJ = $(CXXSRC:%.cc=%.o)

include ../MakeRules

ifneq ($(DODEPEND),no)
$(BINOBJ:%.o=%.d): $(DEPENDINCLUDE)
include $(BINOBJ:%.o=%.d)
endif

