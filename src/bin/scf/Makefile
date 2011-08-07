srcdir = /Users/masaaki/psi4/src/bin/scf
VPATH = /Users/masaaki/psi4/src/bin/scf

include ../MakeVars

PSITARGET = $(top_objdir)/lib/libPSI_scf.a

CXXINC += $(PYTHON_INCLUDE)
CXXSRC = main.cc

BINOBJ = $(CXXSRC:%.cc=%.o)

include ../MakeRules

ifneq ($(DODEPEND),no)
$(BINOBJ:%.o=%.d): $(DEPENDINCLUDE)
include $(BINOBJ:%.o=%.d)
endif

