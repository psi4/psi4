srcdir = /Users/masaaki/psi4/src/bin/cusp
VPATH = /Users/masaaki/psi4/src/bin/cusp

include ../MakeVars

LDLIBS += $(BLAS)

PSILIBS = -lPSI_dpd -lPSI_qt -lPSI_chkpt -lPSI_iwl -lPSI_psio -lPSI_ciomr -lPSI_ipv1

CXXSRC = \
cusp.cc phi.cc delta.cc cache.cc local.cc

BINOBJ = $(CXXSRC:%.cc=%.o)

include ../MakeRules

ifneq ($(DODEPEND),no)
$(BINOBJ:%.o=%.d): $(DEPENDINCLUDE)
include $(BINOBJ:%.o=%.d)
endif

