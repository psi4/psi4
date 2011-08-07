srcdir = /Users/masaaki/psi4/src/bin/adc
VPATH = /Users/masaaki/psi4/src/bin/adc

include ../MakeVars

LDLIBS += $(LAPACK) $(BLAS)
PSITARGET = $(top_objdir)/lib/libPSI_adc.a
#PSILIBS = -lPSI_dpd -lPSI_qt -lPSI_chkpt -lPSI_iwl -lPSI_psio -lPSI_ciomr -lPSI_ipv1 -lPSI_deriv -lPSI_int

CXXSRC = \
adc.cc cache.cc get_moinfo.cc get_params.cc primary_space.cc\
diag_adc.cc diagonalize.cc make_diag.cc mp2.cc amp_write.cc\
orbital_relaxation.cc denom.cc init_ints.cc accelerate.cc\
differential_correlation.cc ccamps.cc ocss.cc sort_pert.cc\
transpert.cc\

BINOBJ = $(CXXSRC:%.cc=%.o)

include ../MakeRules

ifneq ($(DODEPEND),no)
$(BINOBJ:%.o=%.d): $(DEPENDINCLUDE)
include $(BINOBJ:%.o=%.d)
endif

