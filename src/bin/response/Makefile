srcdir = /Users/masaaki/psi4/src/bin/response
VPATH = /Users/masaaki/psi4/src/bin/response

include ../MakeVars

LDLIBS += $(LAPACK) $(BLAS)

PSILIBS = -lPSI_dpd -lPSI_qt -lPSI_chkpt -lPSI_iwl -lPSI_psio -lPSI_ciomr -lPSI_ipv1

CXSRC = \
build_A_RHF.cc  cache.cc       get_params.cc 		response.cc \
build_B_RHF.cc  get_moinfo.cc  invert_RPA_RHF.cc \
polar.cc        optrot.cc      dot.cc                   transpert.cc \
dipquad.cc

BINOBJ = $(CXSRC:%.cc=%.o)

include ../MakeRules

ifneq ($(DODEPEND),no)
$(BINOBJ:%.o=%.d): $(DEPENDINCLUDE)
include $(BINOBJ:%.o=%.d)
endif

