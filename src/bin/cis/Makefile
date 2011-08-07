srcdir = /Users/masaaki/psi4/src/bin/cis
VPATH = /Users/masaaki/psi4/src/bin/cis

include ../MakeVars

LDLIBS += $(LAPACK) $(BLAS)

PSITARGET = $(top_objdir)/lib/libPSI_cis.a
#PSILIBS = -lPSI_dpd -lPSI_qt -lPSI_chkpt -lPSI_iwl -lPSI_psio -lPSI_ciomr -lPSI_ipv1

CXXSRC = \
Fab.cc  U.cc          build_A.cc  d_corr.cc  get_moinfo.cc  mp2.cc \
Fij.cc  Z.cc          cache.cc    denom.cc   get_params.cc  v.cc \
Fkc.cc  amp_write.cc  cis.cc      diag.cc    local.cc \

BINOBJ = $(CXXSRC:%.cc=%.o)

include ../MakeRules

ifneq ($(DODEPEND),no)
$(BINOBJ:%.o=%.d): $(DEPENDINCLUDE)
include $(BINOBJ:%.o=%.d)
endif
