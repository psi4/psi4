
srcdir = /Users/masaaki/psi4/src/bin/cclambda
VPATH = /Users/masaaki/psi4/src/bin/cclambda

include ../MakeVars

LDLIBS += $(LAPACK) $(BLAS)

PSITARGET=$(top_objdir)/lib/libPSI_cclambda.a
#PSILIBS = -lPSI_dpd -lPSI_qt -lPSI_chkpt -lPSI_iwl -lPSI_psio -lPSI_ciomr -lPSI_ipv1

CXXSRC = \
BL2_AO.cc  Lamp_write.cc  cache.cc           cclambda.cc     init_amps.cc \
DL2.cc     Lmag.cc        cc2_Gai.cc         check_ortho.cc  local.cc \
FaeL2.cc   Lnorm.cc       cc2_L1.cc          check_sum.cc    ortho_Rs.cc \
FmiL2.cc   Lsave.cc       cc2_L2.cc          converged.cc    overlap.cc \
G.cc       WabeiL1.cc     cc2_faeL2.cc       denom.cc        overlap_LAMPS.cc \
GL2.cc     WefabL2.cc     cc2_fmiL2.cc       diis.cc         projections.cc \
L1.cc      WejabL2.cc     cc2_hbar_extra.cc  dijabL2.cc      pseudoenergy.cc \
L1FL2.cc   WijmbL2.cc     cc3_l3l1.cc        get_moinfo.cc   sort_amps.cc \
L2.cc      WijmnL2.cc     cc3_l3l2.cc        get_params.cc   spinad_amps.cc \
L3_AAA.cc  WmbejL2.cc     cc3_t3x.cc         halftrans.cc    status.cc \
L3_AAB.cc  c_clean.cc     cc3_t3z.cc         hbar_extra.cc   update.cc

BINOBJ = $(CXXSRC:%.cc=%.o)

include ../MakeRules

ifneq ($(DODEPEND),no)
$(BINOBJ:%.o=%.d): $(DEPENDINCLUDE)
include $(BINOBJ:%.o=%.d)
endif

