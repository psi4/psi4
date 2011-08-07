srcdir = /Users/masaaki/psi4/src/bin/cchbar
VPATH = /Users/masaaki/psi4/src/bin/cchbar

include ../MakeVars

LDLIBS += $(BLAS)

PSITARGET = $(top_objdir)/lib/libPSI_cchbar.a

CXXSRC = \
F.cc               Wabei_RHF.cc        cache.cc       get_params.cc \
Fai.cc             Wabei_RHF_FT2_a.cc  cc2_Wabei.cc   norm_HET1.cc \
HET1_Wabef.cc      Wabei_ROHF.cc       cc2_Wmbej.cc   purge.cc \
Wabei.cc           Wabij.cc            cc2_Wmbij.cc   reference.cc \
Wabei_AAAA_UHF.cc  Wamef.cc            cc2_Zmbej.cc   sort_amps.cc \
Wabei_ABAB_UHF.cc  Wmbej.cc            cc3_HET1.cc    status.cc \
Wabei_BABA_UHF.cc  Wmbij.cc            cchbar.cc      tau.cc \
Wabei_BBBB_UHF.cc  Wmnie.cc            get_moinfo.cc  taut.cc

BINOBJ = $(CXXSRC:%.cc=%.o)

include ../MakeRules

ifneq ($(DODEPEND),no)
$(BINOBJ:%.o=%.d): $(DEPENDINCLUDE)
include $(BINOBJ:%.o=%.d)
endif

