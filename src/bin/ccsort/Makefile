srcdir = /Users/masaaki/psi4/src/bin/ccsort
VPATH = /Users/masaaki/psi4/src/bin/ccsort

include ../MakeVars

LDLIBS += $(LAPACK) $(BLAS)

PSITARGET = $(top_objdir)/lib/libPSI_ccsort.a
#PSILIBS = -lPSI_dpd -lPSI_qt -lPSI_chkpt -lPSI_iwl -lPSI_psio -lPSI_ciomr -lPSI_ipv1

CXXSRC = \
b_sort.cc                distribute.cc            idx_error.cc \
build_B_RHF.cc           domain_print.cc          idx_permute.cc \
build_F_RHF.cc           e_sort.cc                idx_permute_multipass.cc \
build_abcd_packed.cc     e_spinad.cc              idx_permute_presort.cc \
c_sort.cc                f_sort.cc                local.cc \
cache.cc                 f_spinad.cc              local_magnetic.cc \
cc_memcheck.cc           file_build.cc            local_polar.cc \
ccsort.cc                file_build_multipass.cc  scf_check.cc \
classify.cc              file_build_presort.cc    sort_oei.cc \
cphf_B.cc                fock.cc                  sort_pert.cc \
cphf_F.cc                sort_tei.cc \
d_sort.cc                get_moinfo.cc            \
d_spinad.cc              get_params.cc            transpert.cc \
denom.cc

BINOBJ = $(CXXSRC:%.cc=%.o)

include ../MakeRules

ifneq ($(DODEPEND),no)
$(BINOBJ:%.o=%.d): $(DEPENDINCLUDE)
include $(BINOBJ:%.o=%.d)
endif

