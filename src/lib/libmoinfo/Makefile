
srcdir = /Users/masaaki/psi4/src/lib/libmoinfo
VPATH = /Users/masaaki/psi4/src/lib/libmoinfo

include ../MakeVars

PSITARGET = $(top_objdir)/lib/libPSI_moinfo.a

TRUECSRC = 

TRUECXXSRC = moinfo_model_space.cc \
             moinfo_mappings.cc moinfo_slaterdeterminant.cc \
             moinfo_base.cc moinfo_scf.cc moinfo.cc \
             orbital_space.cc \
             model_space.cc model_space_build.cc \
             slater_determinant.cc \

INC = libmoinfo.h moinfo_base.h moinfo_scf.h moinfo.h orbital_space.h model_space.h slater_determinant.h

LIBOBJ = $(TRUECSRC:%.c=%.o) $(TRUECXXSRC:%.cc=%.o)

include ../MakeRules

ifneq ($(DODEPEND),no)
include $(LIBOBJ:%.o=%.d)
endif

install_inc:: $(INC)
	$(MKDIRS) $(DESTDIR)$(incdir)/libmoinfo
	$(INSTALL_INCLUDE) $^ $(DESTDIR)$(incdir)/libmoinfo

