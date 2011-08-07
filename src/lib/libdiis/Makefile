
srcdir = /Users/masaaki/psi4/src/lib/libdiis
VPATH = /Users/masaaki/psi4/src/lib/libdiis

include ../MakeVars

PSITARGET = $(top_objdir)/lib/libPSI_diis.a

TRUECXXSRC = \
diismanager.cc \
diisentry.cc

DEPENDINCLUDE = diisentry.h diismanager.h

LIBOBJ = $(TRUECXXSRC:%.cc=%.o)

include ../MakeRules

install_inc:: $(DEPENDINCLUDE)
	$(MKDIRS) $(incdir)/libdiis
	$(INSTALL_INCLUDE) $^ $(incdir)/libdiis

ifneq ($(DODEPEND),no)
include $(LIBOBJ:%.o=%.d)
endif

