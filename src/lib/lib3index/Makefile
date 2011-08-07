
srcdir = /Users/masaaki/psi4/src/lib/lib3index
VPATH = /Users/masaaki/psi4/src/lib/lib3index

include ../MakeVars

PSITARGET = $(top_objdir)/lib/libPSI_3index.a

TRUECXXSRC = $(notdir $(wildcard $(srcdir)/*.cc))
DEPENDINCLUDE = $(notdir $(wildcard $(srcdir)/*.h))

LIBOBJ = $(TRUECXXSRC:%.cc=%.o)

include ../MakeRules

install_inc:: $(DEPENDINCLUDE)
	$(MKDIRS) $(incdir)/lib3index
	$(INSTALL_INCLUDE) $^ $(incdir)/lib3index


ifneq ($(DODEPEND),no)
include $(LIBOBJ:%.o=%.d)
endif

