
srcdir = /Users/masaaki/psi4/src/lib/libmints
VPATH = /Users/masaaki/psi4/src/lib/libmints

FC_SYMBOL = 2

include ../MakeVars

PSITARGET = $(top_objdir)/lib/libPSI_mints.a

CXXDEF += -DFC_SYMBOL=$(FC_SYMBOL)

CXXINC += $(PYTHON_INCLUDE)

TRUECXXSRC = $(notdir $(wildcard $(srcdir)/*.cc))
DEPENDINCLUDE = $(notdir $(wildcard $(srcdir)/*.h))

LIBOBJ = $(TRUECXXSRC:%.cc=%.o)

include ../MakeRules

install_inc:: $(DEPENDINCLUDE)
	$(MKDIRS) $(DESTDIR)$(incdir)/libmints
	$(INSTALL_INCLUDE) $^ $(DESTDIR)$(incdir)/libmints

ifneq ($(DODEPEND),no)
include $(LIBOBJ:%.o=%.d)
endif
