
srcdir = /Users/masaaki/psi4/src/lib/libsmartptr
VPATH = /Users/masaaki/psi4/src/lib/libsmartptr

FC_SYMBOL = 2

include ../MakeVars

PSITARGET = $(top_objdir)/lib/libPSI_smartptr.a

CXXDEF += -DFC_SYMBOL=$(FC_SYMBOL)

CXXINC += $(PYTHON_INCLUDE)

TRUECXXSRC = $(notdir $(wildcard $(srcdir)/*.cc))
DEPENDINCLUDE = $(notdir $(wildcard $(srcdir)/*.h))

LIBOBJ = $(TRUECXXSRC:%.cc=%.o)

include ../MakeRules

install_inc:: $(DEPENDINCLUDE)
	$(MKDIRS) $(DESTDIR)$(incdir)/libsmartptr
	$(INSTALL_INCLUDE) $^ $(DESTDIR)$(incdir)/libsmartptr

ifneq ($(DODEPEND),no)
include $(LIBOBJ:%.o=%.d)
endif
