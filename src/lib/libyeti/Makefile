
srcdir = /Users/masaaki/psi4/src/lib/libyeti
VPATH = /Users/masaaki/psi4/src/lib/libyeti

FC_SYMBOL = 2

include ../MakeVars

PSITARGET = $(top_objdir)/lib/libPSI_yeti.a

CXXDEF += -DFC_SYMBOL=$(FC_SYMBOL) -DHAVE_PSI=1

CXXINC += $(PYTHON_INCLUDE)

TRUECXXSRC = $(notdir $(wildcard $(srcdir)/*.cc))
DEPENDINCLUDE = $(notdir $(wildcard $(srcdir)/*.h))

LIBOBJ = $(TRUECXXSRC:%.cc=%.o)

include ../MakeRules

install_inc:: $(DEPENDINCLUDE)
	$(MKDIRS) $(DESTDIR)$(incdir)/libyeti
	$(INSTALL_INCLUDE) $^ $(DESTDIR)$(incdir)/libyeti

ifneq ($(DODEPEND),no)
include $(LIBOBJ:%.o=%.d)
endif
