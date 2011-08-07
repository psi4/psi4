
srcdir = /Users/masaaki/psi4/src/lib/libboost_regex
VPATH = /Users/masaaki/psi4/src/lib/libboost_regex

include ../MakeVars

PSITARGET = $(top_objdir)/lib/libboost_regex.a

TRUECXXSRC = $(notdir $(wildcard $(srcdir)/*.cc))

LIBOBJ = $(TRUECXXSRC:%.cc=%.o)

include ../MakeRules

install_inc:: $(DEPENDINCLUDE)
	$(MKDIRS) $(DESTDIR)$(incdir)/libboost_regex
	$(INSTALL_INCLUDE) $^ $(DESTDIR)$(incdir)/libboost_regex

ifneq ($(DODEPEND),no)
include $(LIBOBJ:%.o=%.d)
endif
