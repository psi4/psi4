
srcdir = /Users/masaaki/psi4/src/lib/libboost_filesystem
VPATH = /Users/masaaki/psi4/src/lib/libboost_filesystem

include ../MakeVars

PSITARGET = $(top_objdir)/lib/libboost_filesystem.a

TRUECXXSRC = $(notdir $(wildcard $(srcdir)/*.cc))

LIBOBJ = $(TRUECXXSRC:%.cc=%.o)

include ../MakeRules

install_inc:: $(DEPENDINCLUDE)
	$(MKDIRS) $(DESTDIR)$(incdir)/libboost_filesystem
	$(INSTALL_INCLUDE) $^ $(DESTDIR)$(incdir)/libboost_filesystem

ifneq ($(DODEPEND),no)
include $(LIBOBJ:%.o=%.d)
endif
