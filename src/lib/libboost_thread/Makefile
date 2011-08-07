
srcdir = /Users/masaaki/psi4/src/lib/libboost_thread
VPATH = /Users/masaaki/psi4/src/lib/libboost_thread

include ../MakeVars

PSITARGET = $(top_objdir)/lib/libboost_thread.a

TRUECXXSRC = $(notdir $(wildcard $(srcdir)/*.cc))

LIBOBJ = $(TRUECXXSRC:%.cc=%.o)

include ../MakeRules

install_inc:: $(DEPENDINCLUDE)
	$(MKDIRS) $(DESTDIR)$(incdir)/libboost_thread
	$(INSTALL_INCLUDE) $^ $(DESTDIR)$(incdir)/libboost_thread

ifneq ($(DODEPEND),no)
include $(LIBOBJ:%.o=%.d)
endif
