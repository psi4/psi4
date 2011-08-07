
srcdir = /Users/masaaki/psi4/src/lib/libplugin
VPATH = /Users/masaaki/psi4/src/lib/libplugin

include ../MakeVars

PSITARGET = $(top_objdir)/lib/libPSI_plugin.a
CXXDEF += -DINSTALLEDPSIDATADIR=\"$(pkgdatadir)\"

TRUECXXSRC = plugin.cc load_plugin.cc close_plugin.cc

INC = 

LIBOBJ = $(TRUECXXSRC:%.cc=%.o)

include ../MakeRules

ifneq ($(DODEPEND),no)
include $(LIBOBJ:%.o=%.d)
endif

install_inc:: $(INC)
	$(MKDIRS) $(DESTDIR)$(incdir)/libutil
	$(INSTALL_INCLUDE) $^ $(DESTDIR)$(incdir)/libutil

