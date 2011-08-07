srcdir = /Users/masaaki/psi4/src/bin/geom
VPATH = /Users/masaaki/psi4/src/bin/geom

include ../MakeVars

#PSILIBS = libPSI_ciomr.a libPSI_ipv1.a
PSITARGET = $(top_objdir)/lib/libPSI_geom.o

CXXSRC = geom.cc read11.cc read_geom.cc gprgid.cc

BINOBJ = $(CXXSRC:%.cc=%.o)

include ../MakeRules

ifneq ($(DODEPEND),no)
$(BINOBJ:%.o=%.d): $(DEPENDINCLUDE)
include $(BINOBJ:%.o=%.d)
endif

install_man:: geom.1
	$(MKDIRS) $(DESTDIR)$(mandir)/man1
	$(INSTALL_INCLUDE) $^ $(DESTDIR)$(mandir)/man1
