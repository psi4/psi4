
srcdir = /Users/masaaki/psi4/src/util/tocprint
VPATH = /Users/masaaki/psi4/src/util/tocprint

include ../MakeVars

CXXSRC = tocprint.cc cc2unit.cc
BINOBJ = $(CXXSRC:%.cc=%.o)
PSILIBS = -lPSI_chkpt -lPSI_psio -lPSI_ciomr -lPSI_ipv1
ALLOC =

include ../MakeRules

ifneq ($(DODEPEND),no)
$(BINOBJ:%.o=%.d): $(DEPENDINCLUDE)
include $(BINOBJ:%.o=%.d)
endif

install_man:: tocprint.1
	$(MKDIRS) $(DESTDIR)$(mandir)/man1
	$(INSTALL_PROGRAM) $(srcdir)/tocprint.1 $(DESTDIR)$(mandir)/man1
