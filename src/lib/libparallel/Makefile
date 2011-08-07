
srcdir = /Users/masaaki/psi4/src/lib/libparallel
VPATH = /Users/masaaki/psi4/src/lib/libparallel

include ../MakeVars

PSITARGET = $(top_objdir)/lib/libPSI_parallel.a
CXXDEF += -DINSTALLEDPSIDATADIR=\"$(pkgdatadir)\"

TRUECXXSRC = communicator.cc mpi.cc mad.cc local.cc fprintf.cc process.cc

INC = parallel.h serialize.h

LIBOBJ = $(TRUECXXSRC:%.cc=%.o)

include ../MakeRules

ifneq ($(DODEPEND),no)
include $(LIBOBJ:%.o=%.d)
endif

install_inc:: $(INC)
	$(MKDIRS) $(DESTDIR)$(incdir)/libutil
	$(INSTALL_INCLUDE) $^ $(DESTDIR)$(incdir)/libutil

