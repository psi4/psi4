
srcdir = /Users/masaaki/psi4/src/lib/libutil
VPATH = /Users/masaaki/psi4/src/lib/libutil

include ../MakeVars

PSITARGET = $(top_objdir)/lib/libPSI_util.a

TRUECSRC = 

TRUECXXSRC = stl_string.cc combinations.cc exception.cc memory_manager.cc quad.cc

INC = libutil.h ref.h class_macros.h memory_manager.h quad.h

LIBOBJ = $(TRUECSRC:%.c=%.o) $(TRUECXXSRC:%.cc=%.o)

include ../MakeRules

ifneq ($(DODEPEND),no)
include $(LIBOBJ:%.o=%.d)
endif

install_inc:: $(INC)
	$(MKDIRS) $(DESTDIR)$(incdir)/libutil
	$(INSTALL_INCLUDE) $^ $(DESTDIR)$(incdir)/libutil

