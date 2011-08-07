
srcdir = /Users/masaaki/psi4/src/lib/liboptions
VPATH = /Users/masaaki/psi4/src/lib/liboptions

include ../MakeVars

PSITARGET = $(top_objdir)/lib/libPSI_options.a

TRUECSRC = 

TRUECXXSRC = print.cc

INC = liboptions.h

LIBOBJ = $(TRUECSRC:%.c=%.o) $(TRUECXXSRC:%.cc=%.o)

include ../MakeRules

ifneq ($(DODEPEND),no)
include $(LIBOBJ:%.o=%.d)
endif

install_inc:: $(INC)
	$(MKDIRS) $(DESTDIR)$(incdir)/liboptions
	$(INSTALL_INCLUDE) $^ $(DESTDIR)$(incdir)/liboptions

