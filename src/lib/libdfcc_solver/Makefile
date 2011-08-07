
srcdir = /Users/masaaki/psi4/src/lib/libdfcc_solver
VPATH = /Users/masaaki/psi4/src/lib/libdfcc_solver

include ../MakeVars

PSITARGET = $(top_objdir)/lib/libPSI_dfcc_solver.a

TRUECXXSRC = $(notdir $(wildcard $(srcdir)/*.cc))
DEPENDINCLUDE = $(notdir $(wildcard $(srcdir)/*.h))

DEPENDINCLUDE = 

LIBOBJ = $(TRUECXXSRC:%.cc=%.o)

include ../MakeRules

install_inc:: $(DEPENDINCLUDE)
	$(MKDIRS) $(incdir)/libdfcc_solver
	$(INSTALL_INCLUDE) $^ $(incdir)/libdfcc_solver

ifneq ($(DODEPEND),no)
include $(LIBOBJ:%.o=%.d)
endif

