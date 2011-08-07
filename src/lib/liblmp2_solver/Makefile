
srcdir = /Users/masaaki/psi4/src/lib/liblmp2_solver
VPATH = /Users/masaaki/psi4/src/lib/liblmp2_solver

include ../MakeVars

PSITARGET = $(top_objdir)/lib/libPSI_lmp2_solver.a

TRUECXXSRC = $(notdir $(wildcard $(srcdir)/*.cc))
DEPENDINCLUDE = $(notdir $(wildcard $(srcdir)/*.h))

LIBOBJ = $(TRUECXXSRC:%.cc=%.o)

include ../MakeRules

install_inc:: $(DEPENDINCLUDE)
	$(MKDIRS) $(incdir)/liblmp2_solver
	$(INSTALL_INCLUDE) $^ $(incdir)/liblmp2_solver


ifneq ($(DODEPEND),no)
include $(LIBOBJ:%.o=%.d)
endif

