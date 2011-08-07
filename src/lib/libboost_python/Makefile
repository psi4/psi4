
srcdir = /Users/masaaki/psi4/src/lib/libboost_python
VPATH = /Users/masaaki/psi4/src/lib/libboost_python

include ../MakeVars

PSITARGET = $(top_objdir)/lib/libboost_python.a

CXXINC += $(PYTHON_INCLUDE)

TRUECXXSRC = $(notdir $(wildcard $(srcdir)/*.cc)) 

LIBOBJ = $(TRUECXXSRC:%.cc=%.o)

include ../MakeRules

install_inc:: $(DEPENDINCLUDE)
	$(MKDIRS) $(DESTDIR)$(incdir)/libboost_python
	$(INSTALL_INCLUDE) $^ $(DESTDIR)$(incdir)/libboost_python

ifneq ($(DODEPEND),no)
include $(LIBOBJ:%.o=%.d)
endif
