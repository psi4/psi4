# Location of your PSI4 source
psi_top_srcdir = /home/parrish/psi4/
# Location of your PSI4 install, by default as listed
psi_top_objdir = /home/parrish/psi4/obj-plugin/

#
#
# None of this stuff should need to be modified
#
#

# Start by figuring out whether we're on Linux or Mac (sorry, Mr. Gates)
UNAME := $(shell uname)

include $(psi_top_objdir)/src/bin/MakeVars

PSITARGET = $(shell basename `pwd`).so
PSIINCLUDE = -I$(psi_top_srcdir)/include -I$(psi_top_objdir)/include -I$(psi_top_srcdir)/src/lib -I$(psi_top_objdir)/src/lib
CXXINC += -DINSTALLEDPSIDATADIR=\"$(pkgdatadir)\" $(PYTHON_INCLUDE)
PSILIBS = -L$(psi_top_objdir)/lib -lPSI_plugin

CXXSRC = $(notdir $(wildcard *.cc))
DEPENDINCLUDE = $(notdir $(wildcard *.h))

BINOBJ = $(CXXSRC:%.cc=%.o)

default:: $(PSITARGET)

# Add the flags needed for shared library creation
ifeq ($(UNAME), Linux)
    LDFLAGS = -fPIC -shared
    CXXOTH += -fPIC -g
endif
ifeq ($(UNAME), Darwin)
    LDFLAGS = -fPIC -shared -undefined dynamic_lookup
    CXXOTH += -fPIC -fno-common
endif

# The object files
%.o: %.cc
	$(CXX) $(CXXDBG) $(CXXOPT) $(CXXDEF) $(CXXOTH) $(PSIINCLUDE) -c $< $(OUTPUT_OPTION)

#include $(psi_top_objdir)/src/bin/MakeRules
$(PSITARGET): $(BINOBJ)
	$(CXX) $(LDFLAGS) -o $@ $^ $(CXXDBG) $(PSILIBS)

# Erase all compiled intermediate files
clean:
	rm -f $(BINOBJ) $(PSITARGET) *.d

# Dependency handling
%.d: %.cc
	$(CXXDEPEND) $(CXXDEPENDFLAGS) $(CXXDEF) $(CXXOTH) $(PSIINCLUDE) $< | sed 's/^$*.o/$*.o $*.d/g' > $(@F)

ifneq ($(DODEPEND),no)
$(BINOBJ:%.o=%.d): $(DEPENDINCLUDE)
include $(BINOBJ:%.o=%.d)
endif
