# Location of your PSI4 source
psi_top_srcdir = /Users/masaaki/psi4
# Location of your PSI4 install, by default as listed
psi_top_objdir = ../..

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
PSILIBS =  -L$(psi_top_objdir)/lib -lPSI_plugin
MADFLAGS = -I$(top_objdir)/src/lib -I$(top_srcdir)/madness/include -I$(top_srcdir)/madness/src/lib -I$(top_objdir)/madness/include

CXXSRC = $(notdir $(wildcard *.cc))
DEPENDINCLUDE = $(notdir $(wildcard *.h))

BINOBJ = $(CXXSRC:%.cc=%.o)

default:: $(PSITARGET)

# Add the flags needed for shared library creation
ifeq ($(UNAME), Linux)
    LDFLAGS = -fPIC -shared
    CXXOTH += -fPIC
endif
ifeq ($(UNAME), Darwin)
    LDFLAGS = -fPIC -shared -undefined dynamic_lookup
    CXXOTH += -fPIC -fno-common
endif

# The object files
%.o: %.cc
	$(CXX) $(MADFLAGS) $(CXXDBG) $(CXXOPT) $(CXXDEF) $(CXXOTH) $(PSIINCLUDE) -c $< $(OUTPUT_OPTION)

#include $(psi_top_objdir)/src/bin/MakeRules
$(PSITARGET): $(BINOBJ)
	$(CXX) $(MADFLAGS) $(LDFLAGS) -o $@ $^ $(CXXDBG) $(PSILIBS)

# Erase all compiled intermediate files
clean:
	rm -f $(BINOBJ) $(PSITARGET) *.d

# Dependency handling
%.d: %.cc
	$(CXXDEPEND) $(CXXDEPENDFLAGS) $(CXXDEF) $(CXXOTH) $(PSIINCLUDE) $< | sed 's/^$*.o/$*.o $*.d/g' > $(@F)

$(BINOBJ:%.o=%.d): $(DEPENDINCLUDE)

