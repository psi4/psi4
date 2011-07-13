# Location of your PSI4 source
psi_top_srcdir = /Users/andysim/programming/psi4
# Location of your PSI4 install, by default as listed
psi_top_objdir = /Users/andysim/programming/psi4/objdir_mpi_debug

#
#
# None of this stuff should need to be modified
#
#

# Start by figuring out whether we're on Linux or Mac (sorry, Mr. Gates)
UNAME := $(shell uname)


include $(psi_top_objdir)/src/bin/MakeVars
PSITARGET = $(shell basename `pwd`).so

PSIINCLUDE = -I$(psi_top_srcdir)/madness/include  -I$(psi_top_objdir)/madness/include -I$(psi_top_srcdir)/madness/src/lib -I$(psi_top_srcdir)/include -I$(psi_top_objdir)/include -I$(psi_top_srcdir)/src/lib -I$(psi_top_objdir)/src/lib

LDLIBS += $(LAPACK) $(BLAS) $(PYTHON_LDFLAGS) -lpthread
CXXINC += -DINSTALLEDPSIDATADIR=\"$(pkgdatadir)\" $(PYTHON_INCLUDE)

# NOTE: Keep PSI_parallel last in the PSI library list (but before boost). It contains the fprintf override.
PSILIBS = -L$(psi_top_objdir)/lib -lPSI_plugin -lPSI_yeti -lPSI_smartptr 
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
	$(CXX) $(CXXDBG) $(CXXOPT) $(CXXDEF) $(CXXOTH) $(PSIINCLUDE) -c $< $(OUTPUT_OPTION)

#include $(psi_top_objdir)/src/bin/MakeRules
$(PSITARGET): $(BINOBJ)
	$(CXX) $(LDFLAGS) -o $@ $^ $(CXXDBG) $(LDLIBS) $(PSILIBS)

# Erase all compiled intermediate files
clean:
	rm -f $(BINOBJ) $(PSITARGET) *.d

# Dependency handling
%.d: %.cc
	$(CXXDEPEND) $(CXXDEPENDFLAGS) $(CXXDEF) $(CXXOTH) $(PSIINCLUDE) $< | sed 's/^$*.o/$*.o $*.d/g' > $(@F)

$(BINOBJ:%.o=%.d): $(DEPENDINCLUDE)

include $(BINOBJ:%.o=%.d)
