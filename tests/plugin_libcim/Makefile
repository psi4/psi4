# Location of your PSI4 source
psi_top_srcdir = /home/deprince/psi4/
# Location of your PSI4 install, by default as listed
psi_top_objdir = /home/deprince/psi4/objdir

# Start by figuring out whether we're on Linux or Mac (sorry, Mr. Gates)
UNAME := $(shell uname)

include $(psi_top_objdir)/src/bin/MakeVars

PSITARGET = $(shell basename `pwd`).so
PSIINCLUDE = -I$(psi_top_objdir)/boost/boost_1_48_0 -I$(psi_top_srcdir)/include -I$(psi_top_objdir)/include -I$(psi_top_srcdir)/src/lib -I$(psi_top_objdir)/src/lib -DINSTALLEDPSIDATADIR=\"$(pkgdatadir)\" $(PYTHON_INCLUDE)
PSILIBS = -L$(psi_top_objdir)/lib -lPSI_plugin /home/deprince/psi4/tests/plugin_localcc/plugin_localcc.so

CXXSRC = $(notdir $(wildcard *.cc))
DEPENDINCLUDE = $(notdir $(wildcard *.h))

BINOBJ = $(CXXSRC:%.cc=%.o)

default:: $(PSITARGET)

# Add the flags needed for shared library creation
ifeq ($(UNAME), Linux)
    LDFLAGS = -shared -g
endif
ifeq ($(UNAME), Darwin)
    LDFLAGS = -shared -undefined dynamic_lookup
    CXXOTH += -fno-common
endif

# The object files
%.o: %.cc
	$(CXX) $(CXXDBG) $(CXXOPT) $(CXXDEF) $(CXXOTH) $(PSIINCLUDE) $(PYTHON_INCLUDE) -c $< $(OUTPUT_OPTION)

#include $(psi_top_objdir)/src/bin/MakeRules
$(PSITARGET): $(BINOBJ)
	$(CXX) $(LDFLAGS) -o $@ $^ $(CXXDBG) $(PSILIBS) && touch $(shell basename `pwd`).test

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

