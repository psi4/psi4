# Location of your PSI4 source
psi_top_srcdir = /home/deprince/psi4
# Location of your PSI4 install, by default as listed
psi_top_objdir = /home/deprince/psi4/objdir

# Start by figuring out whether we're on Linux or Mac (sorry, Mr. Gates)
UNAME := $(shell uname)

include $(psi_top_objdir)/src/bin/MakeVars

PSITARGET = $(shell basename `pwd`).so
#PSIINCLUDE = -I$(psi_top_srcdir)/include -I$(psi_top_objdir)/include -I$(psi_top_srcdir)/src/lib -I$(psi_top_objdir)/src/lib
PSIINCLUDE = -I$(psi_top_objdir)/boost/boost_1_48_0 -I$(psi_top_srcdir)/include -I$(psi_top_objdir)/include -I$(psi_top_srcdir)/src/lib -I$(psi_top_objdir)/src/lib -DINSTALLEDPSIDATADIR=\"$(pkgdatadir)\" $(PYTHON_INCLUDE)
CXXINC += -DINSTALLEDPSIDATADIR=\"$(pkgdatadir)\" $(PYTHON_INCLUDE)
PSILIBS = -L$(psi_top_objdir)/lib -lPSI_plugin

CXXSRC = $(notdir $(wildcard *.cc))
CUDASRC = $(notdir $(wildcard *.cu))
DEPENDINCLUDE = $(notdir $(wildcard *.h))

NVCCFLAGS = -arch sm_13 -Xcompiler -fopenmp --compiler-options '-fPIC'

BINOBJ = $(CXXSRC:%.cc=%.o)
CUDAOBJ = $(CUDASRC:%.cu=%.o)

default:: $(PSITARGET)

#GPULIB=/global/u1/d/deprince/psi4/tests/plugin_libfock/gpujk.so

# Add the flags needed for shared library creation
ifeq ($(UNAME), Linux)
    LDFLAGS = -shared -lcublas -L/usr/common/usg/cuda/3.2/lib64 #$(GPULIB) -DGPU
endif
ifeq ($(UNAME), Darwin)
    LDFLAGS = -shared -undefined dynamic_lookup
    CXXOTH += -fno-common
endif

NVCC=nvcc

# The object files
%.o: %.cc
	$(CXX) $(CXXDBG) $(CXXOPT) $(CXXDEF) $(CXXOTH) $(PSIINCLUDE) -c $< $(OUTPUT_OPTION)

%.o: %.cu
	$(NVCC) $(NVCCFLAGS) $(CXXDBG) $(CXXOPT) $(CXXDEF) $(PSIINCLUDE) -c $< $(OUTPUT_OPTION)

#include $(psi_top_objdir)/src/bin/MakeRules
$(PSITARGET): $(BINOBJ) $(CUDAOBJ)
	$(NVCC) $(LDFLAGS) -o $@ $^ $(CXXDBG) $(PSILIBS)
	#$(CXX) $(LDFLAGS) -o $@ $^ $(CXXDBG) $(PSILIBS)

# Erase all compiled intermediate files
clean:
	rm -f $(BINOBJ) $(CUDAOBJ) $(PSITARGET) *.d

# Dependency handling
%.d: %.cc
	$(CXXDEPEND) $(CXXDEPENDFLAGS) $(CXXDEF) $(CXXOTH) $(PSIINCLUDE) $< | sed 's/^$*.o/$*.o $*.d/g' > $(@F)
%.d: %.cu
	$(NVCC) $(CXXDEPENDFLAGS) $(CXXDEF) $(PSIINCLUDE) $< | sed 's/^$*.o/$*.o $*.d/g' > $(@F)

ifneq ($(DODEPEND),no)
$(BINOBJ:%.o=%.d): $(DEPENDINCLUDE)
include $(BINOBJ:%.o=%.d)
$(CUDAOBJ:%.o=%.d): $(DEPENDINCLUDE)
include $(CUDAOBJ:%.o=%.d)
endif
