# Location of your PSI4 source
psi_top_srcdir = /Users/deprince/psi4
# Location of your PSI4 install, by default as listed
psi_top_objdir = /Users/deprince/psi4/objdir

#############################################################
#   set CUDA=anything if you want to use cuda. if it is 
#   left commented out, the .cu file will not be compiled
#############################################################
#CUDA=uncommentifyouwantcuda

ifdef CUDA
HAVECUDA=-DCUDA
NVCC=nvcc
LCUBLAS=-L/usr/local/cuda/lib64 -L/usr/local/cuda/lib -lcublas
%.o: %.cu
	$(NVCC) $(NVCCFLAGS) $(CXXDBG) $(CXXOPT) $(CXXDEF) $(PSIINCLUDE) $(HAVECUDA) -c $< $(OUTPUT_OPTION) 
%.d: %.cu
	$(NVCC) $(CXXDEPENDFLAGS) $(CXXDEF) $(PSIINCLUDE) $< | sed 's/^$*.o/$*.o $*.d/g' > $(@F)
NVCCFLAGS = -arch sm_13 -Xcompiler -fopenmp --compiler-options '-fPIC' $(HAVECUDA)
CUDAOBJ = $(CUDASRC:%.cu=%.o)
else
LCUBLAS=
NVCC=$(CXX)
CUDAOBJ=
NVCCFLAGS=
endif
#############################################################
#     special cuda stuff is done
#############################################################


# Start by figuring out whether we're on Linux or Mac (sorry, Mr. Gates)
UNAME := $(shell uname)

include $(psi_top_objdir)/src/bin/MakeVars

PSITARGET = $(shell basename `pwd`).so
#PSIINCLUDE = -I$(psi_top_objdir)/boost/boost_1_48_0 -I$(psi_top_srcdir)/include -I$(psi_top_objdir)/include -I$(psi_top_srcdir)/src/lib -I$(psi_top_objdir)/src/lib 
PSIINCLUDE = -I$(psi_top_objdir)/boost/boost_1_48_0 -I$(psi_top_srcdir)/include -I$(psi_top_objdir)/include -I$(psi_top_srcdir)/src/lib -I$(psi_top_objdir)/src/lib -DINSTALLEDPSIDATADIR=\"$(pkgdatadir)\" $(PYTHON_INCLUDE)
CXXINC += -DINSTALLEDPSIDATADIR=\"$(pkgdatadir)\" $(PYTHON_INCLUDE) -DCUDA
PSILIBS = -L$(psi_top_objdir)/lib -lPSI_plugin

CXXSRC = $(notdir $(wildcard *.cc))
CUDASRC = $(notdir $(wildcard *.cu))
DEPENDINCLUDE = $(notdir $(wildcard *.h))

BINOBJ = $(CXXSRC:%.cc=%.o)

default:: $(PSITARGET)

#GPULIB=/global/u1/d/deprince/psi4/tests/plugin_libfock/gpujk.so

# Add the flags needed for shared library creation
ifeq ($(UNAME), Linux)
    LDFLAGS = -shared $(LCUBLAS) $(HAVECUDA)
endif
ifeq ($(UNAME), Darwin)
    LDFLAGS = -shared -undefined dynamic_lookup
    CXXOTH += -fno-common
endif

# The object files
%.o: %.cc
	$(CXX) $(CXXDBG) $(CXXOPT) $(CXXDEF) $(CXXOTH) $(PSIINCLUDE) $(HAVECUDA) -c $< $(OUTPUT_OPTION)


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
%.d: %.cpp
	mpic++ $(CXXDEPENDFLAGS) $(CXXDEF) $(PSIINCLUDE) $< | sed 's/^$*.o/$*.o $*.d/g' > $(@F)

ifneq ($(DODEPEND),no)
$(BINOBJ:%.o=%.d): $(DEPENDINCLUDE)
include $(BINOBJ:%.o=%.d)
$(CUDAOBJ:%.o=%.d): $(DEPENDINCLUDE)
include $(CUDAOBJ:%.o=%.d)
endif
