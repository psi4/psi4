# Name of the executable you want to generate
TARGET = df-mp2

# List of source files to be compiled
CXXSRC = main.cc df-mp2.cc df-mp2-formJ.cc

# Location of your PSI4 source
psi_top_srcdir = /Users/andysim/programming/workspace/psi4

# Location of your PSI4 install, by default as listed
psi_top_objdir = /Users/andysim/programming/workspace/psi4_obj_debug

# Include the BLAS and LAPACK linear algebra libraries
BLAS   = -lblas
LAPACK = -llapack

# The name of the C++ compiler
CXX=g++

# C++ Compiler flags, in this case: no optimization and
# enable native GDB debugging, respectively
CXXFLAGS = -O0 -g $(INCLUDES)
#CXXFLAGS = -O2 $(INCLUDES)


#************ None of the following should be modified *********#
# Name of the object intermediates that make up the 
# executable
OBJS = $(CXXSRC:%.cc=%.o)
# Include the PSI4 include files (this needs to be cleaned up)
INCLUDES = -I$(psi_top_srcdir)/include\
-I$(psi_top_objdir)/include\
-I$(psi_top_srcdir)/src/lib\
-I$(psi_top_objdir)/src/lib
# Include the PSI4 libraries
LIBS = -L$(psi_top_objdir)/lib -lPSI_mints -lPSI_ccenergy -lPSI_ccsort\
-lPSI_input3 -lPSI_cscf -lPSI_cints -lPSI_transqt2 -lPSI_psiclean -lPSI_dpd\
-lPSI_chkpt -lPSI_iwl -lPSI_qt -lPSI_psio -lPSI_ciomr -lPSI_ipv1\
-lPSI_options -lPSI_deriv -lPSI_int -lPSI_util -lPSI_parallel $(BLAS) $(LAPACK)
# Execute the compilation for the specified object, linking
# with all the libraries specified
$(TARGET): $(OBJS)
	$(CXX) -o $(TARGET) $(OBJS) $(LIBS)

# Execute the compilation for all objects of the specified
# target
all:
	$(TARGET)

# Erase all compiled intermediate files
clean:
	rm -f $(OBJS) $(TARGET)
