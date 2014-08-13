# CCAFE_HOME= To be set as an Environment variable:
# Eg: export CCAFE_HOME=/msrc/proj/cca/manoj/dccafe/cxx

# Here we must contend with two kinds of .o file --
# those built normally and those built for use in
# dynamic libraries. To that end we delete .o files
# built for DL as soon as we've converted them to .so.

zeroth: first

EG_ROOT=$(CCA_ROOT)/cxx/eg

include $(CCAFE_HOME)/Makefile.Rules
include $(CCAFE_HOME)/MakeIncl.CCAFE_Component
include $(EG_ROOT)/ccaComponent/MakeIncl.EGComponent

#OTHERFLAGS = -D_DBG_ILOOP -D_DBG_PLATE1

ifdef LINUX64_QUADRICS
  EXTRA_LIBS = -lshmem -lelan3 -lelan -lpthread
endif

GA_INC  = -I$(GA_HOME)/include
GA_LIBS = -L$(GA_HOME)/lib/$(TARGET) -L$(MPI_LIB) -ldl -lc -lglobal\
        -lma -llinalg -larmci -ltcgmsg-mpi $(LIBMPI)\
        -lg2c -lnsl $(EXTRA_LIBS) -lm

LOCALFLAGS= $(OTHERFLAGS) $(EG_INC) $(GA_INC) -I$(MPI_INCLUDE) -I.
WRAPPERFLAGS= -w
WRAPPERFLAGS= 

SHARED_FLAGS = -fpic

%.o : %.C
	$(CXX) -c $(SHARED_FLAGS) $(CXXFLAGS)  $<

%.o : %.cxx
	$(CXX) -c $(SHARED_FLAGS) $(CXXFLAGS) $(LOCALFLAGS) $<

%.i : %.cxx
	-$(RM) $@
	$(CXX) $(CPPFLAGS) $(LOCALFLAGS) $< | egrep -v ^$$ |egrep -v '^# ' > $@


%.class : %.java
	$(JAVAC)  $(JAVAFLAGS)  $<
%.h : %.class
	$(JAVAH) $(JAVAHFLAGS) $(<:%.class=%)

JAVASRC = 

# dl stuff that is not ccaffeine dependent
CXXDLSTRICTSRC= \
TestComponent.cxx

# dl stuff that is ccaffeine dependent
CXXDLLOOSESRC= 

# static link stuff
CXXDLNONESRC=

# stuff built with or without a parallel environment.
CXXSRC = $(CXXDLSTRICTSRC) $(CXXDLLOOSESRC) $(CXXDLNONESRC)

NEWSRC = 


CXXOBJ = $(CXXSRC:%.cxx=%.o)
CXXDSSOBJ = $(CXXDLSTRICTSRC:%.cxx=%.so)
CXXDLSOBJ = $(CXXDLLOOSESRC:%.cxx=%.so)
NEWOBJ = $(NEWSRC:%.cxx=%.o)

$(CXXDSSOBJ) : %$(CCAFE_SHLIB_SUFFIX) : %.cxx
	-$(RM) $@ $*.o 
	#$(CCAFE_HOME)/dc/component/genDLWrapperStrict $* > $*_wrapper.cxx
	#$(CCAFE_HOME)/dc/component/genDLIndex $@ create_$* $* > $*.cca
	$(CXX) -c $(CXXFLAGS) $(LOCALFLAGS) $(CCAFE_SHLIB_CFLAGS) $< 
	$(CXX) -c $(CXXFLAGS) $(LOCALFLAGS) $(CCAFE_SHLIB_CFLAGS) \
		$(WRAPPERFLAGS) $*_wrapper.cxx
	$(CCAFE_SHLIB_LD) $(CCAFE_SHLIB_LDFLAGS) -o $@ $*.o $*_wrapper.o \
		$(NEWOBJ) $(GA_LIBS)
	-$(RM) $*.o $(NEWOBJ) $*_wrapper.o

LIB= libcomponent.a

first all: lib dllib
	@echo done $(CXXOBJ) $(JAVACLS) $(JNIHDR)

lib : allobj
	ar $(ARFLAGS) $(LIB) $(CXXOBJ) $(NEWOBJ)
	$(RANLIB) $(LIB)

dllib: $(CXXDSSOBJ) $(CXXDLSOBJ)

allobj : $(CXXOBJ) $(JAVACLS) $(JNIHDR) $(NEWOBJ)

new: $(NEWOBJ)

$(CXXOBJ) : $(CXXSRC)


clean:
	$(RM) -f $(CXXOBJ) $(LIB) $(JAVACLS) $(JNIHDR) $(NEWOBJ) *.so *_wrapper.o *.i *~

depend:
	$(DEPEND) $(CXXFLAGS) $(LOCALFLAGS) $(CXXSRC) $(NEWSRC)

depend-sys:
	$(DEPENDSYS) $(CXXFLAGS) $(LOCALFLAGS) $(CXXSRC) $(NEWSRC)

#include Makefile.depends
