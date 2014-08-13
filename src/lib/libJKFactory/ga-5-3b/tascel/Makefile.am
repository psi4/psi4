ACLOCAL_AMFLAGS = -I m4 -I ../m4

BUILT_SOURCES =
EXTRA_DIST =

AM_FFLAGS =
AM_FFLAGS += $(FFLAG_INT)

AM_CFLAGS =

AM_CXXFLAGS =

AM_CPPFLAGS =
AM_CPPFLAGS += -I$(top_srcdir)/src
AM_CPPFLAGS += $(GA_CPPFLAGS)
AM_CPPFLAGS += $(LAPACK_CPPFLAGS)
AM_CPPFLAGS += $(BLAS_CPPFLAGS)
AM_CPPFLAGS += $(GA_MP_CPPFLAGS)

LDADD =
LDADD += libtascel.la
LDADD += $(GA_LIBS)
LDADD += $(LAPACK_LIBS)
LDADD += $(BLAS_LIBS)
LDADD += $(GA_MP_LIBS)
LDADD += $(FLIBS)

AM_LDFLAGS =
AM_LDFLAGS += $(GA_LDFLAGS)
AM_LDFLAGS += $(LAPACK_LDFLAGS)
AM_LDFLAGS += $(BLAS_LDFLAGS)
AM_LDFLAGS += $(GA_MP_LDFLAGS)

clean-local: clean-config

include_HEADERS =
include_HEADERS += src/AccessMode.h
include_HEADERS += src/Comm.h
include_HEADERS += src/Counter.h
include_HEADERS += src/DataColl.h
include_HEADERS += src/DenseArray.h
include_HEADERS += src/FuncReg.h
include_HEADERS += src/massert.h
include_HEADERS += src/SharedQueue.h
include_HEADERS += src/Sleep.h
include_HEADERS += src/SplitQueue.h
include_HEADERS += src/SplitQueueOpt.h
include_HEADERS += src/StealingStats.h
include_HEADERS += src/TaskCollProps.h
include_HEADERS += src/Timer.h
include_HEADERS += src/TerminationDetector.h
include_HEADERS += src/UniformTaskCollection.h
include_HEADERS += src/UniformTaskCollectionShared.h
include_HEADERS += src/UniformTaskCollectionSplit.h
include_HEADERS += src/UniformTaskCollSplitData.h

lib_LTLIBRARIES =
lib_LTLIBRARIES += libtascel.la

nodist_libtascel_la_SOURCES =
libtascel_la_SOURCES =
libtascel_la_SOURCES += src/Comm.cc
libtascel_la_SOURCES += src/DenseArray.cc
libtascel_la_SOURCES += src/FuncReg.cc
libtascel_la_SOURCES += src/SharedQueue.cc
libtascel_la_SOURCES += src/Sleep.cc
libtascel_la_SOURCES += src/SplitQueue.cc
libtascel_la_SOURCES += src/SplitQueueOpt.cc
libtascel_la_SOURCES += src/StealingStats.cc
libtascel_la_SOURCES += src/TaskCollProps.cc
libtascel_la_SOURCES += src/TerminationDetector.cc
libtascel_la_SOURCES += src/UniformTaskCollection.cc
libtascel_la_SOURCES += src/UniformTaskCollectionShared.cc
libtascel_la_SOURCES += src/UniformTaskCollectionSplit.cc
libtascel_la_SOURCES += src/UniformTaskCollSplitData.cc

BUILT_SOURCES += config.fh
nodist_libtascel_la_SOURCES += config.fh

config.fh: $(top_build_prefix)config.h
	@-rm -f config.fh
	$(AM_V_SED) \
	$(SED) "/^#/!d" $(top_build_prefix)config.h > config.fh

clean-config:
	-rm -f config.fh

check_PROGRAMS =
check_PROGRAMS += test/TestSharedQueue
check_PROGRAMS += test/TestSplitQueue
if HAVE_BLAS
check_PROGRAMS += examples/scf/scf
endif

test_TestSharedQueue_SOURCES = test/TestQueue.cc
test_TestSharedQueue_CPPFLAGS = $(AM_CPPFLAGS) -DSHARED_QUEUE
test_TestSplitQueue_SOURCES = test/TestQueue.cc
test_TestSplitQueue_CPPFLAGS = $(AM_CPPFLAGS) -DSPLIT_QUEUE

examples_scf_scf_SOURCES = 
examples_scf_scf_SOURCES += examples/scf/cscf.h
examples_scf_scf_SOURCES += examples/scf/input.F
examples_scf_scf_SOURCES += examples/scf/integ.F
examples_scf_scf_SOURCES += examples/scf/main.c
examples_scf_scf_SOURCES += examples/scf/output.F
examples_scf_scf_SOURCES += examples/scf/scf.F
examples_scf_scf_SOURCES += examples/scf/timer.F
examples_scf_scf_SOURCES += examples/scf/twoelcpp.cc

EXTRA_DIST += examples/scf/README
EXTRA_DIST += examples/scf/be.inpt
EXTRA_DIST += examples/scf/be16.inpt

.PHONY: checkprogs
checkprogs: $(check_PROGRAMS)

TESTS = $(SERIAL_TESTS) $(PARALLEL_TESTS)
SERIAL_TESTS =
PARALLEL_TESTS =
PARALLEL_TESTS += test/TestSharedQueue$(EXEEXT)
PARALLEL_TESTS += test/TestSplitQueue$(EXEEXT)
if HAVE_BLAS
PARALLEL_TESTS += examples/scf/scf$(EXEEXT)
endif

if CROSS_COMPILING
LOG_COMPILER = \
maybe_mpiexec=`if echo "$(SERIAL_TESTS)" | $(GREP) "$$p" > /dev/null; then echo "$(MPIEXEC)" | $(SED) 's/%NP%/1/'; else echo "$(MPIEXEC)" | $(SED) 's/%NP%/$(NPROCS)/'; fi`; $$maybe_mpiexec
else
LOG_COMPILER = \
maybe_mpiexec=`if echo "$(SERIAL_TESTS)" | $(GREP) "$$p" > /dev/null; then echo ""; else if echo "$(MPIEXEC)" | $(GREP) "%NP%" > /dev/null; then echo "$(MPIEXEC)" | $(SED) 's/%NP%/$(NPROCS)/'; else echo "$(MPIEXEC)"; fi; fi`; $$maybe_mpiexec
endif # CROSS_COMPILING

# support verbose/silent make rules for additional programs
# sed
AM_V_SED = $(am__v_SED_$(V))
am__v_SED_ = $(am__v_SED_$(AM_DEFAULT_VERBOSITY))
am__v_SED_0 = @echo "  SED   " $@;
