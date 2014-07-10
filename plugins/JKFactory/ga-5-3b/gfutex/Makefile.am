ACLOCAL_AMFLAGS = -I m4 -I ../m4

EXTRA_DIST =

AM_FFLAGS =
AM_CFLAGS =
AM_CXXFLAGS =

AM_CPPFLAGS =
AM_CPPFLAGS += -I$(top_srcdir)/src
AM_CPPFLAGS += $(TBB_CPPFLAGS)
AM_CPPFLAGS += $(BOOST_CPPFLAGS)
AM_CPPFLAGS += $(GA_CPPFLAGS)
AM_CPPFLAGS += $(GA_MP_CPPFLAGS)

LDADD =
LDADD += libgfutex.la
LDADD += $(TBB_LIBS)
LDADD += $(GA_LIBS)
LDADD += $(GA_MP_LIBS)
LDADD += $(GA_FLIBS)

AM_LDFLAGS =
AM_LDFLAGS += $(TBB_LDFLAGS)
AM_LDFLAGS += $(GA_LDFLAGS)
AM_LDFLAGS += $(GA_MP_LDFLAGS)
AM_LDFLAGS += $(GA_FLIBS)

include_HEADERS =
include_HEADERS += src/gfutex.h

lib_LTLIBRARIES =
lib_LTLIBRARIES += libgfutex.la

libgfutex_la_SOURCES =
libgfutex_la_SOURCES += src/gf_cache.C
libgfutex_la_SOURCES += src/gfutex.C

.PHONY: checkprogs
checkprogs: $(check_PROGRAMS)

check_PROGRAMS =
check_PROGRAMS += test/test_gfutex
check_PROGRAMS += examples/scf/scf

test_test_gfutex_SOURCES = test/test_gfutex.C
examples_scf_scf_SOURCES = 
examples_scf_scf_SOURCES += examples/scf/cscc.h
examples_scf_scf_SOURCES += examples/scf/input.c
examples_scf_scf_SOURCES += examples/scf/input.h
examples_scf_scf_SOURCES += examples/scf/integ.c
examples_scf_scf_SOURCES += examples/scf/integ.h
examples_scf_scf_SOURCES += examples/scf/oneel.C
examples_scf_scf_SOURCES += examples/scf/oneel.h
examples_scf_scf_SOURCES += examples/scf/output.c
examples_scf_scf_SOURCES += examples/scf/output.h
examples_scf_scf_SOURCES += examples/scf/scf.c
examples_scf_scf_SOURCES += examples/scf/scf.h
examples_scf_scf_SOURCES += examples/scf/sum.c
examples_scf_scf_SOURCES += examples/scf/sum.h
examples_scf_scf_SOURCES += examples/scf/timer.c
examples_scf_scf_SOURCES += examples/scf/timer.h
examples_scf_scf_SOURCES += examples/scf/twoel.C
examples_scf_scf_SOURCES += examples/scf/twoel.h
