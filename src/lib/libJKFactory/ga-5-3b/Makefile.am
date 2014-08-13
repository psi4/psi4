##############################################################################
# Makefile.am for Global Arrays.
#
# Rationale:
# This Makefile.am follows many of the suggestions outlined in the paper
# "Recursive Make Considered Harmful".  We do not use automake's
# 'include' feature (instead preferring a single, large Makefile.am).
# We do recurse into the armci directory, but since armci is a stand-alone
# package this is the correct behavior.
#
# Additional targets:
# Besides the traditional make targets supplied by automake, we have added the
# "checkprogs" target to build example programs and test programs.
#
# Notes:
# In general, each subdirectory has a corresponding section within this
# Makefile.am with the notable exception being the many examples getting
# rolled up into the global/examples section.
#

# The usual aclocal nonsense to get include paths right.
ACLOCAL_AMFLAGS = -I m4

# All public headers, installed programs, test programs, and example programs
# are listed in these variables.  Appended to throughout.  These are the
# automake variables used.
include_HEADERS =
nodist_include_HEADERS =
bin_PROGRAMS =
bin_SCRIPTS =
check_PROGRAMS =
check_LTLIBRARIES =
lib_LTLIBRARIES = libga.la
libga_la_SOURCES =
libga_la_LIBADD =
nodist_libga_la_SOURCES =
EXTRA_DIST = README
BUILT_SOURCES =
MOSTLYCLEANFILES =
CLEANFILES =
DISTCLEANFILES =
MAINTAINERCLEANFILES =
AM_FFLAGS =
AM_CFLAGS =
AM_CXXFLAGS =
AM_CPPFLAGS =
AM_LDFLAGS =
FLDFLAGS =
LDADD =

AM_FFLAGS += $(GA_FOPT)
AM_FFLAGS += $(GA_F_WARN)
AM_FFLAGS += $(FFLAG_INT)
AM_CFLAGS += $(GA_COPT)
AM_CFLAGS += $(GA_C_WARN)
AM_CXXFLAGS += $(GA_CXXOPT)
AM_CXXFLAGS += $(GA_CXX_WARN)

AM_CPPFLAGS += $(SCALAPACK_CPPFLAGS)
AM_CPPFLAGS += $(LAPACK_CPPFLAGS)
AM_CPPFLAGS += $(BLAS_CPPFLAGS)
AM_CPPFLAGS += $(GA_MP_CPPFLAGS)
AM_CPPFLAGS += $(ARMCI_NETWORK_CPPFLAGS)

AM_LDFLAGS += $(SCALAPACK_LDFLAGS)
AM_LDFLAGS += $(LAPACK_LDFLAGS)
AM_LDFLAGS += $(BLAS_LDFLAGS)
AM_LDFLAGS += $(GA_MP_LDFLAGS)
AM_LDFLAGS += $(ARMCI_NETWORK_LDFLAGS)
if ENABLE_F77
AM_LDFLAGS += $(FLD_NOMAIN)
endif

FLDFLAGS += $(SCALAPACK_LDFLAGS)
FLDFLAGS += $(LAPACK_LDFLAGS)
FLDFLAGS += $(BLAS_LDFLAGS)
FLDFLAGS += $(GA_MP_LDFLAGS)
FLDFLAGS += $(ARMCI_NETWORK_LDFLAGS)

LDADD += $(LTLIBOBJS)
LDADD += libga.la
LDADD += $(SCALAPACK_LIBS)
LDADD += $(LAPACK_LIBS)
LDADD += $(BLAS_LIBS)
LDADD += $(GA_MP_LIBS)
LDADD += $(ARMCI_NETWORK_LIBS)

# Certain trickery when turning Fortran support on or off.
if ENABLE_F77
MAYBE_DUMMY_F_LINK = dummy.f
MAYBE_FLIBS = $(FLIBS)
endif

##############################################################################
# compiler and linker flags
#
# Important for external tools wanting to know how to link to GA. This
# functionality caries over from the pre-autotools GA build.
SED_NORMALIZE_WHITESPACE = $(SED) 's/ [ ]*/ /g;s/" /"/g;s/ "/"/g'
.PHONY: flags
flags:
	@echo '# =========================================================================== '
	@echo 'F77="$(F77)"'
	@echo 'CC="$(CC)"'
if CXX_BINDINGS
	@echo 'CXX="$(CXX)"'
endif
	@echo '# Suggested compiler/linker options are as follows.'
	@echo '# GA libraries are installed in $(libdir)' 
	@echo '# GA headers are installed in $(includedir)' 
	@echo '#'
	@echo 'CPPFLAGS="$(SCALAPACK_CPPFLAGS) $(LAPACK_CPPFLAGS) $(BLAS_CPPFLAGS) $(GA_MP_CPPFLAGS) $(ARMCI_NETWORK_CPPFLAGS) -I$(includedir)"' | $(SED_NORMALIZE_WHITESPACE)
	@echo '#'
	@echo 'LDFLAGS="$(SCALAPACK_LDFLAGS) $(LAPACK_LDFLAGS) $(BLAS_LDFLAGS) $(GA_MP_LDFLAGS) $(ARMCI_NETWORK_LDFLAGS) -L$(libdir)"' | $(SED_NORMALIZE_WHITESPACE)
	@echo '#'
	@echo '# For Fortran Programs: '
	@echo 'FFLAGS="$(AM_FFLAGS)"' | $(SED_NORMALIZE_WHITESPACE)
	@echo 'LIBS="-lga $(SCALAPACK_LIBS) $(LAPACK_LIBS) $(BLAS_LIBS) $(GA_MP_LIBS) $(ARMCI_NETWORK_LIBS)"' | $(SED_NORMALIZE_WHITESPACE)
	@echo '#'
	@echo '# For C Programs: '
	@echo 'CFLAGS="$(AM_CFLAGS)"' | $(SED_NORMALIZE_WHITESPACE)
	@echo 'LIBS="-lga $(SCALAPACK_LIBS) $(LAPACK_LIBS) $(BLAS_LIBS) $(GA_MP_LIBS) $(ARMCI_NETWORK_LIBS) $(MAYBE_FLIBS)"' | $(SED_NORMALIZE_WHITESPACE)
if CXX_BINDINGS
	@echo '#'
	@echo '# For C++ Programs: '
	@echo 'CXXFLAGS="$(AM_CXXFLAGS)"' | $(SED_NORMALIZE_WHITESPACE)
	@echo 'LIBS="-lga++ -lga $(SCALAPACK_LIBS) $(LAPACK_LIBS) $(BLAS_LIBS) $(GA_MP_LIBS) $(ARMCI_NETWORK_LIBS) $(MAYBE_FLIBS)"' | $(SED_NORMALIZE_WHITESPACE)
endif
	@echo '# =========================================================================== '

bin_SCRIPTS += tools/ga-config
CLEANFILES += $(bin_SCRIPTS)

##############################################################################
# config
#
# The config.h file is automatically generated from config.h.in at configure
# time.  However, it contains symbols which certain Fortran compilers do not
# understand.  Therefore, we create a config.fh file which is Fortran-specific
# and works with all Fortran compilers.
#
BUILT_SOURCES += config.fh
CLEANFILES += config.fh
nodist_libga_la_SOURCES += config.fh

config.fh: config.h
	@-rm -f $@
	$(SED_V) \
	$(SED) "/^#/!d" config.h > $@

##############################################################################
# compat
#
# Although the compat directory houses replacements for missing or erroneous
# standard C functions and such sources are conditionally compiled based on
# results from configure tests, without the "random" implementation the
# m4-generated tests always fail for scatter and copy_patch.
libga_la_SOURCES += compat/random.c

##############################################################################
# ma
#
libga_la_SOURCES += ma/error.c
libga_la_SOURCES += ma/f2c.c
libga_la_SOURCES += ma/ma.c
libga_la_SOURCES += ma/ma.h
libga_la_SOURCES += ma/memcpy.h
libga_la_SOURCES += ma/string-util.c
libga_la_SOURCES += ma/table.c
if ENABLE_F77
libga_la_SOURCES += ma/maf.F
endif

include_HEADERS += ma/error.h
include_HEADERS += ma/macdecls.h
include_HEADERS += ma/macommon.h
include_HEADERS += ma/maf2c.fh
include_HEADERS += ma/mafdecls.fh
include_HEADERS += ma/matypes.h
include_HEADERS += ma/scope.h
include_HEADERS += ma/string-util.h
include_HEADERS += ma/table.h

AM_CPPFLAGS += -I$(top_build_prefix)ma
AM_CPPFLAGS += -I$(top_srcdir)/ma

check_PROGRAMS += ma/testc
check_PROGRAMS += ma/test-coalesce
check_PROGRAMS += ma/test-inquire
if ENABLE_F77
check_PROGRAMS += ma/testf
endif

MA_SERIAL_TESTS =
MA_SERIAL_TESTS_XFAIL =
MA_PARALLEL_TESTS =
MA_PARALLEL_TESTS_XFAIL =
MA_TESTS = $(MA_SERIAL_TESTS) $(MA_PARALLEL_TESTS)
MA_TESTS_XFAIL = $(MA_SERIAL_TESTS_XFAIL) $(MA_PARALLEL_TESTS_XFAIL)

#MA_SERIAL_TESTS += ma/testc$(EXEEXT) # iteractive prompt
MA_SERIAL_TESTS += ma/test-coalesce$(EXEEXT)
MA_SERIAL_TESTS += ma/test-inquire$(EXEEXT)
if ENABLE_F77
MA_SERIAL_TESTS += ma/testf$(EXEEXT)
MA_SERIAL_TESTS_XFAIL += ma/testf$(EXEEXT)
endif

ma_testf_SOURCES         = ma/testf.F
ma_testc_SOURCES         = ma/testc.c
ma_test_coalesce_SOURCES = ma/test-coalesce.c
ma_test_inquire_SOURCES  = ma/test-inquire.c

ma_testf_LDFLAGS = $(FLDFLAGS)

nodist_EXTRA_ma_testc_SOURCES         = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_ma_test_coalesce_SOURCES = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_ma_test_inquire_SOURCES  = $(MAYBE_DUMMY_F_LINK)

##############################################################################
# LinAlg/lapack+blas
#
# Since we gave all linalg routines a gal_ prefix, we can unconditionally
# compile our internal linalg routines.
#
if ENABLE_F77
libga_la_SOURCES += LinAlg/lapack+blas/dgetf2.f
libga_la_SOURCES += LinAlg/lapack+blas/dgetrf.f
libga_la_SOURCES += LinAlg/lapack+blas/dgetrs.f
libga_la_SOURCES += LinAlg/lapack+blas/dlacpy.f
libga_la_SOURCES += LinAlg/lapack+blas/dlaset.f
libga_la_SOURCES += LinAlg/lapack+blas/dlaswp.f
libga_la_SOURCES += LinAlg/lapack+blas/ilaenv.f
libga_la_SOURCES += LinAlg/lapack+blas/lsame.f
libga_la_SOURCES += LinAlg/lapack+blas/xerbla.f
libga_la_SOURCES += LinAlg/lapack+blas/cgemm.f
libga_la_SOURCES += LinAlg/lapack+blas/dblas.f
libga_la_SOURCES += LinAlg/lapack+blas/sgemm.f
libga_la_SOURCES += LinAlg/lapack+blas/zaxpy.f
libga_la_SOURCES += LinAlg/lapack+blas/zgemm.f
endif
libga_la_SOURCES += LinAlg/lapack+blas/xgemm.c
libga_la_SOURCES += LinAlg/lapack+blas/xgemm.h
libga_la_SOURCES += LinAlg/lapack+blas/galinalg.h
libga_la_SOURCES += LinAlg/lapack+blas/galinalg.fh

AM_CPPFLAGS += -I$(top_srcdir)/LinAlg/lapack+blas

EXTRA_DIST += LinAlg/lapack+blas/README

##############################################################################
# global
#
EXTRA_DIST += global/README

##############################################################################
# tools
#
# The following to support tracing, profiling, etc.
#
if ENABLE_PROFILING
if HAVE_SYS_WEAK_ALIAS_PRAGMA
lib_LTLIBRARIES += libwapi.la
libwapi_la_SOURCES =
libwapi_la_SOURCES += tools/ga-wapi.c
if ENABLE_GPARRAYS
libwapi_la_SOURCES += tools/gp-wapi.c
endif
else
nodist_libga_la_SOURCES += tools/ga-wapi.c
if ENABLE_GPARRAYS
nodist_libga_la_SOURCES += tools/gp-wapi.c
endif
endif
endif

if ENABLE_GPARRAYS
EXTRA_DIST += tools/gpwapigen.py
endif
EXTRA_DIST += tools/wapigen.py
EXTRA_DIST += tools/wapigen_counts.py
EXTRA_DIST += tools/wapigen_trace.py

##############################################################################
# global/src
#
PAPI_H        = $(top_srcdir)/global/src/ga-papi.h

global/src/ga-wapi.h: $(PAPI_H)
	@-rm -f $@
	$(SED_V) \
	$(SED) -e "s/PAPI/WAPI/g" \
	       -e "s/pnga_/wnga_/g" \
	       -e "/typedef/d" \
		   $(PAPI_H) > $@

global/src/ga-wapidefs.h: $(PAPI_H)
	@-rm -f $@
	$(SED_V) \
	$(SED) -e "s/PAPI/WAPIDEFS/g" \
		   -e "s/^.*p\(nga_.*\)(.*$$/#define w\1 p\1/g" \
		   -e "/^ /d" \
		   -e "/^\/\*/d" \
		   -e "/^ \*/d" \
	       -e "/typedef/d" \
		   $(PAPI_H) > $@

BUILT_SOURCES += global/src/ga-wapi.h
BUILT_SOURCES += global/src/ga-wapidefs.h

CLEANFILES += global/src/ga-wapi.h
CLEANFILES += global/src/ga-wapidefs.h

libga_la_SOURCES += global/src/abstract_ops.h
libga_la_SOURCES += global/src/base.c
libga_la_SOURCES += global/src/base.h
libga_la_SOURCES += global/src/capi.c
libga_la_SOURCES += global/src/cnames.h
libga_la_SOURCES += global/src/collect.c
libga_la_SOURCES += global/src/datatypes.c
libga_la_SOURCES += global/src/decomp.c
libga_la_SOURCES += global/src/diag.fh
libga_la_SOURCES += global/src/DP.c
libga_la_SOURCES += global/src/elem_alg.c
libga_la_SOURCES += global/src/fapi.c
libga_la_SOURCES += global/src/ga_ckpt.h
libga_la_SOURCES += global/src/gaconfig.h
libga_la_SOURCES += global/src/ga_diag_seqc.c
libga_la_SOURCES += global/src/ga_malloc.c
libga_la_SOURCES += global/src/ga_profile.h
libga_la_SOURCES += global/src/ga_solve_seq.c
libga_la_SOURCES += global/src/ga_symmetr.c
libga_la_SOURCES += global/src/ga_trace.c
libga_la_SOURCES += global/src/ghosts.c
libga_la_SOURCES += global/src/global.h
libga_la_SOURCES += global/src/global.nalg.c
libga_la_SOURCES += global/src/global.npatch.c
libga_la_SOURCES += global/src/global.periodic.c
libga_la_SOURCES += global/src/globalp.h
libga_la_SOURCES += global/src/global.util.c
libga_la_SOURCES += global/src/hsort.scat.c
libga_la_SOURCES += global/src/matmul.c
libga_la_SOURCES += global/src/matmul.h
libga_la_SOURCES += global/src/matrix.c
libga_la_SOURCES += global/src/nbutil.c
libga_la_SOURCES += global/src/onesided.c
libga_la_SOURCES += global/src/peigstubs.c
libga_la_SOURCES += global/src/scalapack.fh
libga_la_SOURCES += global/src/sclstubs.c
libga_la_SOURCES += global/src/select.c
libga_la_SOURCES += global/src/sparse.c
libga_la_SOURCES += global/src/types.xh
libga_la_SOURCES += global/src/types2.xh
nodist_libga_la_SOURCES += global/src/ga-wapidefs.h
if ENABLE_F77
libga_la_SOURCES += global/src/complex.F
libga_la_SOURCES += global/src/ga_diag_seq.F
libga_la_SOURCES += global/src/rsg.f
if CRAY_XT_NETWORKS
libga_la_SOURCES += global/src/ga_dgemmf.F
endif # CRAY_XT_NETWORKS
if ENABLE_PEIGS
libga_la_SOURCES += global/src/ga_diag.F
endif # ENABLE_PEIGS
if HAVE_SCALAPACK
libga_la_SOURCES += global/src/scalapack.F
endif # HAVE_SCALAPACK
endif # ENABLE_F77
if ENABLE_CHECKPOINT
libga_la_SOURCES += global/src/ga_ckpt.c
endif

include_HEADERS += global/src/gacommon.h
include_HEADERS += global/src/ga.h
include_HEADERS += global/src/global.fh
include_HEADERS += global/src/ga-papi.h
if MSG_COMMS_MPI
include_HEADERS += global/src/ga-mpi.h
include_HEADERS += global/src/ga-mpi.fh
endif
nodist_include_HEADERS += global/src/ga-wapi.h

AM_CPPFLAGS += -I$(top_build_prefix)global/src
AM_CPPFLAGS += -I$(top_srcdir)/global/src

##############################################################################
# global/testing
#
AM_CPPFLAGS += -I$(top_srcdir)/global/testing

check_PROGRAMS += global/testing/big
check_PROGRAMS += global/testing/elempatch
check_PROGRAMS += global/testing/gatscat
check_PROGRAMS += global/testing/getmem
check_PROGRAMS += global/testing/mtest
check_PROGRAMS += global/testing/mulmatpatchc
check_PROGRAMS += global/testing/ntestc
check_PROGRAMS += global/testing/ntestfc
check_PROGRAMS += global/testing/packc
check_PROGRAMS += global/testing/patch_enumc
check_PROGRAMS += global/testing/perf2
check_PROGRAMS += global/testing/print
check_PROGRAMS += global/testing/scan_addc
check_PROGRAMS += global/testing/scan_copyc
check_PROGRAMS += global/testing/testc
check_PROGRAMS += global/testing/testmatmultc
check_PROGRAMS += global/testing/testmult
check_PROGRAMS += global/testing/testmultrect
check_PROGRAMS += global/testing/unpackc
if ENABLE_F77
check_PROGRAMS += global/testing/bin
check_PROGRAMS += global/testing/blktest
check_PROGRAMS += global/testing/d2test
check_PROGRAMS += global/testing/g2test
check_PROGRAMS += global/testing/g3test
check_PROGRAMS += global/testing/ga_lu
check_PROGRAMS += global/testing/ga_shift
check_PROGRAMS += global/testing/ghosts
check_PROGRAMS += global/testing/jacobi
check_PROGRAMS += global/testing/mir_perf2
check_PROGRAMS += global/testing/mmatrix
check_PROGRAMS += global/testing/mulmatpatch
check_PROGRAMS += global/testing/nb2test
check_PROGRAMS += global/testing/ndim
check_PROGRAMS += global/testing/patch
check_PROGRAMS += global/testing/patch2
check_PROGRAMS += global/testing/patch_enumf
check_PROGRAMS += global/testing/perfmod
check_PROGRAMS += global/testing/perform
check_PROGRAMS += global/testing/perf
check_PROGRAMS += global/testing/pg2testmatmult
check_PROGRAMS += global/testing/pg2test
check_PROGRAMS += global/testing/pgtestmatmult
check_PROGRAMS += global/testing/pgtest
check_PROGRAMS += global/testing/random
check_PROGRAMS += global/testing/scan
check_PROGRAMS += global/testing/simple_groups
check_PROGRAMS += global/testing/sparse
check_PROGRAMS += global/testing/sprsmatmult
check_PROGRAMS += global/testing/stride
check_PROGRAMS += global/testing/testeig
check_PROGRAMS += global/testing/testmatmult
check_PROGRAMS += global/testing/testsolve
check_PROGRAMS += global/testing/test
check_PROGRAMS += global/testing/test_mirrored
check_PROGRAMS += global/testing/types_test
check_PROGRAMS += global/testing/field_test
if HAVE_SCALAPACK
check_PROGRAMS += global/testing/testspd
endif # HAVE_SCALAPACK
else # !ENABLE_F77
if HAVE_LAPACK
check_PROGRAMS += global/testing/ga_lu
endif
endif # ENABLE_F77
if MSG_COMMS_MPI
check_PROGRAMS += global/testing/ga-mpi
check_PROGRAMS += global/testing/lock
check_PROGRAMS += global/testing/simple_groups_comm
check_PROGRAMS += global/testing/simple_groups_commc
endif
if SYSV
check_PROGRAMS += global/testing/ipc_clean
endif

# TODO somehow unit tests depend on MPI -- need to fix
if MSG_COMMS_MPI
check_PROGRAMS += global/testing/unit-tests/ga_abs_value
check_PROGRAMS += global/testing/unit-tests/ga_acc
check_PROGRAMS += global/testing/unit-tests/ga_add
check_PROGRAMS += global/testing/unit-tests/ga_add_constant
check_PROGRAMS += global/testing/unit-tests/ga_add_constantpatch
check_PROGRAMS += global/testing/unit-tests/ga_add_diagonal
check_PROGRAMS += global/testing/unit-tests/ga_add_patch
check_PROGRAMS += global/testing/unit-tests/ga_copy2
check_PROGRAMS += global/testing/unit-tests/ga_copy3
check_PROGRAMS += global/testing/unit-tests/ga_copy
check_PROGRAMS += global/testing/unit-tests/ga_copypatch2
check_PROGRAMS += global/testing/unit-tests/ga_copypatch
check_PROGRAMS += global/testing/unit-tests/ga_create1
check_PROGRAMS += global/testing/unit-tests/ga_create2
check_PROGRAMS += global/testing/unit-tests/ga_create3
check_PROGRAMS += global/testing/unit-tests/ga_create
check_PROGRAMS += global/testing/unit-tests/ga_create_handle
check_PROGRAMS += global/testing/unit-tests/ga_create_irreg2
check_PROGRAMS += global/testing/unit-tests/ga_create_irreg3
check_PROGRAMS += global/testing/unit-tests/ga_create_irreg
check_PROGRAMS += global/testing/unit-tests/ga_destroy
check_PROGRAMS += global/testing/unit-tests/ga_dgop
check_PROGRAMS += global/testing/unit-tests/ga_dot
check_PROGRAMS += global/testing/unit-tests/ga_duplicate
check_PROGRAMS += global/testing/unit-tests/ga_elem_divide
check_PROGRAMS += global/testing/unit-tests/ga_elem_dividepatch
check_PROGRAMS += global/testing/unit-tests/ga_elem_maximum
check_PROGRAMS += global/testing/unit-tests/ga_elem_maximumpatch
check_PROGRAMS += global/testing/unit-tests/ga_elem_minimum
check_PROGRAMS += global/testing/unit-tests/ga_elem_minimumpatch
check_PROGRAMS += global/testing/unit-tests/ga_elem_multiply
check_PROGRAMS += global/testing/unit-tests/ga_elem_multiplypatch
check_PROGRAMS += global/testing/unit-tests/ga_fill
check_PROGRAMS += global/testing/unit-tests/ga_fillpatch1
check_PROGRAMS += global/testing/unit-tests/ga_fillpatch
check_PROGRAMS += global/testing/unit-tests/ga_gather2
check_PROGRAMS += global/testing/unit-tests/ga_gather3
check_PROGRAMS += global/testing/unit-tests/ga_gather
check_PROGRAMS += global/testing/unit-tests/ga_get_blockinfo
check_PROGRAMS += global/testing/unit-tests/ga_get
check_PROGRAMS += global/testing/unit-tests/ga_get_diagonal
check_PROGRAMS += global/testing/unit-tests/ga_igop2
check_PROGRAMS += global/testing/unit-tests/ga_igop
check_PROGRAMS += global/testing/unit-tests/ga_inquire
check_PROGRAMS += global/testing/unit-tests/ga_intialize
check_PROGRAMS += global/testing/unit-tests/ga_lgop
check_PROGRAMS += global/testing/unit-tests/ga_median
check_PROGRAMS += global/testing/unit-tests/ga_ndim2
check_PROGRAMS += global/testing/unit-tests/ga_ndim
check_PROGRAMS += global/testing/unit-tests/ga_nnodes
check_PROGRAMS += global/testing/unit-tests/ga_nodeid
check_PROGRAMS += global/testing/unit-tests/ga_pgroup_create2
check_PROGRAMS += global/testing/unit-tests/ga_pgroup_create3
check_PROGRAMS += global/testing/unit-tests/ga_pgroup_create4
check_PROGRAMS += global/testing/unit-tests/ga_pgroup_create5
check_PROGRAMS += global/testing/unit-tests/ga_pgroup_create
check_PROGRAMS += global/testing/unit-tests/ga_pgroup_destroy2
check_PROGRAMS += global/testing/unit-tests/ga_pgroup_destroy
check_PROGRAMS += global/testing/unit-tests/ga_pgroup_nnodes_nodeid
check_PROGRAMS += global/testing/unit-tests/ga_pgroup_setdefault
check_PROGRAMS += global/testing/unit-tests/ga_put2
check_PROGRAMS += global/testing/unit-tests/ga_put
check_PROGRAMS += global/testing/unit-tests/ga_scale2
check_PROGRAMS += global/testing/unit-tests/ga_scale
check_PROGRAMS += global/testing/unit-tests/ga_scale_cols
check_PROGRAMS += global/testing/unit-tests/ga_scale_patch
check_PROGRAMS += global/testing/unit-tests/ga_scale_rows
check_PROGRAMS += global/testing/unit-tests/ga_scatter
check_PROGRAMS += global/testing/unit-tests/ga_set_data
check_PROGRAMS += global/testing/unit-tests/ga_set_diagonal
check_PROGRAMS += global/testing/unit-tests/ga_set_restricted
check_PROGRAMS += global/testing/unit-tests/ga_solve
check_PROGRAMS += global/testing/unit-tests/ga_sync
check_PROGRAMS += global/testing/unit-tests/ga_transpose2
check_PROGRAMS += global/testing/unit-tests/ga_transpose3
check_PROGRAMS += global/testing/unit-tests/ga_transpose
check_PROGRAMS += global/testing/unit-tests/ga_zero
check_PROGRAMS += global/testing/unit-tests/ga_zerodiagonal
check_PROGRAMS += global/testing/unit-tests/ga_zeropatch2
check_PROGRAMS += global/testing/unit-tests/ga_zeropatch
endif

GLOBAL_SERIAL_TESTS =
GLOBAL_SERIAL_TESTS_XFAIL =
GLOBAL_PARALLEL_TESTS =
GLOBAL_PARALLEL_TESTS_XFAIL =
GLOBAL_TESTS = $(GLOBAL_SERIAL_TESTS) $(GLOBAL_PARALLEL_TESTS)
GLOBAL_TESTS_XFAIL = $(GLOBAL_SERIAL_TESTS_XFAIL) $(GLOBAL_PARALLEL_TESTS_XFAIL)

#GLOBAL_PARALLEL_TESTS += global/testing/big$(EXEEXT) # needs lots of memory
GLOBAL_PARALLEL_TESTS += global/testing/elempatch$(EXEEXT)
#GLOBAL_PARALLEL_TESTS += global/testing/gatscat$(EXEEXT) # broken
GLOBAL_PARALLEL_TESTS += global/testing/getmem$(EXEEXT)
GLOBAL_PARALLEL_TESTS += global/testing/mtest$(EXEEXT)
GLOBAL_PARALLEL_TESTS += global/testing/mulmatpatchc$(EXEEXT)
GLOBAL_PARALLEL_TESTS += global/testing/ntestc$(EXEEXT)
GLOBAL_PARALLEL_TESTS += global/testing/ntestfc$(EXEEXT)
GLOBAL_PARALLEL_TESTS += global/testing/packc$(EXEEXT)
GLOBAL_PARALLEL_TESTS += global/testing/patch_enumc$(EXEEXT)
GLOBAL_PARALLEL_TESTS += global/testing/print$(EXEEXT)
GLOBAL_PARALLEL_TESTS += global/testing/scan_addc$(EXEEXT)
GLOBAL_PARALLEL_TESTS += global/testing/scan_copyc$(EXEEXT)
GLOBAL_PARALLEL_TESTS += global/testing/testc$(EXEEXT)
GLOBAL_PARALLEL_TESTS += global/testing/testmatmultc$(EXEEXT)
GLOBAL_PARALLEL_TESTS += global/testing/testmult$(EXEEXT)
GLOBAL_PARALLEL_TESTS += global/testing/testmultrect$(EXEEXT)
GLOBAL_PARALLEL_TESTS += global/testing/unpackc$(EXEEXT)
if ENABLE_F77
GLOBAL_PARALLEL_TESTS += global/testing/bin$(EXEEXT)
GLOBAL_PARALLEL_TESTS += global/testing/blktest$(EXEEXT)
#GLOBAL_PARALLEL_TESTS += global/testing/d2test$(EXEEXT) # needs input file
GLOBAL_PARALLEL_TESTS += global/testing/g2test$(EXEEXT)
GLOBAL_PARALLEL_TESTS += global/testing/g3test$(EXEEXT)
GLOBAL_PARALLEL_TESTS += global/testing/ga_lu$(EXEEXT)
GLOBAL_PARALLEL_TESTS += global/testing/ga_shift$(EXEEXT)
GLOBAL_PARALLEL_TESTS += global/testing/ghosts$(EXEEXT)
GLOBAL_PARALLEL_TESTS += global/testing/jacobi$(EXEEXT)
GLOBAL_PARALLEL_TESTS += global/testing/mir_perf2$(EXEEXT)
GLOBAL_PARALLEL_TESTS += global/testing/mmatrix$(EXEEXT)
GLOBAL_PARALLEL_TESTS += global/testing/mulmatpatch$(EXEEXT)
GLOBAL_PARALLEL_TESTS += global/testing/nb2test$(EXEEXT)
GLOBAL_PARALLEL_TESTS += global/testing/ndim$(EXEEXT)
GLOBAL_PARALLEL_TESTS += global/testing/patch$(EXEEXT)
GLOBAL_PARALLEL_TESTS += global/testing/patch2$(EXEEXT)
GLOBAL_PARALLEL_TESTS += global/testing/patch_enumf$(EXEEXT)
GLOBAL_PARALLEL_TESTS += global/testing/perfmod$(EXEEXT)
GLOBAL_PARALLEL_TESTS += global/testing/perform$(EXEEXT)
GLOBAL_PARALLEL_TESTS += global/testing/perf$(EXEEXT)
GLOBAL_PARALLEL_TESTS += global/testing/perf2$(EXEEXT)
#GLOBAL_PARALLEL_TESTS += global/testing/pg2testmatmult$(EXEEXT) # needs 8 procs exactly
GLOBAL_PARALLEL_TESTS += global/testing/pg2test$(EXEEXT)
#GLOBAL_PARALLEL_TESTS += global/testing/pgtestmatmult$(EXEEXT) # needs 8 procs exactly
GLOBAL_PARALLEL_TESTS += global/testing/pgtest$(EXEEXT)
#GLOBAL_PARALLEL_TESTS += global/testing/random$(EXEEXT) # takes too long
GLOBAL_PARALLEL_TESTS += global/testing/scan$(EXEEXT)
GLOBAL_PARALLEL_TESTS += global/testing/simple_groups$(EXEEXT)
GLOBAL_PARALLEL_TESTS += global/testing/sparse$(EXEEXT)
GLOBAL_PARALLEL_TESTS += global/testing/sprsmatmult$(EXEEXT)
GLOBAL_PARALLEL_TESTS += global/testing/stride$(EXEEXT)
GLOBAL_PARALLEL_TESTS += global/testing/testeig$(EXEEXT)
GLOBAL_PARALLEL_TESTS += global/testing/testmatmult$(EXEEXT)
GLOBAL_PARALLEL_TESTS += global/testing/testsolve$(EXEEXT)
GLOBAL_PARALLEL_TESTS += global/testing/test$(EXEEXT)
if HAVE_SCALAPACK
GLOBAL_PARALLEL_TESTS += global/testing/testspd$(EXEEXT)
endif
else
if HAVE_LAPACK
GLOBAL_PARALLEL_TESTS += global/testing/ga_lu$(EXEEXT)
endif
endif
if MSG_COMMS_MPI
GLOBAL_PARALLEL_TESTS += global/testing/ga-mpi$(EXEEXT)
GLOBAL_PARALLEL_TESTS += global/testing/lock$(EXEEXT)
GLOBAL_PARALLEL_TESTS += global/testing/simple_groups_comm$(EXEEXT)
GLOBAL_PARALLEL_TESTS += global/testing/simple_groups_commc$(EXEEXT)
endif
if SYSV
#GLOBAL_PARALLEL_TESTS += global/testing/ipc_clean$(EXEEXT)
endif

#GLOBAL_PARALLEL_TESTS += global/testing/unit-tests/ga_abs_value$(EXEEXT)
#GLOBAL_PARALLEL_TESTS += global/testing/unit-tests/ga_acc$(EXEEXT)
#GLOBAL_PARALLEL_TESTS += global/testing/unit-tests/ga_add$(EXEEXT)
#GLOBAL_PARALLEL_TESTS += global/testing/unit-tests/ga_add_constant$(EXEEXT)
#GLOBAL_PARALLEL_TESTS += global/testing/unit-tests/ga_add_constantpatch$(EXEEXT)
#GLOBAL_PARALLEL_TESTS += global/testing/unit-tests/ga_add_diagonal$(EXEEXT)
#GLOBAL_PARALLEL_TESTS += global/testing/unit-tests/ga_add_patch$(EXEEXT)
#GLOBAL_PARALLEL_TESTS += global/testing/unit-tests/ga_copy2$(EXEEXT)
#GLOBAL_PARALLEL_TESTS += global/testing/unit-tests/ga_copy3$(EXEEXT)
#GLOBAL_PARALLEL_TESTS += global/testing/unit-tests/ga_copy$(EXEEXT)
#GLOBAL_PARALLEL_TESTS += global/testing/unit-tests/ga_copypatch2$(EXEEXT)
#GLOBAL_PARALLEL_TESTS += global/testing/unit-tests/ga_copypatch$(EXEEXT)
#GLOBAL_PARALLEL_TESTS += global/testing/unit-tests/ga_create1$(EXEEXT)
#GLOBAL_PARALLEL_TESTS += global/testing/unit-tests/ga_create2$(EXEEXT)
#GLOBAL_PARALLEL_TESTS += global/testing/unit-tests/ga_create3$(EXEEXT)
#GLOBAL_PARALLEL_TESTS += global/testing/unit-tests/ga_create$(EXEEXT)
#GLOBAL_PARALLEL_TESTS += global/testing/unit-tests/ga_create_handle$(EXEEXT)
#GLOBAL_PARALLEL_TESTS += global/testing/unit-tests/ga_create_irreg2$(EXEEXT)
#GLOBAL_PARALLEL_TESTS += global/testing/unit-tests/ga_create_irreg3$(EXEEXT)
#GLOBAL_PARALLEL_TESTS += global/testing/unit-tests/ga_create_irreg$(EXEEXT)
#GLOBAL_PARALLEL_TESTS += global/testing/unit-tests/ga_destroy$(EXEEXT)
#GLOBAL_PARALLEL_TESTS += global/testing/unit-tests/ga_dgop$(EXEEXT)
#GLOBAL_PARALLEL_TESTS += global/testing/unit-tests/ga_dot$(EXEEXT)
#GLOBAL_PARALLEL_TESTS += global/testing/unit-tests/ga_duplicate$(EXEEXT)
#GLOBAL_PARALLEL_TESTS += global/testing/unit-tests/ga_elem_divide$(EXEEXT)
#GLOBAL_PARALLEL_TESTS += global/testing/unit-tests/ga_elem_dividepatch$(EXEEXT)
#GLOBAL_PARALLEL_TESTS += global/testing/unit-tests/ga_elem_maximum$(EXEEXT)
#GLOBAL_PARALLEL_TESTS += global/testing/unit-tests/ga_elem_maximumpatch$(EXEEXT)
#GLOBAL_PARALLEL_TESTS += global/testing/unit-tests/ga_elem_minimum$(EXEEXT)
#GLOBAL_PARALLEL_TESTS += global/testing/unit-tests/ga_elem_minimumpatch$(EXEEXT)
#GLOBAL_PARALLEL_TESTS += global/testing/unit-tests/ga_elem_multiply$(EXEEXT)
#GLOBAL_PARALLEL_TESTS += global/testing/unit-tests/ga_elem_multiplypatch$(EXEEXT)
#GLOBAL_PARALLEL_TESTS += global/testing/unit-tests/ga_fill$(EXEEXT)
#GLOBAL_PARALLEL_TESTS += global/testing/unit-tests/ga_fillpatch1$(EXEEXT)
#GLOBAL_PARALLEL_TESTS += global/testing/unit-tests/ga_fillpatch$(EXEEXT)
#GLOBAL_PARALLEL_TESTS += global/testing/unit-tests/ga_gather2$(EXEEXT)
#GLOBAL_PARALLEL_TESTS += global/testing/unit-tests/ga_gather3$(EXEEXT)
#GLOBAL_PARALLEL_TESTS += global/testing/unit-tests/ga_gather$(EXEEXT)
#GLOBAL_PARALLEL_TESTS += global/testing/unit-tests/ga_get_blockinfo$(EXEEXT)
#GLOBAL_PARALLEL_TESTS += global/testing/unit-tests/ga_get$(EXEEXT)
#GLOBAL_PARALLEL_TESTS += global/testing/unit-tests/ga_get_diagonal$(EXEEXT)
#GLOBAL_PARALLEL_TESTS += global/testing/unit-tests/ga_igop2$(EXEEXT)
#GLOBAL_PARALLEL_TESTS += global/testing/unit-tests/ga_igop$(EXEEXT)
#GLOBAL_PARALLEL_TESTS += global/testing/unit-tests/ga_inquire$(EXEEXT)
#GLOBAL_PARALLEL_TESTS += global/testing/unit-tests/ga_intialize$(EXEEXT)
#GLOBAL_PARALLEL_TESTS += global/testing/unit-tests/ga_lgop$(EXEEXT)
#GLOBAL_PARALLEL_TESTS += global/testing/unit-tests/ga_median$(EXEEXT)
#GLOBAL_PARALLEL_TESTS += global/testing/unit-tests/ga_ndim2$(EXEEXT)
#GLOBAL_PARALLEL_TESTS += global/testing/unit-tests/ga_ndim$(EXEEXT)
#GLOBAL_PARALLEL_TESTS += global/testing/unit-tests/ga_nnodes$(EXEEXT)
#GLOBAL_PARALLEL_TESTS += global/testing/unit-tests/ga_nodeid$(EXEEXT)
#GLOBAL_PARALLEL_TESTS += global/testing/unit-tests/ga_pgroup_create2$(EXEEXT)
#GLOBAL_PARALLEL_TESTS += global/testing/unit-tests/ga_pgroup_create3$(EXEEXT)
#GLOBAL_PARALLEL_TESTS += global/testing/unit-tests/ga_pgroup_create4$(EXEEXT)
#GLOBAL_PARALLEL_TESTS += global/testing/unit-tests/ga_pgroup_create5$(EXEEXT)
#GLOBAL_PARALLEL_TESTS += global/testing/unit-tests/ga_pgroup_create$(EXEEXT)
#GLOBAL_PARALLEL_TESTS += global/testing/unit-tests/ga_pgroup_destroy2$(EXEEXT)
#GLOBAL_PARALLEL_TESTS += global/testing/unit-tests/ga_pgroup_destroy$(EXEEXT)
#GLOBAL_PARALLEL_TESTS += global/testing/unit-tests/ga_pgroup_nnodes_nodeid$(EXEEXT)
#GLOBAL_PARALLEL_TESTS += global/testing/unit-tests/ga_pgroup_setdefault$(EXEEXT)
#GLOBAL_PARALLEL_TESTS += global/testing/unit-tests/ga_put2$(EXEEXT)
#GLOBAL_PARALLEL_TESTS += global/testing/unit-tests/ga_put$(EXEEXT)
#GLOBAL_PARALLEL_TESTS += global/testing/unit-tests/ga_scale2$(EXEEXT)
#GLOBAL_PARALLEL_TESTS += global/testing/unit-tests/ga_scale$(EXEEXT)
#GLOBAL_PARALLEL_TESTS += global/testing/unit-tests/ga_scale_cols$(EXEEXT)
#GLOBAL_PARALLEL_TESTS += global/testing/unit-tests/ga_scale_patch$(EXEEXT)
#GLOBAL_PARALLEL_TESTS += global/testing/unit-tests/ga_scale_rows$(EXEEXT)
#GLOBAL_PARALLEL_TESTS += global/testing/unit-tests/ga_scatter$(EXEEXT)
#GLOBAL_PARALLEL_TESTS += global/testing/unit-tests/ga_set_data$(EXEEXT)
#GLOBAL_PARALLEL_TESTS += global/testing/unit-tests/ga_set_diagonal$(EXEEXT)
#GLOBAL_PARALLEL_TESTS += global/testing/unit-tests/ga_set_restricted$(EXEEXT)
#GLOBAL_PARALLEL_TESTS += global/testing/unit-tests/ga_solve$(EXEEXT)
#GLOBAL_PARALLEL_TESTS += global/testing/unit-tests/ga_sync$(EXEEXT)
#GLOBAL_PARALLEL_TESTS += global/testing/unit-tests/ga_transpose2$(EXEEXT)
#GLOBAL_PARALLEL_TESTS += global/testing/unit-tests/ga_transpose3$(EXEEXT)
#GLOBAL_PARALLEL_TESTS += global/testing/unit-tests/ga_transpose$(EXEEXT)
#GLOBAL_PARALLEL_TESTS += global/testing/unit-tests/ga_zero$(EXEEXT)
#GLOBAL_PARALLEL_TESTS += global/testing/unit-tests/ga_zerodiagonal$(EXEEXT)
#GLOBAL_PARALLEL_TESTS += global/testing/unit-tests/ga_zeropatch2$(EXEEXT)
#GLOBAL_PARALLEL_TESTS += global/testing/unit-tests/ga_zeropatch$(EXEEXT)

if HAVE_M4
GLOBAL_TESTING_M4_DEPS =
GLOBAL_TESTING_M4_DEPS += global/testing/ngatest_src/ndim_GA_FILL.src
GLOBAL_TESTING_M4_DEPS += global/testing/ngatest_src/ndim_main.src
GLOBAL_TESTING_M4_DEPS += global/testing/ngatest_src/ndim_NGA_ACC.src
GLOBAL_TESTING_M4_DEPS += global/testing/ngatest_src/ndim_NGA_ADD_PATCH.src
GLOBAL_TESTING_M4_DEPS += global/testing/ngatest_src/ndim_NGA_COPY_PATCH.src
GLOBAL_TESTING_M4_DEPS += global/testing/ngatest_src/ndim_NGA_DOT_PATCH.src
GLOBAL_TESTING_M4_DEPS += global/testing/ngatest_src/ndim_NGA_FILL_PATCH.src
GLOBAL_TESTING_M4_DEPS += global/testing/ngatest_src/ndim_NGA_GATHER.src
GLOBAL_TESTING_M4_DEPS += global/testing/ngatest_src/ndim_NGA_GET.src
GLOBAL_TESTING_M4_DEPS += global/testing/ngatest_src/ndim_NGA_PERIODIC_ACC.src
GLOBAL_TESTING_M4_DEPS += global/testing/ngatest_src/ndim_NGA_PERIODIC_GET.src
GLOBAL_TESTING_M4_DEPS += global/testing/ngatest_src/ndim_NGA_PERIODIC_PUT.src
GLOBAL_TESTING_M4_DEPS += global/testing/ngatest_src/ndim_NGA_PUT.src
GLOBAL_TESTING_M4_DEPS += global/testing/ngatest_src/ndim_NGA_SCALE_PATCH.src
GLOBAL_TESTING_M4_DEPS += global/testing/ngatest_src/ndim_NGA_SCATTER_ACC.src
GLOBAL_TESTING_M4_DEPS += global/testing/ngatest_src/ndim_NGA_SCATTER.src
GLOBAL_TESTING_M4_DEPS += global/testing/ngatest_src/ndim_util_comm.src
GLOBAL_TESTING_M4_DEPS += global/testing/ngatest_src/ndim_util.src
GLOBAL_TESTING_M4_DEPS += global/testing/ngatest_src/ngatest.def
GLOBAL_TESTING_M4_DEPS += global/testing/nga-onesided.m4
GLOBAL_TESTING_M4_DEPS += global/testing/nga-patch.m4
GLOBAL_TESTING_M4_DEPS += global/testing/nga-periodic.m4
GLOBAL_TESTING_M4_DEPS += global/testing/nga-scatter.m4
GLOBAL_TESTING_M4_DEPS += global/testing/ngatest.m4
GLOBAL_TESTING_M4_DEPS += global/testing/nga-util.m4

EXTRA_DIST += $(GLOBAL_TESTING_M4_DEPS)
EXTRA_DIST += global/testing/README

global/testing/nga-onesided.F: $(GLOBAL_TESTING_M4_DEPS)
global/testing/nga-patch.F:    $(GLOBAL_TESTING_M4_DEPS)
global/testing/nga-periodic.F: $(GLOBAL_TESTING_M4_DEPS)
global/testing/nga-scatter.F:  $(GLOBAL_TESTING_M4_DEPS)
global/testing/nga-util.F:     $(GLOBAL_TESTING_M4_DEPS)
global/testing/ngatest.F:      $(GLOBAL_TESTING_M4_DEPS)

CLEANFILES += global/testing/nga-onesided.F
CLEANFILES += global/testing/nga-patch.F
CLEANFILES += global/testing/nga-periodic.F
CLEANFILES += global/testing/nga-scatter.F
CLEANFILES += global/testing/nga-util.F
CLEANFILES += global/testing/ngatest.F

if ENABLE_F77
check_PROGRAMS += global/testing/nga-onesided
check_PROGRAMS += global/testing/nga-patch
check_PROGRAMS += global/testing/nga-periodic
check_PROGRAMS += global/testing/nga-scatter
check_PROGRAMS += global/testing/nga-util
check_PROGRAMS += global/testing/ngatest
endif

if ENABLE_F77
GLOBAL_PARALLEL_TESTS += global/testing/nga-onesided$(EXEEXT)
GLOBAL_PARALLEL_TESTS += global/testing/nga-patch$(EXEEXT)
GLOBAL_PARALLEL_TESTS += global/testing/nga-periodic$(EXEEXT)
GLOBAL_PARALLEL_TESTS += global/testing/nga-scatter$(EXEEXT)
GLOBAL_PARALLEL_TESTS += global/testing/nga-scatter$(EXEEXT)
GLOBAL_PARALLEL_TESTS += global/testing/nga-util$(EXEEXT)
GLOBAL_PARALLEL_TESTS += global/testing/ngatest$(EXEEXT)
endif

.m4.F:
	$(M4_V) \
	$(M4) -I$(top_srcdir)/global/testing $< > $@ 
endif

gtsrc =
gtsrc += global/testing/testutil.fh
gtsrc += global/testing/testutil.h
gtsrc += global/testing/util.c
if ENABLE_F77
gtsrc += global/testing/ffflush.F
gtsrc += global/testing/testblas.F
endif

global_testing_big_SOURCES                 = global/testing/big.c $(gtsrc)
global_testing_bin_SOURCES                 = global/testing/bin.F $(gtsrc)
global_testing_blktest_SOURCES             = global/testing/blktest.F $(gtsrc)
global_testing_d2test_SOURCES              = global/testing/d2test.F $(gtsrc)
global_testing_elempatch_SOURCES           = global/testing/elempatch.c $(gtsrc)
global_testing_field_test_SOURCES          = global/testing/field-test.F $(gtsrc)
global_testing_g2test_SOURCES              = global/testing/g2test.F $(gtsrc)
global_testing_g3test_SOURCES              = global/testing/g3test.F $(gtsrc)
global_testing_ga_lu_SOURCES               = global/testing/ga_lu.c $(gtsrc)
global_testing_ga_mpi_SOURCES              = global/testing/ga-mpi.c $(gtsrc)
global_testing_ga_shift_SOURCES            = global/testing/ga_shift.F $(gtsrc)
global_testing_gatscat_SOURCES             = global/testing/gatscat.c $(gtsrc)
global_testing_getmem_SOURCES              = global/testing/getmem.c $(gtsrc)
global_testing_ghosts_SOURCES              = global/testing/ghosts.F $(gtsrc)
global_testing_ipc_clean_SOURCES           = global/testing/ipc.clean.c $(gtsrc)
global_testing_jacobi_SOURCES              = global/testing/jacobi.F $(gtsrc)
global_testing_lock_SOURCES                = global/testing/lock.c $(gtsrc)
global_testing_mir_perf2_SOURCES           = global/testing/mir_perf2.F $(gtsrc)
global_testing_mmatrix_SOURCES             = global/testing/mmatrix.F $(gtsrc)
global_testing_mtest_SOURCES               = global/testing/mtest.c $(gtsrc)
global_testing_mulmatpatch_SOURCES         = global/testing/mulmatpatch.F $(gtsrc)
global_testing_mulmatpatchc_SOURCES        = global/testing/mulmatpatchc.c $(gtsrc)
global_testing_nb2test_SOURCES             = global/testing/nb2test.F $(gtsrc)
global_testing_ndim_SOURCES                = global/testing/ndim.F $(gtsrc)
global_testing_ntestc_SOURCES              = global/testing/ntestc.c $(gtsrc)
global_testing_ntestfc_SOURCES             = global/testing/ntestfc.c $(gtsrc)
global_testing_packc_SOURCES               = global/testing/packc.c $(gtsrc)
global_testing_patch_SOURCES               = global/testing/patch.F $(gtsrc)
global_testing_patch2_SOURCES              = global/testing/patch2.F $(gtsrc)
global_testing_patch_enumc_SOURCES         = global/testing/patch_enumc.c $(gtsrc)
global_testing_patch_enumf_SOURCES         = global/testing/patch_enumf.F $(gtsrc)
global_testing_perf_SOURCES                = global/testing/perf.F $(gtsrc)
global_testing_perf2_SOURCES               = global/testing/perf2.c $(gtsrc)
global_testing_perfmod_SOURCES             = global/testing/perfmod.F $(gtsrc)
global_testing_perform_SOURCES             = global/testing/perform.F $(gtsrc)
global_testing_pg2test_SOURCES             = global/testing/pg2test.F $(gtsrc)
global_testing_pg2testmatmult_SOURCES      = global/testing/pg2testmatmult.F $(gtsrc)
global_testing_pgtest_SOURCES              = global/testing/pgtest.F $(gtsrc)
global_testing_pgtestmatmult_SOURCES       = global/testing/pgtestmatmult.F $(gtsrc)
global_testing_print_SOURCES               = global/testing/print.c $(gtsrc)
global_testing_random_SOURCES              = global/testing/random.F $(gtsrc)
global_testing_scan_SOURCES                = global/testing/scan.F $(gtsrc)
global_testing_scan_addc_SOURCES           = global/testing/scan_addc.c $(gtsrc)
global_testing_scan_copyc_SOURCES          = global/testing/scan_copyc.c $(gtsrc)
global_testing_simple_groups_SOURCES       = global/testing/simple_groups.F $(gtsrc)
global_testing_simple_groups_comm_SOURCES  = global/testing/simple_groups_comm.F $(gtsrc)
global_testing_simple_groups_commc_SOURCES = global/testing/simple_groups_commc.c $(gtsrc)
global_testing_sparse_SOURCES              = global/testing/sparse.F $(gtsrc)
global_testing_sprsmatmult_SOURCES         = global/testing/sprsmatmult.F $(gtsrc)
global_testing_stride_SOURCES              = global/testing/stride.F $(gtsrc)
global_testing_test_SOURCES                = global/testing/test.F $(gtsrc)
global_testing_test_mirrored_SOURCES       = global/testing/test.F $(gtsrc)
global_testing_testc_SOURCES               = global/testing/testc.c $(gtsrc)
global_testing_testeig_SOURCES             = global/testing/testeig.F $(gtsrc)
global_testing_testmatmult_SOURCES         = global/testing/testmatmult.F $(gtsrc)
global_testing_testmatmultc_SOURCES        = global/testing/testmatmultc.c $(gtsrc)
global_testing_testmult_SOURCES            = global/testing/testmult.c $(gtsrc)
global_testing_testmultrect_SOURCES        = global/testing/testmultrect.c $(gtsrc)
global_testing_testsolve_SOURCES           = global/testing/testsolve.F $(gtsrc)
global_testing_testspd_SOURCES             = global/testing/testspd.F $(gtsrc)
global_testing_types_test_SOURCES          = global/testing/types-test.F $(gtsrc)
global_testing_unpackc_SOURCES             = global/testing/unpackc.c $(gtsrc)
nodist_global_testing_nga_onesided_SOURCES = global/testing/nga-onesided.F $(gtsrc)
nodist_global_testing_nga_patch_SOURCES    = global/testing/nga-patch.F $(gtsrc)
nodist_global_testing_nga_periodic_SOURCES = global/testing/nga-periodic.F $(gtsrc)
nodist_global_testing_nga_scatter_SOURCES  = global/testing/nga-scatter.F $(gtsrc)
nodist_global_testing_nga_util_SOURCES     = global/testing/nga-util.F $(gtsrc)
nodist_global_testing_ngatest_SOURCES      = global/testing/ngatest.F $(gtsrc)

UNIT_TEST_EXTRA_SRC =
UNIT_TEST_EXTRA_SRC += global/testing/unit-tests/mock.c
UNIT_TEST_EXTRA_SRC += global/testing/unit-tests/mock.h

global_testing_unit_tests_ga_abs_value_SOURCES            = global/testing/unit-tests/ga_abs_value.c $(UNIT_TEST_EXTRA_SRC)
global_testing_unit_tests_ga_acc_SOURCES                  = global/testing/unit-tests/ga_acc.c $(UNIT_TEST_EXTRA_SRC)
global_testing_unit_tests_ga_add_SOURCES                  = global/testing/unit-tests/ga_add.c $(UNIT_TEST_EXTRA_SRC)
global_testing_unit_tests_ga_add_constant_SOURCES         = global/testing/unit-tests/ga_add_constant.c $(UNIT_TEST_EXTRA_SRC)
global_testing_unit_tests_ga_add_constantpatch_SOURCES    = global/testing/unit-tests/ga_add_constantpatch.c $(UNIT_TEST_EXTRA_SRC)
global_testing_unit_tests_ga_add_diagonal_SOURCES         = global/testing/unit-tests/ga_add_diagonal.c $(UNIT_TEST_EXTRA_SRC)
global_testing_unit_tests_ga_add_patch_SOURCES            = global/testing/unit-tests/ga_add_patch.c $(UNIT_TEST_EXTRA_SRC)
global_testing_unit_tests_ga_copy2_SOURCES                = global/testing/unit-tests/ga_copy2.c $(UNIT_TEST_EXTRA_SRC)
global_testing_unit_tests_ga_copy3_SOURCES                = global/testing/unit-tests/ga_copy3.c $(UNIT_TEST_EXTRA_SRC)
global_testing_unit_tests_ga_copy_SOURCES                 = global/testing/unit-tests/ga_copy.c $(UNIT_TEST_EXTRA_SRC)
global_testing_unit_tests_ga_copypatch2_SOURCES           = global/testing/unit-tests/ga_copypatch2.c $(UNIT_TEST_EXTRA_SRC)
global_testing_unit_tests_ga_copypatch_SOURCES            = global/testing/unit-tests/ga_copypatch.c $(UNIT_TEST_EXTRA_SRC)
global_testing_unit_tests_ga_create1_SOURCES              = global/testing/unit-tests/ga_create1.c $(UNIT_TEST_EXTRA_SRC)
global_testing_unit_tests_ga_create2_SOURCES              = global/testing/unit-tests/ga_create2.c $(UNIT_TEST_EXTRA_SRC)
global_testing_unit_tests_ga_create3_SOURCES              = global/testing/unit-tests/ga_create3.c $(UNIT_TEST_EXTRA_SRC)
global_testing_unit_tests_ga_create_SOURCES               = global/testing/unit-tests/ga_create.c $(UNIT_TEST_EXTRA_SRC)
global_testing_unit_tests_ga_create_handle_SOURCES        = global/testing/unit-tests/ga_create_handle.c $(UNIT_TEST_EXTRA_SRC)
global_testing_unit_tests_ga_create_irreg2_SOURCES        = global/testing/unit-tests/ga_create_irreg2.c $(UNIT_TEST_EXTRA_SRC)
global_testing_unit_tests_ga_create_irreg3_SOURCES        = global/testing/unit-tests/ga_create_irreg3.c $(UNIT_TEST_EXTRA_SRC)
global_testing_unit_tests_ga_create_irreg_SOURCES         = global/testing/unit-tests/ga_create_irreg.c $(UNIT_TEST_EXTRA_SRC)
global_testing_unit_tests_ga_destroy_SOURCES              = global/testing/unit-tests/ga_destroy.c $(UNIT_TEST_EXTRA_SRC)
global_testing_unit_tests_ga_dgop_SOURCES                 = global/testing/unit-tests/ga_dgop.c $(UNIT_TEST_EXTRA_SRC)
global_testing_unit_tests_ga_dot_SOURCES                  = global/testing/unit-tests/ga_dot.c $(UNIT_TEST_EXTRA_SRC)
global_testing_unit_tests_ga_duplicate_SOURCES            = global/testing/unit-tests/ga_duplicate.c $(UNIT_TEST_EXTRA_SRC)
global_testing_unit_tests_ga_elem_divide_SOURCES          = global/testing/unit-tests/ga_elem_divide.c $(UNIT_TEST_EXTRA_SRC)
global_testing_unit_tests_ga_elem_dividepatch_SOURCES     = global/testing/unit-tests/ga_elem_dividepatch.c $(UNIT_TEST_EXTRA_SRC)
global_testing_unit_tests_ga_elem_maximum_SOURCES         = global/testing/unit-tests/ga_elem_maximum.c $(UNIT_TEST_EXTRA_SRC)
global_testing_unit_tests_ga_elem_maximumpatch_SOURCES    = global/testing/unit-tests/ga_elem_maximumpatch.c $(UNIT_TEST_EXTRA_SRC)
global_testing_unit_tests_ga_elem_minimum_SOURCES         = global/testing/unit-tests/ga_elem_minimum.c $(UNIT_TEST_EXTRA_SRC)
global_testing_unit_tests_ga_elem_minimumpatch_SOURCES    = global/testing/unit-tests/ga_elem_minimumpatch.c $(UNIT_TEST_EXTRA_SRC)
global_testing_unit_tests_ga_elem_multiply_SOURCES        = global/testing/unit-tests/ga_elem_multiply.c $(UNIT_TEST_EXTRA_SRC)
global_testing_unit_tests_ga_elem_multiplypatch_SOURCES   = global/testing/unit-tests/ga_elem_multiplypatch.c $(UNIT_TEST_EXTRA_SRC)
global_testing_unit_tests_ga_fill_SOURCES                 = global/testing/unit-tests/ga_fill.c $(UNIT_TEST_EXTRA_SRC)
global_testing_unit_tests_ga_fillpatch1_SOURCES           = global/testing/unit-tests/ga_fillpatch1.c $(UNIT_TEST_EXTRA_SRC)
global_testing_unit_tests_ga_fillpatch_SOURCES            = global/testing/unit-tests/ga_fillpatch.c $(UNIT_TEST_EXTRA_SRC)
global_testing_unit_tests_ga_gather2_SOURCES              = global/testing/unit-tests/ga_gather2.c $(UNIT_TEST_EXTRA_SRC)
global_testing_unit_tests_ga_gather3_SOURCES              = global/testing/unit-tests/ga_gather3.c $(UNIT_TEST_EXTRA_SRC)
global_testing_unit_tests_ga_gather_SOURCES               = global/testing/unit-tests/ga_gather.c $(UNIT_TEST_EXTRA_SRC)
global_testing_unit_tests_ga_get_blockinfo_SOURCES        = global/testing/unit-tests/ga_get_blockinfo.c $(UNIT_TEST_EXTRA_SRC)
global_testing_unit_tests_ga_get_SOURCES                  = global/testing/unit-tests/ga_get.c $(UNIT_TEST_EXTRA_SRC)
global_testing_unit_tests_ga_get_diagonal_SOURCES         = global/testing/unit-tests/ga_get_diagonal.c $(UNIT_TEST_EXTRA_SRC)
global_testing_unit_tests_ga_igop2_SOURCES                = global/testing/unit-tests/ga_igop2.c $(UNIT_TEST_EXTRA_SRC)
global_testing_unit_tests_ga_igop_SOURCES                 = global/testing/unit-tests/ga_igop.c $(UNIT_TEST_EXTRA_SRC)
global_testing_unit_tests_ga_inquire_SOURCES              = global/testing/unit-tests/ga_inquire.c $(UNIT_TEST_EXTRA_SRC)
global_testing_unit_tests_ga_intialize_SOURCES            = global/testing/unit-tests/ga_intialize.c $(UNIT_TEST_EXTRA_SRC)
global_testing_unit_tests_ga_lgop_SOURCES                 = global/testing/unit-tests/ga_lgop.c $(UNIT_TEST_EXTRA_SRC)
global_testing_unit_tests_ga_median_SOURCES               = global/testing/unit-tests/ga_median.c $(UNIT_TEST_EXTRA_SRC)
global_testing_unit_tests_ga_ndim2_SOURCES                = global/testing/unit-tests/ga_ndim2.c $(UNIT_TEST_EXTRA_SRC)
global_testing_unit_tests_ga_ndim_SOURCES                 = global/testing/unit-tests/ga_ndim.c $(UNIT_TEST_EXTRA_SRC)
global_testing_unit_tests_ga_nnodes_SOURCES               = global/testing/unit-tests/ga_nnodes.c $(UNIT_TEST_EXTRA_SRC)
global_testing_unit_tests_ga_nodeid_SOURCES               = global/testing/unit-tests/ga_nodeid.c $(UNIT_TEST_EXTRA_SRC)
global_testing_unit_tests_ga_pgroup_create2_SOURCES       = global/testing/unit-tests/ga_pgroup_create2.c $(UNIT_TEST_EXTRA_SRC)
global_testing_unit_tests_ga_pgroup_create3_SOURCES       = global/testing/unit-tests/ga_pgroup_create3.c $(UNIT_TEST_EXTRA_SRC)
global_testing_unit_tests_ga_pgroup_create4_SOURCES       = global/testing/unit-tests/ga_pgroup_create4.c $(UNIT_TEST_EXTRA_SRC)
global_testing_unit_tests_ga_pgroup_create5_SOURCES       = global/testing/unit-tests/ga_pgroup_create5.c $(UNIT_TEST_EXTRA_SRC)
global_testing_unit_tests_ga_pgroup_create_SOURCES        = global/testing/unit-tests/ga_pgroup_create.c $(UNIT_TEST_EXTRA_SRC)
global_testing_unit_tests_ga_pgroup_destroy2_SOURCES      = global/testing/unit-tests/ga_pgroup_destroy2.c $(UNIT_TEST_EXTRA_SRC)
global_testing_unit_tests_ga_pgroup_destroy_SOURCES       = global/testing/unit-tests/ga_pgroup_destroy.c $(UNIT_TEST_EXTRA_SRC)
global_testing_unit_tests_ga_pgroup_nnodes_nodeid_SOURCES = global/testing/unit-tests/ga_pgroup_nnodes_nodeid.c $(UNIT_TEST_EXTRA_SRC)
global_testing_unit_tests_ga_pgroup_setdefault_SOURCES    = global/testing/unit-tests/ga_pgroup_setdefault.c $(UNIT_TEST_EXTRA_SRC)
global_testing_unit_tests_ga_put2_SOURCES                 = global/testing/unit-tests/ga_put2.c $(UNIT_TEST_EXTRA_SRC)
global_testing_unit_tests_ga_put_SOURCES                  = global/testing/unit-tests/ga_put.c $(UNIT_TEST_EXTRA_SRC)
global_testing_unit_tests_ga_scale2_SOURCES               = global/testing/unit-tests/ga_scale2.c $(UNIT_TEST_EXTRA_SRC)
global_testing_unit_tests_ga_scale_SOURCES                = global/testing/unit-tests/ga_scale.c $(UNIT_TEST_EXTRA_SRC)
global_testing_unit_tests_ga_scale_cols_SOURCES           = global/testing/unit-tests/ga_scale_cols.c $(UNIT_TEST_EXTRA_SRC)
global_testing_unit_tests_ga_scale_patch_SOURCES          = global/testing/unit-tests/ga_scale_patch.c $(UNIT_TEST_EXTRA_SRC)
global_testing_unit_tests_ga_scale_rows_SOURCES           = global/testing/unit-tests/ga_scale_rows.c $(UNIT_TEST_EXTRA_SRC)
global_testing_unit_tests_ga_scatter_SOURCES              = global/testing/unit-tests/ga_scatter.c $(UNIT_TEST_EXTRA_SRC)
global_testing_unit_tests_ga_set_data_SOURCES             = global/testing/unit-tests/ga_set_data.c $(UNIT_TEST_EXTRA_SRC)
global_testing_unit_tests_ga_set_diagonal_SOURCES         = global/testing/unit-tests/ga_set_diagonal.c $(UNIT_TEST_EXTRA_SRC)
global_testing_unit_tests_ga_set_restricted_SOURCES       = global/testing/unit-tests/ga_set_restricted.c $(UNIT_TEST_EXTRA_SRC)
global_testing_unit_tests_ga_solve_SOURCES                = global/testing/unit-tests/ga_solve.c $(UNIT_TEST_EXTRA_SRC)
global_testing_unit_tests_ga_sync_SOURCES                 = global/testing/unit-tests/ga_sync.c $(UNIT_TEST_EXTRA_SRC)
global_testing_unit_tests_ga_transpose2_SOURCES           = global/testing/unit-tests/ga_transpose2.c $(UNIT_TEST_EXTRA_SRC)
global_testing_unit_tests_ga_transpose3_SOURCES           = global/testing/unit-tests/ga_transpose3.c $(UNIT_TEST_EXTRA_SRC)
global_testing_unit_tests_ga_transpose_SOURCES            = global/testing/unit-tests/ga_transpose.c $(UNIT_TEST_EXTRA_SRC)
global_testing_unit_tests_ga_zero_SOURCES                 = global/testing/unit-tests/ga_zero.c $(UNIT_TEST_EXTRA_SRC)
global_testing_unit_tests_ga_zerodiagonal_SOURCES         = global/testing/unit-tests/ga_zerodiagonal.c $(UNIT_TEST_EXTRA_SRC)
global_testing_unit_tests_ga_zeropatch2_SOURCES           = global/testing/unit-tests/ga_zeropatch2.c $(UNIT_TEST_EXTRA_SRC)
global_testing_unit_tests_ga_zeropatch_SOURCES            = global/testing/unit-tests/ga_zeropatch.c $(UNIT_TEST_EXTRA_SRC)

global_testing_bin_LDFLAGS                = $(FLDFLAGS)
global_testing_blktest_LDFLAGS            = $(FLDFLAGS)
global_testing_d2test_LDFLAGS             = $(FLDFLAGS)
global_testing_field_test_LDFLAGS         = $(FLDFLAGS)
global_testing_g2test_LDFLAGS             = $(FLDFLAGS)
global_testing_g3test_LDFLAGS             = $(FLDFLAGS)
global_testing_ga_shift_LDFLAGS           = $(FLDFLAGS)
global_testing_ghosts_LDFLAGS             = $(FLDFLAGS)
global_testing_jacobi_LDFLAGS             = $(FLDFLAGS)
global_testing_mir_perf2_LDFLAGS          = $(FLDFLAGS)
global_testing_mmatrix_LDFLAGS            = $(FLDFLAGS)
global_testing_mulmatpatch_LDFLAGS        = $(FLDFLAGS)
global_testing_nb2test_LDFLAGS            = $(FLDFLAGS)
global_testing_ndim_LDFLAGS               = $(FLDFLAGS)
global_testing_nga_onesided_LDFLAGS       = $(FLDFLAGS)
global_testing_nga_patch_LDFLAGS          = $(FLDFLAGS)
global_testing_nga_periodic_LDFLAGS       = $(FLDFLAGS)
global_testing_nga_scatter_LDFLAGS        = $(FLDFLAGS)
global_testing_ngatest_LDFLAGS            = $(FLDFLAGS)
global_testing_nga_util_LDFLAGS           = $(FLDFLAGS)
global_testing_patch_enumf_LDFLAGS        = $(FLDFLAGS)
global_testing_patch_LDFLAGS              = $(FLDFLAGS)
global_testing_patch2_LDFLAGS             = $(FLDFLAGS)
global_testing_perf_LDFLAGS               = $(FLDFLAGS)
global_testing_perfmod_LDFLAGS            = $(FLDFLAGS)
global_testing_perform_LDFLAGS            = $(FLDFLAGS)
global_testing_pg2test_LDFLAGS            = $(FLDFLAGS)
global_testing_pg2testmatmult_LDFLAGS     = $(FLDFLAGS)
global_testing_pgtest_LDFLAGS             = $(FLDFLAGS)
global_testing_pgtestmatmult_LDFLAGS      = $(FLDFLAGS)
global_testing_random_LDFLAGS             = $(FLDFLAGS)
global_testing_scan_LDFLAGS               = $(FLDFLAGS)
global_testing_simple_groups_LDFLAGS      = $(FLDFLAGS)
global_testing_simple_groups_comm_LDFLAGS = $(FLDFLAGS)
global_testing_sparse_LDFLAGS             = $(FLDFLAGS)
global_testing_sprsmatmult_LDFLAGS        = $(FLDFLAGS)
global_testing_stride_LDFLAGS             = $(FLDFLAGS)
global_testing_testeig_LDFLAGS            = $(FLDFLAGS)
global_testing_test_LDFLAGS               = $(FLDFLAGS)
global_testing_test_mirrored_LDFLAGS      = $(FLDFLAGS)
global_testing_testmatmult_LDFLAGS        = $(FLDFLAGS)
global_testing_testsolve_LDFLAGS          = $(FLDFLAGS)
global_testing_testspd_LDFLAGS            = $(FLDFLAGS)
global_testing_types_test_LDFLAGS         = $(FLDFLAGS)

nodist_EXTRA_global_testing_big_SOURCES                 = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_elempatch_SOURCES           = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_ga_lu_SOURCES               = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_ga_mpi_SOURCES              = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_gatscat_SOURCES             = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_getmem_SOURCES              = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_ipc_clean_SOURCES           = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_lock_SOURCES                = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_mtest_SOURCES               = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_mulmatpatchc_SOURCES        = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_ntestc_SOURCES              = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_ntestfc_SOURCES             = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_packc_SOURCES               = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_patch_enumc_SOURCES         = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_perf2_SOURCES               = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_print_SOURCES               = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_scan_addc_SOURCES           = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_scan_copyc_SOURCES          = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_simple_groups_commc_SOURCES = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_testc_SOURCES               = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_testmatmultc_SOURCES        = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_testmult_SOURCES            = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_testmultrect_SOURCES        = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_unpackc_SOURCES             = $(MAYBE_DUMMY_F_LINK)

global_testing_test_mirrored_CPPFLAGS = $(AM_CPPFLAGS) -DMIRROR=1

nodist_EXTRA_global_testing_unit_tests_ga_abs_value_SOURCES            = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_unit_tests_ga_acc_SOURCES                  = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_unit_tests_ga_add_SOURCES                  = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_unit_tests_ga_add_constant_SOURCES         = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_unit_tests_ga_add_constantpatch_SOURCES    = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_unit_tests_ga_add_diagonal_SOURCES         = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_unit_tests_ga_add_patch_SOURCES            = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_unit_tests_ga_copy2_SOURCES                = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_unit_tests_ga_copy3_SOURCES                = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_unit_tests_ga_copy_SOURCES                 = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_unit_tests_ga_copypatch2_SOURCES           = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_unit_tests_ga_copypatch_SOURCES            = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_unit_tests_ga_create1_SOURCES              = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_unit_tests_ga_create2_SOURCES              = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_unit_tests_ga_create3_SOURCES              = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_unit_tests_ga_create_SOURCES               = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_unit_tests_ga_create_handle_SOURCES        = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_unit_tests_ga_create_irreg2_SOURCES        = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_unit_tests_ga_create_irreg3_SOURCES        = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_unit_tests_ga_create_irreg_SOURCES         = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_unit_tests_ga_destroy_SOURCES              = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_unit_tests_ga_dgop_SOURCES                 = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_unit_tests_ga_dot_SOURCES                  = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_unit_tests_ga_duplicate_SOURCES            = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_unit_tests_ga_elem_divide_SOURCES          = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_unit_tests_ga_elem_dividepatch_SOURCES     = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_unit_tests_ga_elem_maximum_SOURCES         = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_unit_tests_ga_elem_maximumpatch_SOURCES    = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_unit_tests_ga_elem_minimum_SOURCES         = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_unit_tests_ga_elem_minimumpatch_SOURCES    = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_unit_tests_ga_elem_multiply_SOURCES        = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_unit_tests_ga_elem_multiplypatch_SOURCES   = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_unit_tests_ga_fill_SOURCES                 = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_unit_tests_ga_fillpatch1_SOURCES           = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_unit_tests_ga_fillpatch_SOURCES            = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_unit_tests_ga_gather2_SOURCES              = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_unit_tests_ga_gather3_SOURCES              = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_unit_tests_ga_gather_SOURCES               = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_unit_tests_ga_get_blockinfo_SOURCES        = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_unit_tests_ga_get_SOURCES                  = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_unit_tests_ga_get_diagonal_SOURCES         = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_unit_tests_ga_igop2_SOURCES                = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_unit_tests_ga_igop_SOURCES                 = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_unit_tests_ga_inquire_SOURCES              = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_unit_tests_ga_intialize_SOURCES            = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_unit_tests_ga_lgop_SOURCES                 = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_unit_tests_ga_median_SOURCES               = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_unit_tests_ga_ndim2_SOURCES                = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_unit_tests_ga_ndim_SOURCES                 = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_unit_tests_ga_nnodes_SOURCES               = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_unit_tests_ga_nodeid_SOURCES               = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_unit_tests_ga_pgroup_create2_SOURCES       = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_unit_tests_ga_pgroup_create3_SOURCES       = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_unit_tests_ga_pgroup_create4_SOURCES       = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_unit_tests_ga_pgroup_create5_SOURCES       = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_unit_tests_ga_pgroup_create_SOURCES        = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_unit_tests_ga_pgroup_destroy2_SOURCES      = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_unit_tests_ga_pgroup_destroy_SOURCES       = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_unit_tests_ga_pgroup_nnodes_nodeid_SOURCES = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_unit_tests_ga_pgroup_setdefault_SOURCES    = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_unit_tests_ga_put2_SOURCES                 = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_unit_tests_ga_put_SOURCES                  = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_unit_tests_ga_scale2_SOURCES               = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_unit_tests_ga_scale_SOURCES                = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_unit_tests_ga_scale_cols_SOURCES           = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_unit_tests_ga_scale_patch_SOURCES          = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_unit_tests_ga_scale_rows_SOURCES           = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_unit_tests_ga_scatter_SOURCES              = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_unit_tests_ga_set_data_SOURCES             = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_unit_tests_ga_set_diagonal_SOURCES         = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_unit_tests_ga_set_restricted_SOURCES       = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_unit_tests_ga_solve_SOURCES                = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_unit_tests_ga_sync_SOURCES                 = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_unit_tests_ga_transpose2_SOURCES           = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_unit_tests_ga_transpose3_SOURCES           = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_unit_tests_ga_transpose_SOURCES            = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_unit_tests_ga_zero_SOURCES                 = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_unit_tests_ga_zerodiagonal_SOURCES         = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_unit_tests_ga_zeropatch2_SOURCES           = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_global_testing_unit_tests_ga_zeropatch_SOURCES            = $(MAYBE_DUMMY_F_LINK)

##############################################################################
# global/examples
#
check_PROGRAMS += global/examples/lennard-jones/lennard
if ENABLE_F77
check_PROGRAMS += global/examples/boltzmann/boltz
check_PROGRAMS += global/examples/conjugate_gradient/ga_cg
check_PROGRAMS += global/examples/md_cluster/grp_sim
check_PROGRAMS += global/examples/scf/scf
endif

GLOBAL_PARALLEL_TESTS += global/examples/lennard-jones/lennard$(EXEEXT)
if ENABLE_F77
GLOBAL_PARALLEL_TESTS += global/examples/boltzmann/boltz$(EXEEXT)
#GLOBAL_PARALLEL_TESTS += global/examples/conjugate_gradient/ga_cg$(EXEEXT) # needs input
#GLOBAL_PARALLEL_TESTS += global/examples/md_cluster/grp_sim$(EXEEXT) # needs input
#GLOBAL_PARALLEL_TESTS += global/examples/scf/scf$(EXEEXT) # needs input
endif

global_examples_boltzmann_boltz_SOURCES = \
global/examples/boltzmann/boltzmann.F     \
global/examples/boltzmann/common          \
global/examples/boltzmann/equil.F         \
global/examples/boltzmann/get_patch.F     \
global/examples/boltzmann/initpar.F       \
global/examples/boltzmann/main.F          \
global/examples/boltzmann/printdat.F      \
global/examples/boltzmann/properties.F    \
global/examples/boltzmann/setup.F         \
global/examples/boltzmann/timestep.F      \
global/examples/boltzmann/vorticity.F
global_examples_boltzmann_boltz_LDFLAGS = $(FLDFLAGS)
global_examples_boltzmann_boltz_CPPFLAGS = $(AM_CPPFLAGS) -I$(top_srcdir)/global/examples/boltzmann

global_examples_conjugate_gradient_ga_cg_SOURCES = \
global/examples/conjugate_gradient/computeloops.F  \
global/examples/conjugate_gradient/finclude.h      \
global/examples/conjugate_gradient/ga_cg.c         \
global/examples/conjugate_gradient/matvecmul.c     \
global/examples/conjugate_gradient/read_input.c

global_examples_lennard_jones_lennard_SOURCES = \
global/examples/lennard-jones/lennard.c
nodist_EXTRA_global_examples_lennard_jones_lennard_SOURCES = $(MAYBE_DUMMY_F_LINK)

global_examples_md_cluster_grp_sim_SOURCES = \
global/examples/md_cluster/cell_list.F       \
global/examples/md_cluster/cl_sim.F          \
global/examples/md_cluster/cluster.F         \
global/examples/md_cluster/common.fh         \
global/examples/md_cluster/estats.F          \
global/examples/md_cluster/factor.F          \
global/examples/md_cluster/force.F           \
global/examples/md_cluster/gather.F          \
global/examples/md_cluster/grp_sim.F         \
global/examples/md_cluster/hash.F            \
global/examples/md_cluster/heapsort.F        \
global/examples/md_cluster/kin.F             \
global/examples/md_cluster/mdinit.F          \
global/examples/md_cluster/mdout.F           \
global/examples/md_cluster/mdstep.F          \
global/examples/md_cluster/newcfg.F          \
global/examples/md_cluster/nextc.F           \
global/examples/md_cluster/pairs.F           \
global/examples/md_cluster/ran3.F            \
global/examples/md_cluster/rdpar.F           \
global/examples/md_cluster/scatter.F         \
global/examples/md_cluster/shuffle.F         \
global/examples/md_cluster/sort.F            \
global/examples/md_cluster/tstats.F          \
global/examples/md_cluster/update.F
global_examples_md_cluster_grp_sim_LDFLAGS = $(FLDFLAGS)

global_examples_scf_scf_SOURCES = \
global/examples/scf/cscf.h        \
global/examples/scf/input.F       \
global/examples/scf/integ.F       \
global/examples/scf/output.F      \
global/examples/scf/scf.F         \
global/examples/scf/timer.F
global_examples_scf_scf_LDFLAGS = $(FLDFLAGS)

EXTRA_DIST += global/examples/boltzmann/README
EXTRA_DIST += global/examples/lennard-jones/README
EXTRA_DIST += global/examples/md_cluster/atom.inp
EXTRA_DIST += global/examples/scf/be16.inpt
EXTRA_DIST += global/examples/scf/be.inpt
EXTRA_DIST += global/examples/scf/ft-scf.F
EXTRA_DIST += global/examples/scf/README

##############################################################################
# global/X
#
if ENABLE_XREGION
if NO_X
else
bin_PROGRAMS += global/X/xregion

global_X_xregion_SOURCES =    \
global/X/xregion_colormap.c   \
global/X/xregion_dialog.c     \
global/X/xregion_display.c    \
global/X/xregion_fileio.c     \
global/X/xregion_main.c       \
global/X/xregion_overview.c   \
global/X/xregion_pixregion.c  \
global/X/xregion_scrollbars.c \
global/X/xregion_util.c       \
global/X/xregion_view.c       \
global/X/xregion_xcmap.c

global_X_xregion_CPPFLAGS = $(X_CPPFLAGS)
global_X_xregion_LDFLAGS = $(X_LDFLAGS)
global_X_xregion_LDADD = -lXaw -lXmu -lXt -lX11 -lXau -lXdmcp -lXext

endif
endif

EXTRA_DIST += global/X/README

##############################################################################
# global/trace
#
bin_PROGRAMS += global/trace/adjust
bin_PROGRAMS += global/trace/collisions

global_trace_adjust_SOURCES = global/trace/adjust.c
global_trace_adjust_LDFLAGS =
global_trace_adjust_LDADD =

global_trace_collisions_SOURCES = global/trace/collisions.c
global_trace_collisions_LDFLAGS =
global_trace_collisions_LDADD =

EXTRA_DIST += global/trace/README
EXTRA_DIST += global/trace/test.f

##############################################################################
# global/doc
#
EXTRA_DIST += global/doc/README
EXTRA_DIST += global/doc/Supercomputing94.pdf
EXTRA_DIST += global/doc/user.pdf

##############################################################################
# pario
#
libga_la_SOURCES += pario/dra/buffers.c
libga_la_SOURCES += pario/dra/buffers.h
libga_la_SOURCES += pario/dra/capi.c
libga_la_SOURCES += pario/dra/disk.arrays.c
libga_la_SOURCES += pario/dra/disk.param.c
libga_la_SOURCES += pario/dra/draf2c.h
libga_la_SOURCES += pario/dra/drap.h
libga_la_SOURCES += pario/dra/env.c
libga_la_SOURCES += pario/dra/fortran.c
libga_la_SOURCES += pario/dra/patch.util.c
libga_la_SOURCES += pario/eaf/eaf.c
libga_la_SOURCES += pario/eaf/eaf_f2c.c
libga_la_SOURCES += pario/eaf/eafP.h
libga_la_SOURCES += pario/elio/elio.c
libga_la_SOURCES += pario/elio/eliop.h
libga_la_SOURCES += pario/elio/pablo.h
libga_la_SOURCES += pario/elio/stat.c
libga_la_SOURCES += pario/sf/sf_capi.c
libga_la_SOURCES += pario/sf/sff2c.h
libga_la_SOURCES += pario/sf/sf_fortran.c
libga_la_SOURCES += pario/sf/shared.files.c

include_HEADERS += pario/dra/dra.fh
include_HEADERS += pario/dra/dra.h
include_HEADERS += pario/eaf/eaf.fh
include_HEADERS += pario/eaf/eaf.h
include_HEADERS += pario/elio/chemio.h
include_HEADERS += pario/elio/elio.h
include_HEADERS += pario/sf/coms.h
include_HEADERS += pario/sf/sf.fh
include_HEADERS += pario/sf/sf.h

AM_CPPFLAGS += -I$(top_srcdir)/pario/dra
AM_CPPFLAGS += -I$(top_srcdir)/pario/eaf
AM_CPPFLAGS += -I$(top_srcdir)/pario/elio
AM_CPPFLAGS += -I$(top_srcdir)/pario/sf

EXTRA_DIST += pario/dra/dbg_read.c
EXTRA_DIST += pario/dra/dbg_write.c
EXTRA_DIST += pario/dra/dra2arviz.c
EXTRA_DIST += pario/dra/global.unsup.c

if ENABLE_F77
check_PROGRAMS += pario/dra/big
check_PROGRAMS += pario/dra/dra_mxm
check_PROGRAMS += pario/dra/ntest
check_PROGRAMS += pario/dra/perf
check_PROGRAMS += pario/dra/perf2
check_PROGRAMS += pario/dra/perf3
check_PROGRAMS += pario/dra/test
check_PROGRAMS += pario/dra/test_mxm
check_PROGRAMS += pario/dra/time_mxm
check_PROGRAMS += pario/eaf/test
check_PROGRAMS += pario/sf/test
endif

check_PROGRAMS += pario/dra/bign
check_PROGRAMS += pario/dra/dbg_read
check_PROGRAMS += pario/dra/dbg_write
check_PROGRAMS += pario/dra/dra2arviz
check_PROGRAMS += pario/dra/ntestc
check_PROGRAMS += pario/dra/perfn
check_PROGRAMS += pario/dra/rate

PARIO_SERIAL_TESTS =
PARIO_SERIAL_TESTS_XFAIL =
PARIO_PARALLEL_TESTS =
PARIO_PARALLEL_TESTS_XFAIL =
PARIO_TESTS = $(PARIO_SERIAL_TESTS) $(PARIO_PARALLEL_TESTS)
PARIO_TESTS_XFAIL = $(PARIO_SERIAL_TESTS_XFAIL) $(PARIO_PARALLEL_TESTS_XFAIL)

PARIO_PARALLEL_TESTS += pario/dra/bign$(EXEEXT)
#PARIO_PARALLEL_TESTS += pario/dra/dbg_read$(EXEEXT) # barely compiles, wrong test
#PARIO_PARALLEL_TESTS += pario/dra/dbg_write$(EXEEXT) # barely compiles, wrong test
#PARIO_PARALLEL_TESTS += pario/dra/dra2arviz$(EXEEXT) # not a test?
PARIO_PARALLEL_TESTS += pario/dra/ntestc$(EXEEXT)
PARIO_PARALLEL_TESTS += pario/dra/perfn$(EXEEXT)
PARIO_PARALLEL_TESTS += pario/dra/rate$(EXEEXT)

dtsrc =
dtsrc += pario/dra/ffflush.F
dtsrc += pario/dra/util.c

pario_dra_big_SOURCES       = pario/dra/big.F $(dtsrc)
pario_dra_bign_SOURCES      = pario/dra/bign.c
pario_dra_dbg_read_SOURCES  = pario/dra/dbg_read.c
pario_dra_dbg_write_SOURCES = pario/dra/dbg_write.c
pario_dra_dra2arviz_SOURCES = pario/dra/dra2arviz.c
pario_dra_dra_mxm_SOURCES   = pario/dra/dra_mxm.F $(dtsrc)
pario_dra_ntest_SOURCES     = pario/dra/ntest.F $(dtsrc)
pario_dra_ntestc_SOURCES    = pario/dra/ntestc.c
pario_dra_perf_SOURCES      = pario/dra/perf.F $(dtsrc)
pario_dra_perf2_SOURCES     = pario/dra/perf2.F $(dtsrc)
pario_dra_perf3_SOURCES     = pario/dra/perf3.F $(dtsrc)
pario_dra_perfn_SOURCES     = pario/dra/perfn.c
pario_dra_rate_SOURCES      = pario/dra/rate.c
pario_dra_test_SOURCES      = pario/dra/test.F $(dtsrc)
pario_dra_test_mxm_SOURCES  = pario/dra/test_mxm.F $(dtsrc)
pario_dra_time_mxm_SOURCES  = pario/dra/time_mxm.F $(dtsrc)
pario_eaf_test_SOURCES      = pario/eaf/test.F $(dtsrc)
pario_sf_test_SOURCES       = pario/sf/test.F $(dtsrc)

nodist_EXTRA_pario_dra_bign_SOURCES      = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_pario_dra_dbg_read_SOURCES  = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_pario_dra_dbg_write_SOURCES = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_pario_dra_dra2arviz_SOURCES = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_pario_dra_ntestc_SOURCES    = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_pario_dra_perfn_SOURCES     = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_pario_dra_rate_SOURCES      = $(MAYBE_DUMMY_F_LINK)

pario_dra_big_LDFLAGS       = $(FLDFLAGS)
pario_dra_dra_mxm_LDFLAGS   = $(FLDFLAGS)
pario_dra_ntest_LDFLAGS     = $(FLDFLAGS)
pario_dra_perf_LDFLAGS      = $(FLDFLAGS)
pario_dra_perf2_LDFLAGS     = $(FLDFLAGS)
pario_dra_perf3_LDFLAGS     = $(FLDFLAGS)
pario_dra_test_LDFLAGS      = $(FLDFLAGS)
pario_dra_test_mxm_LDFLAGS  = $(FLDFLAGS)
pario_dra_time_mxm_LDFLAGS  = $(FLDFLAGS)
pario_eaf_test_LDFLAGS      = $(FLDFLAGS)
pario_sf_test_LDFLAGS       = $(FLDFLAGS)

EXTRA_DIST += pario/dra/README

##############################################################################
# ga++/src
#
if CXX_BINDINGS

lib_LTLIBRARIES += libga++.la

libga___la_SOURCES =
libga___la_SOURCES += ga++/src/GAServices.cc
libga___la_SOURCES += ga++/src/GlobalArray.cc
libga___la_SOURCES += ga++/src/PGroup.cc
libga___la_SOURCES += ga++/src/init_term.cc
libga___la_SOURCES += ga++/src/overload.cc
libga___la_SOURCES += ga++/src/services.cc

include_HEADERS += ga++/src/ga++.h
include_HEADERS += ga++/src/GAServices.h
include_HEADERS += ga++/src/GlobalArray.h
include_HEADERS += ga++/src/init_term.h
include_HEADERS += ga++/src/PGroup.h
include_HEADERS += ga++/src/services.h

AM_CPPFLAGS += -I$(top_srcdir)/ga++/src

endif # CXX_BINDINGS

EXTRA_DIST += ga++/README

##############################################################################
# ga++/testing
#
if CXX_BINDINGS

check_PROGRAMS += ga++/testing/elempatch
check_PROGRAMS += ga++/testing/mtest
check_PROGRAMS += ga++/testing/ntestc
check_PROGRAMS += ga++/testing/testc
check_PROGRAMS += ga++/testing/testmult

CXX_SERIAL_TESTS =
CXX_SERIAL_TESTS_XFAIL =
CXX_PARALLEL_TESTS =
CXX_PARALLEL_TESTS_XFAIL =
CXX_TESTS = $(CXX_SERIAL_TESTS) $(CXX_PARALLEL_TESTS)
CXX_TESTS_XFAIL = $(CXX_SERIAL_TESTS_XFAIL) $(CXX_PARALLEL_TESTS_XFAIL)

CXX_PARALLEL_TESTS += ga++/testing/elempatch$(EXEEXT)
CXX_PARALLEL_TESTS += ga++/testing/mtest$(EXEEXT)
CXX_PARALLEL_TESTS += ga++/testing/ntestc$(EXEEXT)
CXX_PARALLEL_TESTS += ga++/testing/testc$(EXEEXT)
CXX_PARALLEL_TESTS += ga++/testing/testmult$(EXEEXT)

ga___testing_elempatch_SOURCES = ga++/testing/elempatch.cc
ga___testing_mtest_SOURCES     = ga++/testing/mtest.cc
ga___testing_ntestc_SOURCES    = ga++/testing/ntestc.cc \
                                 ga++/testing/util.c
ga___testing_testc_SOURCES     = ga++/testing/testc.cc
ga___testing_testmult_SOURCES  = ga++/testing/testmult.cc

ga___testing_elempatch_LDADD   = libga++.la $(LDADD) $(MAYBE_FLIBS)
ga___testing_mtest_LDADD       = libga++.la $(LDADD) $(MAYBE_FLIBS)
ga___testing_ntestc_LDADD      = libga++.la $(LDADD) $(MAYBE_FLIBS)
ga___testing_testc_LDADD       = libga++.la $(LDADD) $(MAYBE_FLIBS)
ga___testing_testmult_LDADD    = libga++.la $(LDADD) $(MAYBE_FLIBS)

ga___testing_elempatch_LDFLAGS = $(FLDFLAGS)
ga___testing_mtest_LDFLAGS     = $(FLDFLAGS)
ga___testing_ntestc_LDFLAGS    = $(FLDFLAGS)
ga___testing_testc_LDFLAGS     = $(FLDFLAGS)
ga___testing_testmult_LDFLAGS  = $(FLDFLAGS)

endif # CXX_BINDINGS

##############################################################################
# python
#
if HAVE_PYTHON

python: install
	(cd $(top_srcdir)/python && $(PYTHON) setup.py build_ext -L$(libdir) -I$(includedir)) || exit 1
	(cd $(top_srcdir)/python && $(PYTHON) setup.py install) || exit 1

endif

##############################################################################
# gparrays
#
if ENABLE_GPARRAYS
GPPAPI_H      = $(top_srcdir)/gparrays/src/gp-papi.h

gparrays/src/gp-wapi.h: $(GPPAPI_H)
	@-rm -f $@
	$(SED_V) \
	$(SED) -e "s/GPPAPI/GPWAPI/g" \
	       -e "s/pgp_/wgp_/g" \
	       -e "/typedef/d" \
		   $(GPPAPI_H) > $@

gparrays/src/gp-wapidefs.h: $(GPPAPI_H)
	@-rm -f $@
	$(SED_V) \
	$(SED) -e "s/GPPAPI/GPWAPIDEFS/g" \
		   -e "s/^.*p\(gp_.*\)(.*$$/#define w\1 p\1/g" \
		   -e "/^ /d" \
		   -e "/^\/\*/d" \
		   -e "/^ \*/d" \
	       -e "/typedef/d" \
	       -e "/include/d" \
		   $(GPPAPI_H) > $@

BUILT_SOURCES += gparrays/src/gp-wapi.h
BUILT_SOURCES += gparrays/src/gp-wapidefs.h

CLEANFILES += gparrays/src/gp-wapi.h
CLEANFILES += gparrays/src/gp-wapidefs.h

libga_la_SOURCES += gparrays/src/gpbase.c
libga_la_SOURCES += gparrays/src/gpbase.h
libga_la_SOURCES += gparrays/src/gpcapi.c
libga_la_SOURCES += gparrays/src/gpfapi.c
libga_la_SOURCES += gparrays/src/gponesided.c
libga_la_SOURCES += gparrays/src/gp-papi.h
nodist_libga_la_SOURCES += gparrays/src/gp-wapidefs.h

include_HEADERS += gparrays/src/gp.h
include_HEADERS += gparrays/src/gp-papi.h
nodist_include_HEADERS += gparrays/src/gp-wapi.h

if ENABLE_F77
check_PROGRAMS += gparrays/testing/test
endif

check_PROGRAMS += gparrays/testing/testc
check_PROGRAMS += gparrays/testing/transpose

gparrays_testing_test_SOURCES      = gparrays/testing/test.F
gparrays_testing_testc_SOURCES     = gparrays/testing/testc.c
gparrays_testing_transpose_SOURCES = gparrays/testing/transpose.c

nodist_EXTRA_gparrays_testing_testc_SOURCES     = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_gparrays_testing_transpose_SOURCES = $(MAYBE_DUMMY_F_LINK)

gparrays_testing_test_LDFLAGS  = $(FLDFLAGS)

AM_CPPFLAGS += -I$(top_build_prefix)gparrays/src
AM_CPPFLAGS += -I$(top_srcdir)/gparrays/src

endif # ENABLE_GPARRAYS

##############################################################################
# gaf2c
#
libga_la_SOURCES += gaf2c/gaf2c.c
libga_la_SOURCES += gaf2c/drand.c
if ENABLE_F77
libga_la_SOURCES += gaf2c/farg.F
endif

include_HEADERS += gaf2c/typesf2c.h
include_HEADERS += gaf2c/farg.h

AM_CPPFLAGS += -I$(top_build_prefix)gaf2c
AM_CPPFLAGS += -I$(top_srcdir)/gaf2c

if ENABLE_F77
check_PROGRAMS += gaf2c/testarg
endif

GAF2C_SERIAL_TESTS =
GAF2C_SERIAL_TESTS_XFAIL =
GAF2C_PARALLEL_TESTS =
GAF2C_PARALLEL_TESTS_XFAIL =
GAF2C_TESTS = $(GAF2C_SERIAL_TESTS) $(GAF2C_PARALLEL_TESTS)
GAF2C_TESTS_XFAIL = $(GAF2C_SERIAL_TESTS_XFAIL) $(GAF2C_PARALLEL_TESTS_XFAIL)

if ENABLE_F77
GAF2C_SERIAL_TESTS += gaf2c/testarg$(EXEEXT)
endif

gaf2c_testarg_SOURCES = gaf2c/testarg.c gaf2c/testargf.f
gaf2c_testarg_LDFLAGS = $(FLDFLAGS)

##############################################################################
# tcgmsg
#
libga_la_SOURCES += tcgmsg/capi.c
libga_la_SOURCES += tcgmsg/fapi.c

include_HEADERS += tcgmsg/tcgmsg.fh
include_HEADERS += tcgmsg/tcgmsg.h
include_HEADERS += tcgmsg/tcgmsg-mpi/msgtypesc.h
include_HEADERS += tcgmsg/tcgmsg-mpi/msgtypesf.h

AM_CPPFLAGS += -I$(top_srcdir)/tcgmsg

##############################################################################
# tcgmsg/examples
#
# NOTE: While we don't build the tcgmsg examples, we include everything in the
# distribution anyway.
#
EXTRA_DIST += tcgmsg/examples/blkdat120lin.f
EXTRA_DIST += tcgmsg/examples/blkdat15.f
EXTRA_DIST += tcgmsg/examples/blkdat240lin.f
EXTRA_DIST += tcgmsg/examples/blkdat30.f
EXTRA_DIST += tcgmsg/examples/blkdat60.f
EXTRA_DIST += tcgmsg/examples/blkdat60lin.f
EXTRA_DIST += tcgmsg/examples/cscf120lin.h
EXTRA_DIST += tcgmsg/examples/cscf15.h
EXTRA_DIST += tcgmsg/examples/cscf240lin.h
EXTRA_DIST += tcgmsg/examples/cscf30.h
EXTRA_DIST += tcgmsg/examples/cscf60.h
EXTRA_DIST += tcgmsg/examples/cscf60lin.h
EXTRA_DIST += tcgmsg/examples/cscf.h
EXTRA_DIST += tcgmsg/examples/daxpy1.s
EXTRA_DIST += tcgmsg/examples/daxpy.f
EXTRA_DIST += tcgmsg/examples/ddot.f
EXTRA_DIST += tcgmsg/examples/demo.proto
EXTRA_DIST += tcgmsg/examples/diagon.f
EXTRA_DIST += tcgmsg/examples/dscal.f
EXTRA_DIST += tcgmsg/examples/fexit.f.proto
EXTRA_DIST += tcgmsg/examples/getmem.c
EXTRA_DIST += tcgmsg/examples/grid.c
EXTRA_DIST += tcgmsg/examples/ieeetrap.c
EXTRA_DIST += tcgmsg/examples/integ.f
EXTRA_DIST += tcgmsg/examples/jacobi.f
EXTRA_DIST += tcgmsg/examples/Makefile.proto
EXTRA_DIST += tcgmsg/examples/mc.f
EXTRA_DIST += tcgmsg/examples/mc.input
EXTRA_DIST += tcgmsg/examples/md.f
EXTRA_DIST += tcgmsg/examples/mxv_daxpy1.f
EXTRA_DIST += tcgmsg/examples/mxv_dgemv.f
EXTRA_DIST += tcgmsg/examples/mxv_fortran.f
EXTRA_DIST += tcgmsg/examples/output.f
EXTRA_DIST += tcgmsg/examples/prtri.f
EXTRA_DIST += tcgmsg/examples/random.c
EXTRA_DIST += tcgmsg/examples/README
EXTRA_DIST += tcgmsg/examples/runit
EXTRA_DIST += tcgmsg/examples/runit.grid
EXTRA_DIST += tcgmsg/examples/scfblas.f
EXTRA_DIST += tcgmsg/examples/scf.f
EXTRA_DIST += tcgmsg/examples/timer.f
EXTRA_DIST += tcgmsg/examples/trace.out
EXTRA_DIST += tcgmsg/examples/xpix.shar

##############################################################################
# tcgmsg/tcgmsg-mpi
#
if MSG_COMMS_MPI
libga_la_SOURCES += tcgmsg/tcgmsg-mpi/checkbyte.c
libga_la_SOURCES += tcgmsg/tcgmsg-mpi/clustercheck.c
libga_la_SOURCES += tcgmsg/tcgmsg-mpi/collect.c
libga_la_SOURCES += tcgmsg/tcgmsg-mpi/drand48.c
libga_la_SOURCES += tcgmsg/tcgmsg-mpi/evlog.c
libga_la_SOURCES += tcgmsg/tcgmsg-mpi/evlog.h
libga_la_SOURCES += tcgmsg/tcgmsg-mpi/evon.c
libga_la_SOURCES += tcgmsg/tcgmsg-mpi/llog.c
libga_la_SOURCES += tcgmsg/tcgmsg-mpi/misc.c
libga_la_SOURCES += tcgmsg/tcgmsg-mpi/niceftn.c
libga_la_SOURCES += tcgmsg/tcgmsg-mpi/nxtval-armci.c
libga_la_SOURCES += tcgmsg/tcgmsg-mpi/p2p.c
libga_la_SOURCES += tcgmsg/tcgmsg-mpi/pbeginf.c
libga_la_SOURCES += tcgmsg/tcgmsg-mpi/pfilecopy.c
libga_la_SOURCES += tcgmsg/tcgmsg-mpi/sizeof.c
libga_la_SOURCES += tcgmsg/tcgmsg-mpi/sndrcv.h
libga_la_SOURCES += tcgmsg/tcgmsg-mpi/srftoc.h
libga_la_SOURCES += tcgmsg/tcgmsg-mpi/tcgmsgP.h

AM_CPPFLAGS += -I$(top_srcdir)/tcgmsg/tcgmsg-mpi

check_PROGRAMS += tcgmsg/tests/test_asyn
check_PROGRAMS += tcgmsg/tests/test
if ENABLE_F77
check_PROGRAMS += tcgmsg/tests/testf
endif

#TCGMSG_PARALLEL_TESTS += tcgmsg/tests/test_asyn$(EXEEXT) # useless test
#TCGMSG_PARALLEL_TESTS += tcgmsg/tests/test$(EXEEXT) # requires user input
#if ENABLE_F77
#TCGMSG_PARALLEL_TESTS += tcgmsg/tests/testf$(EXEEXT)
#endif

tcgmsg_tests_test_asyn_SOURCES = tcgmsg/tests/test_asyn.c
tcgmsg_tests_testf_SOURCES     = tcgmsg/tests/testf.F
tcgmsg_tests_test_SOURCES      = tcgmsg/tests/test.c

tcgmsg_tests_testf_LDFLAGS     = $(FLDFLAGS)

nodist_EXTRA_tcgmsg_tests_test_asyn_SOURCES = $(MAYBE_DUMMY_F_LINK)
nodist_EXTRA_tcgmsg_tests_test_SOURCES      = $(MAYBE_DUMMY_F_LINK)

endif # MSG_COMMS_MPI

EXTRA_DIST += tcgmsg/tcgmsg-mpi/README

##############################################################################
# armci
# if TCGMSG is coming from ARMCI
if ARMCI_NETWORK_ARMCI
else
if MSG_COMMS_TCGMSG4
AM_CPPFLAGS += -I$(top_srcdir)/armci/tcgmsg
endif
if MSG_COMMS_TCGMSG5
AM_CPPFLAGS += -I$(top_srcdir)/armci/tcgmsg
endif
endif

##############################################################################
# armci
if ARMCI_NETWORK_ARMCI
else
if ARMCI_SRC_DIR_COMEX
MAYBE_COMEX = comex
else
MAYBE_ARMCI = armci
endif
endif

SUBDIRS = $(MAYBE_COMEX) $(MAYBE_ARMCI) .

if ARMCI_NETWORK_ARMCI
else
if ARMCI_SRC_DIR_COMEX
comex/libarmci.la:
	($(am__cd) comex && $(MAKE) $(AM_MAKEFLAGS) all) || exit 1;
else
armci/libarmci.la:
	($(am__cd) armci && $(MAKE) $(AM_MAKEFLAGS) all) || exit 1;
endif
endif

if ARMCI_SRC_DIR_COMEX
else
endif

if ARMCI_NETWORK_ARMCI
else
if ARMCI_NETWORK_PORTALS
AM_CPPFLAGS += -I$(top_srcdir)/armci/src-portals
AM_CPPFLAGS += -I$(top_srcdir)/armci/tcgmsg
else
if ARMCI_NETWORK_GEMINI
AM_CPPFLAGS += -I$(top_srcdir)/armci/src-gemini
AM_CPPFLAGS += -I$(top_srcdir)/armci/tcgmsg
else
AM_CPPFLAGS += -I$(top_srcdir)/armci/src/include
AM_CPPFLAGS += -I$(top_srcdir)/armci/tcgmsg
endif
endif
endif

# As of ga-5-2 libarmci is no longer rolled up into libga.
libga_la_LIBADD += $(SCALAPACK_LIBS)
libga_la_LIBADD += $(LAPACK_LIBS)
libga_la_LIBADD += $(BLAS_LIBS)
if ARMCI_NETWORK_ARMCI
libga_la_LIBADD += $(ARMCI_NETWORK_LIBS)
libga_la_LIBADD += $(GA_MP_LIBS)
else
if ARMCI_SRC_DIR_COMEX
libga_la_LIBADD += comex/libarmci.la
else
libga_la_LIBADD += armci/libarmci.la
endif
libga_la_LIBADD += $(GA_MP_LIBS)
libga_la_LIBADD += $(ARMCI_NETWORK_LIBS)
endif

##############################################################################
# tascel
# The tascel directory is not automatically configured nor added to SUBDIRS.
# However, it is still distributed with GA.
EXTRA_DIST += tascel
dist-hook: tascel-dist-hook
.PHONY: tascel-dist-hook
tascel-dist-hook:
	rm -rf $(distdir)/tascel/autom4te.cache

##############################################################################
# gfutex
# The gfutex directory is not automatically configured nor added to SUBDIRS.
# However, it is still distributed with GA.
EXTRA_DIST += gfutex
dist-hook: gfutex-dist-hook
.PHONY: gfutex-dist-hook
gfutex-dist-hook:
	rm -rf $(distdir)/gfutex/autom4te.cache

##############################################################################
# the end
#
.PHONY: checkprogs checkprogs-recursive checkprogs-am
checkprogs: checkprogs-recursive
checkprogs-am: $(check_PROGRAMS)

GA_RECURSIVE_TARGETS = checkprogs-recursive
# COPIED FROM MAKEFILE GENERATED BY AUTOMAKE 1.11.1
# Modified to use our own recursive targets.
$(GA_RECURSIVE_TARGETS): all
	@fail= failcom='exit 1'; \
	for f in x $$MAKEFLAGS; do \
	  case $$f in \
	    *=* | --[!k]*);; \
	    *k*) failcom='fail=yes';; \
	  esac; \
	done; \
	dot_seen=no; \
	target=`echo $@ | sed s/-recursive//`; \
	list='$(SUBDIRS)'; for subdir in $$list; do \
	  echo "Making $$target in $$subdir"; \
	  if test "$$subdir" = "."; then \
	    dot_seen=yes; \
	    local_target="$$target-am"; \
	  else \
	    local_target="$$target"; \
	  fi; \
	  ($(am__cd) $$subdir && $(MAKE) $(AM_MAKEFLAGS) $$local_target) \
	  || eval $$failcom; \
	done; \
	if test "$$dot_seen" = "no"; then \
	  $(MAKE) $(AM_MAKEFLAGS) "$$target-am" || exit 1; \
	fi; test -z "$$fail"
# END COPY

# support verbose/silent make rules for additional programs
# m4
M4_V = $(M4__v_$(V))
M4__v_ = $(M4__v_$(AM_DEFAULT_VERBOSITY))
M4__v_0 = @echo "  M4    " $@;
# sed
SED_V = $(SED__v_$(V))
SED__v_ = $(SED__v_$(AM_DEFAULT_VERBOSITY))
SED__v_0 = @echo "  SED   " $@;
# python
PYTHON_V = $(PYTHON__v_$(V))
PYTHON__v_ = $(PYTHON__v_$(AM_DEFAULT_VERBOSITY))
PYTHON__v_0 = @echo "  PYTHON" $@;
# cp
CP = cp
CP_V = $(CP__v_$(V))
CP__v_ = $(CP__v_$(AM_DEFAULT_VERBOSITY))
CP__v_0 = @echo "  CP    " $@;

.PHONY: pfiles clean-pfiles
if MSG_COMMS_TCGMSG4
pfiles: $(check_PROGRAMS)
	for p in $(check_PROGRAMS); do \
		echo "`whoami` `hostname` $(NPROCS) `pwd`/$$p /tmp" > $$p.p; \
	done
check: pfiles
clean-pfiles:
	find . -name '*.p' -exec rm {} \; ;
clean-local: clean-pfiles
MAYBE_PFILES = pfiles
endif # MSG_COMMS_TCGMSG4

##############################################################################
# test suite
#
# Tests were broken up into groups depending on the piece of the distribution
# which they tested.
#
# Some tests were commented out either because they required an input file or
# they were failing for unknown reasons and we didn't want to further debug.
# Some tests were based on old versions of GA which no longer applied but we
# kept anyway (see pario).
#
# We leave pario tests out of the main set of tests because not all systems
# can handle the disk requirements.
#
SERIAL_TESTS =
SERIAL_TESTS += $(MA_SERIAL_TESTS)
SERIAL_TESTS += $(GLOBAL_SERIAL_TESTS)
SERIAL_TESTS += $(CXX_SERIAL_TESTS)

SERIAL_TESTS_XFAIL =
SERIAL_TESTS_XFAIL += $(MA_SERIAL_TESTS_XFAIL)
SERIAL_TESTS_XFAIL += $(GLOBAL_SERIAL_TESTS_XFAIL)
SERIAL_TESTS_XFAIL += $(CXX_SERIAL_TESTS_XFAIL)

PARALLEL_TESTS =
PARALLEL_TESTS += $(MA_PARALLEL_TESTS)
PARALLEL_TESTS += $(GLOBAL_PARALLEL_TESTS)
PARALLEL_TESTS += $(CXX_PARALLEL_TESTS)

PARALLEL_TESTS_XFAIL =
PARALLEL_TESTS_XFAIL += $(MA_PARALLEL_TESTS_XFAIL)
PARALLEL_TESTS_XFAIL += $(GLOBAL_PARALLEL_TESTS_XFAIL)
PARALLEL_TESTS_XFAIL += $(CXX_PARALLEL_TESTS_XFAIL)

TESTS =
TESTS += $(SERIAL_TESTS)
TESTS += $(PARALLEL_TESTS)

XFAIL_TESTS =
XFAIL_TESTS += $(SERIAL_TESTS_XFAIL)
XFAIL_TESTS += $(PARALLEL_TESTS_XFAIL)

if MSG_COMMS_TCGMSG4
LOG_COMPILER = \
maybe_tcgexec=`if echo "$(SERIAL_TESTS)" | $(GREP) "$$p" > /dev/null; then echo ""; else echo "$(TCGEXEC)"; fi`; eval $$maybe_tcgexec
else
if MSG_COMMS_TCGMSG5
LOG_COMPILER = tst="$$dir$$f -np $(NPROCS)"; eval
else
if CROSS_COMPILING
maybe_mpiexec=`if echo "$(SERIAL_TESTS)" | $(GREP) "$$p" > /dev/null; then echo "$(MPIEXEC)" | $(SED) 's/%NP%/1/'; else echo "$(MPIEXEC)" | $(SED) 's/%NP%/$(NPROCS)/'; fi`; eval $$maybe_mpiexec
else
LOG_COMPILER = \
maybe_mpiexec=`if echo "$(SERIAL_TESTS)" | $(GREP) "$$p" > /dev/null; then echo ""; else if echo "$(MPIEXEC)" | $(GREP) "%NP%" > /dev/null; then echo "$(MPIEXEC)" | $(SED) 's/%NP%/$(NPROCS)/'; else echo "$(MPIEXEC)"; fi; fi`; eval $$maybe_mpiexec
endif # CROSS_COMPILING
endif # MSG_COMMS_TCGMSG5
endif # MSG_COMMS_TCGMSG4

.PHONY: check-ga
check-ga: $(all-am) $(check_LTLIBRARIES) $(check_PROGRAMS) $(MAYBE_PFILES)
	$(MAKE) $(AM_MAKEFLAGS) check-TESTS

.PHONY: check-ma
check-ma: $(all-am) $(check_LTLIBRARIES) $(MA_TESTS) $(MAYBE_PFILES)
	$(MAKE) $(AM_MAKEFLAGS) TESTS="$(MA_TESTS)" XFAIL_TESTS="$(MA_TESTS_XFAIL)" check-TESTS

.PHONY: check-global
check-global: $(all-am) $(check_LTLIBRARIES) $(GLOBAL_TESTS) $(MAYBE_PFILES)
	$(MAKE) $(AM_MAKEFLAGS) TESTS="$(GLOBAL_TESTS)" XFAIL_TESTS="$(GLOBAL_TESTS_XFAIL)" check-TESTS

.PHONY: check-pario
check-pario: $(all-am) $(check_LTLIBRARIES) $(PARIO_TESTS) $(MAYBE_PFILES)
	$(MAKE) $(AM_MAKEFLAGS) TESTS="$(PARIO_TESTS)" XFAIL_TESTS="$(PARIO_TESTS_XFAIL)" check-TESTS

.PHONY: check-ga++
if CXX_BINDINGS
check-ga++: $(all-am) $(check_LTLIBRARIES) $(CXX_TESTS) $(MAYBE_PFILES)
	$(MAKE) $(AM_MAKEFLAGS) TESTS="$(CXX_TESTS)" XFAIL_TESTS="$(CXX_TESTS_XFAIL)" check-TESTS
endif

# Because pario tests are not listed in the TESTS variable, the log targets
# are not generated by automake. Do them manually here.
pario/dra/bign.log: pario/dra/bign$(EXEEXT)
	@p='pario/dra/bign$(EXEEXT)'; $(am__check_pre) $(LOG_COMPILE) "$$tst" $(am__check_post)
pario/dra/ntestc.log: pario/dra/ntestc$(EXEEXT)
	@p='pario/dra/ntestc$(EXEEXT)'; $(am__check_pre) $(LOG_COMPILE) "$$tst" $(am__check_post)
pario/dra/perfn.log: pario/dra/perfn$(EXEEXT)
	@p='pario/dra/perfn$(EXEEXT)'; $(am__check_pre) $(LOG_COMPILE) "$$tst" $(am__check_post)
pario/dra/rate.log: pario/dra/rate$(EXEEXT)
	@p='pario/dra/rate$(EXEEXT)'; $(am__check_pre) $(LOG_COMPILE) "$$tst" $(am__check_post)
