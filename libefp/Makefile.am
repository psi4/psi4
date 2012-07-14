EXTRA_DIST =
CLEANFILES =
ACLOCAL_AMFLAGS = -I m4 ${ACLOCAL_FLAGS}
AM_MAKEFLAGS = --no-print-directory

AM_CPPFLAGS = \
	-include $(top_builddir)/config.h \
	-DSYSCONFDIR=\""$(sysconfdir)"\" \
	-DLIBEXECDIR=\""$(libexecdir)"\" \
	-I${top_srcdir}/src

AM_CFLAGS = ${my_CFLAGS} ${OPENMP_CFLAGS} \
	-fvisibility=hidden

AM_LDFLAGS = \
	-Wl,--gc-sections \
	-Wl,--as-needed

SED_PROCESS = \
	$(AM_V_GEN)$(MKDIR_P) $(dir $@) && $(SED) \
	-e 's,@VERSION\@,$(VERSION),g' \
	-e 's,@prefix\@,$(prefix),g' \
	-e 's,@exec_prefix\@,$(exec_prefix),g' \
	-e 's,@libdir\@,$(libdir),g' \
	-e 's,@includedir\@,$(includedir),g' \
	< $< > $@ || rm $@

%.pc: %.pc.in Makefile
	$(SED_PROCESS)

LIBEFP_CURRENT = 0
LIBEFP_REVISION = 0
LIBEFP_AGE = 0

include_HEADERS = src/efp.h
lib_LTLIBRARIES = src/libefp.la

bin_PROGRAMS = \
	efpmd/efpmd \
	efpmd/cubegen/cubegen

src_libefp_la_SOURCES = \
	src/ai_disp.c \
	src/ai_xr.c \
	src/cblas.h \
	src/chtr.c \
	src/disp.c \
	src/disp.h \
	src/efp.c \
	src/efp_private.h \
	src/elec.c \
	src/elec.h \
	src/elec_grad.c \
	src/int.c \
	src/int.h \
	src/int_shift.h \
	src/math_util.h \
	src/parse.c \
	src/parse.h \
	src/phys_const.h \
	src/pol.c \
	src/terms.h \
	src/xr.c

EXTRA_DIST += \
	Doxyfile \
	LICENSE \
	README.md \
	src/libefp.sym \
	fraglib/makefp.inp \
	fraglib/ch3oh.efp \
	fraglib/h2o.efp \
	fraglib/nh3.efp

src_libefp_la_LDFLAGS = $(AM_LDFLAGS) \
	-version-info $(LIBEFP_CURRENT):$(LIBEFP_REVISION):$(LIBEFP_AGE) \
	-Wl,--version-script=$(top_srcdir)/src/libefp.sym
src_libefp_la_DEPENDENCIES = ${top_srcdir}/src/libefp.sym

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = src/libefp.pc
EXTRA_DIST += src/libefp.pc.in
CLEANFILES += src/libefp.pc

efpmd_efpmd_SOURCES = \
	efpmd/cg.c \
	efpmd/common.c \
	efpmd/common.h \
	efpmd/grad.c \
	efpmd/main.c \
	efpmd/nve.c \
	efpmd/nvt.c \
	efpmd/parse.c \
	efpmd/parse.h \
	efpmd/sim.h \
	efpmd/sp.c

efpmd_efpmd_LDADD = src/libefp.la

efpmd_cubegen_cubegen_SOURCES = \
	efpmd/cubegen/cubegen.c

EXTRA_DIST += \
	efpmd/README.md \
	efpmd/tests/h2o-ch3oh_cg.in \
	efpmd/tests/h2o-nh3_grad.in \
	efpmd/tests/h2o_sp.in

TESTS = \
	tests/elec_1 \
	tests/pol_1 \
	tests/disp_1 \
	tests/xr_1 \
	tests/total_1a \
	tests/elec_2 \
	tests/pol_2 \
	tests/disp_2 \
	tests/xr_2 \
	tests/total_2a \
	tests/elec_3 \
	tests/pol_3 \
	tests/disp_3 \
	tests/xr_3

EXTRA_DIST += \
	tests/geometry_1.h \
	tests/geometry_2.h \
	tests/geometry_3.h

check_PROGRAMS = $(TESTS)

TESTS_COMMON_FILES = tests/test_common.c tests/test_common.h
TESTS_LDADD = src/libefp.la
TESTS_CPPFLAGS = $(AM_CPPFLAGS) -DABS_TOP_SRCDIR=\"$(abs_top_srcdir)\"

tests_elec_1_SOURCES = tests/elec_1.c $(TESTS_COMMON_FILES)
tests_elec_1_LDADD = $(TESTS_LDADD)
tests_elec_1_CPPFLAGS = $(TESTS_CPPFLAGS)

tests_pol_1_SOURCES = tests/pol_1.c $(TESTS_COMMON_FILES)
tests_pol_1_LDADD = $(TESTS_LDADD)
tests_pol_1_CPPFLAGS = $(TESTS_CPPFLAGS)

tests_disp_1_SOURCES = tests/disp_1.c $(TESTS_COMMON_FILES)
tests_disp_1_LDADD = $(TESTS_LDADD)
tests_disp_1_CPPFLAGS = $(TESTS_CPPFLAGS)

tests_xr_1_SOURCES = tests/xr_1.c $(TESTS_COMMON_FILES)
tests_xr_1_LDADD = $(TESTS_LDADD)
tests_xr_1_CPPFLAGS = $(TESTS_CPPFLAGS)

tests_total_1a_SOURCES = tests/total_1a.c $(TESTS_COMMON_FILES)
tests_total_1a_LDADD = $(TESTS_LDADD)
tests_total_1a_CPPFLAGS = $(TESTS_CPPFLAGS)

tests_elec_2_SOURCES = tests/elec_2.c $(TESTS_COMMON_FILES)
tests_elec_2_LDADD = $(TESTS_LDADD)
tests_elec_2_CPPFLAGS = $(TESTS_CPPFLAGS)

tests_pol_2_SOURCES = tests/pol_2.c $(TESTS_COMMON_FILES)
tests_pol_2_LDADD = $(TESTS_LDADD)
tests_pol_2_CPPFLAGS = $(TESTS_CPPFLAGS)

tests_disp_2_SOURCES = tests/disp_2.c $(TESTS_COMMON_FILES)
tests_disp_2_LDADD = $(TESTS_LDADD)
tests_disp_2_CPPFLAGS = $(TESTS_CPPFLAGS)

tests_xr_2_SOURCES = tests/xr_2.c $(TESTS_COMMON_FILES)
tests_xr_2_LDADD = $(TESTS_LDADD)
tests_xr_2_CPPFLAGS = $(TESTS_CPPFLAGS)

tests_total_2a_SOURCES = tests/total_2a.c $(TESTS_COMMON_FILES)
tests_total_2a_LDADD = $(TESTS_LDADD)
tests_total_2a_CPPFLAGS = $(TESTS_CPPFLAGS)

tests_elec_3_SOURCES = tests/elec_3.c $(TESTS_COMMON_FILES)
tests_elec_3_LDADD = $(TESTS_LDADD)
tests_elec_3_CPPFLAGS = $(TESTS_CPPFLAGS)

tests_pol_3_SOURCES = tests/pol_3.c $(TESTS_COMMON_FILES)
tests_pol_3_LDADD = $(TESTS_LDADD)
tests_pol_3_CPPFLAGS = $(TESTS_CPPFLAGS)

tests_disp_3_SOURCES = tests/disp_3.c $(TESTS_COMMON_FILES)
tests_disp_3_LDADD = $(TESTS_LDADD)
tests_disp_3_CPPFLAGS = $(TESTS_CPPFLAGS)

tests_xr_3_SOURCES = tests/xr_3.c $(TESTS_COMMON_FILES)
tests_xr_3_LDADD = $(TESTS_LDADD)
tests_xr_3_CPPFLAGS = $(TESTS_CPPFLAGS)
