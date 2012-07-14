AC_PREREQ([2.62])
AC_INIT([libefp],
	[0.1],
	[ilya.kaliman@gmail.com],
	[libefp],
	[http://libefp.github.com/])
AC_CONFIG_SRCDIR([src/efp.c])
AC_CONFIG_AUX_DIR([build-aux])
AM_INIT_AUTOMAKE([
	foreign 1.11 -Wall -Wno-portability silent-rules tar-pax
	no-dist-gzip dist-xz subdir-objects color-tests
])
AC_PROG_CC_STDC
AC_OPENMP
AC_USE_SYSTEM_EXTENSIONS
AC_SYS_LARGEFILE
AC_CONFIG_MACRO_DIR([m4])
AM_PROG_CC_C_O
AM_SILENT_RULES([yes])
LT_INIT([pic-only])
AC_PREFIX_DEFAULT([/usr])

AC_PROG_SED
AC_PROG_MKDIR_P

developer_CFLAGS="-Wall -Wextra -pedantic \
-Wno-unused-parameter \
-Wmissing-declarations -Wmissing-prototypes \
-Wpointer-arith -Wsign-compare \
-Wstrict-prototypes -Wshadow -Wformat-security"

AC_ARG_ENABLE([developer-flags],
	AS_HELP_STRING([--enable-developer-flags], [use developer flags]),
	[enable_developer_flags=$developer_CFLAGS],
	[enable_developer_flags=""])
AC_SUBST([my_CFLAGS], [$enable_developer_flags])

AC_CONFIG_HEADERS(config.h)
AC_CONFIG_FILES([
	Makefile
])

AC_OUTPUT
AC_MSG_RESULT([
        configuration
        =============

        prefix:             ${prefix}
        sysconfdir:         ${sysconfdir}
        libdir:             ${libdir}
        includedir:         ${includedir}

        compiler:           ${CC}
        cflags:             ${CFLAGS}
        ldflags:            ${LDFLAGS}

        developer cflags:   ${my_CFLAGS}
])
