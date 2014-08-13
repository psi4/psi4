# ===========================================================================
#
# SYNOPSIS
#
#   GA_PARIO
#
# DESCRIPTION
#
#   This macro tries to find out how to compile pario.
#   This was converted from pario/makefile.h and pario/*/GNUmakefile
#   Defines the following precious variables
#   PARIO_CPPFLAGS
#   PARIO_LDFLAGS
#   PARIO_CFLAGS
#   PARIO_FFLAGS
#

AC_DEFUN([GA_PARIO], [

dnl ##########################################################################
dnl FROM pario/makefile.h
dnl ##########################################################################
OSNAME=`uname`
PARIO_CPPFLAGS=
if test x$OSNAME = xAIX ; then
  if /usr/bin/oslevel | awk -F. '{ if ($1 > 5 || ($1 == 5 && $2 > 1)) exit 0 }'
  then
    PARIO_CPPFLAGS="$PARIO_CPPFLAGS -DAIX52"
  fi
  if /usr/sbin/lsdev -C -l aio0 2>&1 | grep Legacy ; then
    PARIO_CPPFLAGS="$PARIO_CPPFLAGS -D_AIO_AIX_SOURCE"
  fi
fi
if test x$LARGE_FILES != x ; then
  if test x$OSNAME = xAIX ; then
    if /usr/bin/oslevel|awk -F. '{ if ($1 > 4 || ($1 == 4 && $2 > 1)) exit 0 }'
    then
      PARIO_CPPFLAGS="$PARIO_CPPFLAGS -D_LARGE_FILES -D_LARGE_FILE_API"
    fi
    if /usr/bin/oslevel|awk -F. '{ if ($1 == 4 && $2 == 2 && $3 <= 0 ) exit 0 }'
    then
      PARIO_CPPFLAGS="$PARIO_CPPFLAGS -DNOAIO"
    fi
  fi
  if test x$TARGET = xSOLARIS ; then
    PARIO_CPPFLAGS="$PARIO_CPPFLAGS -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64"
    PARIO_CFLAGS=`getconf LFS_CFLAGS`
  fi
  if test x$TARGET = xLINUX ; then
    PARIO_CPPFLAGS="$PARIO_CPPFLAGS -D_LARGEFILE64_SOURCE"
    PARIO_CFLAGS=`getconf LFS_CFLAGS`
  fi
  if test x$TARGET = xBGL -o x$TARGET = xBGP ; then
    PARIO_CPPFLAGS="$PARIO_CPPFLAGS -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE"
  fi
  if test x$TARGET = xHPUX ; then
    PARIO_CPPFLAGS="$PARIO_CPPFLAGS -D_LARGEFILE64_SOURCE"
    PARIO_CFLAGS=`getconf XBS5_ILP32_OFFBIG_CFLAGS`
  fi
  if test x$TARGET = xHPUX64 ; then
    PARIO_CPPFLAGS="$PARIO_CPPFLAGS -D_LARGEFILE64_SOURCE"
    PARIO_CFLAGS=`getconf XBS5_LP64_OFF64_CFLAGS`
  fi
  PARIO_CPPFLAGS="$PARIO_CPPFLAGS -DLARGE_FILES"
fi
if test x$TARGET = xDECOSF ; then
  PARIO_LDFLAGS="-laio -lpthreads"
fi
if test x$USE_LINUXAIO != x ; then
  PARIO_CPPFLAGS="$PARIO_CPPFLAGS -DLINUXAIO"
  PARIO_LDFLAGS="$PARIO_LDFLAGS -lrt"
fi

dnl ##########################################################################
dnl FROM pario/elio/GNUmakefile
dnl ##########################################################################

dnl on platforms with Posix AIO you can choose not to use it by defining NOAIO
if test x$NOAIO != x ; then
  PARIO_CPPFLAGS="$PARIO_CPPFLAGS -DNOAIO"
fi

if test x$PABLO != x ; then
  PARIO_CPPFLAGS="$PARIO_CPPFLAGS -DPABLO"
fi

dnl ##########################################################################
dnl FROM pario/eaf/GNUmakefile
dnl ##########################################################################
PARIO_CPPFLAGS="$PARIO_CPPFLAGS -DEAF_STATS"
if test x$ga_armci_network = xPORTALS ; then
  PARIO_CPPFLAGS="$PARIO_CPPFLAGS -DCRAY_XT"
fi

dnl ##########################################################################
dnl FROM pario/dra/GNUmakefile
dnl ##########################################################################
if test x$F77 = xfrt ; then
  PARIO_FFLAGS=-O2
fi

dnl ##########################################################################
dnl FROM pario/sf/GNUmakefile
dnl ##########################################################################
dnl nothing


AC_SUBST([PARIO_CPPFLAGS])
AC_SUBST([PARIO_LDFLAGS])
AC_SUBST([PARIO_CFLAGS])
AC_SUBST([PARIO_FFLAGS])

])dnl GA_PARIO
