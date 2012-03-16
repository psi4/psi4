AC_DEFUN([ACX_MADNESS], [
  HAVE_MPI="0"
  HAVE_MADNESS="0"
  BUILD_MADNESS="0"
  CONFIG_MADNESS="0"
  AC_ARG_WITH([madness],
    [AC_HELP_STRING([--with-madness@<:@=madness-libs@:>@],
      [Enable building PSI4 with MADNESS. If no MADNESS library is specified, the internal MADNESS will be compiled and built. @<:@default=no@:>@]) ],
    [
      case $withval in
        yes)
          # Check for an MPI compiler
          ACX_MPI

          if test $HAVE_MPI = "1"; then
            madlibs='$(top_objdir)/madness/src/lib'
            madinc='-I$(top_srcdir)/madness/src/lib -I$(top_srcdir)/madness/include -I$(top_objdir)/madness/include'
            worldfile="$srcdir/madness/src/lib/world"
            CONFIG_MADNESS="1"
            BUILD_MADNESS="1"
          fi
        ;;
        no)
            HAVE_MPI="0"
            HAVE_MADNESS="0"
            BUILD_MADNESS="0"
            CONFIG_MADNESS="0"
        ;;
        *)
          # Check for an MPI compiler
          ACX_MPI

          if test $HAVE_MPI = "1"; then
             madlibs="$withval/lib"
             madinc="-I$withval/include"
             worldfile="$withval/include/world/world.h"
             CONFIG_MADNESS="0"
             BUILD_MADNESS="0"
          fi
        ;;
      esac
    ],
    [
      HAVE_MPI="0"
      HAVE_MADNESS="0"
      BUILD_MADNESS="0"
      CONFIG_MADNESS="0"
    ]
  )

  if test $HAVE_MPI = "1"; then
    if test $BUILD_MADNESS = "1"; then
      AC_ARG_ENABLE([madness-configure],
        [AC_HELP_STRING([--enable-madness-configure@<:@=yes|no@:>@],
        [Specify whether we need to configure MADNESS. If MADNESS runs configure, it will always recompile.  Disabling MADNESS configure will prevent MADNESS from recompiling. @<:@default=yes@:>@]) ],
        [
          case $enableval in
            yes)
              CONFIG_MADNESS="1"
              ;;
            *)
              CONFIG_MADNESS="0"
              ;;
          esac
        ],
        [
          CONFIG_MADNESS="1"
        ]
      )
    fi

    AC_CHECK_FILE($worldfile, HAVE_MADNESS="1", HAVE_MADNESS="0")

    if test $HAVE_MADNESS = "1"; then
      AC_DEFINE_UNQUOTED(HAVE_MADNESS,1)

      AC_SUBST(madlibs)
      AC_SUBST(madinc)
      AC_SUBST(HAVE_MPI)
      AC_SUBST(HAVE_MADNESS)
      AC_SUBST(BUILD_MADNESS)

      acx_mpi_save_CC="$CC"
      acx_mpi_save_CXX="$CXX"
      CC="$MPICC"
      CXX="$MPICXX"
      AC_SUBST(MPICXX)
      AC_SUBST(MPICC)
    else
      HAVE_MPI="0"
      HAVE_MADNESS="0"
      BUILD_MADNESS="0"
      CONFIG_MADNESS="0"
    fi

  fi


])
