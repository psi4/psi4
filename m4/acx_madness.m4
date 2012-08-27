AC_DEFUN([ACX_MADNESS], [
  HAVE_MADNESS="0"
  BUILD_MADNESS="0"
  CONFIG_MADNESS="0"

  AC_MSG_CHECKING([if we are using MADNESS])
  AC_ARG_WITH([madness],
    [AC_HELP_STRING([--with-madness@<:@=madness-libs@:>@],
      [Enable building PSI4 with MADNESS. If no MADNESS library is specified, the internal MADNESS will be compiled and built. @<:@default=no@:>@]) ],
    [
      case $withval in
        yes)
          if test $HAVE_MPI = "1"; then
            madlibs='$(top_objdir)/madness/src/lib'
            madinc='-I$(top_srcdir)/madness/src/lib -I$(top_srcdir)/madness/include -I$(top_objdir)/madness/include'
            worldfile="$srcdir/madness/src/lib/world/world.h"
            BUILD_MADNESS="1"
            HAVE_MADNESS="1" 
            CONFIG_MADNESS="1"
          fi
          AC_MSG_RESULT([yes])
        ;;
        no)
            HAVE_MADNESS="0"
            BUILD_MADNESS="0"
            CONFIG_MADNESS="0"
            AC_MSG_RESULT([no])
        ;;
        *)
          if test $HAVE_MPI = "1"; then
             madlibs="$withval/lib"
             madinc="-I$withval/include"
             worldfile="$withval/include/world/world.h"
             BUILD_MADNESS="0"
             AC_MSG_RESULT([yes])
          else
            AC_MSG_RESULT([no])
          fi
        ;;
      esac
    ],
    [
      HAVE_MADNESS="0"
      BUILD_MADNESS="0"
      CONFIG_MADNESS="0"
      AC_MSG_RESULT([no])
    ]
  )

  if test $HAVE_MPI = "1"; then
    AC_MSG_CHECKING([if we are using the internal MADNESS])

    if test $BUILD_MADNESS = "1"; then
      AC_MSG_RESULT([yes])

      AC_MSG_CHECKING([if we are configuring MADNESS])
      # if we are building madness, see if madness should be reconfigured.  If it is reconfigured
      # madness will automatically be rebuilt.
      AC_ARG_ENABLE([madness-config],
        [AC_HELP_STRING([--enable-madness-config@<:@=yes|no@:>@],        [Specify whether we need to configure MADNESS. If MADNESS runs configure, it will always recompile MADNESS.  Disabling MADNESS configure will prevent MADNESS from recompiling. @<:@default=yes@:>@]) ],
        [
          case $enableval in
            yes)
              CONFIG_MADNESS="1"
              AC_MSG_RESULT([yes])
              ;;
            *)
              CONFIG_MADNESS="0"
              AC_MSG_RESULT([no])
              ;;
          esac
        ],
        [ 
          CONFIG_MADNESS="1"
          AC_MSG_RESULT([yes])
        ]
      )
    else 
      AC_MSG_RESULT([no])

      # check to see if the user supplied madness correctly
      AC_LANG(C++)
      ac_link='$CXX -o conftest$ac_exeext $CXXFLAGS $CPPFLAGS $LDFLAGS -L$madlibs -lMADworld $madinc conftest.$ac_ext $LIBS >&5'
      AC_CHECK_FILE($worldfile, HAVE_MADNESS="1", HAVE_MADNESS="0")
      AC_MSG_CHECKING([if we can link to MADNESS])
      AC_TRY_LINK([#include <$worldfile>],[ madness::finalize(); ],
        [HAVE_MADNESS="1"; AC_MSG_RESULT(yes)],
        [HAVE_MADNESS="0"; AC_MSG_RESULT(no)])
      ac_link='$CC -o conftest$ac_exeext $CFLAGS $CPPFLAGS $LDFLAGS conftest.$ac_ext $LIBS >&5'
      AC_LANG(C)
    fi


    
    if test $HAVE_MADNESS = "1"; then
      AC_DEFINE_UNQUOTED(HAVE_MADNESS,1)
      AC_DEFINE_UNQUOTED(BUILD_MADNESS,1)
      AC_DEFINE_UNQUOTED(CONFIG_MADNESS,1)
      
      AC_SUBST(madlibs)
      AC_SUBST(madinc)
      AC_SUBST(HAVE_MPI)
      AC_SUBST(HAVE_MADNESS)
      AC_SUBST(BUILD_MADNESS)
      AC_SUBST(CONFIG_MADNESS)
      
    else
      HAVE_MPI="0"
      HAVE_MADNESS="0"
      BUILD_MADNESS="0"
      CONFIG_MADNESS="0"
      AC_DEFINE_UNQUOTED(HAVE_MPI,0)
      AC_DEFINE_UNQUOTED(HAVE_MADNESS,0)
      AC_DEFINE_UNQUOTED(BUILD_MADNESS,0)
      AC_DEFINE_UNQUOTED(CONFIG_MADNESS,0)
    fi
 fi
 AC_SUBST(HAVE_MPI)
 AC_SUBST(HAVE_MADNESS)
 AC_SUBST(BUILD_MADNESS)
 AC_SUBST(CONFIG_MADNESS)
])
