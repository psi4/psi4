# AC_PYTHON([MINIMUM-VERSION],[ACTION-IF-FOUND],[ACTION-IF-NOTFOUND])
# -------------------------------------------------------------------
# Run ACTION-IF-FOUND if a Python interpreter meeting the
# MINIMUM-VERSION requirement is found. Else Run ACTION-IF-NOTFOUND
# if no valid Python interpreter is found. Can be overriden if the
# user sets the PYTHON environmental variable.
AC_DEFUN([AC_PYTHON],
[
  AC_ARG_VAR(PYTHON,
[The `Python' interpreter to use. Defaults to the first program found
out of: `python', `python2', `python2.6', `python2.5', `python2.4'.
Make sure you have installed the developer libraries, too.])
  dnl Find a Python interpreter and its matching config program.
  m4_define_default([_AM_PYTHON_INTERPRETER_LIST],
                    [python python2 python2.6 dnl
python2.5 python2.4])

  # turn back to default flags
  ac_save_CPPFLAGS="$CPPFLAGS"
  ac_save_LIBS="$LIBS"
  
  m4_if([$1],[],[
    dnl No version check is needed
    if test -z "$PYTHON"; then
      AC_PATH_PROGS([PYTHON], _AM_PYTHON_INTERPRETER_LIST, :)
    fi
    am_display_PYTHON=python
  ],[
    dnl A Version check is needed.
    if test -n "$PYTHON"; then
      # If the user set $PYTHON, use it and don't search something else.
      AC_MSG_CHECKING([whether $PYTHON version >= $1])
      AM_PYTHON_CHECK_VERSION([$PYTHON], [$1],
                              [AC_MSG_RESULT(yes)],
                              [AC_MSG_ERROR(too old)])
      am_display_PYTHON=$PYTHON
    else
      # Otherwise, try each interpreter until we find one that satisfies VERSION
      AC_CACHE_CHECK([for a Python interpreter with version >= $1],
        [am_cv_pathless_PYTHON],[
        for am_cv_pathless_PYTHON in _AM_PYTHON_INTERPRETER_LIST none; do
          test "$am_cv_pathless_PYTHON" = none && break
          AM_PYTHON_CHECK_VERSION([$am_cv_pathless_PYTHON], [$1], [break])
        done])
      # Set $PYTHON to the absolute pather of $am_cv_pathless_PYTHON
      if test "$am_cv_pathless_PYTHON" = none; then
        PYTHON=:
      else
        AC_PATH_PROG([PYTHON], [$am_cv_pathless_PYTHON])
      fi
      am_display_PYTHON=$am_cv_pathless_PYTHON
    fi
  ])
 
  if test "$PYTHON" = :; then
    dnl Run any user-specific action, or abort.
    m4_default([$3], [AC_MSG_ERROR([no suitable Python interpreter found])])
  else      
    dnl Query Python for its version number.
    AC_CACHE_CHECK([for $am_display_PYTHON version], [am_cv_python_version],
      [am_cv_python_version=`$PYTHON -c "import sys; sys.stdout.write(sys.version[[:3]])"`])
    AC_SUBST([PYTHON_VERSION], [$am_cv_python_version])

    dnl Okay, if we know the name of the Python interpreter then we know the
    dnl name of the python config program
    AC_CACHE_CHECK([for $am_display_PYTHON include statements], [am_cv_python_include],
      [am_cv_python_include=`${PYTHON}-config --includes`])
    AC_SUBST([PYTHON_INCLUDE], [$am_cv_python_include])

    AC_CACHE_CHECK([for $am_display_PYTHON linker flags], [am_cv_python_ldflags],
      [am_cv_python_ldflags=`${PYTHON}-config --ldflags`])
    AC_SUBST([PYTHON_LDFLAGS], [$am_cv_python_ldflags])

    #
    # final check to see if everything compiles alright
    #
    AC_MSG_CHECKING([consistency of all components of python development environment])
    AC_LANG_PUSH([C])
    # save current global flags
    LIBS="$ac_save_LIBS $PYTHON_LDFLAGS"
    CPPFLAGS="$ac_save_CPPFLAGS $PYTHON_INCLUDE"
    AC_TRY_LINK([
        #include <Python.h>
    ],[
        Py_Initialize();
    ],[pythonexists=yes],[pythonexists=no])

    AC_MSG_RESULT([$pythonexists])

        if test ! "$pythonexists" = "yes"; then
       AC_MSG_ERROR([
  Could not link test program to Python. Maybe the main Python library has been
  installed in some non-standard library path. If so, pass it to configure,
  via the LDFLAGS environment variable.
  Example: ./configure LDFLAGS="-L/usr/non-standard-path/python/lib"
  ============================================================================
   ERROR!
   You probably have to install the development version of the Python package
   for your distribution.  The exact name of this package varies among them.
  ============================================================================
       ])
      PYTHON_VERSION=""
    fi
    AC_LANG_POP
    # turn back to default flags
    CPPFLAGS="$ac_save_CPPFLAGS"
    LIBS="$ac_save_LIBS"

    dnl This is the end of all things
    $2
  fi
])

