dnl
dnl  Local macros.  This are based on files provided with autoconf
dnl so the GNU GPL applies to them.
dnl
dnl
dnl See if the given C program is processed by the C compiler OK.  Does
dnl not necessarily compile.  If you want to compile give a -c to the compiler
dnl arguments.  To compile and link give a -o conftest.
dnl arg1: echo text
dnl arg2: optional code (main is already provided)
dnl arg3: additional compiler arguments
dnl arg4: success action
dnl arg5: fail action
dnl
define(AC_CC_PROCESS_CHECK,
[AC_PROVIDE([$0])dnl
ifelse([$1], , , [AC_MSG_CHECKING([for $1])]
)dnl
cat > conftest.c <<EOF
[$2]
main() {}
EOF
dnl Don't try to run the program, which would prevent cross-configuring.
if eval $CC [$3] conftest.c >/dev/null 2>&1; then
  ifelse([$4], , :, [rm -rf conftest*
  $4
])
ifelse([$5], , , [else
  rm -rf conftest*
  $5
])dnl
fi
rm -f conftest*
AC_MSG_RESULT(OK)]
)dnl
dnl
dnl
dnl See if the given C++ program is processed by the C++ compiler OK.  Does
dnl not necessarily compile.  If you want to compile give a -c to the compiler
dnl arguments.  To compile and link give a -o conftest.
dnl arg1: echo text
dnl arg2: optional code (main is already provided)
dnl arg3: additional compiler arguments
dnl arg4: success action
dnl arg5: fail action
dnl
define(AC_CXX_PROCESS_CHECK,
[AC_PROVIDE([$0])dnl
ifelse([$1], , , [AC_MSG_CHECKING([for $1])]
)dnl
cat > conftest.cc <<EOF
[$2]
int main(int argc, char** argv) {}
EOF
dnl Don't try to run the program, which would prevent cross-configuring.
if eval $CXX $CPPFLAGS [$3] conftest.cc >/dev/null 2>&1; then
  ifelse([$4], , :, [rm -rf conftest*
  $4
  AC_MSG_RESULT(yes)
])
ifelse([$5], , , [else
  rm -rf conftest*
  $5
  AC_MSG_RESULT(no)
])dnl
fi
rm -f conftest*]
)dnl
dnl
dnl
dnl Check for LaTeX
dnl
define(AC_PROG_LATEX,
[AC_PROVIDE([$0])
AC_CHECK_PROG(LATEX, latex, latex)]
)dnl

dnl
dnl Check for dvips
dnl
define(AC_PROG_DVIPS,
[AC_PROVIDE([$0])
AC_CHECK_PROG(DVIPS, dvips, dvips)]
)dnl


dnl
dnl Check for LaTeX2HTML
dnl
define(AC_PROG_LATEX2HTML,
[AC_PROVIDE([$0])
AC_CHECK_PROG(LATEX2HTML, latex2html, latex2html)]
)dnl

define(AC_PROG_BIBTEX,
[AC_PROVIDE([$0])
AC_CHECK_PROG(BIBTEX, bibtex, bibtex)]
)dnl


# Function to look for Apple vecLib on OS X systems
define(AC_CHECK_VECLIB,
[AC_PROVIDE([$0])
  echo "Entering ac_check_veclib"
  SAVE_LIBS=$LIBS
  echo "$target_vendor"
  if eval $CC 
  LIBS=$SAVE_LIBS
  echo "Finished with ac_check_veclib"
]
)
