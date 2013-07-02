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
AC_DEFUN([AC_CXX_PROCESS_CHECK],
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

