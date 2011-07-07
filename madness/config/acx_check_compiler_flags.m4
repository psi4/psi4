#ACX_CHECK_COMPILER_FLAGS(compiler, flag_var, flag, success, fail)
AC_DEFUN([ACX_CHECK_COMPILER_FLAG], [
  AC_LANG_SAVE
  
  AC_LANG([$1])
  acx_check_compiler_flags="no"
  acx_check_compiler_flags_save=$[$2]
  [$2]="$3"
  AC_MSG_CHECKING([whether $1 compiler accepts $3])
  AC_COMPILE_IFELSE([AC_LANG_PROGRAM([],[])],
    [
      acx_check_compiler_flags="yes"
      AC_MSG_RESULT([yes])
    ], [
      AC_MSG_RESULT([no])
    ])
    
  [$2]=$acx_check_compiler_flags_save
  AC_LANG_RESTORE
  
  AS_IF([test $acx_check_compiler_flags != no], [$4], [$5])
])