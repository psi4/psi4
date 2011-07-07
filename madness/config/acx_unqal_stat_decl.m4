AC_DEFUN([ACX_UNQUALIFIED_STATIC_DECL], [
    AC_LANG_PUSH([C++])
    AC_MSG_CHECKING([if unqualified static declarations are considered])
    AC_COMPILE_IFELSE(
      [AC_LANG_PROGRAM(
        [[template <typename T>
          static inline void f(T* a) {};
          template <typename T> void g(T* a) { f(a); }
          template void g(int* a);]],
        [[]])],
      [AC_DEFINE(HAVE_UNQUALIFIED_STATIC_DECL, [1], [Set if compiler will instantiate static templates]) AC_MSG_RESULT([yes])],
      [AC_MSG_RESULT([no])])
    AC_LANG_POP([C++])
])
