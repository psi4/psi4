AC_DEFUN([ACX_GNU_HASHMAP], [
    AC_LANG_PUSH([C++])
    gotgnuhashmap=no
    AC_MSG_CHECKING([for GNU hashmap in '<ext/hash_map>' and namespace __gnu_cxx])
    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[        
#include <ext/hash_map>
using __gnu_cxx::hash_map; ]])],
        [ AC_DEFINE(GNU_HASHMAP_NAMESPACE, __gnu_cxx, [GNU hashmap namespace]) 
          AC_DEFINE(INCLUDE_EXT_HASH_MAP,[1],[If defined header is in ext directory])
          AC_MSG_RESULT([yes])
          gotgnuhashmap=yes ],
        [AC_MSG_RESULT([no])])

    if test $gotgnuhashmap = no; then
        AC_MSG_CHECKING([for GNU hashmap in '<hash_map>' and namespace _SLTP_STD])
        AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[        
#include <hash_map>
using _STLP_STD::hash_map; ]])],
              [ AC_DEFINE(GNU_HASHMAP_NAMESPACE, _STLP_STD, [GNU hashmap namespace]) 
                AC_MSG_RESULT([yes])
                gotgnuhashmap=yes ],
              [AC_MSG_RESULT([no])])
    fi
    AC_LANG_POP([C++])
    if test $gotgnuhashmap = no; then
        AC_MSG_ERROR([Could not find GNU hashmap with any known combination of header+namespace])
    fi

    AC_DEFINE(HAVE_GNU_HASHMAP,[1],[Enable if have GNU hashmap]) 
])




