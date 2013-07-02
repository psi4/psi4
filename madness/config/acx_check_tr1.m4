AC_DEFUN([ACX_CHECK_SHARED_PTR], [
  AC_LANG_SAVE
  AC_LANG([C++])
  
  # Check for shared_ptr in std namespace
  AC_MSG_CHECKING([for shared_ptr])
  acx_shared_ptr=no
  
  # Check for std::shared_ptr in <memory>
  AC_COMPILE_IFELSE(
    [
      AC_LANG_PROGRAM(
        [[#include <memory>]],
        [[std::shared_ptr<int> p;]]
      )
    ],
    [
      AC_DEFINE([MADNESS_USE_MEMORY],[1],[define if MADNESS is using <memory>.])
      AC_DEFINE([MADNESS_HAS_STD_SHARED_PTR],[1],[define if std::shared_ptr is available.])
      acx_shared_ptr=yes
    ]
  )
  
  # Check for std::tr1::shared_ptr in <memory>
  if test "$acx_shared_ptr" = no; then
    AC_COMPILE_IFELSE(
      [
        AC_LANG_PROGRAM(
          [[#include <memory>]],
          [[std::tr1::shared_ptr<int> p;]]
        )
      ],
      [
        AC_DEFINE([MADNESS_USE_MEMORY],[1],[define if MADNESS is using <memory>.])
        AC_DEFINE([MADNESS_HAS_STD_TR1_SHARED_PTR],[1],[define if std::tr1::shared_ptr is available.])
        acx_shared_ptr=yes
      ]
    )
  fi
  
  # Check for std::shared_ptr in <tr1/memory>
  if test "$acx_shared_ptr" = no; then
    AC_COMPILE_IFELSE(
      [
        AC_LANG_PROGRAM(
          [[#include <tr1/memory>]],
          [[std::shared_ptr<int> p;]]
        )
      ],
      [
        AC_DEFINE([MADNESS_USE_TR1_MEMORY],[1],[define if MADNESS is using <tr1/memory>.])
        AC_DEFINE([MADNESS_HAS_STD_SHARED_PTR],[1],[define if std::tr1::shared_ptr is available.])
        acx_shared_ptr=yes
      ]
    )
  fi
  
  # Check for std::tr1::shared_ptr in <tr1/memory>
  if test "$acx_shared_ptr" = no; then
    AC_COMPILE_IFELSE(
      [
        AC_LANG_PROGRAM(
          [[#include <tr1/memory>]],
          [[std::tr1::shared_ptr<int> p;]]
        )
      ],
      [
        AC_DEFINE([MADNESS_USE_TR1_MEMORY],[1],[define if MADNESS is using <tr1/memory>.])
        AC_DEFINE([MADNESS_HAS_STD_TR1_SHARED_PTR],[1],[define if std::tr1::shared_ptr is available.])
        acx_shared_ptr=yes
      ]
    )
  fi
  
  # Check if we should use boost tr1 memory
  if test "$acx_with_boost$acx_shared_ptr" = yesno; then
    AC_DEFINE([MADNESS_USE_BOOST_TR1_MEMORY_HPP],[1],[define if MADNESS is using <boost/tr1/memory.hpp>.])
    AC_DEFINE([MADNESS_HAS_STD_TR1_SHARED_PTR],[1],[define if std::tr1::shared_ptr is available.])
    acx_shared_ptr=yes
  fi
  
  # post shared_ptr results
  AC_MSG_RESULT([$acx_shared_ptr])
  if test "$acx_shared_ptr" = no; then
    AC_MSG_ERROR([std::shared_ptr and std::tr1::shared_ptr are not supported. Reconfigure MADNESS to with the --with-boost option.])
  fi
  
  #Check for std::make_shared and std::allocate_shared
  acx_std_make_shared=no
  AC_MSG_CHECKING([for std::make_shared and std::allocate_shared])
  
  AC_COMPILE_IFELSE(
    [
      AC_LANG_PROGRAM(
        [[#include <memory>]],
        [[using ::std::make_shared; using ::std::allocate_shared;]]
      )
    ],
    [
      AC_DEFINE([MADNESS_HAS_STD_MAKE_SHARED],[1],
        [define if std::make_shared and std::allocate_shared are available.])
      acx_std_make_shared=yes
    ]
  )
    
  # post make_shared results
  AC_MSG_RESULT([$acx_std_make_shared])
  
  AC_LANG_RESTORE
])

AC_DEFUN([ACX_CHECK_TYPE_TRAITS], [
  AC_LANG_SAVE
  AC_LANG([C++])

  # Check for type traits in <type_traits> and std namespace
  AC_MSG_CHECKING([for type traits])
  acx_type_traits=no

  AC_COMPILE_IFELSE(
    [
      AC_LANG_PROGRAM(
        [[#include <type_traits>]],
        [[typedef std::is_same<int, double> sameT;]]
      )
    ],
    [
      AC_DEFINE([MADNESS_USE_TYPE_TRAITS],[1],[define if MADNESS is using <type_traits>.])
      AC_DEFINE([MADNESS_HAS_STD_TYPE_TRAITS],[1],[define if std type traits are available.])
      acx_type_traits=yes
    ]
  )

  if test "$acx_type_traits" = no; then
    AC_COMPILE_IFELSE(
      [
        AC_LANG_PROGRAM(
          [[#include <type_traits>]],
          [[typedef std::tr1::is_same<int, double> sameT;]]
        )
      ],
      [
        AC_DEFINE([MADNESS_USE_TYPE_TRAITS],[1],[define if MADNESS is using <type_traits>.])
        AC_DEFINE([MADNESS_HAS_STD_TR1_TYPE_TRAITS],[1],[define if std::tr1 type traits are available.])
        acx_type_traits=yes
      ]
    )
  fi
  
  if test "$acx_type_traits" = no; then
    AC_COMPILE_IFELSE(
      [
        AC_LANG_PROGRAM(
          [[#include <tr1/type_traits>]],
          [[typedef std::is_same<int, double> sameT;]]
        )
      ],
      [
        AC_DEFINE([MADNESS_USE_TR1_TYPE_TRAITS],[1],[define if MADNESS is using <tr1/type_traits>.])
        AC_DEFINE([MADNESS_HAS_STD_TYPE_TRAITS],[1],[define if std::tr1 type traits are available.])
        acx_type_traits=yes
      ]
    )
  fi
  
  if test "$acx_type_traits" = no; then
    AC_COMPILE_IFELSE(
      [
        AC_LANG_PROGRAM(
          [[#include <tr1/type_traits>]],
          [[typedef std::tr1::is_same<int, double> sameT;]]
        )
      ],
      [
        AC_DEFINE([MADNESS_USE_TR1_TYPE_TRAITS],[1],[define if MADNESS is using <tr1/type_traits>.])
        AC_DEFINE([MADNESS_HAS_STD_TR1_TYPE_TRAITS],[1],[define if std::tr1 type traits are available.])
        acx_type_traits=yes
      ]
    )
  fi
  
  # Check if we should use boost tr1 type_traits
  if test "$acx_with_boost$acx_type_traits" = yesno; then
    AC_DEFINE([MADNESS_USE_BOOST_TR1_TYPE_TRAITS_HPP],[1],[define if MADNESS is using <boost/tr1/type_traits.hpp>.])
    AC_DEFINE([MADNESS_HAS_STD_TR1_TYPE_TRAITS],[1],[define if std::tr1 type traits are available.])
    acx_type_traits=yes
  fi
  
  # post type traits results
  AC_MSG_RESULT([$acx_type_traits])
  if test "$acx_type_traits" = no; then
    AC_MSG_ERROR([std and std::tr1 type traits are not supported. Reconfigure MADNESS with the --with-boost option.])
  fi
  
  AC_LANG_RESTORE
])

AC_DEFUN([ACX_CHECK_ARRAY], [
  AC_LANG_SAVE
  AC_LANG([C++])
  
  # Check for array in std namespace
  AC_MSG_CHECKING([for array])
  acx_array=no
  
  # Check for std::array in <array>
  AC_COMPILE_IFELSE(
    [
      AC_LANG_PROGRAM(
        [[#include <array>]],
        [[std::array<int,10> a;]]
      )
    ],
    [
      AC_DEFINE([MADNESS_USE_ARRAY],[1],[define if MADNESS is using <array>.])
      AC_DEFINE([MADNESS_HAS_STD_ARRAY],[1],[define if std::array is available.])
      AC_DEFINE([MADNESS_ARRAY_HAS_FILL],[1],[define if array has fill member function.])
      acx_array=yes
    ]
  )
  
  # Check for std::tr1::array in <array>
  if test "$acx_array" = no; then
    AC_COMPILE_IFELSE(
      [
        AC_LANG_PROGRAM(
          [[#include <array>]],
          [[std::tr1::array<int,10> a;]]
        )
      ],
      [
        AC_DEFINE([MADNESS_USE_ARRAY],[1],[define if MADNESS is using <array>.])
        AC_DEFINE([MADNESS_HAS_STD_TR1_ARRAY],[1],[define if std::tr1::array is available.])

        # Check to see if array has fill function
        AC_COMPILE_IFELSE(
          [
            AC_LANG_PROGRAM(
              [[#include <array>]],
              [[std::tr1::array<int,10> a; a.fill(0);]]
            )
          ],
          [AC_DEFINE([MADNESS_ARRAY_HAS_FILL],[1],[define if array has fill member function.])]
        )
        acx_array=yes
      ]
    )
  fi
  
  # Check for std::array in <tr1/array>
  if test "$acx_array" = no; then
    AC_COMPILE_IFELSE(
      [
        AC_LANG_PROGRAM(
          [[#include <tr1/array>]],
          [[std::array<int,10> a;]]
        )
      ],
      [
        AC_DEFINE([MADNESS_USE_TR1_ARRAY],[1],[define if MADNESS is using <tr1/array>.])
        AC_DEFINE([MADNESS_HAS_STD_ARRAY],[1],[define if std::tr1::array is available.])
        AC_DEFINE([MADNESS_ARRAY_HAS_FILL],[1],[define if array has fill member function.])
        acx_array=yes
      ]
    )
  fi
  
  # Check for std::tr1::array in <tr1/array>
  if test "$acx_array" = no; then
    AC_COMPILE_IFELSE(
      [
        AC_LANG_PROGRAM(
          [[#include <tr1/array>]],
          [[std::tr1::array<int,10> a;]]
        )
      ],
      [
        AC_DEFINE([MADNESS_USE_TR1_ARRAY],[1],[define if MADNESS is using <tr1/array>.])
        AC_DEFINE([MADNESS_HAS_STD_TR1_ARRAY],[1],[define if std::tr1::array is available.])        
        # Check to see if array has fill function
        AC_COMPILE_IFELSE(
          [
            AC_LANG_PROGRAM(
              [[#include <tr1/array>]],
              [[std::tr1::array<int,10> a; a.fill(0);]]
            )
          ],
          [AC_DEFINE([MADNESS_ARRAY_HAS_FILL],[1],[define if array has fill member function.])]
        )
        acx_array=yes
      ]
    )
  fi
  
  # Check if we should use boost tr1 array
  if test "$acx_with_boost$acx_array" = yesno; then
    AC_DEFINE([MADNESS_USE_BOOST_TR1_ARRAY_HPP],[1],[define if MADNESS is using <boost/tr1/array.hpp>.])
    AC_DEFINE([MADNESS_HAS_STD_TR1_ARRAY],[1],[define if std::tr1::array is available.])
    # Check to see if array has fill function
    AC_COMPILE_IFELSE(
      [
        AC_LANG_PROGRAM(
          [[#include <tr1/array>]],
          [[std::tr1::array<int,10> a; a.fill(0);]]
        )
      ],
      [AC_DEFINE([MADNESS_ARRAY_HAS_FILL],[1],[define if array has fill member function.])]
    )
    acx_array=yes
  fi
  
  # post array results
  AC_MSG_RESULT([$acx_array])
  if test "$acx_array" = no; then
    AC_MSG_ERROR([std::array and std::tr1::array are not supported. Reconfigure MADNESS to use Boost with the --with-boost option.])
  fi

  AC_LANG_RESTORE
])

AC_DEFUN([ACX_CHECK_HASH], [
  AC_LANG_SAVE
  AC_LANG([C++])
  
  # Check for hash in std namespace
  AC_MSG_CHECKING([for hash])
  acx_hash=no
  
  # Check for std::hash in <functional>
  AC_COMPILE_IFELSE(
    [
      AC_LANG_PROGRAM(
        [[#include <functional>]],
        [[std::hash<int> h; h(1);]]
      )
    ],
    [
      AC_DEFINE([MADNESS_USE_FUNCTIONAL],[1],[define if MADNESS is using <functional>.])
      AC_DEFINE([MADNESS_HAS_STD_HASH],[1],[define if std::hash is available.])
      acx_hash=yes
    ]
  )
  
  # Check for std::tr1::hash in <functional>
  if test "$acx_hash" = no; then
    AC_COMPILE_IFELSE(
      [
        AC_LANG_PROGRAM(
          [[#include <functional>]],
          [[std::tr1::hash<int> h; h(1);]]
        )
      ],
      [
        AC_DEFINE([MADNESS_USE_FUNCTIONAL],[1],[define if MADNESS is using <functional>.])
        AC_DEFINE([MADNESS_HAS_STD_TR1_HASH],[1],[define if std::tr1::hash is available.])
        acx_hash=yes
      ]
    )
  fi
  
  # Check for std::hash in <tr1/functional>
  if test "$acx_hash" = no; then
    AC_COMPILE_IFELSE(
      [
        AC_LANG_PROGRAM(
          [[#include <tr1/functional>]],
          [[std::hash<int> h; h(1);]]
        )
      ],
      [
        AC_DEFINE([MADNESS_USE_TR1_FUNCTIONAL],[1],[define if MADNESS is using <tr1/functional>.])
        AC_DEFINE([MADNESS_HAS_STD_HASH],[1],[define if std::tr1::hash is available.])
        acx_hash=yes
      ]
    )
  fi
  
  # Check for std::tr1::hash in <tr1/functional>
  if test "$acx_hash" = no; then
    AC_COMPILE_IFELSE(
      [
        AC_LANG_PROGRAM(
          [[#include <tr1/functional>]],
          [[std::tr1::hash<int> h; h(1);]]
        )
      ],
      [
        AC_DEFINE([MADNESS_USE_TR1_FUNCTIONAL],[1],[define if MADNESS is using <tr1/functional>.])
        AC_DEFINE([MADNESS_HAS_STD_TR1_HASH],[1],[define if std::tr1::hash is available.])        
        acx_hash=yes
      ]
    )
  fi
  
  # Check if we should use boost tr1 hash
  if test "$acx_with_boost$acx_hash" = yesno; then
    AC_DEFINE([MADNESS_USE_BOOST_TR1_FUNCTIONAL_HPP],[1],[define if MADNESS is using <boost/tr1/functional.hpp>.])
    AC_DEFINE([MADNESS_HAS_STD_TR1_HASH],[1],[define if std::tr1::hash is available.])
    acx_hash=yes
  fi
  
  # post hash results
  AC_MSG_RESULT([$acx_hash])
  if test "$acx_hash" = no; then
    AC_MSG_ERROR([std::hash and std::tr1::hash are not supported. Reconfigure MADNESS with the --with-boost option.])
  fi

  AC_LANG_RESTORE
])

AC_DEFUN([ACX_CHECK_RESULT_OF], [
  AC_LANG_SAVE
  AC_LANG([C++])

  # Check for results_of in <functional> and std namespace
  AC_MSG_CHECKING([for result_of])
  acx_result_of=no

  AC_COMPILE_IFELSE(
    [
      AC_LANG_PROGRAM(
        [[#include <functional>]],
        [[using std::result_of;]]
      )
    ],
    [
      AC_DEFINE([MADNESS_USE_FUNCTIONAL],[1],[define if MADNESS is using <functional>.])
      AC_DEFINE([MADNESS_HAS_STD_RESULT_OF],[1],[define if std::result_of are available.])
      acx_result_of=yes
    ]
  )

  if test "$acx_result_of" = no; then
    AC_COMPILE_IFELSE(
      [
        AC_LANG_PROGRAM(
          [[#include <functional>]],
          [[using std::tr1::result_of;]]
        )
      ],
      [
        AC_DEFINE([MADNESS_USE_FUNCTIONAL],[1],[define if MADNESS is using <functional>.])
        AC_DEFINE([MADNESS_HAS_STD_TR1_RESULT_OF],[1],[define if std::tr1::result_of are available.])
        acx_result_of=yes
      ]
    )
  fi
  
  if test "$acx_result_of" = no; then
    AC_COMPILE_IFELSE(
      [
        AC_LANG_PROGRAM(
          [[#include <tr1/functional>]],
          [[using std::result_of;]]
        )
      ],
      [
        AC_DEFINE([MADNESS_USE_TR1_FUNCTIONAL],[1],[define if MADNESS is using <tr1/functional>.])
        AC_DEFINE([MADNESS_HAS_STD_RESULT_OF],[1],[define if std::result_of are available.])
        acx_result_of=yes
      ]
    )
  fi
  
  if test "$acx_result_of" = no; then
    AC_COMPILE_IFELSE(
      [
        AC_LANG_PROGRAM(
          [[#include <tr1/functional>]],
          [[using std::tr1::result_of;]]
        )
      ],
      [
        AC_DEFINE([MADNESS_USE_TR1_FUNCTIONAL],[1],[define if MADNESS is using <tr1/functional>.])
        AC_DEFINE([MADNESS_HAS_STD_TR1_RESULT_OF],[1],[define if std::tr1::result_of are available.])
        acx_result_of=yes
      ]
    )
  fi
  
  # Check if we should use boost tr1 functional
  if test "$acx_with_boost$acx_result_of" = yesno; then
    AC_DEFINE([MADNESS_USE_BOOST_TR1_FUNCTIONAL_HPP],[1],[define if MADNESS is using <boost/tr1/functional.hpp>.])
    AC_DEFINE([MADNESS_HAS_STD_TR1_RESULT_OF],[1],[define if std::tr1 restult_of are available.])
    acx_shared_ptr=yes
  fi
  
  # post result_of results
  AC_MSG_RESULT([$acx_result_of])
  if test "$acx_result_of" = no; then
    AC_MSG_ERROR([std and std::tr1 result_of are not supported. Reconfigure MADNESS to use Boost with the --with-boost option.])
  fi
  
  AC_LANG_RESTORE
])


AC_DEFUN([ACX_CHECK_TR1],
[
  ACX_CHECK_SHARED_PTR
  ACX_CHECK_TYPE_TRAITS
  ACX_CHECK_ARRAY
  ACX_CHECK_HASH
#  ACX_CHECK_RESULT_OF
])
