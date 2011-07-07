AC_DEFUN([ACX_DETECT_CXX], [
    # Sets environment variable CXXVENDOR to one of
    # [GNU,Intel,Portland,Pathscale,IBM,unknown]
    AC_CACHE_CHECK([compiler vendor], [acx_cv_detect_cxx], [
      acx_cv_detect_cxx=unknown
      if test $acx_cv_detect_cxx = unknown; then
          $CXX --version 2>&1 | egrep -q "GCC|GNU|gcc|gnu|g\+\+|Free S"
          if test $? = 0; then
             acx_cv_detect_cxx=GNU
          fi
      fi
      if test $acx_cv_detect_cxx = unknown; then
          $CXX --version 2>&1 | egrep -q "clang"
          if test $? = 0; then
             acx_cv_detect_cxx=clang
          fi
      fi
      if test $acx_cv_detect_cxx = unknown; then
          $CXX --version 2>&1 | grep -q "Intel"
          if test $? = 0; then
             acx_cv_detect_cxx=Intel
          fi
      fi
      if test $acx_cv_detect_cxx = unknown; then
          $CXX --version 2>&1 | grep -q "Portland"
          if test $? = 0; then
             acx_cv_detect_cxx=Portland
          fi
      fi
      if test $acx_cv_detect_cxx = unknown; then
          $CXX -v 2>&1 | grep -q "Pathscale"
          if test $? = 0; then
             acx_cv_detect_cxx=Pathscale
          fi
      fi
      if test $acx_cv_detect_cxx = unknown; then
          $CXX -qversion 2>&1 | grep -q "IBM"
          if test $? = 0; then
             acx_cv_detect_cxx=IBM
          fi
      fi
    ])
    
    CXXVENDOR="$acx_cv_detect_cxx"
])




