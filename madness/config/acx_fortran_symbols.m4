AC_DEFUN([ACX_FORTRAN_SYMBOLS], [
# Dubiously checks for Fortran linking conventions and BLAS+LAPACK at the same time
# mostly to avoid the need for having a fortran compiler installed

# Check for no underscore first since IBM BG ESSL seems to define dgemm with/without underscore
# but dsyev only without underscore ... but AMD ACML also defines both but with different
# interfaces (fortran and c) ... ugh.  Hardwire linking for bgp and restore to original order.

       echo "Checking Fortran-C linking conventions (dgemm -> ?)"
       fsym=no
   
       if test $host = "powerpc-bgp-linux-gnu"; then
          fsym="lc"
          echo "Fortran linking convention is $fsym" 
          AC_DEFINE([FORTRAN_LINKAGE_LC],[1],[Fortran-C linking convention lower case (no underscore)])
       fi
       if test $fsym = no; then
           AC_CHECK_FUNC([dgemm_],[fsym="lcu"])
       fi
       if test $fsym = no; then
           AC_CHECK_FUNC([dgemm],[fsym="lc"])
       fi
       if test $fsym = no; then
           AC_CHECK_FUNC([dgemm__],[fsym="lcuu"])
       fi
       if test $fsym = no; then
           AC_CHECK_FUNC([DGEMM],[fsym="uc"])
       fi
       if test $fsym = no; then
           AC_CHECK_FUNC([DGEMM_],[fsym="ucu"])
       fi
   

    if test "$acx_with_eigen3" != yes; then

       if test $fsym = no; then
           AC_MSG_ERROR([Could not find dgemm with any known linking conventions])
       fi

       echo "Fortran linking convention is $fsym" 

       if test $fsym = lc; then
           AC_DEFINE([FORTRAN_LINKAGE_LC],[1],[Fortran-C linking convention lower case (no underscore)])
           AC_CHECK_FUNC([dsyev],[],AC_MSG_ERROR([Could not find dsyev with selected linking convention]))
       fi
       if test $fsym = lcu; then
           AC_DEFINE([FORTRAN_LINKAGE_LCU],[1],[Fortran-C linking convention lower case with single underscore])
           AC_CHECK_FUNC([dsyev_],[],AC_MSG_ERROR([Could not find dsyev with selected linking convention]))
       fi
       if test $fsym = lcuu; then
           AC_DEFINE([FORTRAN_LINKAGE_LCUU],[1],[Fortran-C linking convention lower case with double underscore])
           AC_CHECK_FUNC([dsyev__],[],AC_MSG_ERROR([Could not find dsyev with selected linking convention]))
       fi
       if test $fsym = uc; then
           AC_DEFINE([FORTRAN_LINKAGE_UC],[1],[Fortran-C linking convention upper case])
           AC_CHECK_FUNC([DSYEV],[],AC_MSG_ERROR([Could not find dsyev with selected linking convention]))
       fi
       if test $fsym = ucu; then
           AC_DEFINE([FORTRAN_LINKAGE_UCU],[1],[Fortran-C linking convention upper case with single underscore])
           AC_CHECK_FUNC([DSYEV_],[],AC_MSG_ERROR([Could not find dsyev with selected linking convention]))
       fi
    fi
])





