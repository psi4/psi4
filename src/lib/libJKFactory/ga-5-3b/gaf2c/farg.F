#if HAVE_CONFIG_H
#   include "config.fh"
#endif
      integer function f2c_iargc()
      FXX_MODULE
      f2c_iargc = F77_IARGC() + 1
      return
      end
      subroutine f2c_getarg( i, s )
      FXX_MODULE
      integer       i, l, ier
      character *(*) s
      F77_GETARG_DECL
      call F77_GETARG(F77_GETARG_ARGS)
      return
      end
