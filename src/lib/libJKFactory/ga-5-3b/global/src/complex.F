#if HAVE_CONFIG_H
#   include "config.fh"
#endif
c All Fortran functions returning a complex or double complex type must
c be defined here with a corresponding C subroutine e.g. (Fortran) complex
c function ga_cdot --> (C) void gai_cdot(..., SingleComplex *ret).
c 
c from http://www.math.utah.edu/software/c-with-fortran.html
c 
c <BEGIN EXCERPT>
c Function return types
c ...
c Unfortunately, returning composite objects that occupy more than a
c single register, or an adjacent register pair, is fraught with peril.
c Older C and C++ compilers did not support this at all, and newer ones
c may do it differently than Fortran compilers do: thus, you should not
c expect to use Fortran functions that return types such as COMPLEX or
c COMPLEX*16. Write a SUBROUTINE interface to your Fortran function
c instead, and then invoke it as a void function from C or C++. 
c <END EXCERPT>
c 
c We attempted to use various calling conventions, argument order,
c hidden first arguments as the return, etc to eliminate the need
c for any Fortran source such as this, but portability was an issue.
c Intel compilers did not behave similarly to GCC or Pathscale.  At a
c later date perhaps we will revisit this.
c
c For the time being, any function returning a complex type must be
c defined here following the established convention herein.

      complex function ga_cdot(g_a,g_b)
      implicit none
      integer g_a, g_b
      external gai_cdot
      ga_cdot = (0.,0.)
      call gai_cdot(g_a,g_b,ga_cdot)
      end

      complex function nga_cdot(g_a,g_b)
      implicit none
      integer g_a, g_b
      external gai_cdot
      nga_cdot = (0.,0.)
      call ngai_cdot(g_a,g_b,nga_cdot)
      end

      double complex function ga_zdot(g_a,g_b)
      implicit none
      integer g_a, g_b
      external gai_zdot
      ga_zdot = (0.,0.)
      call gai_zdot(g_a,g_b,ga_zdot)
      end

      double complex function nga_zdot(g_a,g_b)
      implicit none
      integer g_a, g_b
      external gai_zdot
      nga_zdot = (0.,0.)
      call ngai_zdot(g_a,g_b,nga_zdot)
      end

      complex function ga_cdot_patch(
     $            g_a, t_a, ailo, aihi, ajlo, ajhi,
     $            g_b, t_b, bilo, bihi, bjlo, bjhi)
      implicit none
      integer g_a, ailo, aihi, ajlo, ajhi
      integer g_b, bilo, bihi, bjlo, bjhi
      character*1 t_a, t_b
      external gai_cdot_patch
      ga_cdot_patch = (0.,0.)
      call gai_cdot_patch(
     $            g_a, t_a, ailo, aihi, ajlo, ajhi,
     $            g_b, t_b, bilo, bihi, bjlo, bjhi, ga_cdot_patch)
      end

      complex function nga_cdot_patch(
     $            g_a, t_a, alo, ahi,
     $            g_b, t_b, blo, bhi)
      implicit none
      integer g_a, alo, ahi 
      integer g_b, blo, bhi
      character*1 t_a, t_b
      external ngai_cdot_patch
      nga_cdot_patch = (0.,0.)
      call ngai_cdot_patch(
     $            g_a, t_a, alo, ahi, 
     $            g_b, t_b, blo, bhi, nga_cdot_patch)
      end

      double complex function ga_zdot_patch(
     $            g_a, t_a, ailo, aihi, ajlo, ajhi,
     $            g_b, t_b, bilo, bihi, bjlo, bjhi)
      implicit none
      integer g_a, ailo, aihi, ajlo, ajhi
      integer g_b, bilo, bihi, bjlo, bjhi
      character*1 t_a, t_b
      external gai_zdot_patch
      ga_zdot_patch = (0.,0.)
      call gai_zdot_patch(
     $            g_a, t_a, ailo, aihi, ajlo, ajhi,
     $            g_b, t_b, bilo, bihi, bjlo, bjhi, ga_zdot_patch)
      end

      double complex function nga_zdot_patch(
     $            g_a, t_a, alo, ahi,
     $            g_b, t_b, blo, bhi)
      implicit none
      integer g_a, alo, ahi
      integer g_b, blo, bhi
      character*1 t_a, t_b
      external ngai_zdot_patch
      nga_zdot_patch = (0.,0.)
      call ngai_zdot_patch(
     $            g_a, t_a, alo, ahi,
     $            g_b, t_b, blo, bhi, nga_zdot_patch)
      end
