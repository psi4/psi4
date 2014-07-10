#if HAVE_CONFIG_H
#   include "config.fh"
#endif
      subroutine initpar
#include "common"
c
      double precision q8_ei(18,3)
      data q8_ei
     + /1,-1, 0, 0, 1,-1, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     +  0, 0, 1,-1, 1,-1,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     +  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
      integer i,j
      double precision rdim,b,c2
c
      do i=1, 8
        do j = 1, NDIM
          ei(i,j) = q8_ei(i,j)
        end do
      end do
c
      rdim = 4.0d00
      c2 = cspd**2
      b = 24.0d00
c
      do i =1, 9
        if (i.eq.1) then
          ffb(i) = 0.0d00
          ffc(i) = 0.0d00
          ffd(i) = -1.0d00/c2-4.0d00*rdim/(2.0d00*b*c2)
        elseif (i.le.5) then
          ffb(i) = 4.0d00*rdim/(b*c2)
          ffc(i) = 4.0d00*rdim*(rdim+2.0d00)/(2.0d00*b*c2**2)
          ffd(i) = -4.0d00*rdim/(2.0d00*b*c2)
        else
          ffb(i) = rdim/(b*c2)
          ffc(i) = rdim*(rdim+2.0d00)/(2.0d00*b*c2**2)
          ffd(i) = -rdim/(2.0d00*b*c2)
        endif
      end do
c
c
c   Initialize hash table
c
      hash(0,0)   = 1
      hash(1,0)   = 2
      hash(-1,0)  = 3
      hash(0,1)   = 4
      hash(0,-1)  = 5
      hash(1,1)   = 6
      hash(-1,-1) = 7
      hash(1,-1)  = 8
      hash(-1,1)  = 9
c
c   Initialise inverse hash table
c
      ihash(0,0)   = 1
      ihash(1,0)   = 3
      ihash(-1,0)  = 2
      ihash(0,1)   = 5
      ihash(0,-1)  = 4
      ihash(1,1)   = 7
      ihash(-1,-1) = 6
      ihash(1,-1)  = 9
      ihash(-1,1)  = 8
c
      do i = 1, 50
        tstats(i) = 0.0d00
      end do
c
      return
      end
