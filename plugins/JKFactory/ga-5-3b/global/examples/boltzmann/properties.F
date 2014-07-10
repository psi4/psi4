#if HAVE_CONFIG_H
#   include "config.fh"
#endif
      subroutine properties(g_fg,  fg,  ld_fg1,  ld_fg2,
     +                      g_fld, fld, ld_fld1, ld_fld2,
     +                      g_bc,  bc,  ld_bc1)
#include "common"
c
      integer ld_fg1, ld_fg2, ld_fld1, ld_fld2, ld_bc1
      double precision fg(ld_fg1,ld_fg2, *)
      double precision fld(ld_fld1, ld_fld2, *)
      integer bc(ld_bc1, *)
      integer g_fg, g_fld, g_bc 
c
#include "mafdecls.fh"
#include "global.fh"
c
c  Evaluate properties on grid
c
      double precision cspd2, ex, ey, sum(3)
      double precision rho
      integer i, j, ii, jj
c
c  Zero accumulators and hydrodynamic fields
c
      rtot = 0.0d00
      uxtot = 0.0d00
      uytot = 0.0d00
      cspd2 = cspd/sqrt(2.0d00)
c
      do jj = width(2) + 1, dims(2) - width(2)
        do i = 1, 4
          do ii = width(1) + 1, dims(1) - width(1)
            fld(ii,jj,i) = 0.0d00
          end do
        end do
      end do
c
c   Evaluate density and momentum
c
      do jj = width(2) + 1, dims(2) - width(2)
        do i = 1, 9
          j = i - 1
          if (j.ne.0) then
            ex = cspd2*ei(j,1)
            ey = cspd2*ei(j,2)
          else
            ex = 0.0d00
            ey = 0.0d00
          endif
          do ii = width(1) + 1, dims(1) - width(1)
            fld(ii,jj,1) = fld(ii,jj,1) + fg(ii,jj,i)
            fld(ii,jj,2) = fld(ii,jj,2) + ex*fg(ii,jj,i)
            fld(ii,jj,3) = fld(ii,jj,3) + ey*fg(ii,jj,i)
          end do
        end do
      end do
c
c   Evaluate pressure and accumulate values
c
      do jj = width(2) + 1, dims(2) - width(2)
        do ii = width(1) + 1, dims(1) - width(1)
          fld(ii,jj,2) = fld(ii,jj,2)/fld(ii,jj,1)
          fld(ii,jj,3) = fld(ii,jj,3)/fld(ii,jj,1)
          rtot = rtot + fld(ii,jj,1)
c          if (bc(ii,jj).eq.1) then
c            fld(ii,jj,2) = 0.0d00
c            fld(ii,jj,3) = 0.0d00
c          elseif (bc(ii,jj).eq.2) then
c            fld(ii,jj,1) = rhobc
c            fld(ii,jj,2) = uxbc
c            fld(ii,jj,3) = 0.0d00
c          endif
          uxtot = uxtot + fld(ii,jj,1)*fld(ii,jj,2)
          uytot = uytot + fld(ii,jj,1)*fld(ii,jj,3)
          rho = fld(ii,jj,1)
          fld(ii,jj,4) = rho*rgas*tmprtr0/(1.0d00-b_vdw*rho)
     +                 - a_vdw*rho**2
        end do
      end do
c
c  Sum values across processors
c
      sum(1) = rtot
      sum(2) = uxtot
      sum(3) = uytot
c
      call ga_dgop(MT_DBL,sum,3,'+')
      rtot = sum(1)
      uxtot = sum(2)
      uytot = sum(3)
c
      return
      end
