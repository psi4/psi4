#if HAVE_CONFIG_H
#   include "config.fh"
#endif
      subroutine equil(g_fg,  fg,  ld_fg1,  ld_fg2,
     +                 g_fld, fld, ld_fld1, ld_fld2,
     +                 g_bc,  bc,  ld_bc1)
#include "common"
c
      integer ld_fg1, ld_fg2, ld_fld1, ld_fld2, ld_bc1
      double precision fg(ld_fg1,ld_fg2, *)
      double precision fld(ld_fld1, ld_fld2, *)
      integer bc(ld_bc1, *)
      integer g_fg, g_fld, g_bc 
c
      double precision rdim, b, c2, cspdi, cspd2, ex, ey
      double precision aa0, aa
      integer i, j, ii, jj
c
c   This subroutine evaluates the equilibrium distribution
c
      rdim = 4.0d00
      b = 24.0d00
      c2 = cspd**2
      cspdi = sqrt(2.0d00)/cspd
      cspd2 = 1.0d00/cspdi
c
c   Calculate equilibrium distributions
c
      do jj = width(2) + 1, dims(2) - width(2)
        do i = 10, 18
          j = i-9
          if (j.ne.1) then
            ex = cspd2*ei(j-1,1)
            ey = cspd2*ei(j-1,2)
          else
            ex = 0.0d00
            ey = 0.0d00
          endif
          do ii = width(1) + 1, dims(1) - width(1)
            aa = (rdim/(b*c2))*fld(ii,jj,4)
            aa0 = fld(ii,jj,1) - b*aa
            if (j.eq.1) then
              ffa(j) = aa0 + 4.0d00*aa
            else if (j.le.5) then
              ffa(j) = 4.0d00*aa
            else
              ffa(j) = aa
            endif
c
c   evaluate distribution
c
            fg(ii,jj,i) = ffa(j) + fld(ii,jj,1)*(ffb(j)*(fld(ii,jj,2)*ex
     +                  + fld(ii,jj,3)*ey)
     +                  + ffc(j)*(fld(ii,jj,2)**2*ex**2
     +                  + 2.0d00*fld(ii,jj,2)*fld(ii,jj,3)*ex*ey
     +                  + fld(ii,jj,3)**2*ey**2)
     +                  + ffd(j)*(fld(ii,jj,2)**2+fld(ii,jj,3)**2))
          end do
        end do
      end do
c
      return
      end
