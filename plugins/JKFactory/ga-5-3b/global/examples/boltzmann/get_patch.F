#if HAVE_CONFIG_H
#   include "config.fh"
#endif
      subroutine get_patch(g_fg,  fg,  ld_fg1,  ld_fg2,
     +                     g_fld, fld, ld_fld1, ld_fld2,
     +                     g_bc,  bc,  ld_bc1, ii, jj)
#include "common"
c
      integer ld_fg1, ld_fg2, ld_fld1, ld_fld2, ld_bc1
      double precision fg(ld_fg1,ld_fg2, *)
      double precision fld(ld_fld1, ld_fld2, *)
      integer bc(ld_bc1, *)
      integer g_fg, g_fld, g_bc 
c
c   Subroutine to handle cells at boundary
c
      double precision rho, ux, uy, cspd2, dux, duy, a, b
      integer i, j, ii, jj, k, l, kk, ll
      integer bval
      integer mask(-1:1,-1:1)
c
c   Locate local and absolute indices of cell
c
      i = ii - width(1) - 1 + lo(1)
      j = jj - width(2) - 1 + lo(2)
c
c   Check values of neighboring cells
c
      do l = -1, 1
        do k = -1, 1
          if (i+k.ge.1.and.i+k.le.size(1).and.
     +        j+l.ge.1.and.j+l.le.size(2)) then
            if (bc(ii+k,jj+l).ne.0) then
              mask(k,l) = 2
            else
              mask(k,l) = 0
            endif
          else
            mask(k,l) = 2
          endif
        end do
      end do
c
c   Determine if cells in mask represent interior cells or are
c   on the boundary
c
      do l = -1, 1
        do k = -1, 1
          if (mask(k,l).ne.0.and.(k.ne.0.or.l.ne.0)) then
            do ll = max(l-1,-1), min(l+1,1)
              do kk = max(k-1,-1), min(k+1,1)
                if (mask(kk,ll).eq.0) mask(k,l) = 1
              end do
            end do
          endif
        end do
      end do
c
c  Evaluate distribution in boundary patch
c
      bval = bc(ii,jj)
c
c  Apply simple bounce back condition
c
      if (bval.eq.1) then
        do l = -1, 1
          do k = -1, 1
            if (mask(k,l).eq.2) then
              fgp(k,l,ihash(k,l)) = fg(ii,jj,hash(k,l)+18)
            else
              fgp(k,l,ihash(k,l)) = fg(ii+k,jj+l,ihash(k,l)+18)
            endif
          end do
        end do
c
c   Apply constant velocity boundary condition to flat upper boundary
c
      else if (bval.eq.2) then
        cspd2 = cspd/sqrt(2.0d00)
        rho = 0.0d00
        ux = 0.0d00
        uy = 0.0d00
        do l = -1, 1
          do k = -1, 1
            if (mask(k,l).eq.2) then
              fgp(k,l,ihash(k,l)) = fg(ii,jj,hash(k,l)+18)
            else
              fgp(k,l,ihash(k,l)) = fg(ii+k,jj+l,ihash(k,l)+18)
            endif
            rho = rho + fgp(k,l,ihash(k,l))
            ux = ux - cspd2*dble(k)*fgp(k,l,ihash(k,l))
            uy = uy - cspd2*dble(l)*fgp(k,l,ihash(k,l))
          end do
        end do
        ux = ux/rho
        uy = uy/rho
        ux = uxbc - ux
c
c   Add corrections needed to adjust for velocity mismatch
c
        fgp(1,1,ihash(1,1)) = fgp(1,1,ihash(1,1))
     +                      - 0.5d00*rho*ux/cspd2
        fgp(-1,1,ihash(-1,1)) = fgp(-1,1,ihash(-1,1))
     +                        + 0.5d00*rho*ux/cspd2
cc        fgp(1,1,ihash(1,1)) = fgp(1,1,ihash(1,1))
cc     +                      - 0.5d00*rho*ux/cspd2
cc     +                      - rho*uy/(3.0d00*cspd2)
cc        fgp(0,1,ihash(0,1)) = fgp(0,1,ihash(0,1))
cc     +                      - rho*uy/(3.0d00*cspd2)
cc        fgp(-1,1,ihash(-1,1)) = fgp(-1,1,ihash(-1,1))
cc     +                        + 0.5d00*rho*ux/cspd2
cc     +                        - rho*uy/(3.0d00*cspd2)
c      ux = uxbc
c      uy = 0.0d00
c      rho = (cspd2/(uy+cspd2))*(fg(ii,jj,ihash(0,0)+18)
c     +    + fg(ii+1,jj,ihash(1,0)+18)+fg(ii-1,jj,ihash(-1,0)+18)
c     +    + 2.0d00*(fg(ii-1,jj-1,ihash(-1,-1)+18)
c     +    + fg(ii,jj-1,ihash(0,-1)+18)+fg(ii+1,jj-1,ihash(1,-1)+18)))
c      fgp(0,1,ihash(0,1)) = fg(ii,jj-1,ihash(0,-1)+18)
c      a = (rho*ux-cspd2*(-fg(ii-1,jj,ihash(1,0)+18)
c     +  + fg(ii+1,jj,ihash(1,0)+18)+fg(ii+1,jj+1,ihash(1,1)+18)
c     +  - fg(ii-1,jj+1,ihash(-1,1)+18)))/cspd2
c      b = fg(ii+1,jj+1,ihash(1,1)+18)+fg(ii+1,jj-1,ihash(1,-1)+18)
c      fgp(1,1,ihash(1,1)) = -0.5d00*(a+b)
c      fgp(1,1,ihash(1,1)) = 0.5d00*(a-b)
c
c   Apply constant velocity boundary condition
c
      else if (bval.eq.3) then
        cspd2 = sqrt(2.0d00)*cspd
        rho = 0.0d00
        ux = 0.0d00
        uy = 0.0d00
        do l = -1, 1
          do k = -1, 1
            if (mask(k,l).eq.2.and.(k.ne.0.or.l.ne.0)) then
              fgp(k,l,ihash(k,l)) = fg(ii,jj,hash(k,l))
              rho = rho + fg(ii,jj,hash(k,l))
              ux = ux + cspd2*dble(k)*fg(ii,jj,hash(k,l))
              uy = uy + cspd2*dble(l)*fg(ii,jj,hash(k,l))
            else
              fgp(k,l,ihash(k,l)) = fg(ii+k,jj+l,ihash(i,l))
              rho = rho + fg(ii+k,jj+l,ihash(k,l))
              ux = ux + cspd2*dble(k)*fg(ii+k,jj+l,ihash(k,l))
              uy = uy + cspd2*dble(l)*fg(ii+k,jj+l,ihash(k,l))
            endif
          end do
        end do
c
c   Determine value of correction needed to get specified final
c   velocity
c
        dux = bcpar(2,1) - ux
        duy = bcpar(2,2) - uy
        
      endif
c
      return
      end
