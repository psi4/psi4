#if HAVE_CONFIG_H
#   include "config.fh"
#endif
      subroutine vorticity(g_fg,  fg,  ld_fg1,  ld_fg2,
     +                     g_fld, fld, ld_fld1, ld_fld2,
     +                     g_bc,  bc,  ld_bc1, ic, jc)
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
c   Subroutine to evaluate the vorticity of velocity field
c
      double precision drho(2), dux(2), duy(2)
      double precision dr,dx,dy,adxi,adyi,drxi,dryi,dri,drni,dxi,dyi
      double precision rhot,uxt,uyt
      integer ic, jc, i, j, ii, jj
c
      call ga_update_ghosts(g_fld)
      do jj = width(2) + 1, dims(2) - width(2)
        do ii = width(1) + 1, dims(1) - width(1)
          drni = 0.0d00
          drxi = 0.0d00
          dryi = 0.0d00
          do j = 1, 2
            drho(j) = 0.0d00
            dux(j) = 0.0d00
            duy(j) = 0.0d00
          end do
          do j = -1, 1
            dy = dble(j)*delta_x
            if (dy.ne.0.0d00) then
              dyi = 1.0d00/dy
            else
              dyi = 0.0d00
            endif
            do i = -1, 1
              if (bc(ii+i,jj+j).eq.0.and.(i.ne.0.or.j.ne.0)) then
                dx = dble(i)*delta_x
                if (dx.ne.0.0d00) then
                  dxi = 1.0d00/dx
                else
                  dxi = 0.0d00
                endif
                dr = sqrt(dx**2+dy**2)
                dri = 1.0d00/dr
                adxi = abs(dx*dri**2)
                adyi = abs(dy*dri**2)
c
                rhot = fld(ii+i,jj+j,1)
                uxt = fld(ii+i,jj+j,2)
                uyt = fld(ii+i,jj+j,3)
                drho(1) = drho(1) + (rhot-fld(ii,jj,1))*dxi*adxi
                drho(2) = drho(2) + (rhot-fld(ii,jj,1))*dyi*adyi
                dux(1) = dux(1) + (uxt-fld(ii,jj,2))*dxi*adxi
                dux(2) = dux(2) + (uxt-fld(ii,jj,2))*dyi*adyi
                duy(1) = duy(1) + (uyt-fld(ii,jj,3))*dxi*adxi
                duy(2) = duy(2) + (uyt-fld(ii,jj,3))*dyi*adyi
                drxi = drxi+adxi
                dryi = dryi+adyi
              endif
            end do
          end do
          drho(1) = drho(1)/drxi
          drho(2) = drho(2)/dryi
          dux(1) = dux(1)/drxi
          dux(2) = dux(2)/dryi
          duy(1) = duy(1)/drxi
          duy(2) = duy(2)/dryi
          fld(ii,jj,6) = duy(1)-dux(2)
        end do
      end do
      return
      end
