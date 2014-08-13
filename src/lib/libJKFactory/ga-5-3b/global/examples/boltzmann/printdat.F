#if HAVE_CONFIG_H
#   include "config.fh"
#endif
      subroutine printdat(g_fg,  fg,  ld_fg1,  ld_fg2,
     +                    g_fld, fld, ld_fld1, ld_fld2,
     +                    g_bc,  bc,  ld_bc1)
#include "common"
c
      integer ld_fg1, ld_fg2, ld_fld1, ld_fld2, ld_bc1
      double precision fg(ld_fg1, ld_fg2, *)
      double precision fld(ld_fld1, ld_fld2, *)
      integer bc(ld_bc1, *)
      integer g_fg, g_fld, g_bc 
c
#include "mafdecls.fh"
#include "global.fh"
c
      double precision buffer(NSIZE,10,3)
      double precision dx, dy
      integer glo(NDIM+1), ghi(NDIM+1), bld(NDIM)
      integer i,j,k,imax,icnt1,icnt2,inc1,inc2,jcnt
      integer maxelem
c
c   Subroutine to print current value of fields to a file
c
      imax = size(2)/10 + 1
      glo(1) = 1
      ghi(1) = size(1)
      glo(3) = 1
      ghi(3) = 3
      bld(1) = NSIZE
      bld(2) = 10
c
      maxelem = 256
c
c    Check dimensions to see if size needs to be reduced
c
      if (size(1).gt.maxelem) then
        inc1 = size(1)/maxelem
        icnt1 = 0
        do i = 1, size(1), inc1
          icnt1 = icnt1 + 1
        end do
      else
        inc1 = 1
        icnt1= size(1)
      endif
      if (size(2).gt.maxelem) then
        inc2 = size(2)/maxelem
        icnt2 = 0
        do i = 1, size(2), inc2
          icnt2 = icnt2 + 1
        end do
      else
        inc2 = 1
        icnt2= size(2)
      endif
c
      if (ga_nodeid().eq.0) then
        open(unit=7,file='bltz.gmv',status='unknown')
c
c   Write out header type information
c
        write(7,100)
        write(7,200) icnt1,icnt2,1
        dx = dble(inc1)*xmax/dble(size(1)-1)
        dy = dble(inc1)*xmax/dble(size(1)-1)
        write(7,300) (dble(i-1)*dx,i=1,icnt1)
        write(7,300) (dble(i-1)*dy,i=1,icnt2)
        write(7,300) 0.0
        write(7,1000)
        write(7,400)
        write(7,500)
        jcnt = 0
        do i = 1, imax
          glo(2) = (i-1)*10 + 1
          ghi(2) = i*10
          if (ghi(2).gt.size(2)) ghi(2) = size(2)
          if (ghi(2).ge.glo(2)) then
            call nga_get(g_fld, glo, ghi, buffer, bld)
            do k = 1, ghi(2)-glo(2)+1
              jcnt = jcnt + 1
              if (mod((jcnt-1),inc2).eq.0) then
                write(7,300) (buffer(j,k,1),j=1,size(1),inc1)
              endif
            end do
          endif
        end do
        write(7,600)
        jcnt = 0
        do i = 1, imax
          glo(2) = (i-1)*10 + 1
          ghi(2) = i*10
          if (ghi(2).gt.size(2)) ghi(2) = size(2)
          if (ghi(2).ge.glo(2)) then
            call nga_get(g_fld, glo, ghi, buffer, bld)
            do k = 1, ghi(2)-glo(2)+1
              jcnt = jcnt + 1
              if (mod((jcnt-1),inc2).eq.0) then
                write(7,300) (buffer(j,k,2),j=1,size(1),inc1)
              endif
            end do
          endif
        end do
        write(7,700)
        do i = 1, imax
          glo(2) = (i-1)*10 + 1
          ghi(2) = i*10
          if (ghi(2).gt.size(2)) ghi(2) = size(2)
          if (ghi(2).ge.glo(2)) then
            call nga_get(g_fld, glo, ghi, buffer, bld)
            do k = 1, ghi(2)-glo(2)+1
              jcnt = jcnt + 1
              if (mod((jcnt-1),inc2).eq.0) then
                write(7,300) (buffer(j,k,3),j=1,size(1),inc1)
              endif
            end do
          endif
        end do
        write(7,1100)
        glo(3) = 4
        ghi(3) = 6
        do i = 1, imax
          glo(2) = (i-1)*10 + 1
          ghi(2) = i*10
          if (ghi(2).gt.size(2)) ghi(2) = size(2)
          if (ghi(2).ge.glo(2)) then
            call nga_get(g_fld, glo, ghi, buffer, bld)
            do k = 1, ghi(2)-glo(2)+1
              jcnt = jcnt + 1
              if (mod((jcnt-1),inc2).eq.0) then
                write(7,300) (buffer(j,k,3),j=1,size(1),inc1)
              endif
            end do
          endif
        end do
        write(7,800)
        write(7,900)
        close(7)
      endif
c
  100 format('gmvinput ascii')
  200 format('nodes       -1',3i10)
  300 format('    ',5f12.4)
  400 format('variable')
  500 format('rho                           1')
  600 format('ux                            1')
  700 format('uy                            1')
  800 format('endvars')
  900 format('endgmv')
 1000 format('cells     0')
 1100 format('vorticity                     1')
      return
      end
