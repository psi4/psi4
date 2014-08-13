#if HAVE_CONFIG_H
#   include "config.fh"
#endif
      subroutine heapsort(iflg)
#include "common.fh"
c
c
c  This routine sorts either the locally owned particles (iflg=0) or
c  the particles stored in the buffer arrays (iflg=1) so that they are
c  arranged in order of increasing index
c
      double precision rr1(3,6),rmss
      double precision tbeg,wraptime
c     
      integer l, ir, i,j
      integer iflg,jj,kk,iat,idx
c     
      if (iflg.eq.1) go to 1000
      if ( antot .eq. 0 ) return
      if ( antot .eq. 1 ) return
c
      tbeg = wraptime()
c     
      l = antot/2 + 1
      ir = antot
c     
 10   continue
      if ( l .gt. 1 ) then
        l = l-1
        do jj = 1, 3
          do kk = 1, 6
            rr1(jj,kk) = ra(l,jj,kk)
          end do
        end do
        rmss = mass(l)
        iat = at(l)
        idx = aidx(l)
      else
        do jj = 1, 3
          do kk = 1, 6
            rr1(jj,kk) = ra(ir,jj,kk)
          end do
        end do
        rmss = mass(ir)
        iat = at(ir)
        idx = aidx(ir)
c
        do jj = 1, 3
          do kk = 1, 6
            ra(ir,jj,kk) = ra(1,jj,kk)
          end do
        end do
        mass(ir) = mass(1)
        at(ir) = at(1)
        aidx(ir) = aidx(1)
        ir = ir - 1
        if ( ir .eq. 1 ) then
          do jj = 1, 3
            do kk = 1, 6
              ra(1,jj,kk) = rr1(jj,kk)
            end do
          end do
          mass(1) = rmss
          at(1) = iat
          aidx(1) = idx
          tmstat(9) = tmstat(9) + wraptime() - tbeg
          return
        endif
      endif
c     
      if ( ir .eq. 1 ) then
        tmstat(9) = tmstat(9) + wraptime() - tbeg
        return
      endif
c     
      i = l
      j = l + l
 20   continue
      if ( j .le. ir ) then
        if ( j .lt. ir ) then
          if ( aidx(j) .lt. aidx(j+1)) j = j+1
        endif
        if ( idx .lt. aidx(j)) then
          do jj = 1, 3
            do kk = 1, 6
              ra(i,jj,kk) = ra(j,jj,kk)
            end do
          end do
          mass(i) = mass(j)
          at(i) = at(j)
          aidx(i) = aidx(j)
          i = j
          j = j + j
        else
          j = ir + 1
        endif
        go to 20
      endif
      do jj = 1, 3
        do kk = 1, 6
          ra(i,jj,kk) = rr1(jj,kk)
        end do
      end do
      mass(i) = rmss
      at(i) = iat
      aidx(i) = idx
      goto 10
c
 1000 if ( btot .eq. 0 ) return
      if ( btot .eq. 1 ) return
c     
      tbeg = wraptime()
c
      l = btot/2 + 1
      ir = btot
c     
 30   continue
      if ( l .gt. 1 ) then
        l = l-1
        rr1(1,1) = xcrd(l)
        rr1(2,1) = ycrd(l)
        rr1(3,1) = zcrd(l)
        rr1(1,2) = xfrc(l)
        rr1(2,2) = yfrc(l)
        rr1(3,2) = zfrc(l)
        rr1(3,3) = zacc(l)
        rmss = mbuf(l)
        iat = bat(l)
        idx = bidx(l)
      else
        rr1(1,1) = xcrd(ir)
        rr1(2,1) = ycrd(ir)
        rr1(3,1) = zcrd(ir)
        rr1(1,2) = xfrc(ir)
        rr1(2,2) = yfrc(ir)
        rr1(3,2) = zfrc(ir)
        rr1(3,3) = zacc(ir)
        rmss = mbuf(ir)
        iat = bat(ir)
        idx = bidx(ir)
c
        xcrd(ir) = xcrd(1)
        ycrd(ir) = ycrd(1)
        zcrd(ir) = zcrd(1)
        xfrc(ir) = xfrc(1)
        yfrc(ir) = yfrc(1)
        zfrc(ir) = zfrc(1)
        zacc(ir) = zacc(1)
        mbuf(ir) = mbuf(1)
        bat(ir) = bat(1)
        bidx(ir) = bidx(1)
        ir = ir - 1
        if ( ir .eq. 1 ) then
          xcrd(1) = rr1(1,1)
          ycrd(1) = rr1(2,1)
          zcrd(1) = rr1(3,1)
          xfrc(1) = rr1(1,2)
          yfrc(1) = rr1(2,2)
          zfrc(1) = rr1(3,2)
          zacc(1) = rr1(3,3)
          mbuf(1) = rmss
          bat(1) = iat
          bidx(1) = idx
          tmstat(10) = tmstat(10) + wraptime() - tbeg
          return
        endif
      endif
c     
      if ( ir .eq. 1 ) then
        tmstat(10) = tmstat(10) + wraptime() - tbeg
        return
      endif
c     
      i = l
      j = l + l
 40   continue
      if ( j .le. ir ) then
        if ( j .lt. ir ) then
          if ( bidx(j) .lt. bidx(j+1)) j = j+1
        endif
        if ( idx .lt. bidx(j)) then
          xcrd(i) = xcrd(j)
          ycrd(i) = ycrd(j)
          zcrd(i) = zcrd(j)
          xfrc(i) = xfrc(j)
          yfrc(i) = yfrc(j)
          zfrc(i) = zfrc(j)
          zacc(i) = zacc(j)
          mbuf(i) = mbuf(j)
          bat(i) = bat(j)
          bidx(i) = bidx(j)
          i = j
          j = j + j
        else
          j = ir + 1
        endif
        go to 40
      endif
      xcrd(i) = rr1(1,1)
      ycrd(i) = rr1(2,1)
      zcrd(i) = rr1(3,1)
      xfrc(i) = rr1(1,2)
      yfrc(i) = rr1(2,2)
      zfrc(i) = rr1(3,2)
      zacc(i) = rr1(3,3)
      mbuf(i) = rmss
      bat(i) = iat
      bidx(i) = idx
      goto 30
c
      end
