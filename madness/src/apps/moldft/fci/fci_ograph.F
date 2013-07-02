*
* $Id: fci_ograph.F,v 1.2 1999/07/28 00:23:41 d3e129 Exp $
*
      subroutine fci_ograph
      implicit none
#include "fciinfo.fh"
c     
c     generate orbital configuration graphs for nelec electron systems
c     
      integer iorb, irow, k, krow, kdead, isum, ir,
     $     iwalk, iocc(nactmx), numwalko
      logical oprint
c     
      oprint = .false.
c     
      call ifill(0,karc,1,3*(nelmax+1)*(nactmx+1))
      call ifill(0,karcwt,1,3*(nelmax+1)*(nactmx+1))
      call ifill(0,kwt,1,(nelmax+1)*(nactmx+1))
c     
c     
c     generate graph including some incorrect tails
c     
      kwt(1,1) = 1
      do 10 iorb = 1,nactiv
         do 20 irow = 1,min(2*(iorb-1),nelec)+1
            if (kwt(irow,iorb).gt.0) then
               do 30 k=1,3
                  krow = irow + (k-1)
                  if (krow.le.nelec+1) then
                     karc(k,irow,iorb) = krow
                     kwt(krow,iorb+1) = kwt(krow,iorb+1) +
     &                    kwt(irow,iorb)
                  endif
 30            continue
            endif
 20      continue
 10   continue
      nwalko = kwt(nelec+1,nactiv+1)
c     
c     remove incorrect tails
c     
      do 40 irow = 1,nelec
         kwt(irow,nactiv+1) = 0
 40   continue
      do 50 iorb = nactiv,1,-1
         do 60 irow = 1,min(2*(iorb-1),nelec)+1
            kdead = 0
            do 70 k = 1,3
               krow = irow + (k-1)
               if (krow.gt.nelec+1) then
                  kdead = kdead + 1
               else
                  if (kwt(krow,iorb+1).eq.0) then
                     karc(k,irow,iorb) = 0
                     kdead = kdead + 1
                  endif
               endif
 70         continue
            if (kdead.eq.3) kwt(irow,iorb) = 0
 60      continue
 50   continue
c     
c     generate arc weight array
c     
      do 80 iorb = nactiv+1,2,-1
         do 90 irow = 1,min(2*(iorb-1),nelec)+1
            if(kwt(irow,iorb).gt.0) then
               isum = 0
               do 100 k = 1,3
                  krow = irow - (k-1)
                  if (krow.gt.0) then
                     if (karc(k,krow,iorb-1).gt.0) then
                        karcwt(k,krow,iorb-1) = isum
                        isum = isum + kwt(krow,iorb-1)
                     endif
                  endif
 100           continue
            endif
 90      continue
 80   continue
c     
      if (oprint) then
         write(6,3) nelec,nactiv,nwalko
 3       format(1x,' orbital configuration graph '//
     &        ' no. of electrons ',i6/
     &        ' no. of orbitals  ',i6/
     &        ' no. of walks     ',i6/)
         do 110 iorb = 1,nactiv+1
            ir = min(2*(iorb-1),nelec) + 1
            write(6,1) iorb,(kwt(irow,iorb),irow=1,ir)
 1          format(1x,i4,8(5x,i4,5x):/5x,8(5x,i4,5x))
            write(6,2) ((karc(k,irow,iorb),k=1,3),irow=1,ir)
 2          format(5x,8(1x,3i4,1x):/5x,8(1x,3i4,1x))
            write(6,2) ((karcwt(k,irow,iorb),k=1,3),irow=1,ir)
            write(6,*)
 110     continue
c     
         do iwalk = 1, numwalko
            call fci_ofind_occ(iwalk, iocc)
            write(6,222) iwalk, (iocc(k),k=1,nactiv)
 222        format(1x,i8,2x,20i2)
         enddo
      endif
c     
      end
