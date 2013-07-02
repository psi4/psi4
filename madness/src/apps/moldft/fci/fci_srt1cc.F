*
* $Id: fci_srt1cc.F,v 1.2 1999/07/28 00:23:42 d3e129 Exp $
*
      subroutine fci_srt1cc(ionel,nonel,kbase,kdim,kpt,mkpt)
      implicit none
#include "fciinfo.fh"
c     
c     generate an index vector into the one electron coupling coefficients
c     according to the bra determinant and also subsorted by ket determinant
c     
c     ionel ... one electron coupling coefficients from onelcc
c     kdim ... no. of k (bra) determinants
c     kpt ... index vector to be
c     mkpt ... pointer into index vector indexed by k
c     
      integer nonel, kbase, kdim
      integer ionel(5,nonel),kpt(nonel),mkpt(kdim+1)
c
      integer iel, k, ksum, mmkpt, jdim, koff
*      integer ji, jj, jel, iiel, jjel  ! used in slow, but correct, sort
c     
c     first generate index according to k
c     
      call ifill(0,mkpt,1,kdim+1)
      do 10 iel = 1,nonel
         k = ionel(3,iel) - kbase
         mkpt(k) = mkpt(k) + 1
 10   continue
      ksum = 0
      do 20 k = 1,kdim
         mmkpt = mkpt(k)
         mkpt(k) = ksum
         ksum = ksum + mmkpt
 20   continue
      if (ksum.ne.nonel) call fci_err(' screwed ksum',ksum)
      do 30 iel = 1,nonel
         k = ionel(3,iel) - kbase
         mkpt(k) = mkpt(k) + 1
         kpt(mkpt(k)) = iel
 30   continue
      if (mkpt(kdim).ne.nonel) call fci_err(' screwed mkpt ',mkpt(kdim))
      do 40 k = kdim + 1,2,-1
         mkpt(k) = mkpt(k-1)
 40   continue
      mkpt(1) = 0
c     
c     now naively sort according to j for given k
c     might have to get clever here for large cases
c     
      do 50 k = 1,kdim
         jdim = mkpt(k+1)-mkpt(k)
         koff = mkpt(k)
c
         call fci_hsort(kpt(koff+1), ionel, jdim)
c$$$
c$$$         do 60 iiel = 1,jdim
c$$$            iel = kpt(iiel+koff)
c$$$            ji = ionel(4,iel)
c$$$            do 70 jjel = iiel+1,jdim
c$$$               jel = kpt(jjel+koff)
c$$$               jj = ionel(4,jel)
c$$$               if (jj.lt.ji) then
c$$$                  kpt(iiel+koff) = jel
c$$$                  kpt(jjel+koff) = iel
c$$$                  iel = jel
c$$$                  ji = jj
c$$$               endif
c$$$ 70         continue
c$$$ 60      continue
 50   continue
c     
      return
      end
