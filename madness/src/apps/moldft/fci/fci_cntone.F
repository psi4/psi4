*
* $Id: fci_cntone.F,v 1.2 1999/07/28 00:23:39 d3e129 Exp $
*
      subroutine fci_cntone(nonelt,nonemx)
      implicit none
#include "fciinfo.fh"
c     
c     1) count all non-zero one-electron coupling coefficents -> nonelt
c     2) find the maximum involving any given orbital configuration -> nonemx
c
c     This new version just estimates the information.
c     nonemx is about 2x too big, nonelt is about 50% too big
c     
      integer nonelt            ! [output]
      integer nonemx            ! [output]
c
      integer kocc(nactmx),indbar(nactmx)
      logical oprint
c
      integer ns, k, nspt, mkdim, ntest,kwalk, kdet, nd
c
      oprint =.false.
      nonelt = 0
      nonemx = 0
      
      call fci_tstart(11,.false.)
c     find occupancy, compute no. of open shells and indbar vector
      do 5 kwalk = 1,nwalko
         call fci_ofind_occ(kwalk,kocc)
         kdet = indexo(kwalk)
         ns = 0
         nd = 0
         do 10 k = 1,nactiv
            if (kocc(k).eq.2) ns = ns + 1
            if (kocc(k).eq.3) nd = nd + 1
            indbar(k) = ns
 10      continue
         nspt = ns/2 + 1
         mkdim = mnsdet(ns)
         if (mkdim.eq.0) goto 5
         ntest = nd*(1+2*(nactiv+1-ns-nd)) + ns*(nactiv+1)
         ntest = ntest*mkdim
         if (nonemx .lt. ntest) nonemx = ntest
         nonelt = nonelt + ntest
 5    continue
      call fci_tfinish(11,.false.)
c
      end
