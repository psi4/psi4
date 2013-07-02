*
* $Id: fci_load.F,v 1.2 1999/07/28 00:23:40 d3e129 Exp $
*
      subroutine fci_load(h,g,ecore)
      implicit none
#include "fciinfo.fh"
c
c     Read the 1e/2e integrals from a simple ASCII list ... this
c     for debug only from testcase integrals.
c
      double precision h(nacttr) ! [output]
      double precision g(nacttr,nacttr) ! [output]
      double precision ecore    ! [output]
c
      integer i, j, k, l, ij, kl, ind
      double precision gg
c      
      ind(i,j) = max(i,j)*(max(i,j)-1)/2 + min(i,j)
c
      call dfill(0.0d0,g,1,nacttr*nacttr)
      call dfill(0.0d0,h,1,nacttr)
c
      ecore = 0.0d0
c
10    read (5,*,end=20) gg,i,j,k,l
      if (i.gt.nactiv .or. j.gt.nactiv .or. k.gt.nactiv .or.
     $     l.gt.nactiv) goto 10
      if (k.gt.0) then
        ij = ind(i,j)
        kl = ind(k,l)
        g(ij,kl) = gg
        g(kl,ij) = gg
      else if (i.gt.0) then
         if (abs(gg) .lt. 1d-6) gg = 0.0d0
        ij = ind(i,j)
        h(ij) = gg
      else
        ecore = gg
      endif
      goto 10
c
20    continue
c
      end
