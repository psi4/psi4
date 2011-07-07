*
* $Id: fci_sigma.F,v 1.2 1999/07/28 00:23:42 d3e129 Exp $
*
      subroutine fci_sigma(c,s,ecore,h,g,ionel,i1dim,nonemx,
     $     dkl,bij,kblmax)
      implicit none
#include "fciinfo.fh"
c     
c     s = hamiltonian * c (MUST PRESENTLY BE IN DETERMINANTS)
c     
      double precision c(ndets) ! [input] CI vector
      double precision s(ndets) ! [output] H*c
      double precision ecore    ! [input] Nuclear repulsion+frozen core
      double precision h(nacttr)! [input] Modified 1-e integrals
      double precision g(nacttr,nacttr) ! [input] 2-e integrals
      double precision dkl(*), bij(*) ! [scratch] Workspace
      integer i1dim             ! [input] last dimension of ionel
      integer nonemx            ! [input] Max # coupling coeff from onelcc
      integer ionel(5,i1dim)    ! [scratch] workspace to hold coupling coeff
      integer kblmax            ! [input] max size of any intermdediate state
c      
      integer kblock, koff, num1cc, kwalk, kdim, kel,
     $     k, l, kl, kaddr, jdet, idet, nonel
      integer kocc(nactmx), kbase, ksym, kns
      double precision ddot
      integer idamax
      external ddot, idamax
c
      call fci_tstart(8,.false.)
      call dfill(0.0d0,s,1,ndets)
      call dscal(nacttr,2.0d0,h,1)
c     loop thru intermediate states ... process when have all can hold
      kblock = 0
      num1cc = 0
      kwalk  = 1
c
 10   if (kwalk.gt.nwalko) goto 20
      call fci_owalk_info(kwalk,kocc,kns,kbase,kdim,ksym)
      if (kblock+kdim.le.kblmax .and. num1cc+nonemx.le.i1dim) then
         koff = kblock
         call fci_onelcc(kocc,koff,kdim,ksym,
     $        ionel(1,num1cc+1),nonel,nonemx)
c
c     Do the 1-e bit while we know the symmetry of K and real kbase
c
         if (ksym .eq. statesym) then
            do kel = num1cc+1,num1cc+nonel
               k = max(ionel(1,kel),ionel(2,kel))
               l = min(ionel(1,kel),ionel(2,kel))
               kl = iky(k) + l
               idet = ionel(3,kel) + kbase - koff
               jdet = ionel(4,kel)
               s(idet) = s(idet) + h(kl)*c(jdet)*ionel(5,kel)
            enddo
         endif
c
         num1cc = num1cc + nonel
         kwalk = kwalk + 1
         kblock = kblock + kdim
         goto 10
      endif
c
 20   continue
c
      if (kblock .gt. 0) then
c     vector compute required addresses
         do 30 kel = 1,num1cc
            k = max(ionel(1,kel),ionel(2,kel))
            l = min(ionel(1,kel),ionel(2,kel))
            kl = iky(k) + l
            ionel(1,kel) = ionel(3,kel) + (kl - 1)*kblock
 30      continue
c     make dkl
         call dfill(0.0d0,dkl,1,nacttr*kblock)
         do 40 kel = 1,num1cc
            kaddr = ionel(1,kel)
            jdet = ionel(4,kel)
            dkl(kaddr) = dkl(kaddr) + ionel(5,kel) * c(jdet)
 40      continue
c     make bij
         call dgemm('n','n',kblock,nacttr,nacttr,1.0d0,dkl,kblock,
     $        g,nacttr,0.0d0,bij,kblock)
*         call fci_axb(kblock,nacttr,nacttr,dkl,g,bij)
c     add in 2-e sigma component 
         do 50 kel = 1,num1cc
            kaddr = ionel(1,kel)
            idet = ionel(4,kel)
            s(idet) = s(idet) + ionel(5,kel) * bij(kaddr) 
 50      continue
      else if (kwalk.le.nwalko) then
         write(6,*) ' kblock, kdim, kblmax, num1cc, nonemx, i1dim ',
     $        kblock, kdim, kblmax, num1cc, nonemx, i1dim 
         call fci_err('fci_sigma: kblock is zero ', kblock)
      endif
c     goto next kwalk if any left
      kblock = 0
      num1cc = 0
      if (kwalk.le.nwalko) goto 10
c     
      call dscal(ndets,0.5d0,s,1)
      call dscal(nacttr,0.5d0,h,1)
      call daxpy(ndets, ecore, c, 1, s, 1)
c
      call fci_tfinish(8,.false.)
c
      end
