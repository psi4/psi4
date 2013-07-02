*
* $Id: fci_hami.F,v 1.2 1999/07/28 00:23:40 d3e129 Exp $
*
      subroutine fci_hami(hamil, ecore, h, g, work, lenwork)
      implicit none
#include "fciinfo.fh"
c
c     Simple wrapper routine to build the full triangular Hamiltonian
c     from the unmodified one- and two-electron integrals
c
      double precision hamil(*) ! [output] Returns Hamiltonian
      double precision ecore    ! [input] Frozen core + nuclear replulsion
      double precision h(nacttr) ! [input] Triangle of unmodified 1-e integrals
      double precision g(nacttr,nacttr) ! [input] 2-e integrals
      integer lenwork           ! [input]
      double precision work(lenwork) ! [scratch]
c
      integer ih, io, need, nonelt, nonemx, ik, im, i
c
      call fci_cntone(nonelt, nonemx)
      ih = 1
      io = ih + nacttr
      ik = io + (nonemx*5 + 1) / nav
      im = ik + (nonemx + 1) / nav
      need = im + ndetmx
      if (need .gt. lenwork) call fci_err('hami: need ', need)
c
c     Modify the 1-e integrals
c
      call fci_hmod(h, g, work(ih))
c
      call dfill(0.0d0, hamil, 1, ndets*(ndets+1)/2)
      call fci_hamil2(hamil,work(ih),g,work(io),nonemx, 
     $     work(ik),work(im))
c
      do i = 1, ndets
         hamil(i*(i+1)/2) = hamil(i*(i+1)/2) + ecore
      enddo
c
      end
      subroutine fci_hamil2(hamil,h,g,ionel,nonemx,kpt,mkpt)
      implicit none
#include "fciinfo.fh"
c     
c     construct the entire hamiltonian over determinants
c     
      integer nonemx
      double precision hamil(*),h(nacttr),g(nacttr,nacttr)
      integer ionel(5,nonemx),kpt(nonemx),mkpt(*)
c     
c     loop thru intermediate walks, compute all non-zero 1 e- coeffs
c     and also provide index into the coeffs based on intermediate
c     determinant
c     
      integer kwalk, nonel, kdet, kelbas, keltop,
     $     jprev, kkel, kel, k, l, kl, jdet, jjdet, ieltop, iel, iiel,
     $     idet, ijdet, i, j, ij
      integer kocc(nactmx), kns, kbase, kdim, ksym
      double precision ekl, eijkl

      call fci_tstart(13,.false.)
      call dfill(0.0d0,hamil,1,ndets*(ndets+1)/2)
      do 10 kwalk = 1,nwalko
         call fci_owalk_info(kwalk,kocc,kns,kbase,kdim,ksym)
         call fci_onelcc(kocc,kbase,kdim,ksym,ionel,nonel,nonemx)
         if (nonel.gt.0) then
            call fci_srt1cc(ionel,nonel,kbase,kdim,kpt,mkpt)
c     loop thru non-zero <k|ekl|j>
            do 20 kdet = 1,kdim
               kelbas = mkpt(kdet) + 1
               keltop = mkpt(kdet+1)
               jprev = -1
               do 30 kkel = kelbas,keltop
                  kel = kpt(kkel)
                  k = max(ionel(1,kel),ionel(2,kel))
                  l = min(ionel(1,kel),ionel(2,kel))
                  kl = k*(k-1)/2 + l
                  jdet = ionel(4,kel)
                  jjdet = jdet*(jdet-1)/2
                  ekl = float(ionel(5,kel))
c     do the one electron term
                  if (ksym .eq. statesym) then
                     if (jdet .ge. kdet+kbase) then
                        hamil(jjdet+kdet+kbase) = 
     $                       hamil(jjdet+kdet+kbase) + ekl*h(kl)
                     endif
                  endif
                  ekl = ekl * 0.5d0
                  if (jdet.ne.jprev) then
c     look for last element with this jdet and kdet ... will be nearby
                     jprev = jdet
                     ieltop = kkel
                     do 40 iiel = kkel+1,keltop
                        if (ionel(4,kpt(iiel)).ne.jdet) goto 50
                        ieltop = ieltop+1
 40                  continue
                  endif
c     loop thru non-zero <i|eij|k> for idet<=jdet
 50               do 60 iiel = kelbas,ieltop
                     iel = kpt(iiel)
                     idet = ionel(4,iel)
                     ijdet = jjdet + idet
                     i = max(ionel(1,iel),ionel(2,iel))
                     j = min(ionel(1,iel),ionel(2,iel))
                     ij = i*(i-1)/2 + j
                     eijkl = dfloat(ionel(5,iel)) * ekl
                     hamil(ijdet) = hamil(ijdet) + eijkl*g(ij,kl)
 60               continue
 30            continue
 20         continue
         endif
 10   continue
c     
      call fci_tfinish(13,.false.)
c
      end
