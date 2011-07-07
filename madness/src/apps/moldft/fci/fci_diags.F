*
* $Id: fci_diags.F,v 1.2 1999/07/28 00:23:39 d3e129 Exp $
*
      subroutine fci_diags(hd,ecore,h,g,ionel,nonemx,kpt,mkpt)
      implicit none
#include "fciinfo.fh"
c     
c     construct the diagonals of the hamiltonian 
c     
      double precision hd(ndets) ! [output] The diagonals
      double precision ecore    ! [input] Nuclear repulsion + frozen core
      double precision h(nacttr) ! [input] 1-e integrals
      double precision g(nacttr,nacttr) ! [input] 2-e integrals
      integer nonemx            ! [input] dimension of ionel
      integer ionel(5,nonemx)   ! [scratch] workspace for coupling coeff
      integer kpt(nonemx)       ! [scratch] workspace
      integer mkpt(ndetmx+1)    ! [scratch] workspace
c
      integer kwalk, kelbas, keltop, jprev, kkel,
     $     kel, k, l, kl, jdet, iello, ielhi, iel, idet, i, j, ij,
     $     nonel, iiel
      integer kocc(nactmx), kdim, kbase, ksym, kns, kdet
      double precision ekl, eijkl
c     
c     loop thru intermediate walks, compute all non-zero 1 e- coeffs
c     and also provide index into the coeffs based on intermediate
c     determinant
c     
      call fci_tstart(5,.false.)
      call dfill(ecore,hd,1,ndets)
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
                  ekl = ionel(5,kel)
c     do the one electron term
                  if (ksym .eq. statesym) then
                     if (jdet .eq. kdet+kbase) then
                        hd(jdet) = hd(jdet) +  ekl * h(kl)
                     endif
                  endif
                  ekl = ekl * 0.5d0
                  if (jdet.ne.jprev) then
c     look for last element with this jdet and kdet ... will be nearby
                     jprev = jdet
                     iello = kkel
                     ielhi = kkel
                     do 40 iiel = kkel+1,keltop
                        if (ionel(4,kpt(iiel)).ne.jdet) goto 50
                        ielhi = ielhi+1
 40                  continue
                  endif
c     loop thru non-zero <i|eij|k> for idet=jdet
 50               do 60 iiel = iello,ielhi
                     iel = kpt(iiel)
                     idet = ionel(4,iel)
                     i = max(ionel(1,iel),ionel(2,iel))
                     j = min(ionel(1,iel),ionel(2,iel))
                     ij = i*(i-1)/2 + j
                     eijkl = ionel(5,iel) * ekl
                     hd(jdet) = hd(jdet) + eijkl*g(ij,kl)
 60               continue
 30            continue
 20         continue
         endif
 10   continue
*      write(6,*) ' diagonals from diags '
*      write(6,1) hd
* 1    format(1x,6f12.6)
c     
      call fci_tfinish(5,.false.)
c
      end
