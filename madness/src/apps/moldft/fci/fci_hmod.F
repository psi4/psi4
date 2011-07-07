*
* $Id: fci_hmod.F,v 1.2 1999/07/28 00:23:40 d3e129 Exp $
*
      subroutine fci_hmod(h, g, hmod)
      implicit none
#include "fciinfo.fh"
c
c     hmod = h - 1/2 sum(k) (ik|kj) to accomodate use of
c     Eij,kl = EijEkl - delta(j,k)Eil
c
      double precision h(nacttr) ! [input] 1e Hamiltonian
      double precision g(nacttr,nacttr) ! [input] 2-e integrals
      double precision hmod(nacttr) ! [output] Modified integrals
c     
      integer i,j,k,ij,ind
      ind(i,j) = max(i,j)*(max(i,j)-1)/2 + min(i,j)
c
      do i = 1, nactiv
         do j = 1, i
            ij = ind(i,j)
            hmod(ij) = h(ij)
            do k = 1, nactiv
               hmod(ij) = hmod(ij) - 0.5d0*g(ind(i,k),ind(k,j))
            enddo
         enddo
      enddo
c
*      write(6,*) ' HMOD '
*      call fci_prntri(hmod,nactiv)
c
      end
