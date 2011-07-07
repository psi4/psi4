*
* $Id: fci_intchk.F,v 1.2 1999/07/28 00:23:40 d3e129 Exp $
*
      subroutine fci_intchk(h, g)
      implicit none
#include "fciinfo.fh"
      double precision h(nacttr) ! [input/output]
      double precision g(nacttr,nacttr) ! [input/output]
c
c     Ensure that the integrals have the correct symmetry
c     and force small integrals that should be zero to a hard zero.
c
c     This is just an interim routine until we have full symmetry packing
c
      integer i, j, k, l, ij, kl, symij, symkl
#include "bitops.fh"
c
      do i = 1, nactiv
         do j = 1, i
            ij = iky(i)+j
            symij = ieor(orbsym(i),orbsym(j))
            if (symij .ne. 0) then
               if (abs(h(ij)) .gt. 1d-6) then 
                  write(6,*) ' h ', i, j, h(ij)
                  call fci_err('fci_intchk:should be zero',ij)
               endif
               h(ij) = 0.0d0
            endif
            do k = 1, nactiv
               do l = 1, k
                  kl = iky(k)+l
                  symkl = ieor(orbsym(k),orbsym(l))
                  if (ieor(symij,symkl).ne.0) then
                     if (abs(g(ij,kl)) .gt. 1d-6) then 
                        write(6,*) ' g ', i, j, k, l, g(ij,kl)
                        call fci_err('fci_intchk:should be zero',0)
                     endif
                     g(ij,kl) = 0.0d0
                  endif
               enddo
            enddo
         enddo
      enddo
c
c     While we are it it check the ij -- kl symmetry
c
      do ij = 1, (nactiv*(nactiv+1))/2
         do kl = 1, ij-1
            if (abs(g(ij,kl)-g(kl,ij)) .gt. 1d-10) 
     $           call fci_err('fci_intchk: g is not symmetric',0)
            g(ij,kl) = 0.5d0*(g(ij,kl)+g(kl,ij))
            g(kl,ij) = g(ij,kl)
         enddo
      enddo
c
      end
