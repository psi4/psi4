*
* $Id: fci_hv.F,v 1.2 1999/07/28 00:23:40 d3e129 Exp $
*
      subroutine fci_hv(ecore, h, g, work, lenwrk, v, s)
      implicit none
#include "fciinfo.fh"
c
c     Compute the Hamiltonian matrix-vector product s = Hv 
c
      double precision ecore    ! [input] Frozen-core+nuclear replusion energy
      double precision h(nacttr) ! [input] Effective 1-e integrals Core+T+V
      double precision g(nacttr,nacttr) ! [input] Two-electron integrals
      integer lenwrk            ! [input] Size of workspace
      double precision work(lenwrk) ! [scratch] Workspace - the more the better
      double precision v(ndets) ! [input] CI vector
      double precision s(ndets) ! [output] s = H*v
c
      integer nonelt, nonemx, kdimmx, kblmax, i1dim, numav
      integer ihmod, ionel, idkl, ibij, itop, ns
      double precision avgdim   ! Average no. of dets per orbital occ
      double precision avgone   ! Average no. of coupling coeff per orbital occ
c
c     Figure out max number of coupling coeffs and intermediate states
c     
      call fci_cntone(nonelt, nonemx)
      if (nelec .le. nactiv) then ! Max number of open shell in this calc
         ns = nelec
      else
         ns = 2*nactiv - nelec
      endif
      kdimmx = mnsdet(ns) ! Max no. of dets for an intermediate state
c
c     Allocate memory.  Attempt to hold numav average size intermediate
c     orbital occupations per pass, noting must be able to hold at least kdimmx
c
      avgdim = dble(ndets)/dble(nwalkosym)
      avgone = dble(nonelt)/dble(nwalko)
c
      numav = min(nwalko, 256)
c
 10   kblmax = max(kdimmx, nint(numav*avgdim)) ! Max dets per block
      i1dim  = max(nonemx, nint(numav*avgone))
c
      ihmod = 1                 ! Modified h ... h + 1/2 sum(k) (ik|kj)
      ionel = ihmod + nacttr    ! 1-e coupling coefficients
      idkl  = ionel + (5*i1dim)/nav + 1 ! Intermediate matrix D
      ibij  = idkl  + nacttr*kblmax ! Intermediate matrix B
      itop  = ibij  + nacttr*kblmax ! Total usage
      if (itop .gt. lenwrk) then
         if (numav .le. 1)call fci_err('fci_hv: alloc error',itop)
         numav = numav / 2
         goto 10
      endif
c
c     Modify the input 1-e integrals 
c
      call fci_hmod(h,g,work(ihmod))
c
c     Form the product
c
      call fci_sigma(v,s,ecore,work(ihmod),g,work(ionel),
     &     i1dim,nonemx,work(idkl),work(ibij),kblmax)
c
      end


