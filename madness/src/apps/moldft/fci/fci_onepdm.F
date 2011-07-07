*
* $Id: fci_onepdm.F,v 1.2 1999/07/28 00:23:41 d3e129 Exp $
*
      subroutine fci_onepdm(vec, nvec, weight, onepdm)
      implicit none
#include "fciinfo.fh"
c     
c     Compute the AVERAGE 1-particle density matrix for the input
c     determinant CI vectors and weights.  Force the weights to
c     add up to unity.
c     
c     onepdm(p,q) returns sum(ivec) weight(ivec) <V(ivec)|Epq|V(ivec)>
c     
      integer nvec              ! [input] No. of vectors
      double precision vec(ndets, nvec) ! [input] The CI vectors
      double precision weight(nvec) ! [input] The weights
      double precision onepdm(nactiv,nactiv) ! [output] 1 PDM
c     
      integer nvecmx            ! Maximum value of nvec
      parameter (nvecmx=10) 
      double precision epq(ndetmx*nvecmx*nactmx) ! Should dynamically allocate
      integer numq, qlist(nactmx)
      integer iwalk, idim, iocc(nactmx), nsi, ioff, p, q, qq
      integer ivec, ivq, w, sym
      double precision sum
c     
      if (nvec.gt.nvecmx) call fci_err('fci_onepdm: nvecmx? ', nvec)
c
c     Normalize the weights
c
      sum = 0.0d0
      do ivec = 1, nvec
         sum = sum + weight(ivec)
      enddo
      if (sum .eq. 0.0d0) call fci_err('fci_onepdm: zero weight sum',0)
      do ivec = 1, nvec
         weight(ivec) = weight(ivec) / sum
      enddo
      write(6,11) (weight(ivec),ivec = 1, nvec)
 11   format(/' Weights for averaging density matrix ', 20f10.6)
      write(6,*)
c     
      call dfill(0.0d0, onepdm, 1, nactiv*nactiv)
c     
c     Loop thru orbital walks I and get information
c     
      do 10 iwalk = 1, nwalko
         call fci_owalk_info(iwalk, iocc, nsi, ioff, idim, sym)
         if (sym .ne. statesym) goto 10
c     
c     Loop thru orbitals p
c     
         do p = 1, nactiv
c     
c     Form  <V|Epq|Iw> for all q
c     
            call fci_couple1(vec, nvec, p, iocc, idim, 
     $           numq, qlist, epq)
c     
c     Multiply with the CI vectors and add into the matrices
c     
            do qq = 1, numq
               q = qlist(qq)
               do ivec = 1, nvec
                  ivq = idim*(ivec-1 + nvec*(qq-1)) ! epq(idim,nvec,nactiv)
                  sum = 0.0d0
                  do w = 1, idim
                     sum = sum + epq(ivq+w)*vec(ioff+w,ivec)
                  enddo
                  onepdm(p,q) = onepdm(p,q) + sum*weight(ivec)
               enddo
            enddo
         enddo
 10   continue                  ! iwalk loop
c     
c     Done.  Check the symmetry and then force exact symmetry.
c     
      do p = 1, nactiv
         do q = 1, nactiv
            if (abs(onepdm(p,q)-onepdm(q,p)) .gt. 1d-10) then
               write(6,*) ivec, p, q, onepdm(p,q), onepdm(q,p)
               call fci_err('fci_onepdm: not symmetric ', 100*p+q)
            endif
            sum = (onepdm(p,q)+onepdm(q,p))*0.5d0
            onepdm(p,q) = sum
            onepdm(q,p) = sum
         enddo
      enddo
c     
c     Check the trace
c     
      sum = 0.0d0
      do p = 1, nactiv
         sum = sum + onepdm(p,p)
      enddo
      if (abs(sum-dfloat(nelec)).gt.1d-10) call fci_err
     $     ('fci_onepdm: incorrect trace ', ivec)
c
      write(6,*) ' Average one particle density matrix '
      write(6,*)
      call fci_output(onepdm(1,1),1,nactiv,1,nactiv,nactiv,
     $     nactiv,1)
      write(6,*)
c     
      end
