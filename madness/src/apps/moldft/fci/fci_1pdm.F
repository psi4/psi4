*
* $Id: fci_1pdm.F,v 1.2 1999/07/28 00:23:39 d3e129 Exp $
*
      subroutine fci_1pdm(vec, nvec, onepdm)
      implicit none
#include "fciinfo.fh"
c
c     Compute the 1-particle density matrices connecting the nvec 
c     CI vectors.  This routine aids in testing the fci_couple1 routine.
c
c     onepdm(p,q,ivec,jvec) returns <V(ivec)|Epq|V(jvec)>
c
      integer nvec              ! [input] No. of vectors
      double precision vec(ndets, nvec) ! [input] The CI vectors
      double precision onepdm(nactiv,nactiv,nvec,nvec) ! [output] 1 PDMs
c
      integer nvecmx           ! Maximum value of nvec
      parameter (nvecmx=3) 
      double precision epq(ndetmx*nvecmx*nactmx)
      integer numq, qlist(nactmx)
      integer iwalk, idim, iocc(nactmx), nsi, ioff, p, q, qq
      integer ivec, jvec, ivq, w, sym
      double precision sum
c
c     Workspace for diagonalization of 1-pdm for debugging
c
      double precision pdm(nactmx*nactmx), evecs(nactmx,nactmx),
     $     evals(nactmx), iwork(nactmx)
c
      if (nvec.gt.nvecmx) call fci_err('fci_1pdm: nvecmx? ', nvec)
c
      call dfill(0.0d0, onepdm, 1, nactiv*nactiv*nvec*nvec)
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
                  do jvec = 1, nvec
                     sum = 0.0d0
                     do w = 1, idim
                        sum = sum + epq(ivq+w)*vec(ioff+w,jvec)
                     enddo
                     onepdm(p,q,ivec,jvec) = onepdm(p,q,ivec,jvec) + sum
                  enddo
               enddo
            enddo
         enddo
 10      continue               ! iwalk loop
c
c     Done.  
c
c     For debugging check the <ivec|Epq|jvec> = <jvec|Eqp|ivec>
c
      do ivec = 1, nvec
         do jvec = 1, nvec
            write(6,*) ' PDM ', ivec, jvec
            call fci_output(onepdm(1,1,ivec,jvec),1,nactiv,
     $           1,nactiv,nactiv,nactiv,1)
            do p = 1, nactiv
               do q = 1, nactiv
                  if (abs(onepdm(p,q,ivec,jvec)-onepdm(q,p,jvec,ivec))
     $                 .gt. 1d-10) call fci_err
     $                 ('fci_1pdm: not symmetric ', 0)
               enddo
            enddo
            write(6,*) ' 1pdm ', ivec,jvec,' has correct symmetry'
         enddo
      enddo
c
c     For debugging check that sum(p) <ivec|Epp|ivec> = nelec
c
      do ivec = 1, nvec
         sum = 0.0d0
         do p = 1, nactiv
            sum = sum + onepdm(p,p,ivec,ivec)
         enddo
         write(6,*) ' 1pdm ', ivec,ivec,' has trace', sum
         if (abs(sum-dfloat(nelec)).gt.1d-10) call fci_err
     $        ('fci_1pdm: incorrect trace ', ivec)
      enddo
c
c     For debugging diagonalize <ivec|Epp|ivec> and ensure it
c     it is positive definite with all eigenvalues < 2
c
      do ivec = 1, nvec
         do p = 1, nactiv
            do q = 1, p
               pdm(p*(p-1)/2+q) = onepdm(p,q,ivec,ivec)
            enddo
         enddo
         call yacobi(pdm,evecs,nactiv,evals,iwork,2)
         write(6,*) ' Natural orbital occupations for ', ivec
         call fci_output(evals,1,nactiv,1,1,nactiv,1,1)
         write(6,*) ' Natural orbitals for ', ivec
         call fci_output(evecs,1,nactiv,1,nactiv,nactiv,nactiv,1)
         do p = 1, nactiv
            if (evals(p).lt.-1d-10 .or. (evals(p)-2.0d0).gt.1d-10)
     $           call fci_err('fci_1pdm: invalid 1pdm eigenvalue',p)
         enddo
      enddo
c
      end

