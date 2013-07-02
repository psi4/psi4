*
* $Id: fci_2pdm.F,v 1.2 1999/07/28 00:23:39 d3e129 Exp $
*
      subroutine fci_2pdm(vec, nvec, twopdm, onepdm, ecore, h, g)
      implicit none
#include "fciinfo.fh"
c     
c     Compute the 2-particle density matrices connecting the nvec 
c     CI vectors.  This routine aids in testing the fci_couple2 routine.
c     It is more efficient to compute the actual 2-pdm with another
c     factorization <Epq,rs> = <EpaErb|Kx><Kx|EbsEaq> (or variations
c     thereof).
c     
c     The 1-pdm and integrals are passed in here solely to test 
c     the 2-pdm trace and energies
c     
c     twopdm(p,q,r,s,ivec,jvec) returns <V(ivec)|Epq,rs|V(jvec)>
c     
      integer nvec              ! [input] No. of vectors
      double precision vec(ndets, nvec) ! [input] The CI vectors
      double precision twopdm(nactiv,nactiv,nactiv,nactiv,nvec,nvec) ! [output]
      double precision onepdm(nactiv,nactiv,nvec,nvec)
      double precision ecore
      double precision h(nacttr)
      double precision g(nacttr,nacttr)
c     
      integer nvecmx            ! Maximum value of nvec
      parameter (nvecmx=3) 
      double precision epqrs(ndetmx*nvecmx*nactmx*nactmx) ! Workspace
      integer numqs, qslist(2,nactmx*nactmx)
      integer iwalk, idim, iocc(nactmx), nsi, ioff, p, q, r, s, qs
      integer ivec, jvec, ivqs, w, shi, sym
      double precision sum, diff, test, count
c     
      integer ind
      ind(p,q) = max(p,q)*(max(p,q)-1)/2 + min(p,q)
c     
      test = 1d-4
      count = 0.0d0
c     
      if (nvec.gt.nvecmx) call fci_err('fci_2pdm: nvecmx? ', nvec)
c     
      call dfill(0.0d0,twopdm,1,nactiv*nactiv*nactiv*nactiv*nvec*nvec)
c     
c     Loop thru orbital walks I and get information
c     
      do 10 iwalk = 1, nwalko
         call fci_owalk_info(iwalk, iocc, nsi, ioff, idim, sym)
         if (sym .ne. statesym) goto 10
*     write(6,*) ' iwalk ', iwalk
c     
c     Loop thru orbitals p and r
c     
         do p = 1, nactiv
            do r = 1, p
c     
c     Form  <V|Epq,rs|Iw> for all q,s
c     
               call fci_couple2(vec, nvec, p, r, iocc, idim, 
     $              numqs, qslist, epqrs)
               count = count + numqs*idim*nvec
c     
c     Multiply with the CI vectors and add into the matrices
c     
               do qs = 1, numqs
                  q = qslist(1,qs)
                  s = qslist(2,qs)
                  do ivec = 1, nvec
                     ivqs = idim*(ivec-1 + nvec*(qs-1))
                     do jvec = 1, nvec
                        sum = 0.0d0
                        do w = 1, idim
                           sum = sum + epqrs(ivqs+w)*vec(ioff+w,jvec)
                        enddo
                        twopdm(p,q,r,s,ivec,jvec) = 
     $                       twopdm(p,q,r,s,ivec,jvec) + sum
                     enddo
                  enddo
               enddo
            enddo
         enddo
 10   continue                  ! iwalk loop
c     
c     Symmetrize (pq) <-> (rs)
c     
      do ivec = 1, nvec
         do jvec = 1, nvec
            do p = 1, nactiv
               do r = 1, p
                  do q = 1, nactiv
                     shi = nactiv
                     if (p .eq. r) shi = q
                     do s = 1, shi
                        twopdm(r,s,p,q,ivec,jvec) = 
     $                       twopdm(p,q,r,s,ivec,jvec)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
c     
      write(6,31) count
 31   format(/' No. of non-zero two-body coupling coeffs ', 1p, d9.2/)
c     
c     Done.  
c     
c     For debugging check : 
c     <ivec|Epq,rs|jvec> = <jvec|Eqp,sr|ivec>
c     <ivec|Epq,rs|jvec> = <ivec|Ers,pq|jvec> ... now that we use this
c     symmetry in the computation this test is nearly useless.
c     
      do ivec = 1, nvec
         do jvec = 1, nvec
c$$$  write(6,*) ' TWOPDM ', ivec, jvec
c$$$  do p = 1, nactiv
c$$$  do q = 1, nactiv
c$$$  do r = 1, nactiv
c$$$  do s = 1, nactiv
c$$$  write(6,112) p,q,r,s,twopdm(p,q,r,s,ivec,jvec)
c$$$  112                    format(4i5,2x,f10.6)
c$$$  enddo
c$$$  enddo
c$$$  enddo
c$$$  enddo
            do p = 1, nactiv
               do q = 1, nactiv
                  do r = 1, nactiv
                     do s = 1, nactiv
c     
                        diff = twopdm(p,q,r,s,ivec,jvec) - 
     $                       twopdm(q,p,s,r,jvec,ivec)
                        if (abs(diff) .gt. test) call fci_err
     $                       ('fci_2pdm: symmetry 1 fails',
     $                       1000*p+100*q+10*r+s)
c     
                        diff = twopdm(p,q,r,s,ivec,jvec) - 
     $                       twopdm(r,s,p,q,ivec,jvec)
                        if (abs(diff) .gt. test) call fci_err
     $                       ('fci_2pdm: symmetry 2 fails',
     $                       1000*p+100*q+10*r+s)
                     enddo
                  enddo
               enddo
            enddo
            write(6,*) ' 2pdm ', ivec,jvec,' has correct symmetry'
         enddo
      enddo
c     
c     For debugging check that sum(s) Fpq,ss = (nelec-1)*Ppq
c     
      do ivec = 1, nvec
         do jvec = 1, nvec
            do p = 1, nactiv
               do q = 1, nactiv
                  sum = 0.0d0
                  do s = 1, nactiv
                     sum = sum + twopdm(p,q,s,s,ivec,jvec)
                  enddo
                  sum = sum / dfloat(nelec-1)
                  if (abs(sum-onepdm(p,q,ivec,jvec)) .gt. test) 
     $                 call fci_err('2pdm - trace is incorrect', 0)
               enddo
            enddo
            write(6,*) ' 2pdm ', ivec, jvec, ' has correct trace'
         enddo
      enddo
c     
c     Compute the expectation value of the Hamiltonian from
c     the density matrices.
c     
      do ivec = 1, nvec
         do jvec = 1, nvec
            sum = 0.0d0
            do p = 1, nactiv
               do q = 1, nactiv
                  sum = sum + onepdm(p,q,ivec,jvec)*h(ind(p,q))
c     
                  do r = 1, nactiv
                     do s = 1,nactiv
                        sum = sum + 0.5d0*twopdm(p,q,r,s,ivec,jvec)*
     $                       g(ind(p,q),ind(r,s))
                     enddo
                  enddo
               enddo
            enddo
            if (ivec.eq.jvec) sum = sum + ecore
            write(6,77) ivec, jvec, sum
 77         format(' H(',i2,',',i2,') = ', f20.8)
         enddo
      enddo
c     
      end

