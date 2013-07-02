*
* $Id: fci_3pdm.F,v 1.2 1999/07/28 00:23:39 d3e129 Exp $
*
      subroutine fci_3pdm(vec, nvec, work, twopdm)
      implicit none
#include "fciinfo.fh"
c     
      integer nvec              ! [input] No. of vectors
      double precision vec(ndets, nvec) ! [input] The CI vectors
      double precision work(nactiv,nactiv,nactiv,nactiv,nvec,nvec)
      double precision twopdm(nactiv,nactiv,nactiv,nactiv,nvec,nvec) ! [output]
      integer nvecmx            ! Maximum value of nvec
      parameter (nvecmx=3) 
      double precision epqrstu(ndetmx*nvecmx*nactmx*nactmx) ! Workspace
      integer numsu, sulist(2,nactmx*nactmx)
      integer iwalk, idim, iocc(nactmx), nsi, ioff, p, q, r, s, t, u, su
      integer pq, rs, tu, sym
      integer ivec, jvec, ivsu, w, uhi
      double precision sum, diff, test, count
c     
c     Test code for the 3-pdm coupling coeff generation ... there are
c     much more efficient ways to compute the actual density matrix.
c     
      integer indsq
      indsq(p,q) = (p-1)*nactiv + q
c     
      test = 1d-4
      count = 0.0d0
c     
      if (nvec.gt.nvecmx) call fci_err('fci_2pdm: nvecmx? ', nvec)
c     
c     Generate all of the coupling coefficients test the code
c     but only accumulate the sum(p) Epp,rs,tu for comparison
c     with the two-pdm
c
      call dfill(0.0d0, work, 1, (nactiv**4)*(nvec**2))
      count = 0.0d0
      do p = 1, nactiv
         do q = 1, nactiv
            pq = indsq(p,q)
            do 10 iwalk = 1, nwalko
               call fci_owalk_info(iwalk, iocc, nsi, ioff, idim, sym)
               if (sym.ne.statesym .or. iocc(q).eq.1) goto 10
               do r = 1, p
                  do t = 1, r
                     call fci_couple3(vec, nvec, p, q, r, t, 
     $                    iocc, idim, numsu, sulist, epqrstu)
                     count = count + numsu*idim*nvec
c     
                     do su = 1, numsu
                        s = sulist(1,su)
                        u = sulist(2,su)
                        rs = indsq(r,s)
                        tu = indsq(t,u)
                        do ivec = 1, nvec
                           ivsu = idim*(ivec-1 + nvec*(su-1))
                           do jvec = 1, nvec
                              sum = 0.0d0
                              do w = 1, idim
                                 sum = sum + epqrstu(ivsu+w)*
     $                                vec(ioff+w,jvec)
                              enddo
c
                              if (p.eq.q) 
     $                             work(r,s,t,u,ivec,jvec) = 
     $                                work(r,s,t,u,ivec,jvec) + sum
                              if (pq.ne.rs .and. r.eq.s)
     $                             work(p,q,t,u,ivec,jvec) = 
     $                             work(p,q,t,u,ivec,jvec) + sum
                              if (pq.ne.tu .and. rs.ne.tu .and. t.eq.u) 
     $                             work(r,s,p,q,ivec,jvec) = 
     $                             work(r,s,p,q,ivec,jvec) + sum
                           enddo
                        enddo
                     enddo
                  enddo
               enddo
 10         continue            ! iwalk loop
         enddo
      enddo
c
c     Must still symmetrize off-diagonal over rs,tu
c
      do ivec = 1, nvec
         do jvec = 1, nvec
            do r = 1, nactiv
               do s = 1, nactiv
                  do t = 1, r
                     uhi = nactiv
                     if (t.eq.r) uhi = s - 1 ! Off diagonal only
                     do u = 1, uhi
                        sum = work(r,s,t,u,ivec,jvec) + 
     $                       work(t,u,r,s,ivec,jvec)
                        work(r,s,t,u,ivec,jvec) = sum
                        work(t,u,r,s,ivec,jvec) = sum
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
c     
      write(6,1) count
 1    format(/' No. of 3-particle coupling coeffs = ', 1p,d9.2)
c     
c     Test the PQ trace
c     
      do ivec = 1, nvec
         do jvec = 1, nvec
            do u = 1, nactiv
               do t = 1, nactiv
                  do s = 1, nactiv
                     do r = 1, nactiv
                        diff = twopdm(r,s,t,u,ivec,jvec) - 
     $                       work(r,s,t,u,ivec,jvec)/dfloat(nelec-2)
                        if (abs(diff) .gt. test) call fci_err
     $                       ('3pdm: bad trace ',1000*r+100*s+10*t+u)
                     enddo
                  enddo
               enddo
            enddo
            write(6,*) ' 3-pdm ', ivec, jvec, ' has good pq trace'
         enddo
      enddo
c     
      end
