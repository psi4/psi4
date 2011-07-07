*
* $Id: fci_make_eaj.F,v 1.2 1999/07/28 00:23:40 d3e129 Exp $
*
      subroutine fci_make_all_eaj
      implicit none
#include "fciinfo.fh"
c
      integer ns, ioff, len, mdim
c
      call fci_tstart(6,.false.)
      ioff = 1
      do ns = multi-1,nsmax,2
         mdim = mnsdet(ns)
         len = 3*(ns+1)*mdim
         if (ioff+len-1 .gt. eajmax) 
     $        call fci_err('eijmax too small', ioff+len-1)
         call fci_make_eaj(ns,mdim,eajdata(ioff),eajphase(1,1,ns))
         eajptr(ns) = ioff
         ioff = ioff + len
      enddo
      call fci_tfinish(6,.false.)
c
      end
      subroutine fci_make_eaj(ns,mdim,eaj,eajphase)
      implicit none
#include "fciparam.fh"
      integer ns                ! [input]
      integer mdim              ! [input]
      integer eaj(mdim,3,ns+1)  ! [output]
      double precision eajphase(3,*) ! [output]
c
c     Form the transformation of primtive spin functions under
c     the action of Eaj, i.e., compute the matrix elements
c     <Ku|Eaj|Iv>.  The three cases are explained below.
c     Only the action of the permutation on the spin function is
c     accounted for.  The phase information is stored in eajphase
c
      integer j, v, u, i
      integer vocc(nsmax), uocc(nsmax)
c
c     Type 1: n_j(I) = 1   <Ku|Eaj|Iv>
c
      do v = 1, mdim
         call fci_mfind_occ(v,vocc,ns)
         do j = 1, ns
c
c     |u> = (s...j) |v>
c
            do i = 1,j-1
               uocc(i) = vocc(i)
            enddo
            do i = j,ns-1
               uocc(i) = vocc(i+1)
            enddo
            uocc(ns) = vocc(j)
c
            call fci_mfind_weight(u,uocc,ns)
c
            eaj(v,1,j) = u
            eajphase(1,j) = dfloat(mod(ns-j+1,2)*2-1)
         enddo
      enddo
c
      if (ns+2 .gt. nsmax) return
c
c     Types 2 & 3: n_j(I) = 2 <Ku|Eaj|Iv>
c
      do v = 1, mdim
         call fci_mfind_occ(v,vocc,ns)
         do j = 1, ns+1
c
c     Type 2:   |u> = (j...1)(s+2...1) |+-v>  (alpha->a, beta->j)
c     .             = (j...1)|-v+> * -1^(s+1)
c
            do i = 1,j-1
               uocc(i) = vocc(i)
            enddo
            uocc(j) = 2         ! 2 = - step
            do i = j+1,ns+1
               uocc(i) = vocc(i-1)
            enddo
            uocc(ns+2) = 1      ! 1 = + step
c
            call fci_mfind_weight(u,uocc,ns+2)
            eaj(v,2,j) = u
            eajphase(2,j) = dfloat(mod(ns+j+1,2)*2-1)
c
c     Type 3:   |u> = (j...1)(s+2...1)(12) |+-v>  (beta->a, alpha->j)
c     .             = (j...1)|+v-> * -1^(s)
c
            uocc(j) = 1
            uocc(ns+2) = 2
            call fci_mfind_weight(u,uocc,ns+2)
            eaj(v,3,j) = u
            eajphase(3,j) = -eajphase(2,j)
         enddo
      enddo
c
      end
      subroutine fci_apply_eaj(jbar, nj, jdim, kdim, eaj, eajphase,
     $     ounitj, jperm, jval, kperm, kval)
      implicit none
c     
c     Apply the operator Eaj to a vector.  The input vectors (jperm, jval)
c     are permuted by the action of the operator and the phase information
c     is accumulated resulting in the output vectors (kperm, kval).
c     The input vectors are NOT modified.
c     
c     I.e, given X(w) form  Y(v) = <Kv|Eaj|Jw> * X(w) where X=jval and
c     Y=kval, and accumulating the action of the permutation in kperm
c     by kperm(v)=jperm(w).
c     
c     eaj should be passed in as eajdata(eajptr(ns))     [fciinfo.fh]
c     eahphase should be passed in as eajphase(1,1,ns)   [fciinfo.fh]
c     
c     if (ounitj) then
c     .  the program assumes that jperm(j)=j and jval(j)=1.0d0
c     .  and does NOT reference jperm or jval
c     else
c     .  the program operates on jperm/jval as described above
c     endif
c     
      integer jbar              ! [input] Posn of j in SOCC orbs in J (1...ns)
      integer nj                ! [input] Occupancy of j in J
      integer jdim              ! [input] No. of determinants for ns(J)
      integer kdim              ! [input] No. of determinants for ns(K)
      integer eaj(jdim,3,*)     ! [input] Eaj transformation matrix for ns(J)
      double precision eajphase(3,*) ! [input] Eaj phase information
      logical ounitj            ! [input] If true jperm/jval are assumed unit
      integer jperm(jdim)       ! [input] Input permutation
      double precision jval(jdim) ! [input] Input values X(w)
      integer kperm(kdim)       ! [output] Output permutation kperm(v)=jperm(w)
      double precision kval(kdim) ! [output] Output values <Kv|Eaj|Jw> * X(w)
c     
      integer i, case, type, w, v
      double precision phase
c     
      if (nj .eq. 2) then
         do i = 1, kdim
            kval(i)  = 0.0d0
            kperm(i) = 1        ! Valid value so that don't need to test on zero
         enddo
      endif
c     
      do case = 1, nj
         type = case + nj - 1
         phase = eajphase(type,jbar)
*         write(6,*) ' DOING type ', type, phase
c     
         if (ounitj) then
            do w = 1, jdim
               v = eaj(w,type,jbar)
               kperm(v) = w
               kval(v)  = phase
*               write(6,*) ' w, v, phase ', w, v, phase
            enddo
         else
            do w = 1, jdim
               v = eaj(w,type,jbar)
               kperm(v) = jperm(w)
               kval(v)  = jval(w)*phase
            enddo
         endif
      enddo
c     
      end
      subroutine fci_apply_eaj3(
     $     pbar, np, eap, eapphase,
     $     rbar, nr, ebr, ebrphase,
     $     tbar, nt, ect, ectphase,
     $     jdim, ndim, hdim, odim,
     $     perm, val, hperm, hval)
      implicit none
c     
c     This bundles 3-calls to eaj into one to reduce subroutine call
c     overhead inside the 3-particle coupling coeff generation routine.
c
c     . -p   -r   -t
c     J -> N -> H -> O
c
c     or
c
c     <Ny|Eap|Jv><Hg|Ebr|Ny><Oz|Ect|Hg>
c
c     or
c
c     identity -> jperm/jval -> nperm/nval -> hperm/hval
c
c     See eaj for more details on the algorithm and parameters
c     
      integer pbar, rbar, tbar  ! [input] 
      integer np, nr, nt        ! [input] 
      integer jdim, ndim, hdim, odim ! [input]
      integer eap(jdim,3,*), ebr(ndim,3,*), ect(hdim,3,*) ! [input]
      double precision eapphase(3,*), ebrphase(3,*), ectphase(3,*) ! [input]
      integer perm(*)           ! [scratch] at least hdim
      double precision val(*)   ! [scratch] at least hdim
      integer hperm(odim)       ! [output]
      double precision hval(odim) ! [output]
c     
      integer i, case, type, w, v
      double precision phase
c   
c     <N|Eap|J>
c
      if (np .eq. 2) then
         do i = 1, ndim
            hval(i)  = 0.0d0
            hperm(i) = 1        ! Valid value so that don't need to test on zero
         enddo
      endif
      do case = 1, np
         type = case + np - 1
         phase = eapphase(type,pbar)
         do w = 1, jdim
            v = eap(w,type,pbar)
            hperm(v) = w
            hval(v)  = phase
         enddo
      enddo
c
c     <H|Ebr|N>
c
      if (nr .eq. 2) then
         do i = 1, hdim
            val(i) = 0.0d0
            perm(i) = 1
         enddo
      endif
      do case = 1, nr
         type = case + nr - 1
         phase = ebrphase(type,rbar)
         do w = 1, ndim
            v = ebr(w,type,rbar)
            perm(v) = hperm(w)
            val(v)  = hval(w)*phase
         enddo
      enddo
c
c     <O|Ect|H>
c
      if (nt .eq. 2) then
         do i = 1, odim
            hval(i) = 0.0d0
            hperm(i) = 1
         enddo
      endif
      do case = 1, nt
         type = case + nt - 1
         phase = ectphase(type,tbar)
         do w = 1, hdim
            v = ect(w,type,tbar)
            hperm(v) = perm(w)
            hval(v)  = val(w)*phase
         enddo
      enddo
c     
      end
      
