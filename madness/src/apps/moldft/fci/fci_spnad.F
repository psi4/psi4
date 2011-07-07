*
* $Id: fci_spnad.F,v 1.2 1999/07/28 00:23:42 d3e129 Exp $
*
      subroutine fci_spin_adapt_nvec
     $     (ns, direction, det, ldd, csf, ldc, nvec)
      implicit none
#include "fciinfo.fh"
c     
c     Given NVEC sets of coefficients of either determinants or CSF
c     corresponding to a single orbital occupation with ns open shells,
c     then do the following
c     
c     if (direction == 1)
c     ... project from determinants to CSF with multiplicity multi (fci_init)
c     ... if there is spin contamination norm will not be conserved
c     else if (direction == 2)
c     ... project from CSF with multiplicity multi to determinants
c     ... norm should be conserved
c     else
c     ... invalid value for direction
c     endif
c     
      integer ns                ! [input] No. of open shells
      integer direction         ! [input] 1=det to csf, 2=csf to det
      integer ldd, ldc          ! [input] Leading dimensions
      integer nvec              ! [input] No. of sets of coefficients
      double precision det(ldd,nvec) ! [in/out] coefficients of determinants
      double precision csf(ldc,nvec) ! [in/out] coefficients of csf
c     
      integer iocc(nsmax)       ! CSF walk
      integer mocc(nsmax)       ! determinant walk
      double precision csum(nsmax) ! CSF intermediate spin values
      double precision dsum     ! determinant intermediate Sz values
c
      integer mdet, mcsf, iwalk, icsf, is, mwalk, ivec
      double precision sum, rnum, rden, rs, rm, coeff
c     
      mdet = mnsdet(ns)
      mcsf = mnscsf(ns)
c
      if (direction .eq. 1) then
         call dfill(0.0d0, csf, 1, mcsf)
      else
         call dfill(0.0d0, det, 1, mdet)
      endif
c     
c     Loop thru CSF by looping thru walks on the determinant graph.  
c     Any that only have positive intermediate spin values can be used 
c     as a CSF walk and they will occur in the correct lexical order. 
c     
      iwalk = 0
      do 30 icsf=1,mdet
         call fci_mfind_occ(icsf,iocc,ns)
         sum = 0.0d0
         do 20 is = 1,ns
            sum = sum + 1.5d0 - dfloat(iocc(is))
            if (sum .lt. -0.0001d0) goto 30 ! Not a valid CSF
            csum(is) = sum
 20      continue
         iwalk = iwalk + 1      ! Lexical index of this CSF
c     
         do 40 mwalk = 1,mdet
            call fci_mfind_occ(mwalk,mocc,ns)
            rnum = 1.0d0
            rden = 1.0d0
            dsum = 0.0d0
            do 50 is = 1,ns
               dsum = dsum + 1.5d0 - dfloat(mocc(is))
c
               rs = csum(is)
               rm = dsum
               if (iocc(is).eq.1) then
                  if (mocc(is).eq.1) then
                     rnum = rnum*(rs+rm) ! s=+ m=+
                     rden = rden*2.0d0*rs
                  else
                     rnum = rnum*(rs-rm) ! s=+ m=-
                     rden = rden*2.0d0*rs
                  endif
               else
                  if (mocc(is).eq.1) then
                     rnum = -rnum*(rs-rm+1) ! s=- m=+
                     rden = rden*(2.0d0*rs+2.0d0)
                  else
                     rnum = rnum*(rs+rm+1.0d0) ! s=- m=-
                     rden = rden*(2.0d0*rs+2.0d0)
                  endif
               endif
               if (rnum .eq. 0.0d0) goto 40
 50         continue
c     
c     Coefficient of determinant mwalk in CSF iwalk
c     
            coeff = sign(sqrt(abs(rnum/rden)),rnum)
c     
            if (direction .eq. 1) then
               do ivec = 1, nvec
                  csf(iwalk,ivec) = csf(iwalk,ivec) + 
     $                 coeff*det(mwalk,ivec)
               enddo
            else
               do ivec = 1, nvec
                  det(mwalk,ivec) = det(mwalk,ivec) + 
     $                 coeff*csf(iwalk,ivec)
               enddo
            endif
c
 40      continue
 30   continue
      if (iwalk .ne. mcsf) call fci_err('spnad: bad #csf ', iwalk)
c     
      end
      subroutine fci_spin_adapt(ns, direction, det, csf)
      implicit none
#include "fciinfo.fh"
c     
c     Given a set of coefficients of either determinants or CSF
c     corresponding to a single orbital occupation with ns open shells,
c     then do the following
c     
c     if (direction == 1)
c     ... project from determinants to CSF with multiplicity multi (fci_init)
c     ... if there is spin contamination norm will not be conserved
c     else if (direction == 2)
c     ... project from CSF with multiplicity multi to determinants
c     ... norm should be conserved
c     else
c     ... invalid value for direction
c     endif
c     
      integer ns                ! [input] No. of open shells
      integer direction         ! [input] 1=det to csf, 2=csf to det
      double precision det(*)   ! [in/out] coefficients of determinants
      double precision csf(*)   ! [in/out] coefficients of csf
c     
      integer iocc(nsmax)       ! CSF walk
      integer mocc(nsmax)       ! determinant walk
      double precision csum(nsmax) ! CSF intermediate spin values
      double precision dsum     ! determinant intermediate Sz values
c
      integer mdet, mcsf, iwalk, icsf, is, mwalk
      double precision sum, rnum, rden, rs, rm, coeff
c     
      mdet = mnsdet(ns)
      mcsf = mnscsf(ns)
c
      if (direction .eq. 1) then
         call dfill(0.0d0, csf, 1, mcsf)
      else
         call dfill(0.0d0, det, 1, mdet)
      endif
c     
c     Loop thru CSF by looping thru walks on the determinant graph.  
c     Any that only have positive intermediate spin values can be used 
c     as a CSF walk and they will occur in the correct lexical order. 
c     
      iwalk = 0
      do 30 icsf=1,mdet
         call fci_mfind_occ(icsf,iocc,ns)
         sum = 0.0d0
         do 20 is = 1,ns
            sum = sum + 1.5d0 - dfloat(iocc(is))
            if (sum .lt. -0.0001d0) goto 30 ! Not a valid CSF
            csum(is) = sum
 20      continue
         iwalk = iwalk + 1      ! Lexical index of this CSF
c     
         do 40 mwalk = 1,mdet
            call fci_mfind_occ(mwalk,mocc,ns)
            rnum = 1.0d0
            rden = 1.0d0
            dsum = 0.0d0
            do 50 is = 1,ns
               dsum = dsum + 1.5d0 - dfloat(mocc(is))
c
               rs = csum(is)
               rm = dsum
               if (iocc(is).eq.1) then
                  if (mocc(is).eq.1) then
                     rnum = rnum*(rs+rm) ! s=+ m=+
                     rden = rden*2.0d0*rs
                  else
                     rnum = rnum*(rs-rm) ! s=+ m=-
                     rden = rden*2.0d0*rs
                  endif
               else
                  if (mocc(is).eq.1) then
                     rnum = -rnum*(rs-rm+1) ! s=- m=+
                     rden = rden*(2.0d0*rs+2.0d0)
                  else
                     rnum = rnum*(rs+rm+1.0d0) ! s=- m=-
                     rden = rden*(2.0d0*rs+2.0d0)
                  endif
               endif
               if (rnum .eq. 0.0d0) goto 40
 50         continue
c     
c     Coefficient of determinant mwalk in CSF iwalk
c     
            coeff = sign(sqrt(abs(rnum/rden)),rnum)
c     
            if (direction .eq. 1) then
               csf(iwalk) = csf(iwalk) + coeff*det(mwalk)
            else
               det(mwalk) = det(mwalk) + coeff*csf(iwalk)
            endif
c
 40      continue
 30   continue
      if (iwalk .ne. mcsf) call fci_err('spnad: bad #csf ', iwalk)
c     
      end
