!
module f90_kind
  implicit none
  integer, parameter :: single = kind(0.0)
  integer, parameter :: double = kind(0.0d0)
end module f90_kind
!
!
!
module dkh_main 
!
!-----------------------------------------------------------------------
      USE f90_kind,                           ONLY: double
!-----------------------------------------------------------------------
      IMPLICIT NONE 
      PUBLIC                                 :: DKH
!-----------------------------------------------------------------------
      CONTAINS
!-----------------------------------------------------------------------
!     NAME                                                             *
!     dkh_main                                                         *
!                                                                      *
!     FUNCTION                                                         *
!     2nd to 4th order DKH calculation                                 *
!                                                                      *
!     THIS SOURCE HAS BEEN DOCKED TO: (as of 03/2012)                  *
!     at35 (original f77 DKH1-DKH5 development environment 2002)       *
!     MOLCAS/MOLPRO (2003, later replaced by arbitrary-order DKH 2006) *
!     ORCA (given to Frank Neese, end of 2003)                         *
!     cp2k (original f77 code changed to f90 by                        *
!                  Jens Thar/Barbara Kirchner in 2006/2007)            *
!     ACES2/ACES3 (given to student of Peter Schwerdtfeger/Rod Bartlett*
!                  in about 2009)                                      *
!                                                                      *
!     AUTHORS                                                          *
!     written by:                                                      *
!     Alexander Wolf, Markus Reiher, Bernd Hess (Erlangen, 2002)       *
!                                                                      *
!     2nd to 4th order ported to F90 for cp2k by:                      *
!     Jens Thar, Barbara Kirchner: Uni Bonn (09/2006)                  *
!     final polishing by:                                              *                                            
!     Markus Reiher: ETH Zurich (03/2012)                              *
!                                                                      *
!     REFERENCES:                                                      *
!  M. Reiher, A. Wolf, J. Chem. Phys. 121 (2004) 10944-10956           *
!  A. Wolf, M. Reiher, B. A. Hess, J. Chem. Phys. 117 (2002) 9215-9226 *
!                                                                      *
!-----------------------------------------------------------------------
!
     SUBROUTINE DKH(s,v,h,pVp,n,dkh_order) 
!
!-----------------------------------------------------------------------
!                                                                      *
!  INPUT:                                                              *
!    n          Number of primitive gaussians                          *
!    s    (:,:) overlap matrix                                         *
!    pVp  (:,:) pVp matrix                                             *
!    dkh_order  order of scalar-relativistic DKH calc \in [1,5]        *
!                                                                      *
!  IN_OUT:                                                             *
!    v    (:,:) input: nonrelativistic potential energy matrix         *
!               output: (ev1+ev2)                                      *
!    h    (:,:) input: kinetic energy matrix                           *
!               output: kinetic+ext-potential part of Hamiltonian in   *
!                       position space                                 *
!                                                                      *
!  INTERNAL                                                            *
!    sinv (:,:) inverted, orthogonalized overlap matrix                *
!    ev0t (:)   DKH-even0 matrix in T-basis                            *
!    e    (:)   e=SQRT(p^2c^2+c^4)                                     *
!    eig  (:,:) eigenvectors of sinv' h sinv                           *
!    tt   (:)   eigenvalues of sinv' h sinv                            *
!    revt (:,:) reverse transformation matrix T-basis -> position space*
!    aa   (:)   kinematical factors f. DKH DSQRT((c^2+e(i))/(2.0*e(i)))*
!    rr   (:)   kinematical factors f. DKH c/(c^2+e(i))                *
!    vt   (:,:) non relativistic potential matrix in T-basis           *
!    pvpt (:,:) pvp integral matrix in T-basis                         * 
!    ev1t (:,:) DKH-even1 matrix in T-basis                            *
!    evt2 (:,:) DKH-even2 matrix in T-basis                            *
!    ev1  (:,:) DKH-even1 matrix in position space                     *
!    ev2  (:,:) DKH-even2 matrix in position space                     *
!    ove (:,:) scratch                                                 *
!    aux (:,:) scratch                                                 *
!    velit  velocity of light 137 a.u.                                 *
!    prea   prefactor, 1/137^2                                         *
!    con2   prefactor, 2/137^2                                         *
!    con    prefactor, 137^2                                           *
!-----------------------------------------------------------------------
!      
      INTEGER,INTENT(IN)                        :: n
      REAL(KIND=double),DIMENSION(n,n),INTENT(INOUT)::  s,v,pVp,h 
      REAL(KIND=double),DIMENSION(:),ALLOCATABLE    ::  ev0t,e,aa,rr,&
                                                    tt
      Real(KIND=double),DIMENSION(:,:), &
      ALLOCATABLE                               ::  eig,sinv,revt,aux,&
                                                    ove,ev1t,pev1tp,&
                                                    vt,pVpt,ev1,ev2,ev2t,&
                                                    ev3,ev4,ev3t,ev4t  
      REAL(KIND=double)                             ::  velit
!CAW  pp: p^2-values (in momentum-space), stored as matrix!!
      INTEGER                                   :: i,j,k,dkh_order
!
!-----------------------------------------------------------------------
!     Define velocity of light
!-----------------------------------------------------------------------
!
      velit = 137.0359895d0
!
!-----------------------------------------------------------------------
!     Allocate some matrices
!-----------------------------------------------------------------------
!
      ALLOCATE(eig(n,n))
      ALLOCATE(sinv(n,n))
      ALLOCATE(revt(n,n))
      ALLOCATE(aux(n,n))
      ALLOCATE(ove(n,n)) 
      ALLOCATE(ev0t(n))
      ALLOCATE(e(n))
      ALLOCATE(aa(n))
      ALLOCATE(rr(n))
      ALLOCATE(tt(n))
      ALLOCATE(ev1t(n,n))
      ALLOCATE(ev2t(n,n))
      ALLOCATE(ev3t(n,n))
      ALLOCATE(ev4t(n,n))
      ALLOCATE(vt(n,n))
      ALLOCATE(pVpt(n,n))
      ALLOCATE(pev1tp(n,n)) 
      ALLOCATE(ev1(n,n))
      ALLOCATE(ev2(n,n))
      ALLOCATE(ev3(n,n))
      ALLOCATE(ev4(n,n))
!
!-----------------------------------------------------------------------
!     Check (and set) DKH order
!-----------------------------------------------------------------------
!
      if (dkh_order.ge.5.or.dkh_order.lt.2) dkh_order=4  
!
!-----------------------------------------------------------------------
!     Schmidt-orthogonalize overlap matrix
!-----------------------------------------------------------------------
!
      CALL sog (n,s,sinv)      
!
!-----------------------------------------------------------------------
!     Calculate matrix representation from nonrelativistic T matrix
!-----------------------------------------------------------------------
!
      CALL diag ( h,n,eig,tt,sinv,aux,0 )
!
!-----------------------------------------------------------------------
!     Calculate kinetic part of Hamiltonian in T-basis 
!-----------------------------------------------------------------------
!
      CALL kintegral (n,ev0t,tt,e,velit)
!
!-----------------------------------------------------------------------
!     Calculate reverse transformation matrix revt
!-----------------------------------------------------------------------
!
      CALL dgemm("N","N",n,n,n,1.0d0,sinv,n,eig,n,0.0d0,aux,n)
      CALL dgemm("N","N",n,n,n,1.0d0,s,n,aux,n,0.0d0,revt,n)      
!
!-----------------------------------------------------------------------
!     Calculate kinetic part of the Hamiltonian  
!-----------------------------------------------------------------------
!
      h = 0.0d0
      DO i=1,n
        DO j=1,i
          DO k=1,n
            h(i,j)=h(i,j)+revt(i,k)*revt(j,k)*ev0t(k)
            h(j,i)=h(i,j)
          END DO
        END DO
      END DO
!
!-----------------------------------------------------------------------
!     Calculate kinematical factors for DKH
!-----------------------------------------------------------------------
!
      DO i=1,n
        aa(i)=DSQRT((velit*velit+e(i)) / (2.0d0*e(i)))            
        rr(i)=DSQRT(velit*velit)/(velit*velit+e(i))                        
      END DO
!
!----------------------------------------------------------------------- 
!     Transform v integrals to T-basis (v -> vt)
!-----------------------------------------------------------------------
!
      CALL trsm(v,sinv,ove,n,aux)
      CALL trsm(ove,eig,vt,n,aux)
!
!-----------------------------------------------------------------------
!     Transform pVp integrals to T-basis (pVp -> pVpt)
!-----------------------------------------------------------------------

      CALL trsm(pVp,sinv,ove,n,aux)
      CALL trsm(ove,eig,pVpt,n,aux)
!
!-----------------------------------------------------------------------
!     Calculate even1 in T-basis
!-----------------------------------------------------------------------
! 
      CALL even1(n,ev1t,vt,pvpt,aa,rr)
!
!----------------------------------------------------------------------
!     Transform even1 back to position space
!----------------------------------------------------------------------
!
      CALL dgemm("N","N",n,n,n,1.0d0,revt,n,ev1t,n,0.0d0,aux,n)
      CALL dgemm("N","T",n,n,n,1.0d0,aux,n,revt,n,0.0d0,ev1,n)  
!
!-----------------------------------------------------------------------
!     Calculate even2 in T-basis
!-----------------------------------------------------------------------
!
      CALL even2c (n,ev2t,vt,pvpt,aa,rr,tt,e)
!
!-----------------------------------------------------------------------
!     Transform even2 back to position space
!-----------------------------------------------------------------------
      aux=0.0d0
      CALL dgemm("N","N",n,n,n,1.0d0,revt,n,ev2t,n,0.0d0,aux,n)
      CALL dgemm("N","T",n,n,n,1.0d0,aux,n,revt,n,0.0d0,ev2,n)
!
!-----------------------------------------------------------------------
!     Calculate even3 in T-basis, only if requested
!-----------------------------------------------------------------------
!
      IF (dkh_order.ge.3) THEN       
        CALL peven1p(n,pev1tp,vt,pvpt,aa,rr,tt)
        CALL even3b(n,ev3t,ev1t,pev1tp,vt,pvpt,aa,rr,tt,e)
!
!-----------------------------------------------------------------------
!     Transform even3 back to position space
!-----------------------------------------------------------------------
!
        aux=0.0d0
        CALL dgemm("N","N",n,n,n,1.0d0,revt,n,ev3t,n,0.0d0,aux,n)
        CALL dgemm("N","T",n,n,n,1.0d0,aux,n,revt,n,0.0d0,ev3,n)
!
!-----------------------------------------------------------------------
!     Calculate even4 in T-basis, only if requested
!-----------------------------------------------------------------------
!
        IF (dkh_order.ge.4) THEN
          CALL even4a(n,ev4t,ev1t,pev1tp,vt,pvpt,aa,rr,tt,e)
!
!-----------------------------------------------------------------------
!     Transform even4 back to position space
!-----------------------------------------------------------------------
!
          aux=0.0d0
          CALL dgemm("N","N",n,n,n,1.0d0,revt,n,ev4t,n,0.0d0,aux,n)
          CALL dgemm("N","T",n,n,n,1.0d0,aux,n,revt,n,0.0d0,ev4,n)
        END IF
      END IF
!          
!-----------------------------------------------------------------------
!     Calculate v in position space
!-----------------------------------------------------------------------
!
      CALL mat_add(v,1.0d0,ev1,1.0d0,ev2,n)
      IF(dkh_order.ge.3) THEN
        CALL mat_add(v,1.0d0,v,1.0d0,ev3,n)
        IF(dkh_order.ge.4) THEN
          CALL mat_add(v,1.0d0,v,1.0d0,ev4,n)
        END IF
      END IF
!
!-----------------------------------------------------------------------
!
      DEALLOCATE(eig,sinv,revt,ove,aux,vt,pVpt,ev1,ev2,ev3,ev4,ev1t,ev2t,ev3t,ev4t,pev1tp) 
      DEALLOCATE(ev0t,e,aa,rr,tt)
!
      RETURN
      END SUBROUTINE DKH
!
!-----------------------------------------------------------------------
!
      SUBROUTINE write_vec (v)
      IMPLICIT NONE
      REAL(KIND=double),DIMENSION(:),INTENT(IN)    :: v
      INTEGER                                      :: i
!
      DO i=1,12
        write(*,999)v(i)
      END DO
      DO i=1,12
        write(*,999)v(i+33)
      END DO
999   format(12es12.4)
      RETURN
      END SUBROUTINE write_vec
!
!-----------------------------------------------------------------------
!
      SUBROUTINE write_mat (m)
      IMPLICIT NONE
      REAL(KIND=double),DIMENSION(:,:),INTENT(IN)  :: m
      INTEGER                                      :: i
!
      DO i=1,12
        write(*,999)m(i,1:12)
      END DO
      DO i=1,12
        write(*,999)m(i+33,34:45)
      END DO
999   format(12es12.4)
      RETURN
      END SUBROUTINE write_mat
!
!-----------------------------------------------------------------------
!

      SUBROUTINE kintegral (n,ev0t,tt,e,velit)
      IMPLICIT NONE
      INTEGER,INTENT(IN)                              :: n
      REAL(KIND =double),DIMENSION(:),INTENT(OUT)         :: ev0t
      REAL(KIND =double),DIMENSION(:),INTENT(IN)          :: tt
      REAL(KIND =double),DIMENSION(:),INTENT(OUT)         :: e
      REAL(KIND =double),INTENT(IN)                       :: velit
      REAL(KIND =double)           :: ratio,tv1,tv2,tv3,tv4,prea,con,con2
      INTEGER                  :: i
!
      DO i=1,n
        IF (tt(i).LT.0.0d0) THEN
          write(*,*) ' dkh_main.F | tt(',i,') = ',tt(i)
        END IF
        prea=1/(velit*velit)
        con2=prea+prea
        con=1.0d0/prea
!BAH 2000:
!       If T is sufficiently small, use series expansion to avoid
!       cancellation, otherwise calculate SQRT directly
!
        ev0t(i)=tt(i)
        ratio=tt(i)/velit
        IF (ratio.LE.0.02d0) THEN
          tv1=tt(i)
          tv2=-tv1*tt(i)*prea/2.0d0
          tv3=-tv2*tt(i)*prea
          tv4=-tv3*tt(i)*prea*1.25d0
          ev0t(i)=tv1+tv2+tv3+tv4
        ELSE
          ev0t(i)=con*(DSQRT(1.0d0+con2*tt(i))-1.0d0)   
        END IF
        e(i)=ev0t(i)+con           
      END DO
!
      RETURN
      END SUBROUTINE kintegral
!
!-----------------------------------------------------------------------
!
      SUBROUTINE even1(n,ev1t,vt,pvpt,aa,rr)
!
!-----------------------------------------------------------------------
!                                                                      -
!     1st order DKH-approximation                                      -
!                                                                      -
!     n    in   dimension of matrices                                  -
!     ev1t out  even1 output matrix                                    -
!     vt   in   potential matrix v in T-space                          -
!     pvpt in   pvp matrix in T-space                                  -
!     aa   in   A-factors (diagonal)                                   -
!     rr   in   R-factors (diagonal)                                   -
!                                                                      -
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
      INTEGER,INTENT(IN)                          :: n
      REAL(KIND=double),DIMENSION(:,:),INTENT(IN) :: vt,pvpt
      REAL(KIND=double),DIMENSION(:,:),INTENT(OUT):: ev1t
      REAL(KIND=double),DIMENSION(:),INTENT(IN)   :: aa,rr
      INTEGER                                     :: i,j
!
!-----------------------------------------------------------------------
!
      DO i=1,n
        DO j=1,i
          ev1t(i,j)=vt(i,j)*aa(i)*aa(j)+pVpt(i,j)*aa(i)*rr(i)*aa(j)*rr(j)
          ev1t(j,i)=ev1t(i,j)
        END DO
      END DO
!
      RETURN
      END SUBROUTINE even1
!
!-----------------------------------------------------------------------
!
      SUBROUTINE peven1p(n,pev1tp,vt,pvpt,aa,rr,tt)
!
!-----------------------------------------------------------------------
!                                                                      -
!     1st order DKH-approximation                                      -
!                                                                      -
!     n      in   dimension of matrices                                -
!     pev1tp out  peven1p output matrix                                -
!     vt     in   potential matrix v in T-space                        -
!     pvpt   in   pvp matrix in T-space                                -
!     aa     in   A-factors (diagonal)                                 -
!     rr     in   R-factors (diagonal)                                 -
!                                                                      -
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
      INTEGER,INTENT(IN)                          :: n
      REAL(KIND=double),DIMENSION(:,:),INTENT(IN) :: vt,pvpt
      REAL(KIND=double),DIMENSION(:,:),INTENT(OUT):: pev1tp
      REAL(KIND=double),DIMENSION(:),INTENT(IN)   :: aa,rr,tt
      INTEGER                                     :: i,j
!
!-----------------------------------------------------------------------
!
      DO i=1,n
        DO j=1,i
          pev1tp(i,j)=4.0d0*vt(i,j)*aa(i)*aa(j)*rr(i)*rr(i)*rr(j)*rr(j)*tt(i)*tt(j)+&
          pVpt(i,j)*aa(i)*rr(i)*aa(j)*rr(j)
          pev1tp(j,i)=pev1tp(i,j)
        END DO
      END DO

      RETURN
      END SUBROUTINE peven1p
!
!-----------------------------------------------------------------------
!
      SUBROUTINE even2c (n,ev2,vv,gg,aa,rr,tt,e)
!
!***********************************************************************
!                                                                      *
!     Alexander Wolf, last modified: 20.02.2002 - DKH2                 *
!                                                                      *
!     2nd order DK-approximation ( original DK-transformation with     *
!                                       U = SQRT(1+W^2) + W        )   *
!                                                                      *
!     Version: 1.1  (20.2.2002) :  Usage of SR mat_add included        *
!              1.0  (6.2.2002)                                         *
!     Modification history:                                            *
!     30.09.2006 Jens Thar: deleted obsolete F77 memory manager        *
!                                                                      *
!     ev2 = 1/2 [W1,O1]                                                *
!                                                                      *
!         ----  Meaning of Parameters  ----                            *
!                                                                      *
!     n       in   Dimension of matrices                               *
!     ev2     out  even2 output matrix = final result                  *
!     vv      in   potential v                                         *
!     gg      in   pvp                                                 *
!     aa      in   A-Factors (DIAGONAL)                                *
!     rr      in   R-Factors (DIAGONAL)                                *
!     tt      in   Nonrel. kinetic Energy (DIAGONAL)                   *
!     e       in   Rel. Energy = SQRT(p^2*c^2 + c^4)  (DIAGONAL)       *
!                                                                      *
!***********************************************************************
!
      IMPLICIT NONE 
      INTEGER,INTENT(IN)                       :: n
      INTEGER                                  :: i,j
      REAL(KIND=double),DIMENSION(:),INTENT(IN)    :: aa,rr,tt,e
      REAL(KIND=double),DIMENSION(:,:),INTENT(OUT) :: ev2
      REAL(KIND=double),DIMENSION(:,:),INTENT(IN)  :: vv,gg 
      REAL(KIND=double),DIMENSION(:,:),ALLOCATABLE ::v,pvp,vh,pvph,&
                                                 w1o1,o1w1
!
!-----------------------------------------------------------------------
!     1.   General Structures and Patterns for DKH2                      
!-----------------------------------------------------------------------
!
      ALLOCATE(v(n,n))
      ALLOCATE(pVp(n,n))
      ALLOCATE(vh(n,n))
      ALLOCATE(pVph(n,n))
      v=0.0d0
      pVp=0.0d0
      vh=0.0d0
      pVph=0.0d0
      v(1:n,1:n)=vv(1:n,1:n)
      vh(1:n,1:n)=vv(1:n,1:n)
      pvp(1:n,1:n)=gg(1:n,1:n)
      pvph(1:n,1:n)=gg(1:n,1:n)
      ev2=0.0d0
!  Calculate  v = A V A:
     CALL mat_axa(v,n,aa)
!  Calculate  pvp = A P V P A:
     CALL mat_arxra(pvp,n,aa,rr)
!  Calculate  vh = A V~ A:
     CALL mat_1_over_h(vh,n,e)
     CALL mat_axa(vh,n,aa)
!  Calculate  pvph = A P V~ P A:
     CALL mat_1_over_h(pvph,n,e)
     CALL mat_arxra(pvph,n,aa,rr)
!  Create/initialize necessary matrices:
     ALLOCATE(w1o1(n,n))
     ALLOCATE(o1w1(n,n))
     w1o1=0.0d0
     o1w1=0.0d0
!  Calculate w1o1:
     CALL dgemm("N","N",n,n,n,-1.0d0,pvph,n,v,n,0.0d0,w1o1,n)
     CALL mat_muld(w1o1,pvph,pvp,n,  1.0d0,1.0d0,tt,rr)
     CALL mat_mulm(w1o1,vh,  v,n,    1.0d0,1.0d0,tt,rr)
     CALL dgemm("N","N",n,n,n,-1.0d0,vh,n,pvp,n,1.0d0,w1o1,n)
!  Calculate o1w1:
     CALL dgemm("N","N",n,n,n,1.0d0,pvp,n,vh,n,0.0d0,o1w1,n)
     CALL mat_muld(o1w1,pvp,pvph,n,  -1.0d0,1.0d0,tt,rr)
     CALL mat_mulm(o1w1,v,  vh,n,    -1.0d0,1.0d0,tt,rr)
     CALL dgemm("N","N",n,n,n,1.0d0,v,n,pvph,n,1.0d0,o1w1,n)
!  Calculate in symmetric pakets
!-----------------------------------------------------------------------
!     2.   1/2 [W1,O1] = 1/2 W1O1 -  1/2 O1W1  
!-----------------------------------------------------------------------
!
      CALL mat_add (ev2,0.5d0,w1o1,-0.5d0,o1w1,n)
!
!-----------------------------------------------------------------------
!     3.   Finish up 
!-----------------------------------------------------------------------
!
      DEALLOCATE(v,vh,pvp,pvph,w1o1,o1w1)
!
      RETURN
      END SUBROUTINE even2c
!
!-----------------------------------------------------------------------
!
      SUBROUTINE even3b (n,ev3,e1,pe1p,vv,gg,aa,rr,tt,e)
!
!***********************************************************************
!                                                                      *
!     Alexander Wolf, last modified:  20.2.2002 - DKH3                 *
!                                                                      *
!     3rd order DK-approximation (generalised DK-transformation)       *
!                                                                      *
!     Version: 1.1  (20.2.2002) :  Usage of SR mat_add included        *
!              1.0  (7.2.2002)                                         *
!                                                                      *
!     ev3 = 1/2 [W1,[W1,E1]]                                           *
!                                                                      *
!     Modification history:                                            *
!     30.09.2006 Jens Thar: deleted obsolete F77 memory manager        *
!                                                                      *
!         ----  Meaning of Parameters  ----                            *
!                                                                      *
!     n       in   Dimension of matrices                               *
!     ev3     out  even3 output matrix = final result                  *
!     e1      in   E1 = even1-operator                                 *
!     pe1p    in   pE1p                                                *
!     vv      in   potential v                                         *
!     gg      in   pvp                                                 *
!     aa      in   A-Factors (DIAGONAL)                                *
!     rr      in   R-Factors (DIAGONAL)                                *
!     tt      in   Nonrel. kinetic Energy (DIAGONAL)                   *
!     e       in   Rel. Energy = SQRT(p^2*c^2 + c^4)  (DIAGONAL)       *
!                                                                      *
!***********************************************************************
!
      IMPLICIT NONE
      INTEGER,INTENT(IN)                       :: n
      INTEGER                                  :: i,j
      REAL(KIND=double),DIMENSION(:),INTENT(IN)    :: aa,rr,tt,e
      REAL(KIND=double),DIMENSION(:,:),INTENT(OUT) :: ev3
      REAL(KIND=double),DIMENSION(:,:),INTENT(IN)  :: vv,gg,e1,pe1p
      REAL(KIND=double),DIMENSION(:,:),ALLOCATABLE ::vh,pvph,&
                                                 w1w1,w1e1w1,scr_1,scr_2
!
!-----------------------------------------------------------------------
!     1.   General Structures and Patterns for DKH3
!-----------------------------------------------------------------------
!
      ALLOCATE(vh(n,n))
      ALLOCATE(pVph(n,n))
      vh=0.0d0
      pVph=0.0d0
      vh(1:n,1:n)=vv(1:n,1:n)
      pvph(1:n,1:n)=gg(1:n,1:n)
      ev3=0.0d0
!  Calculate  vh = A V~ A:
     CALL mat_1_over_h(vh,n,e)
     CALL mat_axa(vh,n,aa)
!  Calculate  pvph = A P V~ P A:
     CALL mat_1_over_h(pvph,n,e)
     CALL mat_arxra(pvph,n,aa,rr)
!  Create/Initialize necessary matrices:
      ALLOCATE(w1w1(n,n))
      ALLOCATE(w1e1w1(n,n))
      ALLOCATE(scr_1(n,n))
      ALLOCATE(scr_2(n,n))
      w1w1=0.0d0
      w1e1w1=0.0d0
      scr_1=0.0d0
      scr_2=0.0d0
!  Calculate w1w1:
      CALL dgemm("N","N",n,n,n,1.0d0,pvph,n,vh,n,0.0d0,w1w1,n)
      CALL mat_muld(w1w1,pvph,pvph,n,-1.0d0,1.0d0,tt,rr)
      CALL mat_mulm(w1w1,vh,  vh,n,  -1.0d0,1.0d0,tt,rr)
      CALL dgemm("N","N",n,n,n,1.0d0,vh,n,pvph,n,1.0d0,w1w1,n)
!  Calculate w1e1w1:
      CALL mat_muld(scr_1 ,pvph ,pe1p,n, 1.0d0,0.0d0,tt,rr)
      CALL dgemm("N","N",n,n,n,1.0d0,vh,n,pe1p,n,0.0d0,scr_2,n)
      CALL dgemm("N","N",n,n,n,1.0d0,scr_1,n,vh,n,0.0d0,w1e1w1,n)
      CALL mat_muld(w1e1w1,scr_1,pvph,n, -1.0d0,1.0d0,tt,rr)
      CALL dgemm("N","N",n,n,n,-1.0d0,scr_2,n,vh,n,1.0d0,w1e1w1,n)
      CALL mat_muld(w1e1w1,scr_2,pvph,n, 1.0d0,1.0d0,tt,rr)
!-----------------------------------------------------------------------
!     2.   ev3 = 1/2 (W1^2)E1 + 1/2 E1(W1^2) - W1E1W1
!-----------------------------------------------------------------------
      CALL dgemm("N","N",n,n,n,0.5d0,w1w1,n,e1,n,0.0d0,ev3,n)
      CALL dgemm("N","N",n,n,n,0.5d0,e1,n,w1w1,n,1.0d0,ev3,n)
      CALL mat_add (ev3,1.0d0,ev3,-1.0d0,w1e1w1,n)
!-----------------------------------------------------------------------
!     3.   Finish up the stuff!!
!-----------------------------------------------------------------------
      DEALLOCATE(vh,pvph,w1w1,w1e1w1,scr_1,scr_2)
!
      RETURN
      END SUBROUTINE even3b
!
!-----------------------------------------------------------------------
!
      SUBROUTINE even4a (n,ev4,e1,pe1p,vv,gg,aa,rr,tt,e)
!
!***********************************************************************
!                                                                      *
!     Alexander Wolf,   last modified: 25.02.2002   --   DKH4          *
!                                                                      *
!     4th order DK-approximation (scalar = spin-free)                  *
!                                                                      *
!     Version: 1.2  (25.2.2002) :  Elegant (short) way of calculation  *
!                                  included                            *
!              1.1  (20.2.2002) :  Usage of SR mat_add included        *
!              1.0  (8.2.2002)                                         *
!                                                                      *
!     ev4  =  1/2 [W2,[W1,E1]] + 1/8 [W1,[W1,[W1,O1]]]  =              *
!                                                                      *
!          =      sum_1        +         sum_2                         *
!                                                                      *
!                                                                      *
!     Modification history:                                            *
!     30.09.2006 Jens Thar: deleted obsolete F77 memory manager        *
!                                                                      *
!         ----  Meaning of Parameters  ----                            *
!                                                                      *
!     n       in   Dimension of matrices                               *
!     ev4     out  even4 output matrix = final result                  *
!     e1     in   E1                                                   *
!     pe1p   in   p(E1)p                                               *
!     vv      in   potential v                                         *
!     gg      in   pvp                                                 *
!     aa      in   A-Factors (DIAGONAL)                                *
!     rr      in   R-Factors (DIAGONAL)                                *
!     tt      in   Nonrel. kinetic Energy (DIAGONAL)                   *
!     e       in   Rel. Energy = SQRT(p^2*c^2 + c^4)  (DIAGONAL)       *
!                                                                      *
!***********************************************************************
!
      IMPLICIT NONE 
      INTEGER,INTENT(IN)                       :: n
      INTEGER                                  :: i,j
      REAL(KIND=double),DIMENSION(:),INTENT(IN)    :: aa,rr,tt,e
      REAL(KIND=double),DIMENSION(:,:),INTENT(OUT) :: ev4
      REAL(KIND=double),DIMENSION(:,:),INTENT(IN)  :: vv,gg,e1,pe1p
      REAL(KIND=double),DIMENSION(:,:),ALLOCATABLE ::v,pvp,vh,pvph,&
                                                 w1w1,w1o1,o1w1,scr_1,scr_2,&
                                                 scr_3,scr_4,scrh_1,&
                                                 scrh_2,scrh_3,scrh_4,&
                                                 sum_1,sum_2
!
!-----------------------------------------------------------------------
!     1.   General Structures and Patterns for DKH4
!-----------------------------------------------------------------------
      ALLOCATE(v(n,n))
      ALLOCATE(pVp(n,n))
      ALLOCATE(vh(n,n))
      ALLOCATE(pVph(n,n))
      v=0.0d0
      pVp=0.0d0
      vh=0.0d0
      pVph=0.0d0
      v(1:n,1:n)=vv(1:n,1:n)
      vh(1:n,1:n)=vv(1:n,1:n)
      pvp(1:n,1:n)=gg(1:n,1:n)
      pvph(1:n,1:n)=gg(1:n,1:n)
      ev4=0.0d0
!  Calculate  v = A V A:
     CALL mat_axa(v,n,aa)
!  Calculate  pvp = A P V P A:
     CALL mat_arxra(pvp,n,aa,rr)
!  Calculate  vh = A V~ A:
     CALL mat_1_over_h(vh,n,e)
     CALL mat_axa(vh,n,aa)
!  Calculate  pvph = A P V~ P A:
     CALL mat_1_over_h(pvph,n,e)
     CALL mat_arxra(pvph,n,aa,rr)
!  Create/initialize necessary matrices:
      ALLOCATE(w1w1(n,n))
      w1w1 = 0.0d0
      ALLOCATE(w1o1(n,n))
      w1o1 = 0.0d0
      ALLOCATE(o1w1(n,n))
      o1w1 = 0.0d0
      ALLOCATE(sum_1(n,n))
      sum_1 = 0.0d0
      ALLOCATE(sum_2(n,n))
      sum_2 = 0.0d0
      ALLOCATE(scr_1(n,n))
      scr_1 = 0.0d0
      ALLOCATE(scr_2(n,n))
      scr_2 = 0.0d0
      ALLOCATE(scr_3(n,n))
      scr_3 = 0.0d0
      ALLOCATE(scr_4(n,n))
      scr_4 = 0.0d0
      ALLOCATE(scrh_1(n,n))
      scrh_1 = 0.0d0
      ALLOCATE(scrh_2(n,n))
      scrh_2 = 0.0d0
      ALLOCATE(scrh_3(n,n))
      scrh_3 = 0.0d0
      ALLOCATE(scrh_4(n,n))
      scrh_4 = 0.0d0
!  Calculate w1w1:
      CALL dgemm("N","N",n,n,n,1.0d0,pvph,n,vh,n,0.0d0,w1w1,n)
      CALL mat_muld(w1w1,pvph,pvph,n, -1.0d0,1.0d0,tt,rr)
      CALL mat_mulm(w1w1,vh,  vh,n,   -1.0d0,1.0d0,tt,rr)
      CALL dgemm("N","N",n,n,n,1.0d0,vh,n,pvph,n,1.0d0,w1w1,n)
!  Calculate w1o1:
      CALL dgemm("N","N",n,n,n,-1.0d0,pvph,n,v,n,0.0d0,w1o1,n)
      CALL mat_muld(w1o1,pvph,pvp,n,  1.0d0,1.0d0,tt,rr)
      CALL mat_mulm(w1o1,vh,  v,n,    1.0d0,1.0d0,tt,rr)
      CALL dgemm("N","N",n,n,n,-1.0d0,vh,n,pvp,n,1.0d0,w1o1,n)
!  Calculate o1w1:
      CALL dgemm("N","N",n,n,n,1.0d0,pvp,n,vh,n,0.0d0,o1w1,n)
      CALL mat_muld(o1w1,pvp,pvph,n,  -1.0d0,1.0d0,tt,rr)
      CALL mat_mulm(o1w1,v,  vh,n,    -1.0d0,1.0d0,tt,rr)
      CALL dgemm("N","N",n,n,n,1.0d0,v,n,pvph,n,1.0d0,o1w1,n)
!-----------------------------------------------------------------------
!   2. sum_1 = 1/2 [W2,[W1,E1]] = 1/2 (W2W1E1 - W2E1W1 - W1E1W2 + E1W1W2)
!-----------------------------------------------------------------------
      CALL dgemm("N","N",n,n,n,1.0d0,vh,n,e1,n,0.0d0,scr_1,n)
      CALL dgemm("N","N",n,n,n,1.0d0,pvph,n,e1,n,0.0d0,scr_2,n)
      CALL dgemm("N","N",n,n,n,1.0d0,pe1p,n,vh,n,0.0d0,scr_3,n)
      CALL mat_muld(scr_4, pe1p,pvph,n,1.0d0,0.0d0,tt,rr)
      CALL mat_muld(scrh_1,pvph,pe1p,n,1.0d0,0.0d0,tt,rr)
      CALL mat_1_over_h(scrh_1,n,e)
      CALL dgemm("N","N",n,n,n,1.0d0,vh,n,pe1p,n,0.0d0,scrh_2,n)
      CALL mat_1_over_h(scrh_2,n,e)
      CALL dgemm("N","N",n,n,n,1.0d0,e1,n,pvph,n,0.0d0,scrh_3,n)
      CALL mat_1_over_h(scrh_3,n,e)
      CALL dgemm("N","N",n,n,n,1.0d0,e1,n,vh,n,0.0d0,scrh_4,n)
      CALL mat_1_over_h(scrh_4,n,e)
      CALL dgemm("N","N",n,n,n,0.5d0,scrh_1,n,scr_1,n,0.0d0,sum_1,n)
      CALL mat_muld(sum_1,scrh_1,scr_2,n,-0.5d0,1.0d0,tt,rr)
      CALL dgemm("N","N",n,n,n,-0.5d0,scrh_2,n,scr_1,n,1.0d0,sum_1,n)
      CALL mat_muld(sum_1,scrh_2,scr_2,n, 0.5d0,1.0d0,tt,rr)
      CALL dgemm("N","N",n,n,n,-0.5d0,scrh_3,n,scr_1,n,1.0d0,sum_1,n)
      CALL mat_muld(sum_1,scrh_3,scr_2,n, 0.5d0,1.0d0,tt,rr)
      CALL mat_mulm(sum_1,scrh_4,scr_1,n, 0.5d0,1.0d0,tt,rr)
      CALL dgemm("N","N",n,n,n,-0.5d0,scrh_4,n,scr_2,n,1.0d0,sum_1,n)
      CALL mat_muld(sum_1,scrh_1,scr_3,n,-0.5d0,1.0d0,tt,rr)
      CALL mat_muld(sum_1,scrh_1,scr_4,n, 0.5d0,1.0d0,tt,rr)
      CALL mat_muld(sum_1,scrh_2,scr_3,n, 0.5d0,1.0d0,tt,rr)
      CALL mat_muld(sum_1,scrh_2,scr_4,n,-0.5d0,1.0d0,tt,rr)
      CALL mat_muld(sum_1,scrh_3,scr_3,n, 0.5d0,1.0d0,tt,rr)
      CALL mat_muld(sum_1,scrh_3,scr_4,n,-0.5d0,1.0d0,tt,rr)
      CALL dgemm("N","N",n,n,n,-0.5d0,scrh_4,n,scr_3,n,1.0d0,sum_1,n)
      CALL dgemm("N","N",n,n,n,0.5d0,scrh_4,n,scr_4,n,1.0d0,sum_1,n)
      CALL mat_muld(scr_1, pvph,pe1p,n,1.0d0,0.0d0,tt,rr)
      CALL dgemm("N","N",n,n,n,1.0d0,vh,n,pe1p,n,0.0d0,scr_2,n)
      CALL dgemm("N","N",n,n,n,1.0d0,e1,n,pvph,n,0.0d0,scr_3,n)
      CALL dgemm("N","N",n,n,n,1.0d0,e1,n,vh,n,0.0d0,scr_4,n)
      CALL dgemm("N","N",n,n,n,1.0d0,vh,n,e1,n,0.0d0,scrh_1,n)
      CALL mat_1_over_h(scrh_1,n,e)
      CALL dgemm("N","N",n,n,n,1.0d0,pvph,n,e1,n,0.0d0,scrh_2,n)
      CALL mat_1_over_h(scrh_2,n,e)
      CALL dgemm("N","N",n,n,n,1.0d0,pe1p,n,vh,n,0.0d0,scr_3,n)
      CALL mat_1_over_h(scrh_3,n,e)
      CALL mat_muld(scrh_4,pe1p,pvph,n,1.0d0,0.0d0,tt,rr)
      CALL mat_1_over_h(scrh_4,n,e)
      CALL dgemm("N","N",n,n,n,0.5d0,scr_1,n,scrh_1,n,0.0d0,sum_1,n)
      CALL mat_muld(sum_1,scr_1,scrh_2,n,-0.5d0,1.0d0,tt,rr) 
      CALL dgemm("N","N",n,n,n,-0.5d0,scr_2,n,scrh_1,n,1.0d0,sum_1,n)
      CALL mat_muld(sum_1,scr_2,scrh_2,n, 0.5d0,1.0d0,tt,rr)
      CALL mat_muld(sum_1,scr_1,scrh_3,n,-0.5d0,1.0d0,tt,rr)
      CALL mat_muld(sum_1,scr_1,scrh_4,n, 0.5d0,1.0d0,tt,rr)
      CALL mat_muld(sum_1,scr_2,scrh_3,n, 0.5d0,1.0d0,tt,rr)
      CALL mat_muld(sum_1,scr_2,scrh_4,n,-0.5d0,1.0d0,tt,rr)
      CALL dgemm("N","N",n,n,n,-0.5d0,scr_3,n,scrh_1,n,0.0d0,sum_1,n)
      CALL mat_muld(sum_1,scr_3,scrh_2,n, 0.5d0,1.0d0,tt,rr)
      CALL mat_mulm(sum_1,scr_4,scrh_1,n, 0.5d0,1.0d0,tt,rr)
      CALL dgemm("N","N",n,n,n,-0.5d0,scr_4,n,scrh_2,n,1.0d0,sum_1,n)
      CALL mat_muld(sum_1,scr_3,scrh_3,n, 0.5d0,1.0d0,tt,rr)
      CALL mat_muld(sum_1,scr_3,scrh_4,n,-0.5d0,1.0d0,tt,rr)
      CALL dgemm("N","N",n,n,n,-0.5d0,scr_4,n,scrh_3,n,1.0d0,sum_1,n)
      CALL dgemm("N","N",n,n,n,0.5d0,scr_4,n,scrh_4,n,1.0d0,sum_1,n)
!-----------------------------------------------------------------------
!   3.  sum_2 = 1/8 [W1,[W1,[W1,O1]]] =
!
!             = 1/8 ( (W1^3)O1 - 3(W1^2)O1W1 + 3 W1O1(W1^2) - O1(W1^3) )
!-----------------------------------------------------------------------
      CALL dgemm("N","N",n,n,n,0.125d0,w1w1,n,w1o1,n,0.0d0,sum_2,n)
      CALL dgemm("N","N",n,n,n,-0.375d0,w1w1,n,o1w1,n,1.0d0,sum_2,n)
      CALL dgemm("N","N",n,n,n,0.375d0,w1o1,n,w1w1,n,1.0d0,sum_2,n)
      CALL dgemm("N","N",n,n,n,-0.125d0,o1w1,n,w1w1,n,1.0d0,sum_2,n)
!-----------------------------------------------------------------------
!   4.  result = sum_1 + sum_2
!-----------------------------------------------------------------------
      CALL mat_add(ev4,1.0d0,sum_1,1.0d0,sum_2,n)
!-----------------------------------------------------------------------
!   5. Finish up the stuff!!
!-----------------------------------------------------------------------
      DEALLOCATE(v,pvp,vh,pvph,w1w1,w1o1,o1w1,sum_1,sum_2)
      DEALLOCATE(scr_1,scr_2,scr_3,scr_4,scrh_1,scrh_2,scrh_3,scrh_4)
!
      RETURN
      END SUBROUTINE even4a
!
!
!
!-----------------------------------------------------------------------
!
!
!
!-----------------------------------------------------------------------
!                                                                      -
      SUBROUTINE mat_1_over_h (p,n,e)
!
!***********************************************************************
!                                                                      *
!   2. SR mat_1_over_h: Transform matrix p into matrix p/(e(i)+e(j))   *
!                                                                      *
!   p    in  REAL(:,:) :   matrix p                                    *
!   e    in  REAL(:)   :   rel. energy (diagonal)                      *
!   n    in  INTEGER                                                   *
!                                                                      *
!***********************************************************************
!
      IMPLICIT none
      REAL(KIND=double),DIMENSION(:),INTENT(IN)         ::  e
      REAL(KIND=double),DIMENSION(:,:),INTENT(INOUT)    ::  p
      INTEGER,INTENT(IN)                                ::  n
      INTEGER                                           :: i,j
!
      DO i=1,n
        DO j=1,n
          p(i,j)=p(i,j)/(e(i)+e(j))
        ENDDO
      ENDDO
!
      RETURN
      END SUBROUTINE mat_1_over_h
!
!-----------------------------------------------------------------------
!
      SUBROUTINE mat_axa (p,n,a)
!***********************************************************************
!                                                                      *
!   3. SR mat_axa: Transform matrix p into matrix  a*p*a               *
!                                                                      *
!   p    in  REAL(:,:):   matrix p                                     *
!   a    in  REAL(:)  :   A-factors (diagonal)                         *
!JT n    in  INTEGER  :   dimension of matrix p                        *  
!                                                                      *
!***********************************************************************
!
      IMPLICIT none
      REAL(KIND=double),DIMENSION(:),INTENT(IN)        :: a
      REAL(KIND=double),DIMENSION(:,:),INTENT(INOUT)   :: p
      INTEGER,INTENT(IN)                               :: n
      INTEGER                                          :: i,j
!
      DO i=1,n
        DO j=1,n
           p(i,j)=p(i,j)*a(i)*a(j)
        ENDDO
      ENDDO
!
      RETURN
      END SUBROUTINE mat_axa 
!
!-----------------------------------------------------------------------
!
      SUBROUTINE mat_arxra (p,n,a,r)
!
!***********************************************************************
!                                                                      *
!   4. SR mat_arxra: Transform matrix p into matrix  a*r*p*r*a         *
!                                                                      *
!   p    in  REAL(:,:) :   matrix p                                    *
!   a    in  REAL(:)   :   A-factors (diagonal)                        *
!   r    in  REAL(:)   :   R-factors (diagonal)                        *
!   n    in  INTEGER   :   dimension of matrix p                       *
!                                                                      *
!***********************************************************************
!
      IMPLICIT none
      REAL(KIND=double),DIMENSION(:),INTENT(IN)       :: a,r
      REAL(KIND=double),DIMENSION(:,:),INTENT(INOUT)  :: p
      INTEGER,INTENT(IN)                              :: n
      INTEGER                                         :: i,j
!
      DO i=1,n
        DO j=1,n
           p(i,j)=p(i,j)*a(i)*a(j)*r(i)*r(j)   
        ENDDO
      ENDDO
!
      RETURN
      END SUBROUTINE mat_arxra 
!
!-----------------------------------------------------------------------
!
      SUBROUTINE mat_mulm (p,q,r,n,alpha,beta,t,rr)
!
!***********************************************************************
!                                                                      *
!   5. SR mat_mulm:  Multiply matrices according to:                   *
!                                                                      *
!                      p = alpha*q*(..P^2..)*r + beta*p                *
!                                                                      *
!   p      out  REAL(:,:):   matrix p                                  *
!   q      in   REAL(:,:):   matrix q                                  *
!   r      in   REAL(:,.):   matrix r                                  *
!   n      in   INTEGER  :   dimension n of matrices                   *
!   alpha  in   REAL(double) :                                             *
!   beta   in   REAL(double) :                                             *
!   t      in   REAL(:)  :   non-rel. kinetic energy  (diagonal)       *
!   rr     in   REAL(:)  :   R-factors  (diagonal)                     *
!                                                                      *
!***********************************************************************
!
      IMPLICIT none
      REAL(KIND=double),DIMENSION(:,:),INTENT(IN)      :: r
      REAL(KIND=double),DIMENSION(:,:),INTENT(IN)      :: q,p
      REAL(KIND=double),DIMENSION(:,:),ALLOCATABLE     :: qtemp
      REAL(KIND=double),DIMENSION(:),INTENT(IN)        :: t,rr
      REAL(KIND=double),INTENT(IN)                     :: alpha,beta
      INTEGER,INTENT(IN)                               :: n
      INTEGER                                          :: i,j
!
      ALLOCATE(qtemp(n,n))
      DO i=1,n
        DO j=1,n
          qtemp(i,j)=q(i,j)*2.0d0*t(j)*rr(j)*rr(j)
        ENDDO
      ENDDO
      CALL dgemm("N","N",n,n,n,alpha,qtemp,n,r,n,beta,p,n)
      DEALLOCATE(qtemp)
!
      RETURN
      END SUBROUTINE mat_mulm 
!
!-----------------------------------------------------------------------
!
      SUBROUTINE mat_muld (p,q,r,n,alpha,beta,t,rr)
!
!***********************************************************************
!                                                                      *
!   16. SR mat_muld:  Multiply matrices according to:                  *
!                                                                      *
!                      p = alpha*q*(..1/P^2..)*r + beta*p              *
!                                                                      *
!   p      out  REAL(:,:):   matrix p                                  *
!   q      in   REAL(:,:):   matrix q                                  *
!   r      in   REAL(:,:):   matrix r                                  *
!   n      in   INTEGER  :   Dimension of all matrices                 *
!   alpha  in   REAL(double) :                                         *
!   beta   in   REAL(double) :                                         *
!   t      in   REAL(:)  :   non-rel. kinetic energy  (diagonal)       *
!   rr     in   REAL(:)  :   R-factors  (diagonal)                     *
!                                                                      *
!***********************************************************************
!
      IMPLICIT none
      REAL(KIND=double),DIMENSION(:,:),INTENT(IN)      :: r
      REAL(KIND=double),DIMENSION(:,:),INTENT(IN)      :: q,p
      REAL(KIND=double),DIMENSION(:,:),ALLOCATABLE     :: qtemp
      REAL(KIND=double),DIMENSION(:),INTENT(IN)        :: t,rr
      REAL(KIND=double),INTENT(IN)                     :: alpha,beta
      INTEGER,INTENT(IN)                               :: n
      INTEGER                                          :: i,j
!
      ALLOCATE(qtemp(n,n)) 
      DO i=1,n
        DO j=1,n
          qtemp(i,j)=q(i,j)*0.5d0/(t(j)*rr(j)*rr(j))
        ENDDO
      ENDDO
      CALL dgemm("N","N",n,n,n,alpha,qtemp,n,r,n,beta,p,n)
      DEALLOCATE(qtemp) 
!       
      RETURN
      END SUBROUTINE mat_muld 
!       
!-----------------------------------------------------------------------
!       
      SUBROUTINE mat_add (p,alpha,q,beta,r,n)
!
!***********************************************************************
!                                                                      *
!   19. SR mat_add:  Add two matrices of the same size according to:   *
!                                                                      *
!                            p = alpha*q + beta*r                      *
!                                                                      *
!   p      out  REAL(:,:)  :   matrix p                                *
!   q      in   REAL(:,:)  :   matrix q                                *
!   r      in   REAL(:,:)  :   matrix r                                *
!   alpha  in   REAL(double)                                           *
!   beta   in   REAL(double)                                           *
!                                                                      *
!   Matrix p must already exist before calling this SR!!               *
!                                                                      *
!  [written by: Alexander Wolf,  20.2.2002,  v1.0]                     *
!                                                                      *
!***********************************************************************
!
      IMPLICIT none
      REAL(KIND=double),DIMENSION(:,:),INTENT(IN)      :: q,r
      REAL(KIND=double),DIMENSION(:,:),INTENT(OUT)     :: p
      REAL(KIND=double),INTENT(IN)                     :: alpha,beta
      INTEGER,INTENT(IN)                               :: n
      INTEGER                                          :: i,j
!
      DO i=1,n
        DO j=1,n
          p(i,j) = alpha*q(i,j) + beta*r(i,j)
        ENDDO
      ENDDO
!
      RETURN
      END SUBROUTINE mat_add
!
!---------------------------------------------------------------------
!
      SUBROUTINE TRSM ( W,B,C,N,H)
!
      IMPLICIT NONE
      REAL(KIND=double),DIMENSION(:,:) :: B,C,H,W
      INTEGER                     :: N,I,J,K,L,IJ
!
      IJ=0
      DO I=1,N
        DO J=1,I
          IJ=IJ+1
          C(I,J)=0.0d0
          C(J,I)=0.0d0
          H(I,J)=0.0d0
          H(J,I)=0.0d0
        END DO
      END DO
      DO I=1,N
        DO L=1,N
          DO K=1,N
            H(I,L)=B(K,I)*W(K,L)+H(I,L)
          END DO
        END DO
      END DO 
      IJ=0
      DO I=1,N
        DO J=1,I
          IJ=IJ+1
          DO L=1,N
            C(I,J)=H(I,L)*B(L,J)+C(I,J)
            C(J,I)=C(I,J)
          END DO
        END DO
      END DO
!
      RETURN
      END SUBROUTINE TRSM
!
!----------------------------------------------------------------------
!
      SUBROUTINE DIAG (matrix_t_pgf,n,eig,ew,matrix_sinv_pgf,aux,ic)
!
      IMPLICIT NONE
      REAL(KIND=double), DIMENSION(:,:),INTENT(INOUT)    :: matrix_t_pgf,eig,aux
      REAL(KIND=double), DIMENSION(:,:),INTENT(IN)       :: matrix_sinv_pgf  
      REAL(KIND=double), DIMENSION(:),INTENT(INOUT)      :: ew
      INTEGER,INTENT(IN)                                 :: n
      INTEGER                                            :: i,j,k,ic
!
      eig = 0.0d0
      aux = 0.0d0
      CALL dgemm("N","N",n,n,n,1.0d0,matrix_t_pgf,n,matrix_sinv_pgf,n,0.0d0,eig,n)
      aux = 0.0d0
      CALL dgemm("T","N",n,n,n,1.0d0,matrix_sinv_pgf,n,eig,n,0.0d0,aux,n)
      CALL JACOB2 ( AUX,EIG,EW,N,IC )
!
      RETURN
      END SUBROUTINE DIAG
!
!---------------------------------------------------------------------
!
      SUBROUTINE JACOB2 ( sogt,eigv,eigw,n,ic )
!
      IMPLICIT NONE
      INTEGER, INTENT(IN)                              :: n
      REAL (KIND=double),DIMENSION(n,n),INTENT(OUT)    :: eigv        
      REAL (KIND=double),DIMENSION(n),INTENT(OUT)      :: eigw       
      REAL (KIND=double),DIMENSION(n,n),INTENT(INOUT)  :: sogt      
      REAL (KIND=double)                               :: tol 
      INTEGER, INTENT(IN)                              :: ic
      INTEGER                      :: i,j,ind,m,l,im,mm,il,ll,k
      REAL (KIND=double)           :: ext_norm,thr_min,u1,thr,x,y,xy,&
                                      sint,sint2,cost,cost2,sincs
!
      tol=1.0E-15
      ext_norm=0.0d0
      u1=real(n)
      DO i=1,n
        eigv(i,i)=1.0d0
        eigw(i)=sogt(i,i)
        DO j=1,i
          IF(i.ne.j) THEN
            eigv(i,j)=0.0d0
            eigv(j,i)=0.0d0
            ext_norm=ext_norm+sogt(i,j)*sogt(i,j)
          END IF
        END DO
      END DO
      IF (ext_norm.gt.0.0d0) THEN 
        ext_norm=DSQRT(2.0d0*ext_norm)
        thr_min=ext_norm*tol/u1                                  
        ind=0
        thr=ext_norm
        DO
          thr=thr/u1                                            
          DO
            l=1                       
            DO                                                 
              m=l+1
              DO                                              
                IF ((DABS(sogt(m,l))-thr).ge.0.0d0) THEN    
                  ind=1
                  x=0.5d0*(eigw(l)-eigw(m))                       
                  y=-sogt(m,l)/DSQRT(sogt(m,l)*sogt(m,l)+x*x)   
                  IF (x.lt.0.0d0) y=-y                     
                                                           
                  IF (y.gt.1.0d0) y=1.0d0               
                  IF (y.lt.-1.0d0) y=-1.0d0                      
                  xy=1.0d0-y*y                                    
                  sint=y/DSQRT(2.0d0*(1.0d0+DSQRT(xy)))
                  sint2=sint*sint
                  cost2=1.0d0-sint2
                  cost=DSQRT(cost2)
                  sincs=sint*cost
                  DO i=1,n                              
                    IF((i-m).ne.0) THEN
                      IF ((i-m).lt.0) THEN
                        im=m
                        mm=i
                      ELSE
                        im=i
                        mm=m
                      END IF
                      IF ((i-l).ne.0) THEN
                        IF ((i-l).lt.0) THEN
                          il=l
                          ll=i
                        ELSE
                          il=i
                          ll=l
                        END IF
                        x=sogt(il,ll)*cost-sogt(im,mm)*sint  
                        sogt(im,mm)=sogt(il,ll)*sint+sogt(im,mm)*cost 
                        sogt(il,ll)=x  
                      END IF
                    END IF
                    x=eigv(i,l)*cost-eigv(i,m)*sint                  
                    eigv(i,m)=eigv(i,l)*sint+eigv(i,m)*cost
                    eigv(i,l)=x
                  END DO
                  x=2.0d0*sogt(m,l)*sincs                          
                  y=eigw(l)*cost2+eigw(m)*sint2-x
                  x=eigw(l)*sint2+eigw(m)*cost2+x
                  sogt(m,l)=(eigw(l)-eigw(m))*sincs+sogt(m,l)*(cost2-sint2)
                  eigw(l)=y
                  eigw(m)=x
                END IF
                IF ((m-n).eq.0) EXIT
                m=m+1
              END DO
              IF ((l-m+1).eq.0) EXIT
              l=l+1
            END DO
            IF((ind-1).ne.0.0d0) EXIT
            ind=0
          END DO
         IF ((thr-thr_min).le.0.0d0) EXIT
        END DO
      END IF
      IF (ic.ne.0) THEN
        DO i=1,n
          DO j=1,n
            IF ((eigw(i)-eigw(j)).gt.0.0d0) THEN
              x=eigw(i)
              eigw(i)=eigw(j)
              eigw(j)=x
              DO k=1,n
                y=eigv(k,i)
                eigv(k,i)=eigv(k,j)
                eigv(k,j)=y
              END DO
            END IF
          END DO
        END DO 

      END IF
!
      RETURN
      END SUBROUTINE JACOB2
!
!---------------------------------------------------------------------
!
      SUBROUTINE SOG (n,matrix_s_pgf,matrix_sinv_pgf)
!
      IMPLICIT NONE
      REAL(KIND=double), DIMENSION(:,:),INTENT(IN)   :: matrix_s_pgf
      REAL(KIND=double), DIMENSION(:,:),INTENT(INOUT):: matrix_sinv_pgf
      REAL(KIND=double), DIMENSION(:),ALLOCATABLE    :: a
      REAL(KIND=double), DIMENSION(:,:),ALLOCATABLE  :: g
      INTEGER                                        :: n,i,j,jn,j2,k
      REAL(KIND=double)                              :: diag_s,scalar,row_sum
!
      ALLOCATE(a(n))
      ALLOCATE(g(n,n))
      DO jn=1,n                                                    
        diag_s = matrix_s_pgf(jn,jn)
        g(jn,jn)=1.0d0
        IF(jn.ne.1) THEN                                          
          DO j=1,jn-1                                            
            scalar=0.0d0
            DO i=1,j                                            
              scalar=scalar+matrix_s_pgf(i,jn)*g(i,j)          
            END DO
            diag_s=diag_s-scalar*scalar                       
            a(j) = scalar
          END DO
          DO j=1,jn-1                                        
            row_sum=0.0d0
            DO k=j,jn-1                                     
              row_sum=row_sum+a(k)*g(j,k)
            END DO
            g(j,jn)=-row_sum                               
          END DO   
        END IF                                            
        diag_s=1.0d0/DSQRT(diag_s)                      
        DO i=1,jn
          g(i,jn)=g(i,jn)*diag_s                        
        END DO
      END DO
      DO j=1,n
        DO i=1,j
          matrix_sinv_pgf(j,i)=0.0d0 
          matrix_sinv_pgf(i,j)=g(i,j)
        END DO
      END DO      
      DEALLOCATE(a,g)
!
      RETURN
      END SUBROUTINE SOG
!
!--------------------------------------------------------------------
!
END module dkh_main
!
!--------------------------------------------------------------------
