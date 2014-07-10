C  Copyright (c) 2003-2010 University of Florida
C
C  This program is free software; you can redistribute it and/or modify
C  it under the terms of the GNU General Public License as published by
C  the Free Software Foundation; either version 2 of the License, or
C  (at your option) any later version.

C  This program is distributed in the hope that it will be useful,
C  but WITHOUT ANY WARRANTY; without even the implied warranty of
C  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C  GNU General Public License for more details.

C  The GNU General Public License is included in this distribution
C  in the file COPYRIGHT.
         SUBROUTINE  OED__NAI_DERV_3CEN_PCGTO_BLOCK
     +
     +                    ( NBATCH,
     +                      NINT1DX,NINT1DY,NINT1DZ,
     +                      ATOMIC,
     +                      MIJ,NCEN,MIJCEN,
     +                      NIJ,NIJBEG,NIJEND,
     +                      NGQP,NMOM,NGQSCR,MGIJCEN,
     +                      NPGTOA,NPGTOB,
     +                      NXYZA,NXYZB,
     +                      SHELLA,SHELLB,SHELLP,
     +                      XA,YA,ZA,XB,YB,ZB,
     +                      ABX,ABY,ABZ,
     +                      NUCLEI,
     +                      XN,YN,ZN,NCHARGE,
     +                      NUCCEN,
     +                      DERAX,DERAY,DERAZ,
     +                      DERBX,DERBY,DERBZ,
     +                      DERCX,DERCY,DERCZ,
     +                      DIFFX,DIFFY,DIFFZ,
     +                      DIFFA,DIFFB,DIFFC,
     +                      IXAEQB,
     +                      ALPHAA,ALPHAB,
     +                      FTABLE,MGRID,NGRID,TMAX,TSTEP,TVSTEP,
     +                      PRIMA,PRIMB,
     +                      NORMA,NORMB,
     +                      RHOAB,
     +                      PX,PY,PZ,PAX,PAY,PAZ,PINVHF,SCALE,
     +                      RTS,WTS,GQSCR,TVAL,
     +                      R1X,R1Y,R1Z,R2,
     +                      EXP2A,EXP2B,EXP2AB,
     +                      INT1DX,INT1DY,INT1DZ,
     +
     +                                BATCH )
     +
C------------------------------------------------------------------------
C  OPERATION   : OED__NAI_DERV_3CEN_PCGTO_BLOCK
C  MODULE      : ONE ELECTRON INTEGRALS DIRECT
C  MODULE-ID   : OED
C  SUBROUTINES : OED__RYS_ROOTS_WEIGHTS
C                OED__NAI_1D_COEFFICIENTS
C                OED__NAI_1D_AB_INTEGRALS
C                OED__NAI_1D_CENDERV_INTEGRALS
C                OED__NAI_1D_SHDERV_INTEGRALS
C                OED__NAI_DERV_INT1D_TO_00
C                OED__NAI_DERV_INT1D_TO_A0
C                OED__NAI_DERV_INT1D_TO_AB
C  DESCRIPTION : This operation calculates a batch of derivated unnormed
C                3-center nuclear attraction integrals between cartesian
C                gaussians
C
C                              sum  [A|C|B]
C                               C          ij
C
C                for a block of ij exponent pairs and a set of nuclear
C                attraction centers C different from the shell centers
C                A and B. Differentiation is possible on both A and B
C                shell centers and eventually on a specific nuclear
C                attraction center C. Note, that the atomic case A = B
C                is allowed here and is still considered to be a
C                3-center case.
C
C                After differentiation the total number of nuclear
C                attraction integrals generated here is thus given by
C                the total number of cartesian monomials NXYZA * NXYZB
C                times the total number of exponent pairs MIJ in the
C                present block.
C
C                On exit, the batch elements will be stored as:
C
C                             batch (ij,nxyza*nxyzb)
C
C
C                  Input:
C
C                    NBATCH       =  size of the array that will hold
C                                    the final primitive cartesian
C                                    derivative integral batch as well
C                                    as intermediate differentiated
C                                    1D integrals
C                    NINT1Dx      =  space needed for each of the 1D
C                                    x = X,Y,Z integral arrays (they
C                                    might be different due to different
C                                    orders of differentiation for
C                                    each cartesian component)
C                    ATOMIC       =  indicates, if purely atomic
C                                    integrals will be evaluated
C                    MIJ          =  current # of ij primitive index
C                                    pairs corresponding to the
C                                    contracted shell pairs A,B
C                    NCEN         =  # of nuclear attraction centers
C                    MIJCEN       =  # of ij primitive index pairs
C                                    times # of nuclear attraction
C                                    centers
C                    NIJ          =  total # of ij primitive index
C                                    pairs for the contracted shell
C                                    pair A,B
C                    NIJBEG(END)  =  first(last) ij primitive index
C                                    defining the ij block
C                    NGQP         =  # of gaussian quadrature points
C                                    (roots)
C                    NMOM         =  # of necessary moment integrals
C                                    to calculate the quadrature roots
C                    NGQSCR       =  size of gaussian quadrature
C                                    scratch space needed to calculate
C                                    all the quadrature roots
C                    MGIJCEN      =  # of roots times # of ij primitive
C                                    index pairs times # of nuclear
C                                    attraction centers
C                    NPGTOx       =  # of primitives per contraction
C                                    for contraction shells x = A,B
C                    NXYZx        =  # of cartesian monomials for
C                                    each contraction shell x = A,B
C                    SHELLx       =  the shell type for contraction
C                                    shells x = A,B and P=A+B
C                    Xx,Yx,Zx     =  the x,y,z-coordinates for centers
C                                    x = A,B
C                    ABm          =  the m=x,y,z-coordinate differences
C                                    between centers A and B
C                    NUCLEI       =  # of nuclear attraction centers
C                    XN,YN,ZN     =  the x,y,z-coordinates for all
C                                    nuclear attraction centers
C                    NCHARGE      =  the nuclear charges for all
C                                    nuclear attraction centers
C                    NUCCEN       =  contains those index labels of
C                                    the nuclear attraction centers
C                                    that survived the screening process
C                    DERyp        =  the order of differentiation on
C                                    centers y = A,B and possibly a
C                                    nuclear attraction center C with
C                                    respect to the p = x,y,z
C                                    coordinates
C                    DIFFy        =  is true, if differentiation will
C                                    be performed on centers y = A and B
C                                    and at least on one of the nuclear
C                                    attraction centers y = C different
C                                    from A and B and if differentiation
C                                    will be performed at least once
C                                    one on each of the y = x,y,z
C                                    coordinates
C                    IXAEQB       =  is an index = 0 or 1, depending
C                                    on if center A is different or
C                                    equal from center B, respectively
C                    ALPHAx       =  the primitive exponents for
C                                    contraction shells x = A,B
C                    FTABLE       =  Fm (T) table for interpolation
C                                    in low T region
C                    MGRID        =  maximum m in Fm (T) table
C                    NGRID        =  # of T's for which Fm (T) table
C                                    was set up
C                    TMAX         =  maximum T in Fm (T) table
C                    TSTEP        =  difference between two consecutive
C                                    T's in Fm (T) table
C                    TVSTEP       =  Inverse of TSTEP
C                    PRIMx        =  i,j labels of primitives for the
C                                    respective contraction shells
C                                    x = A,B
C                    NORMx        =  the normalization factors due to
C                                    the primitive exponents for the
C                                    contraction shells x = A,B
C                    RHOAB        =  the complete set of NIJ exponential
C                                    prefactors between contraction
C                                    shells A and B
C                    Px           =  will hold current MIJ coordinates
C                                    x=X,Y,Z for the gaussian product
C                                    centers P=A+B
C                    PAx          =  will hold current MIJ coordinate
C                                    x=X,Y,Z differences P-A between
C                                    centers P and A
C                    PINVHF       =  will hold current MIJ values of
C                                    1/(2*P), where P are the exponent
C                                    sums for contraction shells A
C                                    and B
C                    SCALE        =  will hold current distinct MIJCEN
C                                    (expanded to MGIJCEN) values of
C                                    scaling factors
C                    RTS          =  will hold all current MGIJCEN
C                                    quadrature roots
C                    WTS          =  will hold all current MGIJCEN
C                                    quadrature weights
C                    GQSCR        =  will be used as scratch space
C                                    for determining the quadrature
C                                    roots and weights
C                    TVAL         =  will hold current MIJCEN values
C                                    of T-exponents defining the Rys
C                                    weight functions
C                    R1x          =  will hold the current MGIJCEN
C                                    VRR R1-coefficients (individual
C                                    cartesian components x=X,Y,Z) for
C                                    shell expansion on center P
C                    R2           =  will hold the current MGIJCEN
C                                    coordinate independent VRR
C                                    R2-coefficients
C                    EXP2x        =  will hold current double primitive
C                                    exponent values in MIJ order
C                                    (expanded to MGIJCEN) for each
C                                    contracted shell x = A,B
C                    EXP2AB       =  will hold current double primitive
C                                    exponent sum values in MIJ order
C                                    between both contracted shells
C                                    A and B
C                    INT1Dx       =  will hold all current derivated
C                                    1D integrals for each cartesian
C                                    component (x = X,Y,Z) during all
C                                    stages of differentiation
C
C                  Output:
C
C                    BATCH        =  current batch of primitive
C                                    derivative cartesian nuclear
C                                    attraction [A|B] integrals
C
C
C  AUTHOR      : Norbert Flocke
C------------------------------------------------------------------------
C
C
C             ...include files and declare variables.
C
C
         IMPLICIT    NONE

         LOGICAL     ATOMIC
         LOGICAL     DIFFA,DIFFB,DIFFC
         LOGICAL     DIFFX,DIFFY,DIFFZ
         LOGICAL     PROCEED

         INTEGER     CENTER
         INTEGER     DERA,DERB
         INTEGER     DERAX,DERAY,DERAZ
         INTEGER     DERBX,DERBY,DERBZ
         INTEGER     DERCX,DERCY,DERCZ
         INTEGER     G000,G010,G020,G030,G040,G050,G060
         INTEGER     I,J,L,M,N
         INTEGER     IJ
         INTEGER     IXC,IXAEQB
         INTEGER     MGRID,NGRID
         INTEGER     MIJ,MIJCEN,MGIJCEN
         INTEGER     NBATCH
         INTEGER     NC,NCEN
         INTEGER     NG,NGQP,NGQPCEN,NMOM,NGQSCR
         INTEGER     NIJ,NIJBEG,NIJEND
         INTEGER     NINT,NINT1DX,NINT1DY,NINT1DZ
         INTEGER     NPGTOA,NPGTOB
         INTEGER     NUCLEI
         INTEGER     NXYZA,NXYZB
         INTEGER     SHELLA,SHELLB,SHELLP

         INTEGER     NUCCEN (1:NCEN)

         INTEGER     PRIMA (1:MIJ)
         INTEGER     PRIMB (1:MIJ)

         DOUBLE PRECISION  ABX,ABY,ABZ
         DOUBLE PRECISION  EXPA,EXPB
         DOUBLE PRECISION  P,PINV,PVAL
         DOUBLE PRECISION  PCX,PCY,PCZ
         DOUBLE PRECISION  PXVAL,PYVAL,PZVAL
         DOUBLE PRECISION  RNPCSQ
         DOUBLE PRECISION  SCALEM,SCALEN
         DOUBLE PRECISION  TEMP
         DOUBLE PRECISION  TMAX,TSTEP,TVSTEP
         DOUBLE PRECISION  XA,YA,ZA,XB,YB,ZB,XC,YC,ZC
         DOUBLE PRECISION  ZERO,HALF,ONE

         DOUBLE PRECISION  XN      (1:NUCLEI)
         DOUBLE PRECISION  YN      (1:NUCLEI)
         DOUBLE PRECISION  ZN      (1:NUCLEI)
         DOUBLE PRECISION  NCHARGE (1:NUCLEI)

         DOUBLE PRECISION  ALPHAA  (1:NPGTOA)
         DOUBLE PRECISION  ALPHAB  (1:NPGTOB)

         DOUBLE PRECISION  NORMA   (1:NPGTOA)
         DOUBLE PRECISION  NORMB   (1:NPGTOB)

         DOUBLE PRECISION  RHOAB   (1:NIJ)

         DOUBLE PRECISION  BATCH   (1:NBATCH)

         DOUBLE PRECISION  EXP2AB  (1:MIJ)
         DOUBLE PRECISION  PX      (1:MIJ)
         DOUBLE PRECISION  PY      (1:MIJ)
         DOUBLE PRECISION  PZ      (1:MIJ)
         DOUBLE PRECISION  PAX     (1:MIJ)
         DOUBLE PRECISION  PAY     (1:MIJ)
         DOUBLE PRECISION  PAZ     (1:MIJ)
         DOUBLE PRECISION  PINVHF  (1:MIJ)

         DOUBLE PRECISION  GQSCR   (1:NGQSCR)
         DOUBLE PRECISION  TVAL    (1:MIJCEN)

         DOUBLE PRECISION  R1X     (1:MGIJCEN)
         DOUBLE PRECISION  R1Y     (1:MGIJCEN)
         DOUBLE PRECISION  R1Z     (1:MGIJCEN)
         DOUBLE PRECISION  R2      (1:MGIJCEN)
         DOUBLE PRECISION  EXP2A   (1:MGIJCEN)
         DOUBLE PRECISION  EXP2B   (1:MGIJCEN)
         DOUBLE PRECISION  RTS     (1:MGIJCEN)
         DOUBLE PRECISION  SCALE   (1:MGIJCEN)
         DOUBLE PRECISION  WTS     (1:MGIJCEN)

         DOUBLE PRECISION  INT1DX  (1:NINT1DX)
         DOUBLE PRECISION  INT1DY  (1:NINT1DY)
         DOUBLE PRECISION  INT1DZ  (1:NINT1DZ)

         DOUBLE PRECISION  FTABLE  (0:MGRID,0:NGRID)

         PARAMETER  (ZERO = 0.D0)
         PARAMETER  (HALF = 0.5D0)
         PARAMETER  (ONE  = 1.D0)
C
C
C------------------------------------------------------------------------
C
C
C             ...calculate the quantities needed to establish the
C                1D nuclear attraction integrals.
C
C
         IF (ATOMIC) THEN
             L = 0
             M = 0
             N = 0
             DO IJ = NIJBEG,NIJEND
                M = M + 1
                I = PRIMA (M)
                J = PRIMB (M)
                P = ALPHAA (I) + ALPHAB (J)
                PX (M) = XA
                PY (M) = YA
                PZ (M) = ZA
                PINV = ONE / P
                PINVHF (M) = HALF * PINV
                SCALEM = PINV * NORMA (I) * NORMB (J)
                DO NC = 1,NCEN
                   N = N + 1
                   IXC = NUCCEN (NC)
                   PCX = XA - XN (IXC)
                   PCY = YA - YN (IXC)
                   PCZ = ZA - ZN (IXC)
                   RNPCSQ = PCX * PCX + PCY * PCY + PCZ * PCZ
                   TVAL (N) = P * RNPCSQ
                   SCALEN = SCALEM * NCHARGE (IXC)
                   DO NG = 1,NGQP
                      L = L + 1
                      SCALE (L) = SCALEN
                   END DO
                END DO
             END DO
         ELSE
             L = 0
             M = 0
             N = 0
             DO IJ = NIJBEG,NIJEND
                M = M + 1
                I = PRIMA (M)
                J = PRIMB (M)
                EXPA = ALPHAA (I)
                EXPB = ALPHAB (J)
                P = EXPA + EXPB
                PINV = ONE / P
                PVAL = - EXPB * PINV
                PAX (M) = PVAL * ABX
                PAY (M) = PVAL * ABY
                PAZ (M) = PVAL * ABZ
                PXVAL = PAX (M) + XA
                PYVAL = PAY (M) + YA
                PZVAL = PAZ (M) + ZA
                PX (M) = PXVAL
                PY (M) = PYVAL
                PZ (M) = PZVAL
                PINVHF (M) = HALF * PINV
                SCALEM = PINV * NORMA (I) * NORMB (J) * RHOAB (IJ)
                DO NC = 1,NCEN
                   N = N + 1
                   IXC = NUCCEN (NC)
                   PCX = PXVAL - XN (IXC)
                   PCY = PYVAL - YN (IXC)
                   PCZ = PZVAL - ZN (IXC)
                   RNPCSQ = PCX * PCX + PCY * PCY + PCZ * PCZ
                   TVAL (N) = P * RNPCSQ
                   SCALEN = SCALEM * NCHARGE (IXC)
                   DO NG = 1,NGQP
                      L = L + 1
                      SCALE (L) = SCALEN
                   END DO
                END DO
             END DO
         END IF
C
C
C             ...determine memory allocation offsets for the scratch
C                arrays used to calculate the quadrature roots +
C                weights:
C
C                   G000 = offset for A coefficients (Jacobi/Laguerre)
C                   G010 = offset for B coefficients (Jacobi/Laguerre)
C                   G020 = offset for moments (Jacobi/Laguerre)
C                   G030 = offset for diagonals of symmetric termat
C                   G040 = offset for offdiagonals of symmetric termat
C                   G050 = offset for first row intermediates during
C                          evaluation of symmetric termat
C                   G060 = offset for second row intermediates during
C                          evaluation of symmetric termat
C
C
         G000 = 1
         G010 = G000 + NMOM
         G020 = G010 + NMOM - 1
         G030 = G020 + NMOM
         G040 = G030 + NGQP
         G050 = G040 + NGQP
         G060 = G050 + NMOM
C
C
C             ...calculate all roots and weights for all ij pairs and
C                nuclear centers. Array R2 is passed as a scratch array.
C
C
         CALL    OED__RYS_ROOTS_WEIGHTS
     +
     +                ( MIJCEN,MGIJCEN,
     +                  NGQP,NMOM,
     +                  TVAL,R2,
     +                  FTABLE,MGRID,NGRID,TMAX,TSTEP,TVSTEP,
     +                  GQSCR(G000),GQSCR(G010),
     +                  GQSCR(G020),
     +                  GQSCR(G030),GQSCR(G040),
     +                  GQSCR(G050),GQSCR(G060),
     +
     +                           RTS,
     +                           WTS )
     +
     +
C
C
C             ...calculate the derivative 1D integral coefficents
C                (if needed). These are the double exponent values
C                for differentiation on centers A and/or B and the
C                double exponent sum values between centers A and B.
C
C
         NGQPCEN = NGQP * NCEN

         IF (DIFFA) THEN
             M = 0
             DO IJ = 1,MIJ
                I = PRIMA (IJ)
                EXPA = ALPHAA (I) + ALPHAA (I)
                DO N = 1,NGQPCEN
                   EXP2A (M+N) = EXPA
                END DO
                M = M + NGQPCEN
             END DO
         END IF

         IF (DIFFB) THEN
             M = 0
             DO IJ = 1,MIJ
                J = PRIMB (IJ)
                EXPB = ALPHAB (J) + ALPHAB (J)
                DO N = 1,NGQPCEN
                   EXP2B (M+N) = EXPB
                END DO
                M = M + NGQPCEN
             END DO
         END IF

         IF (DIFFC) THEN
             DO IJ = 1,MIJ
                I = PRIMA (IJ)
                J = PRIMB (IJ)
                P = ALPHAA (I) + ALPHAB (J)
                EXP2AB (IJ) = P + P
             END DO
         END IF
C
C
C             ...generate all VRR coefficients.
C
C
         CALL    OED__NAI_1D_COEFFICIENTS
     +
     +                ( NGQP,NCEN,
     +                  MIJ,MGIJCEN,
     +                  ATOMIC,
     +                  NUCLEI,
     +                  XN,YN,ZN,
     +                  NUCCEN,
     +                  PX,PY,PZ,
     +                  PAX,PAY,PAZ,
     +                  PINVHF,
     +                  RTS,
     +
     +                            R1X,R1Y,R1Z,
     +                            R2 )
     +
     +
C
C
C             ...start assembling the 1D AB nuclear attraction
C                integrals and their differentiation. Assemble first
C                the 1DX AB integrals and perform differentiation
C                sequence, if necessary.
C
C
C         WRITE (*,*) ' ATOMIC = ',ATOMIC
C         WRITE (*,*) ' NCEN = ',NCEN
C         WRITE (*,*) ' NGQP = ',NGQP
C         WRITE (*,*) ' IXAEQB = ',IXAEQB
C         WRITE (*,*) ' DIFFA,DIFFB,DIFFC = ',DIFFA,DIFFB,DIFFC
C         WRITE (*,*) ' DIFFX,DIFFY,DIFFZ = ',DIFFX,DIFFY,DIFFZ
C         WRITE (*,*) ' DERAX,DERAY,DERAZ = ',DERAX,DERAY,DERAZ
C         WRITE (*,*) ' DERBX,DERBY,DERBZ = ',DERBX,DERBY,DERBZ
C         WRITE (*,*) ' DERCX,DERCY,DERCZ = ',DERCX,DERCY,DERCZ

         IF (DIFFX) THEN

             CALL    OED__NAI_1D_AB_INTEGRALS
     +
     +                    ( SHELLP+DERAX+DERBX,
     +                      MIN (SHELLA+DERAX,SHELLB+DERBX),
     +                      SHELLA+DERAX,SHELLB+DERBX,
     +                      MGIJCEN,
     +                      WTS,
     +                      R1X,R2,
     +                      ABX,
     +                      .TRUE.,
     +                      BATCH,
     +
     +                                INT1DX )
     +
     +
C
C
C             ...apply (if any) differentiation of the 1DX AB integrals
C                on a nuclear attraction center different from shell
C                centers A and B.
C
C
             PROCEED = DIFFC .AND. (DERCX.GT.0)

             IF (PROCEED) THEN
                 IXC = NUCCEN (1)
                 XC = XN (IXC)
                 NINT = MGIJCEN * (SHELLA+DERAX+1) * (SHELLB+DERBX+1)
                 DO N = 1,DERCX
                    IF (N.GT.1) THEN
                        DO M = 1,NINT
                           TEMP = BATCH (M)
                           BATCH (M) = INT1DX (M)
                           INT1DX (M) = TEMP
                        END DO
                    ELSE
                        DO M = 1,NINT
                           BATCH (M) = INT1DX (M)
                        END DO
                    END IF

                    CALL    OED__NAI_1D_CENDERV_INTEGRALS
     +
     +                           ( MIJ,MGIJCEN,
     +                             NGQP,
     +                             SHELLA+DERAX,SHELLB+DERBX,
     +                             EXP2AB,
     +                             PX,XC,
     +                             N,
     +                             RTS,
     +                             R1X,
     +                             BATCH,
     +
     +                                      INT1DX )
     +
     +
                 END DO
             END IF
C
C
C             ...perform differentiation (if any) of the 1DX AB
C                integrals on shell center A (and simultaneously on
C                shell center B, if atomic).
C
C
             PROCEED = DIFFA .AND. (DERAX.GT.0)

             IF (PROCEED) THEN
                 DERA = DERAX
                 DERB = DERBX
                 DO N = 1,DERAX
                    NINT = MGIJCEN * (SHELLA+DERA+1) * (SHELLB+DERB+1)
                    DO M = 1,NINT
                       BATCH (M) = INT1DX (M)
                    END DO

                    DERA = DERA - 1
                    DERB = DERB - IXAEQB

                    IF (DERA.LT.0 .OR. DERB.LT.0) THEN
                        WRITE (*,*) ' Problems forming x-derivatives! '
                        WRITE (*,*) ' DERA,DERB = ',DERA,DERB
                        WRITE (*,*) ' oed__nai_derv_3cen_pcgto_block '
                        STOP
                    END IF

                    CALL    OED__NAI_1D_SHDERV_INTEGRALS
     +
     +                           ( MGIJCEN,
     +                             SHELLA+DERA,SHELLB+DERB,
     +                             1,IXAEQB,
     +                             EXP2A,EXP2B,
     +                             BATCH,
     +
     +                                      INT1DX )
     +
     +
                 END DO
             END IF
C
C
C             ...if not atomic, perform differentiation (if any) of
C                the 1DX AB integrals on the remaining shell center B.
C
C
             PROCEED = (.NOT.ATOMIC) .AND. DIFFB .AND. (DERBX.GT.0)

             IF (PROCEED) THEN
                 DERB = DERBX
                 DO N = 1,DERBX
                    NINT = MGIJCEN * (SHELLA+1) * (SHELLB+DERB+1)
                    DO M = 1,NINT
                       BATCH (M) = INT1DX (M)
                    END DO

                    DERB = DERB - 1

                    IF (DERB.LT.0) THEN
                        WRITE (*,*) ' Problems forming x-derivatives! '
                        WRITE (*,*) ' DERB = ',DERB
                        WRITE (*,*) ' oed__nai_derv_3cen_pcgto_block '
                        STOP
                    END IF

                    CALL    OED__NAI_1D_SHDERV_INTEGRALS
     +
     +                           ( MGIJCEN,
     +                             SHELLA,SHELLB+DERB,
     +                             0,1,
     +                             EXP2A,EXP2B,
     +                             BATCH,
     +
     +                                      INT1DX )
     +
     +
                 END DO
             END IF

         ELSE

             CALL    OED__NAI_1D_AB_INTEGRALS
     +
     +                    ( SHELLP,MIN (SHELLA,SHELLB),
     +                      SHELLA,SHELLB,
     +                      MGIJCEN,
     +                      WTS,
     +                      R1X,R2,
     +                      ABX,
     +                      .TRUE.,
     +                      BATCH,
     +
     +                                INT1DX )
     +
     +
         END IF
C
C
C             ...assemble next the 1DY AB integrals and perform
C                differentiation sequence, if necessary.
C
C
         IF (DIFFY) THEN

             CALL    OED__NAI_1D_AB_INTEGRALS
     +
     +                    ( SHELLP+DERAY+DERBY,
     +                      MIN (SHELLA+DERAY,SHELLB+DERBY),
     +                      SHELLA+DERAY,SHELLB+DERBY,
     +                      MGIJCEN,
     +                      WTS,
     +                      R1Y,R2,
     +                      ABY,
     +                      .FALSE.,
     +                      BATCH,
     +
     +                                INT1DY )
     +
     +
C
C
C             ...apply (if any) differentiation of the 1DY AB integrals
C                on a nuclear attraction center different from shell
C                centers A and B.
C
C
             PROCEED = DIFFC .AND. (DERCY.GT.0)

             IF (PROCEED) THEN
                 IXC = NUCCEN (1)
                 YC = YN (IXC)
                 NINT = MGIJCEN * (SHELLA+DERAY+1) * (SHELLB+DERBY+1)
                 DO N = 1,DERCY
                    IF (N.GT.1) THEN
                        DO M = 1,NINT
                           TEMP = BATCH (M)
                           BATCH (M) = INT1DY (M)
                           INT1DY (M) = TEMP
                        END DO
                    ELSE
                        DO M = 1,NINT
                           BATCH (M) = INT1DY (M)
                        END DO
                    END IF

                    CALL    OED__NAI_1D_CENDERV_INTEGRALS
     +
     +                           ( MIJ,MGIJCEN,
     +                             NGQP,
     +                             SHELLA+DERAY,SHELLB+DERBY,
     +                             EXP2AB,
     +                             PY,YC,
     +                             N,
     +                             RTS,
     +                             R1Y,
     +                             BATCH,
     +
     +                                      INT1DY )
     +
     +
                 END DO
             END IF
C
C
C             ...perform differentiation (if any) of the 1DY AB
C                integrals on shell center A (and simultaneously on
C                shell center B, if atomic).
C
C
             PROCEED = DIFFA .AND. (DERAY.GT.0)

             IF (PROCEED) THEN
                 DERA = DERAY
                 DERB = DERBY
                 DO N = 1,DERAY
                    NINT = MGIJCEN * (SHELLA+DERA+1) * (SHELLB+DERB+1)
                    DO M = 1,NINT
                       BATCH (M) = INT1DY (M)
                    END DO

                    DERA = DERA - 1
                    DERB = DERB - IXAEQB

                    IF (DERA.LT.0 .OR. DERB.LT.0) THEN
                        WRITE (*,*) ' Problems forming y-derivatives! '
                        WRITE (*,*) ' DERA,DERB = ',DERA,DERB
                        WRITE (*,*) ' oed__nai_derv_3cen_pcgto_block '
                        STOP
                    END IF

                    CALL    OED__NAI_1D_SHDERV_INTEGRALS
     +
     +                           ( MGIJCEN,
     +                             SHELLA+DERA,SHELLB+DERB,
     +                             1,IXAEQB,
     +                             EXP2A,EXP2B,
     +                             BATCH,
     +
     +                                      INT1DY )
     +
     +
                 END DO
             END IF
C
C
C             ...if not atomic, perform differentiation (if any) of
C                the 1DY AB integrals on the remaining shell center B.
C
C
             PROCEED = (.NOT.ATOMIC) .AND. DIFFB .AND. (DERBY.GT.0)

             IF (PROCEED) THEN
                 DERB = DERBY
                 DO N = 1,DERBY
                    NINT = MGIJCEN * (SHELLA+1) * (SHELLB+DERB+1)
                    DO M = 1,NINT
                       BATCH (M) = INT1DY (M)
                    END DO

                    DERB = DERB - 1

                    IF (DERB.LT.0) THEN
                        WRITE (*,*) ' Problems forming y-derivatives! '
                        WRITE (*,*) ' DERB = ',DERB
                        WRITE (*,*) ' oed__nai_derv_3cen_pcgto_block '
                        STOP
                    END IF

                    CALL    OED__NAI_1D_SHDERV_INTEGRALS
     +
     +                           ( MGIJCEN,
     +                             SHELLA,SHELLB+DERB,
     +                             0,1,
     +                             EXP2A,EXP2B,
     +                             BATCH,
     +
     +                                      INT1DY )
     +
     +
                 END DO
             END IF

         ELSE

             CALL    OED__NAI_1D_AB_INTEGRALS
     +
     +                    ( SHELLP,MIN (SHELLA,SHELLB),
     +                      SHELLA,SHELLB,
     +                      MGIJCEN,
     +                      WTS,
     +                      R1Y,R2,
     +                      ABY,
     +                      .FALSE.,
     +                      BATCH,
     +
     +                                INT1DY )
     +
     +
         END IF
C
C
C             ...assemble finally the 1DZ AB integrals and perform
C                differentiation sequence, if necessary.
C
C
         IF (DIFFZ) THEN

             CALL    OED__NAI_1D_AB_INTEGRALS
     +
     +                    ( SHELLP+DERAZ+DERBZ,
     +                      MIN (SHELLA+DERAZ,SHELLB+DERBZ),
     +                      SHELLA+DERAZ,SHELLB+DERBZ,
     +                      MGIJCEN,
     +                      WTS,
     +                      R1Z,R2,
     +                      ABZ,
     +                      .FALSE.,
     +                      BATCH,
     +
     +                                INT1DZ )
     +
     +
C
C
C             ...apply (if any) differentiation of the 1DZ AB integrals
C                on a nuclear attraction center different from shell
C                centers A and B.
C
C
             PROCEED = DIFFC .AND. (DERCZ.GT.0)

             IF (PROCEED) THEN
                 IXC = NUCCEN (1)
                 ZC = ZN (IXC)
                 NINT = MGIJCEN * (SHELLA+DERAZ+1) * (SHELLB+DERBZ+1)
                 DO N = 1,DERCZ
                    IF (N.GT.1) THEN
                        DO M = 1,NINT
                           TEMP = BATCH (M)
                           BATCH (M) = INT1DZ (M)
                           INT1DZ (M) = TEMP
                        END DO
                    ELSE
                        DO M = 1,NINT
                           BATCH (M) = INT1DZ (M)
                        END DO
                    END IF

                    CALL    OED__NAI_1D_CENDERV_INTEGRALS
     +
     +                           ( MIJ,MGIJCEN,
     +                             NGQP,
     +                             SHELLA+DERAZ,SHELLB+DERBZ,
     +                             EXP2AB,
     +                             PZ,ZC,
     +                             N,
     +                             RTS,
     +                             R1Z,
     +                             BATCH,
     +
     +                                      INT1DZ )
     +
     +
                 END DO
             END IF
C
C
C             ...perform differentiation (if any) of the 1DZ AB
C                integrals on shell center A (and simultaneously on
C                shell center B, if atomic).
C
C
             PROCEED = DIFFA .AND. (DERAZ.GT.0)

             IF (PROCEED) THEN
                 DERA = DERAZ
                 DERB = DERBZ
                 DO N = 1,DERAZ
                    NINT = MGIJCEN * (SHELLA+DERA+1) * (SHELLB+DERB+1)
                    DO M = 1,NINT
                       BATCH (M) = INT1DZ (M)
                    END DO

                    DERA = DERA - 1
                    DERB = DERB - IXAEQB

                    IF (DERA.LT.0 .OR. DERB.LT.0) THEN
                        WRITE (*,*) ' Problems forming z-derivatives! '
                        WRITE (*,*) ' DERA,DERB = ',DERA,DERB
                        WRITE (*,*) ' oed__nai_derv_3cen_pcgto_block '
                        STOP
                    END IF

                    CALL    OED__NAI_1D_SHDERV_INTEGRALS
     +
     +                           ( MGIJCEN,
     +                             SHELLA+DERA,SHELLB+DERB,
     +                             1,IXAEQB,
     +                             EXP2A,EXP2B,
     +                             BATCH,
     +
     +                                      INT1DZ )
     +
     +
                 END DO
             END IF
C
C
C             ...if not atomic, perform differentiation (if any) of
C                the 1DZ AB integrals on the remaining shell center B.
C
C
             PROCEED = (.NOT.ATOMIC) .AND. DIFFB .AND. (DERBZ.GT.0)

             IF (PROCEED) THEN
                 DERB = DERBZ
                 DO N = 1,DERBZ
                    NINT = MGIJCEN * (SHELLA+1) * (SHELLB+DERB+1)
                    DO M = 1,NINT
                       BATCH (M) = INT1DZ (M)
                    END DO

                    DERB = DERB - 1

                    IF (DERB.LT.0) THEN
                        WRITE (*,*) ' Problems forming z-derivatives! '
                        WRITE (*,*) ' DERB = ',DERB
                        WRITE (*,*) ' oed__nai_derv_3cen_pcgto_block '
                        STOP
                    END IF

                    CALL    OED__NAI_1D_SHDERV_INTEGRALS
     +
     +                           ( MGIJCEN,
     +                             SHELLA,SHELLB+DERB,
     +                             0,1,
     +                             EXP2A,EXP2B,
     +                             BATCH,
     +
     +                                      INT1DZ )
     +
     +
                 END DO
             END IF

         ELSE

             CALL    OED__NAI_1D_AB_INTEGRALS
     +
     +                    ( SHELLP,MIN (SHELLA,SHELLB),
     +                      SHELLA,SHELLB,
     +                      MGIJCEN,
     +                      WTS,
     +                      R1Z,R2,
     +                      ABZ,
     +                      .FALSE.,
     +                      BATCH,
     +
     +                                INT1DZ )
     +
     +
         END IF
C
C
C             ...assemble the 1D AB integrals to the [A|B] batch.
C
C
         IF (SHELLP.EQ.0) THEN

             CALL    OED__NAI_DERV_INT1D_TO_00
     +
     +                    ( MIJ,NGQPCEN,MGIJCEN,
     +                      INT1DX,INT1DY,INT1DZ,
     +                      DIFFY,DIFFZ,
     +                      R1X,
     +                      SCALE,
     +
     +                                BATCH )
     +
     +
         ELSE IF (SHELLB.EQ.0) THEN

             CALL    OED__NAI_DERV_INT1D_TO_A0
     +
     +                    ( SHELLA,
     +                      MIJ,NGQPCEN,MGIJCEN,
     +                      NXYZA,
     +                      INT1DX,INT1DY,INT1DZ,
     +                      DIFFY,DIFFZ,
     +                      R1X,R1Y,
     +                      SCALE,
     +
     +                                BATCH )
     +
     +
         ELSE IF (SHELLA.EQ.0) THEN

             CALL    OED__NAI_DERV_INT1D_TO_A0
     +
     +                    ( SHELLB,
     +                      MIJ,NGQPCEN,MGIJCEN,
     +                      NXYZB,
     +                      INT1DX,INT1DY,INT1DZ,
     +                      DIFFY,DIFFZ,
     +                      R1X,R1Y,
     +                      SCALE,
     +
     +                                BATCH )
     +
     +
         ELSE

             CALL    OED__NAI_DERV_INT1D_TO_AB
     +
     +                    ( SHELLA,SHELLB,
     +                      MIJ,NGQPCEN,MGIJCEN,
     +                      NXYZA,NXYZB,
     +                      INT1DX,INT1DY,INT1DZ,
     +                      DIFFY,DIFFZ,
     +                      R1X,R1Y,
     +                      SCALE,
     +
     +                                BATCH )
     +
     +
         END IF
C
C
C             ...ready!
C
C
         RETURN
         END
