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
         SUBROUTINE  OED__NAI_E0_PCGTO_BLOCK
     +
     +                    ( NBATCH,NINT1D,
     +                      ATOMIC,
     +                      MIJ,NCEN,MIJCEN,
     +                      NIJ,NIJBEG,NIJEND,
     +                      NGQP,NMOM,NGQSCR,MGIJCEN,
     +                      NPGTOA,NPGTOB,
     +                      NXYZET,NXYZP,
     +                      SHELLA,SHELLP,
     +                      XA,YA,ZA,XB,YB,ZB,
     +                      ABX,ABY,ABZ,
     +                      NUCLEI,
     +                      XN,YN,ZN,NCHARGE,
     +                      NUCCEN,
     +                      ALPHAA,ALPHAB,
     +                      FTABLE,MGRID,NGRID,TMAX,TSTEP,TVSTEP,
     +                      PRIMA,PRIMB,
     +                      NORMA,NORMB,
     +                      RHOAB,
     +                      PX,PY,PZ,PAX,PAY,PAZ,PINVHF,SCALE,
     +                      RTS,WTS,GQSCR,TVAL,
     +                      R1X,R1Y,R1Z,R2,
     +                      INT1DX,INT1DY,INT1DZ,
     +
     +                                BATCH )
     +
C------------------------------------------------------------------------
C  OPERATION   : OED__NAI_E0_PCGTO_BLOCK
C  MODULE      : ONE ELECTRON INTEGRALS DIRECT
C  MODULE-ID   : OED
C  SUBROUTINES : OED__RYS_ROOTS_WEIGHTS
C                OED__NAI_1D_COEFFICIENTS
C                OED__NAI_1D_P_INTEGRALS
C                OED__NAI_INT1D_TO_E0
C  DESCRIPTION : This operation calculates a batch of unnormed nuclear
C                attraction integrals between primitive cartesian
C                gaussians for the shell range:
C
C                              [E|0]   , E = A to T = A + B
C                                   ij
C
C                and the block of ij exponent pairs. The total number
C                of nuclear attraction integrals generated here is thus
C                given by the total number of cartesian monomials NXYZET
C                times the total number of exponent pairs MIJ in the
C                present block.
C
C                On exit, the batch elements will be stored as:
C
C                             batch (ij,nxyzet)
C
C
C                  Input:
C
C                    NBATCH       =  size of the primitive cartesian
C                                    integral batch
C                    NINT1D       =  space needed for each of the 1D
C                                    X,Y,Z integral arrays
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
C                    NXYZET       =  sum of # of cartesian monomials
C                                    for all shells in the range
C                                    E = A,...,P=A+B
C                    NXYZP        =  # of cartesian monomials for
C                                    the P=A+B shell
C                    SHELLx       =  the shell type for contraction
C                                    shells x = A and P=A+B
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
C                    INT1Dx       =  will hold all current 1D integrals
C                                    for each cartesian component
C                                    (x = X,Y,Z)
C
C                  Output:
C
C                    BATCH        =  current batch of primitive
C                                    cartesian nuclear attraction
C                                    [E|0] integrals
C
C
C  AUTHOR      : Norbert Flocke
C------------------------------------------------------------------------
C
C             ...include files and declare variables.
C
C
         IMPLICIT    NONE

         LOGICAL     ATOMIC

         INTEGER     G000,G010,G020,G030,G040,G050,G060
         INTEGER     I,J,L,M,N
         INTEGER     IXC
         INTEGER     IJ
         INTEGER     MGRID,NGRID
         INTEGER     MIJ,MIJCEN
         INTEGER     NBATCH,NINT1D
         INTEGER     NC,NCEN
         INTEGER     NG,NGQP,NGQPCEN,NMOM,NGQSCR,MGIJCEN
         INTEGER     NIJ,NIJBEG,NIJEND
         INTEGER     NPGTOA,NPGTOB
         INTEGER     NUCLEI
         INTEGER     NXYZET,NXYZP
         INTEGER     SHELLA,SHELLP

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
         DOUBLE PRECISION  SUM
         DOUBLE PRECISION  TMAX,TSTEP,TVSTEP
         DOUBLE PRECISION  XA,YA,ZA,XB,YB,ZB
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
         DOUBLE PRECISION  RTS     (1:MGIJCEN)
         DOUBLE PRECISION  SCALE   (1:MGIJCEN)
         DOUBLE PRECISION  WTS     (1:MGIJCEN)

         DOUBLE PRECISION  INT1DX  (1:NINT1D)
         DOUBLE PRECISION  INT1DY  (1:NINT1D)
         DOUBLE PRECISION  INT1DZ  (1:NINT1D)

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
C             ...perform the following steps:
C
C                1) construct all 1D x,y,z VRR coefficients for all
C                   ij pairs, nuclear centers and quadrature points.
C
C                2) construct all 1D x,y,z nuclear attraction integrals
C                   for all ij pairs, nuclear centers and quadrature
C                   points.
C
C                3) assemble the complete [E|0] overlap batch for all
C                   ij pairs using the 1D integrals. Arrays R1X and R1Y
C                   are passed as scratch arrays.
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
         CALL    OED__NAI_1D_P_INTEGRALS
     +
     +                ( SHELLP,
     +                  MGIJCEN,
     +                  WTS,
     +                  R1X,R1Y,R1Z,
     +                  R2,
     +
     +                            INT1DX,
     +                            INT1DY,
     +                            INT1DZ )
     +
     +
         IF (SHELLP.EQ.0) THEN

             N = 0
             DO J = 1,MIJ
                SUM = ZERO
                DO I = 1,NCEN
                   SUM = SUM + SCALE (N+I) * INT1DX (N+I)
                END DO
                BATCH (J) = SUM
                N = N + NCEN
             END DO

         ELSE
             NGQPCEN = NGQP * NCEN

             CALL    OED__NAI_INT1D_TO_E0
     +
     +                    ( SHELLA,SHELLP,
     +                      MIJ,NGQPCEN,MGIJCEN,
     +                      NXYZET,NXYZP,
     +                      INT1DX,INT1DY,INT1DZ,
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
