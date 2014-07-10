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
         SUBROUTINE  OED__NAI_DERV_2CEN_PCGTO_BLOCK
     +
     +                    ( NBATCH,
     +                      NINT1DX,NINT1DY,NINT1DZ,
     +                      MIJ,
     +                      NIJ,NIJBEG,NIJEND,
     +                      NGQP,NMOM,NGQSCR,MGQPIJ,
     +                      NPGTOA,NPGTOB,
     +                      NXYZA,NXYZB,
     +                      SHELLA,SHELLB,SHELLP,
     +                      XA,YA,ZA,XB,YB,ZB,
     +                      ABX,ABY,ABZ,
     +                      CENTER,CHARGE,
     +                      NDERX,NDERY,NDERZ,
     +                      DERAX,DERAY,DERAZ,
     +                      DERBX,DERBY,DERBZ,
     +                      DIFFX,DIFFY,DIFFZ,
     +                      ALPHAA,ALPHAB,
     +                      FTABLE,MGRID,NGRID,TMAX,TSTEP,TVSTEP,
     +                      PRIMA,PRIMB,
     +                      NORMA,NORMB,
     +                      RHOAB,
     +                      PX,PY,PZ,PAX,PAY,PAZ,PINVHF,SCALE,
     +                      RTS,WTS,GQSCR,TVAL,
     +                      R1X,R1Y,R1Z,R2,
     +                      EXP2A,EXP2B,
     +                      INT1DX,INT1DY,INT1DZ,
     +
     +                                BATCH )
     +
C------------------------------------------------------------------------
C  OPERATION   : OED__NAI_DERV_2CEN_PCGTO_BLOCK
C  MODULE      : ONE ELECTRON INTEGRALS DIRECT
C  MODULE-ID   : OED
C  SUBROUTINES : OED__RYS_ROOTS_WEIGHTS
C                OED__NAI_1D_COEFFICIENTS
C                OED__NAI_1D_AB_INTEGRALS
C                OED__NAI_1D_SHDERV_INTEGRALS
C                OED__NAI_DERV_INT1D_TO_00
C                OED__NAI_DERV_INT1D_TO_A0
C                OED__NAI_DERV_INT1D_TO_AB
C  DESCRIPTION : This operation calculates a batch of derivated unnormed
C                2-center nuclear attraction integrals between cartesian
C                gaussians
C
C                                   [A|C = A or B|B]
C                                                   ij
C
C                for a block of ij exponent pairs. Differentiation is
C                performed entirely either on the A or B shell and the
C                nuclear attraction center C coincides with the center
C                on either the B or A shell, respectively. Note, that
C                the atomic case A = B is not possible here, since
C                this would lead to a purely atomic nuclear attraction
C                integral whose differentiation is always zero to all
C                orders.
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
C                    MIJ          =  current # of ij primitive index
C                                    pairs corresponding to the
C                                    contracted shell pairs A,B
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
C                    MGQPIJ       =  # of roots times # of ij primitive
C                                    index pairs
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
C                    CENTER       =  location of the nuclear attraction
C                                    center. Can be only 'A' or 'B'.
C                    CHARGE       =  the nuclear charge for the center
C                                    specified in variable CENTER
C                    NDERp        =  the total order of differentiation
C                                    with respect to the p = x,y,z
C                                    coordinates
C                    DERyp        =  the order of differentiation on
C                                    shell centers y = A and B
C                                    with respect to the p = x,y,z
C                                    coordinates
C                    DIFFy        =  is true, if differentiation will
C                                    be performed at least once one on
C                                    each of the y = x,y,z coordinates
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
C                    SCALE        =  will hold current distinct MIJ
C                                    (expanded to MGQPIJ) values of
C                                    scaling factors
C                    RTS          =  will hold all current MGQPIJ
C                                    quadrature roots
C                    WTS          =  will hold all current MGQPIJ
C                                    quadrature weights
C                    GQSCR        =  will be used as scratch space
C                                    for determining the quadrature
C                                    roots and weights
C                    TVAL         =  will hold current MIJ values
C                                    of T-exponents defining the Rys
C                                    weight functions
C                    R1x          =  will hold the current MGQPIJ
C                                    VRR R1-coefficients (individual
C                                    cartesian components x=X,Y,Z) for
C                                    shell expansion on center P
C                    R2           =  will hold the current MGQPIJ
C                                    coordinate independent VRR
C                                    R2-coefficients
C                    EXP2x        =  will hold current double primitive
C                                    exponent values in MIJ order
C                                    (expanded to MGQPIJ) for the
C                                    contracted shell x = A or B
C                    INT1Dx       =  will hold all current derivated
C                                    1D integrals for each cartesian
C                                    component (x = X,Y,Z) during all
C                                    stages of differentiation
C
C                  Output:
C
C                    BATCH        =  current batch of primitive case I
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

         CHARACTER*1 CENTER

         LOGICAL     ATOMIC
         LOGICAL     DIFFX,DIFFY,DIFFZ

         INTEGER     DERA,DERB
         INTEGER     DERAX,DERAY,DERAZ
         INTEGER     DERBX,DERBY,DERBZ
         INTEGER     DEROPA,DEROPB
         INTEGER     G000,G010,G020,G030,G040,G050,G060
         INTEGER     I,J,L,M,N
         INTEGER     IJ
         INTEGER     MGRID,NGRID
         INTEGER     MIJ,MGQPIJ
         INTEGER     NBATCH
         INTEGER     NDERX,NDERY,NDERZ
         INTEGER     NG,NGQP,NMOM,NGQSCR
         INTEGER     NIJ,NIJBEG,NIJEND
         INTEGER     NINT,NINT1DX,NINT1DY,NINT1DZ
         INTEGER     NPGTOA,NPGTOB
         INTEGER     NXYZA,NXYZB
         INTEGER     SHELLA,SHELLB,SHELLP

         INTEGER     PRIMA (1:MIJ)
         INTEGER     PRIMB (1:MIJ)

         DOUBLE PRECISION  ABX,ABY,ABZ
         DOUBLE PRECISION  CHARGE
         DOUBLE PRECISION  EXPA,EXPB
         DOUBLE PRECISION  FACTOR
         DOUBLE PRECISION  P,PINV,PVAL
         DOUBLE PRECISION  PARITY
         DOUBLE PRECISION  PCX,PCY,PCZ
         DOUBLE PRECISION  PXVAL,PYVAL,PZVAL
         DOUBLE PRECISION  RNPCSQ
         DOUBLE PRECISION  SCALEM
         DOUBLE PRECISION  TMAX,TSTEP,TVSTEP
         DOUBLE PRECISION  XA,YA,ZA,XB,YB,ZB,XC,YC,ZC
         DOUBLE PRECISION  ZERO,HALF,ONE

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
         DOUBLE PRECISION  TVAL    (1:MIJ)

         DOUBLE PRECISION  R1X     (1:MGQPIJ)
         DOUBLE PRECISION  R1Y     (1:MGQPIJ)
         DOUBLE PRECISION  R1Z     (1:MGQPIJ)
         DOUBLE PRECISION  R2      (1:MGQPIJ)
         DOUBLE PRECISION  EXP2A   (1:MGQPIJ)
         DOUBLE PRECISION  EXP2B   (1:MGQPIJ)
         DOUBLE PRECISION  RTS     (1:MGQPIJ)
         DOUBLE PRECISION  SCALE   (1:MGQPIJ)
         DOUBLE PRECISION  WTS     (1:MGQPIJ)

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
C             ...calculate center dependent stuff.
C
C
         IF (CENTER .EQ. 'A') THEN

             XC = XA
             YC = YA
             ZC = ZA
             DEROPA = 0
             DEROPB = 1

             PARITY = (-ONE) ** (DERAX + DERAY + DERAZ)

             M = 0
             DO IJ = 1,MIJ
                J = PRIMB (IJ)
                EXPB = ALPHAB (J) + ALPHAB (J)
                DO N = 1,NGQP
                   EXP2B (M+N) = EXPB
                END DO
                M = M + NGQP
             END DO

         ELSE IF (CENTER .EQ. 'B') THEN

             XC = XB
             YC = YB
             ZC = ZB
             DEROPA = 1
             DEROPB = 0
             PARITY = (-ONE) ** (DERBX + DERBY + DERBZ)

             M = 0
             DO IJ = 1,MIJ
                I = PRIMA (IJ)
                EXPA = ALPHAA (I) + ALPHAA (I)
                DO N = 1,NGQP
                   EXP2A (M+N) = EXPA
                END DO
                M = M + NGQP
             END DO

         ELSE

             WRITE (*,*) ' Unrecognizable center! '
             WRITE (*,*) ' CENTER = ',CENTER
             WRITE (*,*) ' oed__nai_derv_2cen_pcgto_block '
             WRITE (1,*) ' Unrecognizable center! '
             WRITE (1,*) ' CENTER = ',CENTER
             WRITE (1,*) ' oed__nai_derv_2cen_pcgto_block '
             STOP

         END IF
C
C
C             ...calculate the quantities needed to establish the
C                1D nuclear attraction integrals.
C
C
         FACTOR = PARITY * CHARGE

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
            PCX = PXVAL - XC
            PCY = PYVAL - YC
            PCZ = PZVAL - ZC
            PINVHF (M) = HALF * PINV
            SCALEM = PINV * FACTOR * NORMA (I) * NORMB (J) * RHOAB (IJ)
            RNPCSQ = PCX * PCX + PCY * PCY + PCZ * PCZ
            TVAL (M) = P * RNPCSQ
            DO NG = 1,NGQP
               N = N + 1
               SCALE (N) = SCALEM
            END DO
         END DO
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
     +                ( MIJ,MGQPIJ,
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
C             ...generate all VRR coefficients.
C
C
         ATOMIC = .FALSE.

         CALL    OED__NAI_1D_COEFFICIENTS
     +
     +                ( NGQP,1,
     +                  MIJ,MGQPIJ,
     +                  ATOMIC,
     +                  1,
     +                  XC,YC,ZC,
     +                  1,
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
C             ...start assembling the 2-center 1D AB nuclear attraction
C                integrals and their differentiation. Assemble first the
C                2-center 1DX AB integrals and perform differentiation
C                sequence, if necessary.
C
C
         IF (DIFFX) THEN

             DERA = DEROPA * NDERX
             DERB = DEROPB * NDERX

             CALL    OED__NAI_1D_AB_INTEGRALS
     +
     +                    ( SHELLP+DERA+DERB,
     +                      MIN (SHELLA+DERA,SHELLB+DERB),
     +                      SHELLA+DERA,SHELLB+DERB,
     +                      MGQPIJ,
     +                      WTS,
     +                      R1X,R2,
     +                      ABX,
     +                      .TRUE.,
     +                      BATCH,
     +
     +                                INT1DX )
     +
     +
             DO N = 1,NDERX

                NINT = MGQPIJ * (SHELLA+DERA+1) * (SHELLB+DERB+1)

                DO M = 1,NINT
                   BATCH (M) = INT1DX (M)
                END DO

                DERA = DERA - DEROPA
                DERB = DERB - DEROPB

                CALL    OED__NAI_1D_SHDERV_INTEGRALS
     +
     +                       ( MGQPIJ,
     +                         SHELLA+DERA,SHELLB+DERB,
     +                         DEROPA,DEROPB,
     +                         EXP2A,EXP2B,
     +                         BATCH,
     +
     +                                  INT1DX )
     +
     +
             END DO

         ELSE

             CALL    OED__NAI_1D_AB_INTEGRALS
     +
     +                    ( SHELLP,MIN (SHELLA,SHELLB),
     +                      SHELLA,SHELLB,
     +                      MGQPIJ,
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
C             ...assemble the 2-center 1DY AB integrals and perform
C                differentiation sequence, if necessary.
C
C
         IF (DIFFY) THEN

             DERA = DEROPA * NDERY
             DERB = DEROPB * NDERY

             CALL    OED__NAI_1D_AB_INTEGRALS
     +
     +                    ( SHELLP+DERA+DERB,
     +                      MIN (SHELLA+DERA,SHELLB+DERB),
     +                      SHELLA+DERA,SHELLB+DERB,
     +                      MGQPIJ,
     +                      WTS,
     +                      R1Y,R2,
     +                      ABY,
     +                      .FALSE.,
     +                      BATCH,
     +
     +                                INT1DY )
     +
     +
             DO N = 1,NDERY

                NINT = MGQPIJ * (SHELLA+DERA+1) * (SHELLB+DERB+1)

                DO M = 1,NINT
                   BATCH (M) = INT1DY (M)
                END DO

                DERA = DERA - DEROPA
                DERB = DERB - DEROPB

                CALL    OED__NAI_1D_SHDERV_INTEGRALS
     +
     +                       ( MGQPIJ,
     +                         SHELLA+DERA,SHELLB+DERB,
     +                         DEROPA,DEROPB,
     +                         EXP2A,EXP2B,
     +                         BATCH,
     +
     +                                  INT1DY )
     +
     +
             END DO

         ELSE

             CALL    OED__NAI_1D_AB_INTEGRALS
     +
     +                    ( SHELLP,MIN (SHELLA,SHELLB),
     +                      SHELLA,SHELLB,
     +                      MGQPIJ,
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
C             ...assemble the 2-center 1DZ AB integrals and perform
C                differentiation sequence, if necessary.
C
C
         IF (DIFFZ) THEN

             DERA = DEROPA * NDERZ
             DERB = DEROPB * NDERZ

             CALL    OED__NAI_1D_AB_INTEGRALS
     +
     +                    ( SHELLP+DERA+DERB,
     +                      MIN (SHELLA+DERA,SHELLB+DERB),
     +                      SHELLA+DERA,SHELLB+DERB,
     +                      MGQPIJ,
     +                      WTS,
     +                      R1Z,R2,
     +                      ABZ,
     +                      .FALSE.,
     +                      BATCH,
     +
     +                                INT1DZ )
     +
     +
             DO N = 1,NDERZ

                NINT = MGQPIJ * (SHELLA+DERA+1) * (SHELLB+DERB+1)

                DO M = 1,NINT
                   BATCH (M) = INT1DZ (M)
                END DO

                DERA = DERA - DEROPA
                DERB = DERB - DEROPB

                CALL    OED__NAI_1D_SHDERV_INTEGRALS
     +
     +                       ( MGQPIJ,
     +                         SHELLA+DERA,SHELLB+DERB,
     +                         DEROPA,DEROPB,
     +                         EXP2A,EXP2B,
     +                         BATCH,
     +
     +                                  INT1DZ )
     +
     +
             END DO

         ELSE

             CALL    OED__NAI_1D_AB_INTEGRALS
     +
     +                    ( SHELLP,MIN (SHELLA,SHELLB),
     +                      SHELLA,SHELLB,
     +                      MGQPIJ,
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
C             ...assemble the 2-center 1D AB integrals to the [A|B]
C                batch.
C
C
         IF (SHELLP.EQ.0) THEN

             CALL    OED__NAI_DERV_INT1D_TO_00
     +
     +                    ( MIJ,NGQP,MGQPIJ,
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
     +                      MIJ,NGQP,MGQPIJ,
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
     +                      MIJ,NGQP,MGQPIJ,
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
     +                      MIJ,NGQP,MGQPIJ,
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
