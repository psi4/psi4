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
         SUBROUTINE  OED__KIN_DERV_PCGTO_BLOCK
     +
     +                    ( NBATCH,
     +                      NKIN1D,
     +                      NOVL1DX,NOVL1DY,NOVL1DZ,
     +                      MIJ,NIJ,NIJBEG,NIJEND,
     +                      NPGTOA,NPGTOB,
     +                      NXYZA,NXYZB,
     +                      SHELLA,SHELLB,SHELLP,
     +                      XA,YA,ZA,XB,YB,ZB,
     +                      ABX,ABY,ABZ,
     +                      NDERX,NDERY,NDERZ,
     +                      DERAX,DERAY,DERAZ,
     +                      DERBX,DERBY,DERBZ,
     +                      DIFFX,DIFFY,DIFFZ,
     +                      DIFFA,DIFFB,
     +                      ALPHAA,ALPHAB,
     +                      CENSQX,CENSQY,CENSQZ,
     +                      PRIMA,PRIMB,
     +                      NORMA,NORMB,
     +                      RHOAB,
     +                      EA,EB,E2A,E2B,E2AB,
     +                      PAX,PAY,PAZ,PINVHF,SCALE,
     +                      KIN1DX,KIN1DY,KIN1DZ,
     +                      OVL1DX,OVL1DY,OVL1DZ,
     +
     +                                BATCH )
     +
C------------------------------------------------------------------------
C  OPERATION   : OED__KIN_DERV_PCGTO_BLOCK
C  MODULE      : ONE ELECTRON INTEGRALS DIRECT
C  MODULE-ID   : OED
C  SUBROUTINES : OED__OVL_1D_AB_INTEGRALS
C                OED__OVL_1D_DERV_INTEGRALS
C                OED__KIN_1D_DERV_INTEGRALS
C                OED__KIN_DERV_INT1D_TO_00
C                OED__KIN_DERV_INT1D_TO_A0
C                OED__KIN_DERV_INT1D_TO_AB
C  DESCRIPTION : This operation calculates a batch of derivated unnormed
C                kinetic integrals between primitive cartesian gaussians
C
C                                  [A|B]
C                                       ij
C
C                for a block of ij exponent pairs. The total number
C                of kinetic integrals generated here is thus given by
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
C                    NKIN1D       =  space needed for each of the
C                                    kinetic 1D X,Y,Z integral arrays
C                    NOVL1Dx      =  space needed for each of the 1D
C                                    overlap x = X,Y,Z integral arrays
C                                    (they might be different due to
C                                    different orders of differentiation
C                                    for each cartesian component)
C                    MIJ          =  current # of ij primitive index
C                                    pairs corresponding to the
C                                    contracted shell pairs A,B
C                    NIJ          =  total # of ij primitive index
C                                    pairs for the contracted shell
C                                    pair A,B
C                    NIJBEG(END)  =  first(last) ij primitive index
C                                    defining the ij block
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
C                    NDERp        =  the total order of differentiation
C                                    with respect to the p = x,y,z
C                                    coordinates
C                    DERxp        =  the order of differentiation on
C                                    centers x = A,B with respect to
C                                    the p = x,y,z coordinates
C                    DIFFp        =  is true, if differentiation will
C                                    be performed on centers p = A,B
C                                    involving the p = x,y,z coordinates
C                    ALPHAx       =  the primitive exponents for
C                                    contraction shells x = A,B
C                    CENSQp       =  derivative center sequence array
C                                    for the p = x,y,z coordinates
C                    PRIMx        =  i,j labels of primitives for the
C                                    respective contraction shells
C                                    x = A,B
C                    NORMx        =  the normalization factors due to
C                                    the primitive exponents for the
C                                    contraction shells x = A,B
C                    RHOAB        =  the complete set of NIJ exponential
C                                    prefactors between contraction
C                                    shells A and B
C                    Ex           =  will hold current MIJ exponents
C                                    from centers x = A,B
C                    E2x          =  will hold current MIJ double
C                                    exponent values from centers
C                                    x = A,B
C                    E2AB         =  will hold current MIJ double
C                                    exponent products between centers
C                                    A and B
C                    PAx          =  will hold current MIJ coordinate
C                                    x=X,Y,Z differences P-A between
C                                    centers P and A
C                    PINVHF       =  will hold current MIJ values of
C                                    1/(2*P), where P are the exponent
C                                    sums for contraction shells A
C                                    and B
C                    SCALE        =  will hold current MIJ values of
C                                    scaling factors
C                    EXP2x        =  will hold current double primitive
C                                    exponent values in MIJ order for
C                                    each contracted shell x = A,B
C                    KIN1Dx       =  will hold all current derivated
C                                    kinetic 1D integrals for each
C                                    cartesian component (x = X,Y,Z)
C                    OVL1Dx       =  will hold all current derivated
C                                    overlap 1D integrals for each
C                                    cartesian component (x = X,Y,Z)
C                                    during all stages of derivation
C
C                  Output:
C
C                    BATCH        =  current batch of primitive
C                                    derivative cartesian kinetic
C                                    [A|B] integrals
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

         LOGICAL     DIFFA,DIFFB
         LOGICAL     DIFFX,DIFFY,DIFFZ

         INTEGER     CENTER
         INTEGER     DERA,DERB
         INTEGER     DERAX,DERAY,DERAZ
         INTEGER     DERBX,DERBY,DERBZ
         INTEGER     DEROPA,DEROPB
         INTEGER     I,J,M,N
         INTEGER     IJ
         INTEGER     MIJ
         INTEGER     NBATCH
         INTEGER     NDERX,NDERY,NDERZ
         INTEGER     NIJ,NIJBEG,NIJEND
         INTEGER     NINT
         INTEGER     NKIN1D
         INTEGER     NOVL1DX,NOVL1DY,NOVL1DZ
         INTEGER     NPGTOA,NPGTOB
         INTEGER     NXYZA,NXYZB
         INTEGER     SHELLA,SHELLB,SHELLP

         INTEGER     CENSQX (0:NDERX)
         INTEGER     CENSQY (0:NDERY)
         INTEGER     CENSQZ (0:NDERZ)

         INTEGER     PRIMA (1:MIJ)
         INTEGER     PRIMB (1:MIJ)

         DOUBLE PRECISION  ABX,ABY,ABZ
         DOUBLE PRECISION  EXPA,EXPB
         DOUBLE PRECISION  PINV,PVAL
         DOUBLE PRECISION  XA,YA,ZA,XB,YB,ZB
         DOUBLE PRECISION  ZERO,HALF,ONE,ONEP5,TWO

         DOUBLE PRECISION  ALPHAA  (1:NPGTOA)
         DOUBLE PRECISION  ALPHAB  (1:NPGTOB)

         DOUBLE PRECISION  NORMA   (1:NPGTOA)
         DOUBLE PRECISION  NORMB   (1:NPGTOB)

         DOUBLE PRECISION  RHOAB   (1:NIJ)

         DOUBLE PRECISION  BATCH   (1:NBATCH)

         DOUBLE PRECISION  EA      (1:MIJ)
         DOUBLE PRECISION  EB      (1:MIJ)
         DOUBLE PRECISION  E2A     (1:MIJ)
         DOUBLE PRECISION  E2B     (1:MIJ)
         DOUBLE PRECISION  E2AB    (1:MIJ)
         DOUBLE PRECISION  PAX     (1:MIJ)
         DOUBLE PRECISION  PAY     (1:MIJ)
         DOUBLE PRECISION  PAZ     (1:MIJ)
         DOUBLE PRECISION  PINVHF  (1:MIJ)
         DOUBLE PRECISION  SCALE   (1:MIJ)

         DOUBLE PRECISION  KIN1DX  (1:NKIN1D)
         DOUBLE PRECISION  KIN1DY  (1:NKIN1D)
         DOUBLE PRECISION  KIN1DZ  (1:NKIN1D)
         DOUBLE PRECISION  OVL1DX  (1:NOVL1DX)
         DOUBLE PRECISION  OVL1DY  (1:NOVL1DY)
         DOUBLE PRECISION  OVL1DZ  (1:NOVL1DZ)

         PARAMETER  (ZERO  = 0.D0)
         PARAMETER  (HALF  = 0.5D0)
         PARAMETER  (ONE   = 1.D0)
         PARAMETER  (ONEP5 = 1.5D0)
         PARAMETER  (TWO   = 2.D0)
C
C
C------------------------------------------------------------------------
C
C
C             ...calculate the quantities needed to establish the
C                1D kinetic integrals.
C
C
         M = 0
         DO IJ = NIJBEG,NIJEND
            M = M + 1
            I = PRIMA (M)
            J = PRIMB (M)
            EXPA = ALPHAA (I)
            EXPB = ALPHAB (J)
            EA (M) = EXPA
            EB (M) = EXPB
            E2AB (M) = TWO * EXPA * EXPB
            PINV = ONE / (EXPA + EXPB)
            PVAL = - EXPB * PINV
            PAX (M) = PVAL * ABX
            PAY (M) = PVAL * ABY
            PAZ (M) = PVAL * ABZ
            PINVHF (M) = HALF * PINV
            SCALE  (M) = (PINV ** ONEP5)
     +                   * NORMA (I) * NORMB (J) * RHOAB (IJ)
         END DO
C
C
C             ...calculate the derivative 1D overlap integral
C                coefficents (if needed). These are the double
C                exponent values for differentiation on centers
C                A and/or B.
C
C
         IF (DIFFA) THEN
             DO IJ = 1,MIJ
                I = PRIMA (IJ)
                E2A (IJ) = ALPHAA (I) + ALPHAA (I)
             END DO
         END IF

         IF (DIFFB) THEN
             DO IJ = 1,MIJ
                J = PRIMB (IJ)
                E2B (IJ) = ALPHAB (J) + ALPHAB (J)
             END DO
         END IF
C
C
C             ...start assembling the 1D A+1,B+1 derivative overlap
C                integrals. Assemble first the starting 1DX integrals
C                and perform differentiation sequence, if necessary.
C
C
         IF (DIFFX) THEN

             DERA = DERAX
             DERB = DERBX

             CALL    OED__OVL_1D_AB_INTEGRALS
     +
     +                    ( SHELLP+DERA+DERB+2,
     +                      MIN (SHELLA+DERA+1,SHELLB+DERB+1),
     +                      SHELLA+DERA+1,SHELLB+DERB+1,
     +                      MIJ,
     +                      PAX,
     +                      PINVHF,
     +                      ABX,
     +                      BATCH,
     +
     +                                OVL1DX )
     +
     +
             DO N = 1,NDERX

                NINT = MIJ * (SHELLA+DERA+2) * (SHELLB+DERB+2)

                DO M = 1,NINT
                   BATCH (M) = OVL1DX (M)
                END DO

                CENTER = CENSQX (N)

                IF (CENTER.EQ.1) THEN
                    DEROPA = 1
                    DEROPB = 0
                ELSE
                    DEROPA = 0
                    DEROPB = 1
                END IF

                DERA = DERA - DEROPA
                DERB = DERB - DEROPB

                IF (DERA.LT.0 .OR. DERB.LT.0) THEN
                    WRITE (*,*) ' Problems forming ovlp x-derivatives! '
                    WRITE (*,*) ' DERA,DERB = ',DERA,DERB
                    WRITE (*,*) ' oed__kin_derv_pcgto_block '
                    STOP
                END IF

                CALL    OED__OVL_1D_DERV_INTEGRALS
     +
     +                       ( SHELLA+DERA+1,SHELLB+DERB+1,
     +                         DEROPA,DEROPB,
     +                         MIJ,
     +                         E2A,E2B,
     +                         BATCH,
     +
     +                                  OVL1DX )
     +
     +
             END DO

         ELSE

             CALL    OED__OVL_1D_AB_INTEGRALS
     +
     +                    ( SHELLP+2,MIN (SHELLA+1,SHELLB+1),
     +                      SHELLA+1,SHELLB+1,
     +                      MIJ,
     +                      PAX,
     +                      PINVHF,
     +                      ABX,
     +                      BATCH,
     +
     +                                OVL1DX )
     +
     +
         END IF
C
C
C             ...assemble the starting 1DY AB integrals and perform
C                differentiation sequence, if necessary.
C
C
         IF (DIFFY) THEN

             DERA = DERAY
             DERB = DERBY

             CALL    OED__OVL_1D_AB_INTEGRALS
     +
     +                    ( SHELLP+DERA+DERB+2,
     +                      MIN (SHELLA+DERA+1,SHELLB+DERB+1),
     +                      SHELLA+DERA+1,SHELLB+DERB+1,
     +                      MIJ,
     +                      PAY,
     +                      PINVHF,
     +                      ABY,
     +                      BATCH,
     +
     +                                OVL1DY )
     +
     +
             DO N = 1,NDERY

                NINT = MIJ * (SHELLA+DERA+2) * (SHELLB+DERB+2)

                DO M = 1,NINT
                   BATCH (M) = OVL1DY (M)
                END DO

                CENTER = CENSQY (N)

                IF (CENTER.EQ.1) THEN
                    DEROPA = 1
                    DEROPB = 0
                ELSE
                    DEROPA = 0
                    DEROPB = 1
                END IF

                DERA = DERA - DEROPA
                DERB = DERB - DEROPB

                IF (DERA.LT.0 .OR. DERB.LT.0) THEN
                    WRITE (*,*) ' Problems forming ovlp y-derivatives! '
                    WRITE (*,*) ' DERA,DERB = ',DERA,DERB
                    WRITE (*,*) ' oed__kin_derv_pcgto_block '
                    STOP
                END IF

                CALL    OED__OVL_1D_DERV_INTEGRALS
     +
     +                       ( SHELLA+DERA+1,SHELLB+DERB+1,
     +                         DEROPA,DEROPB,
     +                         MIJ,
     +                         E2A,E2B,
     +                         BATCH,
     +
     +                                  OVL1DY )
     +
     +
             END DO

         ELSE

             CALL    OED__OVL_1D_AB_INTEGRALS
     +
     +                    ( SHELLP+2,MIN (SHELLA+1,SHELLB+1),
     +                      SHELLA+1,SHELLB+1,
     +                      MIJ,
     +                      PAY,
     +                      PINVHF,
     +                      ABY,
     +                      BATCH,
     +
     +                                OVL1DY )
     +
     +
         END IF
C
C
C             ...assemble the starting 1DZ AB integrals and perform
C                differentiation sequence, if necessary.
C
C
         IF (DIFFZ) THEN

             DERA = DERAZ
             DERB = DERBZ

             CALL    OED__OVL_1D_AB_INTEGRALS
     +
     +                    ( SHELLP+DERA+DERB+2,
     +                      MIN (SHELLA+DERA+1,SHELLB+DERB+1),
     +                      SHELLA+DERA+1,SHELLB+DERB+1,
     +                      MIJ,
     +                      PAZ,
     +                      PINVHF,
     +                      ABZ,
     +                      BATCH,
     +
     +                                OVL1DZ )
     +
     +
             DO N = 1,NDERZ

                NINT = MIJ * (SHELLA+DERA+2) * (SHELLB+DERB+2)

                DO M = 1,NINT
                   BATCH (M) = OVL1DZ (M)
                END DO

                CENTER = CENSQZ (N)

                IF (CENTER.EQ.1) THEN
                    DEROPA = 1
                    DEROPB = 0
                ELSE
                    DEROPA = 0
                    DEROPB = 1
                END IF

                DERA = DERA - DEROPA
                DERB = DERB - DEROPB

                IF (DERA.LT.0 .OR. DERB.LT.0) THEN
                    WRITE (*,*) ' Problems forming ovlp z-derivatives! '
                    WRITE (*,*) ' DERA,DERB = ',DERA,DERB
                    WRITE (*,*) ' oed__kin_derv_pcgto_block '
                    STOP
                END IF

                CALL    OED__OVL_1D_DERV_INTEGRALS
     +
     +                       ( SHELLA+DERA+1,SHELLB+DERB+1,
     +                         DEROPA,DEROPB,
     +                         MIJ,
     +                         E2A,E2B,
     +                         BATCH,
     +
     +                                  OVL1DZ )
     +
     +
             END DO

         ELSE

             CALL    OED__OVL_1D_AB_INTEGRALS
     +
     +                    ( SHELLP+2,MIN (SHELLA+1,SHELLB+1),
     +                      SHELLA+1,SHELLB+1,
     +                      MIJ,
     +                      PAZ,
     +                      PINVHF,
     +                      ABZ,
     +                      BATCH,
     +
     +                                OVL1DZ )
     +
     +
         END IF
C
C
C             ...combine the 1D X,Y,Z A+1,B+1 derivative overlap
C                integrals to the 1D X,Y,Z AB kinetic integrals.
C
C
         CALL    OED__KIN_1D_DERV_INTEGRALS
     +
     +                ( SHELLA,SHELLB,
     +                  MIJ,
     +                  EA,EB,E2AB,
     +                  OVL1DX,OVL1DY,OVL1DZ,
     +
     +                            KIN1DX,
     +                            KIN1DY,
     +                            KIN1DZ )
     +
     +
C
C
C             ...assemble the 1D AB derivative overlap integrals
C                and the 1D AB derivative kinetic integrals to the
C                [A|B] batch.
C
C
         IF (SHELLP.EQ.0) THEN

             CALL    OED__KIN_DERV_INT1D_TO_00
     +
     +                    ( MIJ,
     +                      KIN1DX,KIN1DY,KIN1DZ,
     +                      OVL1DX,OVL1DY,OVL1DZ,
     +                      DIFFX,DIFFY,DIFFZ,
     +                      EA,EB,
     +                      SCALE,
     +
     +                                BATCH )
     +
     +
         ELSE IF (SHELLB.EQ.0) THEN

             CALL    OED__KIN_DERV_INT1D_TO_A0
     +
     +                    ( SHELLA,
     +                      MIJ,
     +                      NXYZA,
     +                      KIN1DX,KIN1DY,KIN1DZ,
     +                      OVL1DX,OVL1DY,OVL1DZ,
     +                      DIFFX,DIFFY,DIFFZ,
     +                      EA,EB,
     +                      SCALE,
     +
     +                                BATCH )
     +
     +
         ELSE
             CALL    OED__KIN_DERV_INT1D_TO_AB
     +
     +                    ( SHELLA,SHELLB,
     +                      MIJ,
     +                      NXYZA,NXYZB,
     +                      KIN1DX,KIN1DY,KIN1DZ,
     +                      OVL1DX,OVL1DY,OVL1DZ,
     +                      DIFFX,DIFFY,DIFFZ,
     +                      EA,EB,
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
