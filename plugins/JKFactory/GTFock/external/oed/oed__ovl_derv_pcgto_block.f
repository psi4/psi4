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
         SUBROUTINE  OED__OVL_DERV_PCGTO_BLOCK
     +
     +                    ( NBATCH,
     +                      NINT1DX,NINT1DY,NINT1DZ,
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
     +                      PAX,PAY,PAZ,PINVHF,SCALE,
     +                      EXP2A,EXP2B,
     +                      INT1DX,INT1DY,INT1DZ,
     +
     +                                BATCH )
     +
C------------------------------------------------------------------------
C  OPERATION   : OED__OVL_E0_PCGTO_BLOCK
C  MODULE      : ONE ELECTRON INTEGRALS DIRECT
C  MODULE-ID   : OED
C  SUBROUTINES : OED__OVL_1D_AB_INTEGRALS
C                OED__OVL_1D_DERV_INTEGRALS
C                OED__OVL_DERV_INT1D_TO_00
C                OED__OVL_DERV_INT1D_TO_A0
C                OED__OVL_DERV_INT1D_TO_AB
C  DESCRIPTION : This operation calculates a batch of derivated unnormed
C                overlap integrals between primitive cartesian gaussians
C
C                                  [A|B]
C                                       ij
C
C                for a block of ij exponent pairs. The total number
C                of overlap integrals generated here is thus given by
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
C                    INT1Dx       =  will hold all current derivated
C                                    1D integrals for each cartesian
C                                    component (x = X,Y,Z) during all
C                                    stages of differentiation
C
C                  Output:
C
C                    BATCH        =  current batch of primitive
C                                    derivative cartesian overlap
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
         INTEGER     NINT,NINT1DX,NINT1DY,NINT1DZ
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
         DOUBLE PRECISION  ZERO,HALF,ONE,ONEP5

         DOUBLE PRECISION  ALPHAA  (1:NPGTOA)
         DOUBLE PRECISION  ALPHAB  (1:NPGTOB)

         DOUBLE PRECISION  NORMA   (1:NPGTOA)
         DOUBLE PRECISION  NORMB   (1:NPGTOB)

         DOUBLE PRECISION  RHOAB   (1:NIJ)

         DOUBLE PRECISION  BATCH   (1:NBATCH)

         DOUBLE PRECISION  PAX     (1:MIJ)
         DOUBLE PRECISION  PAY     (1:MIJ)
         DOUBLE PRECISION  PAZ     (1:MIJ)
         DOUBLE PRECISION  PINVHF  (1:MIJ)
         DOUBLE PRECISION  SCALE   (1:MIJ)
         DOUBLE PRECISION  EXP2A   (1:MIJ)
         DOUBLE PRECISION  EXP2B   (1:MIJ)

         DOUBLE PRECISION  INT1DX  (1:NINT1DX)
         DOUBLE PRECISION  INT1DY  (1:NINT1DY)
         DOUBLE PRECISION  INT1DZ  (1:NINT1DZ)

         PARAMETER  (ZERO  = 0.D0)
         PARAMETER  (HALF  = 0.5D0)
         PARAMETER  (ONE   = 1.D0)
         PARAMETER  (ONEP5 = 1.5D0)
C
C
C------------------------------------------------------------------------
C
C
C             ...calculate the quantities needed to establish the
C                1D overlap integrals.
C
C
         M = 0
         DO IJ = NIJBEG,NIJEND
            M = M + 1
            I = PRIMA (M)
            J = PRIMB (M)
            EXPA = ALPHAA (I)
            EXPB = ALPHAB (J)
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
C             ...calculate the derivative 1D integral coefficents
C                (if needed). These are the double exponent values
C                for differentiation on centers A and/or B.
C
C
         IF (DIFFA) THEN
             DO IJ = 1,MIJ
                I = PRIMA (IJ)
                EXP2A (IJ) = ALPHAA (I) + ALPHAA (I)
             END DO
         END IF

         IF (DIFFB) THEN
             DO IJ = 1,MIJ
                J = PRIMB (IJ)
                EXP2B (IJ) = ALPHAB (J) + ALPHAB (J)
             END DO
         END IF
C
C
C             ...start assembling the 1D AB overlap integrals and
C                their differentiation. Assemble first the 1DX AB
C                integrals and perform differentiation sequence,
C                if necessary.
C
C
         IF (DIFFX) THEN

             DERA = DERAX
             DERB = DERBX

             CALL    OED__OVL_1D_AB_INTEGRALS
     +
     +                    ( SHELLP+DERA+DERB,
     +                      MIN (SHELLA+DERA,SHELLB+DERB),
     +                      SHELLA+DERA,SHELLB+DERB,
     +                      MIJ,
     +                      PAX,
     +                      PINVHF,
     +                      ABX,
     +                      BATCH,
     +
     +                                INT1DX )
     +
     +
             DO N = 1,NDERX

                NINT = MIJ * (SHELLA+DERA+1) * (SHELLB+DERB+1)

                DO M = 1,NINT
                   BATCH (M) = INT1DX (M)
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
                    WRITE (*,*) ' Problems forming x-derivatives! '
                    WRITE (*,*) ' DERA,DERB = ',DERA,DERB
                    WRITE (*,*) ' oed__ovl_derv_pcgto_block '
                    STOP
                END IF

                CALL    OED__OVL_1D_DERV_INTEGRALS
     +
     +                       ( SHELLA+DERA,SHELLB+DERB,
     +                         DEROPA,DEROPB,
     +                         MIJ,
     +                         EXP2A,EXP2B,
     +                         BATCH,
     +
     +                                  INT1DX )
     +
     +
             END DO

         ELSE

             CALL    OED__OVL_1D_AB_INTEGRALS
     +
     +                    ( SHELLP,MIN (SHELLA,SHELLB),
     +                      SHELLA,SHELLB,
     +                      MIJ,
     +                      PAX,
     +                      PINVHF,
     +                      ABX,
     +                      BATCH,
     +
     +                                INT1DX )
     +
     +
         END IF
C
C
C             ...assemble the 1DY AB integrals and perform
C                differentiation sequence, if necessary.
C
C
         IF (DIFFY) THEN

             DERA = DERAY
             DERB = DERBY

             CALL    OED__OVL_1D_AB_INTEGRALS
     +
     +                    ( SHELLP+DERA+DERB,
     +                      MIN (SHELLA+DERA,SHELLB+DERB),
     +                      SHELLA+DERA,SHELLB+DERB,
     +                      MIJ,
     +                      PAY,
     +                      PINVHF,
     +                      ABY,
     +                      BATCH,
     +
     +                                INT1DY )
     +
     +
             DO N = 1,NDERY

                NINT = MIJ * (SHELLA+DERA+1) * (SHELLB+DERB+1)

                DO M = 1,NINT
                   BATCH (M) = INT1DY (M)
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
                    WRITE (*,*) ' Problems forming y-derivatives! '
                    WRITE (*,*) ' DERA,DERB = ',DERA,DERB
                    WRITE (*,*) ' oed__ovl_derv_pcgto_block '
                    STOP
                END IF

                CALL    OED__OVL_1D_DERV_INTEGRALS
     +
     +                       ( SHELLA+DERA,SHELLB+DERB,
     +                         DEROPA,DEROPB,
     +                         MIJ,
     +                         EXP2A,EXP2B,
     +                         BATCH,
     +
     +                                  INT1DY )
     +
     +
             END DO

         ELSE

             CALL    OED__OVL_1D_AB_INTEGRALS
     +
     +                    ( SHELLP,MIN (SHELLA,SHELLB),
     +                      SHELLA,SHELLB,
     +                      MIJ,
     +                      PAY,
     +                      PINVHF,
     +                      ABY,
     +                      BATCH,
     +
     +                                INT1DY )
     +
     +
         END IF
C
C
C             ...assemble the 1DZ AB integrals and perform
C                differentiation sequence, if necessary.
C
C
         IF (DIFFZ) THEN

             DERA = DERAZ
             DERB = DERBZ

             CALL    OED__OVL_1D_AB_INTEGRALS
     +
     +                    ( SHELLP+DERA+DERB,
     +                      MIN (SHELLA+DERA,SHELLB+DERB),
     +                      SHELLA+DERA,SHELLB+DERB,
     +                      MIJ,
     +                      PAZ,
     +                      PINVHF,
     +                      ABZ,
     +                      BATCH,
     +
     +                                INT1DZ )
     +
     +
             DO N = 1,NDERZ

                NINT = MIJ * (SHELLA+DERA+1) * (SHELLB+DERB+1)

                DO M = 1,NINT
                   BATCH (M) = INT1DZ (M)
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
                    WRITE (*,*) ' Problems forming z-derivatives! '
                    WRITE (*,*) ' DERA,DERB = ',DERA,DERB
                    WRITE (*,*) ' oed__ovl_derv_pcgto_block '
                    STOP
                END IF

                CALL    OED__OVL_1D_DERV_INTEGRALS
     +
     +                       ( SHELLA+DERA,SHELLB+DERB,
     +                         DEROPA,DEROPB,
     +                         MIJ,
     +                         EXP2A,EXP2B,
     +                         BATCH,
     +
     +                                  INT1DZ )
     +
     +
             END DO

         ELSE

             CALL    OED__OVL_1D_AB_INTEGRALS
     +
     +                    ( SHELLP,MIN (SHELLA,SHELLB),
     +                      SHELLA,SHELLB,
     +                      MIJ,
     +                      PAZ,
     +                      PINVHF,
     +                      ABZ,
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

             CALL    OED__OVL_DERV_INT1D_TO_00
     +
     +                    ( MIJ,
     +                      INT1DX,INT1DY,INT1DZ,
     +                      DIFFX,DIFFY,DIFFZ,
     +                      SCALE,
     +
     +                                BATCH )
     +
     +
         ELSE IF (SHELLB.EQ.0) THEN

             CALL    OED__OVL_DERV_INT1D_TO_A0
     +
     +                    ( SHELLA,
     +                      MIJ,
     +                      NXYZA,
     +                      INT1DX,INT1DY,INT1DZ,
     +                      DIFFX,DIFFY,DIFFZ,
     +                      PAX,PAY,
     +                      SCALE,
     +
     +                                BATCH )
     +
     +
         ELSE IF (SHELLA.EQ.0) THEN

             CALL    OED__OVL_DERV_INT1D_TO_A0
     +
     +                    ( SHELLB,
     +                      MIJ,
     +                      NXYZB,
     +                      INT1DX,INT1DY,INT1DZ,
     +                      DIFFX,DIFFY,DIFFZ,
     +                      PAX,PAY,
     +                      SCALE,
     +
     +                                BATCH )
     +
     +
         ELSE

             CALL    OED__OVL_DERV_INT1D_TO_AB
     +
     +                    ( SHELLA,SHELLB,
     +                      MIJ,
     +                      NXYZA,NXYZB,
     +                      INT1DX,INT1DY,INT1DZ,
     +                      DIFFX,DIFFY,DIFFZ,
     +                      PAX,PAY,
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
