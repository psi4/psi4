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
         SUBROUTINE  OED__KIN_AB_PCGTO_BLOCK
     +
     +                    ( NBATCH,NKIN1D,NOVL1D,
     +                      ATOMIC,
     +                      MIJ,NIJ,NIJBEG,NIJEND,
     +                      NPGTOA,NPGTOB,
     +                      NXYZA,NXYZB,
     +                      SHELLA,SHELLB,SHELLP,
     +                      XA,YA,ZA,XB,YB,ZB,
     +                      ABX,ABY,ABZ,
     +                      ALPHAA,ALPHAB,
     +                      PRIMA,PRIMB,
     +                      NORMA,NORMB,
     +                      RHOAB,
     +                      EA,EB,E2AB,
     +                      PAX,PAY,PAZ,PINVHF,SCALE,
     +                      KIN1DX,KIN1DY,KIN1DZ,
     +                      OVL1DX,OVL1DY,OVL1DZ,
     +
     +                                BATCH )
     +
C------------------------------------------------------------------------
C  OPERATION   : OED__KIN_AB_PCGTO_BLOCK
C  MODULE      : ONE ELECTRON INTEGRALS DIRECT
C  MODULE-ID   : OED
C  SUBROUTINES : OED__OVL_1D_INTEGRALS
C                OED__KIN_1D_INTEGRALS
C                OED__KIN_INT1D_TO_A0
C                OED__KIN_INT1D_TO_AB
C  DESCRIPTION : This operation calculates a batch of unnormed kinetic
C                integrals between primitive cartesian gaussians for
C                shells A and B:
C
C                                   [A|B]
C                                        ij
C
C                and the block of ij exponent pairs. The total number
C                of kinetic integrals generated here is thus given by
C                the total number of cartesian monomials NXYZA*NXYZB
C                times the total number of exponent pairs MIJ in the
C                present block.
C
C                On exit, the batch elements will be stored as:
C
C                             batch (ij,nxyza,nxyzb)
C
C
C                  Input:
C
C                    NBATCH       =  size of the primitive cartesian
C                                    integral batch
C                    NKIN1D       =  space needed for each of the
C                                    kinetic 1D X,Y,Z integral arrays
C                    NOVL1D       =  space needed for each of the
C                                    kinetic 1D X,Y,Z integral arrays
C                    ATOMIC       =  indicates, if purely atomic
C                                    integrals will be evaluated
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
C                                    contraction shells x = A,B
C                    SHELLx       =  the shell type for contraction
C                                    shells x = A,B and P=A+B
C                    Xx,Yx,Zx     =  the x,y,z-coordinates for centers
C                                    x = A,B
C                    ABm          =  the m=x,y,z-coordinate differences
C                                    between centers A and B
C                    ALPHAx       =  the primitive exponents for
C                                    contraction shells x = A,B
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
C                    KIN1Dx       =  will hold all current kinetic
C                                    1D integrals for each cartesian
C                                    component (x = X,Y,Z)
C                    OVL1Dx       =  will hold all current overlap
C                                    1D integrals for each cartesian
C                                    component (x = X,Y,Z)
C
C                  Output:
C
C                    BATCH        =  current batch of primitive
C                                    cartesian kinetic [A|B] integrals
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

         INTEGER     I,J,M
         INTEGER     IJ
         INTEGER     MIJ
         INTEGER     NBATCH,NKIN1D,NOVL1D
         INTEGER     NIJ,NIJBEG,NIJEND
         INTEGER     NPGTOA,NPGTOB
         INTEGER     NXYZA,NXYZB
         INTEGER     SHELLA,SHELLB,SHELLP

         INTEGER     PRIMA (1:MIJ)
         INTEGER     PRIMB (1:MIJ)

         DOUBLE PRECISION  ABX,ABY,ABZ
         DOUBLE PRECISION  EXPA,EXPB
         DOUBLE PRECISION  PINV,PVAL
         DOUBLE PRECISION  X,Y,Z
         DOUBLE PRECISION  XA,YA,ZA,XB,YB,ZB
         DOUBLE PRECISION  HALF,ONE,ONEP5,TWO

         DOUBLE PRECISION  ALPHAA  (1:NPGTOA)
         DOUBLE PRECISION  ALPHAB  (1:NPGTOB)

         DOUBLE PRECISION  NORMA   (1:NPGTOA)
         DOUBLE PRECISION  NORMB   (1:NPGTOB)

         DOUBLE PRECISION  RHOAB   (1:NIJ)

         DOUBLE PRECISION  BATCH   (1:NBATCH)

         DOUBLE PRECISION  EA      (1:MIJ)
         DOUBLE PRECISION  EB      (1:MIJ)
         DOUBLE PRECISION  E2AB    (1:MIJ)
         DOUBLE PRECISION  PAX     (1:MIJ)
         DOUBLE PRECISION  PAY     (1:MIJ)
         DOUBLE PRECISION  PAZ     (1:MIJ)
         DOUBLE PRECISION  PINVHF  (1:MIJ)
         DOUBLE PRECISION  SCALE   (1:MIJ)

         DOUBLE PRECISION  KIN1DX  (1:NKIN1D)
         DOUBLE PRECISION  KIN1DY  (1:NKIN1D)
         DOUBLE PRECISION  KIN1DZ  (1:NKIN1D)
         DOUBLE PRECISION  OVL1DX  (1:NOVL1D)
         DOUBLE PRECISION  OVL1DY  (1:NOVL1D)
         DOUBLE PRECISION  OVL1DZ  (1:NOVL1D)

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
C                1D overlap and kinetic integrals.
C
C
         IF (ATOMIC) THEN
             M = 0
             DO 100 IJ = NIJBEG,NIJEND
                M = M + 1
                I = PRIMA (M)
                J = PRIMB (M)
                EXPA = ALPHAA (I)
                EXPB = ALPHAB (J)
                EA (M) = EXPA
                EB (M) = EXPB
                E2AB (M) = TWO * EXPA * EXPB
                PINV = ONE / (EXPA + EXPB)
                PINVHF (M) = HALF * PINV
                SCALE  (M) = (PINV ** ONEP5) * NORMA (I) * NORMB (J)
  100        CONTINUE
         ELSE
             M = 0
             DO 110 IJ = NIJBEG,NIJEND
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
     +                       * NORMA (I) * NORMB (J) * RHOAB (IJ)
  110        CONTINUE
         END IF
C
C
C             ...perform the following steps:
C
C                1) construct all 1D x,y,z overlap integrals for
C                   all ij pairs for up to shells SHELLA+1 and
C                   SHELLB+1.
C
C                2) using the 1D x,y,z overlap integrals construct
C                   all 1D x,y,z kinetic integrals for all ij pairs.
C
C                3) assemble the complete [A|B] kinetic batch for all
C                   ij pairs using the 1D kinetic and overlap integrals.
C                   Arrays EA and EB are passed as scratch arrays.
C
C
         CALL    OED__OVL_1D_INTEGRALS
     +
     +                ( SHELLP+2,SHELLB+1,
     +                  ATOMIC,
     +                  MIJ,
     +                  PAX,PAY,PAZ,
     +                  PINVHF,
     +                  ABX,ABY,ABZ,
     +
     +                            OVL1DX,
     +                            OVL1DY,
     +                            OVL1DZ )
     +
     +
         CALL    OED__KIN_1D_INTEGRALS
     +
     +                ( SHELLA,SHELLB,SHELLP,
     +                  ATOMIC,
     +                  MIJ,
     +                  EA,EB,E2AB,
     +                  OVL1DX,OVL1DY,OVL1DZ,
     +
     +                            KIN1DX,
     +                            KIN1DY,
     +                            KIN1DZ )
     +
     +
         IF (SHELLP.EQ.0) THEN

             DO M = 1,MIJ
                X = KIN1DX (M)
                Y = KIN1DY (M)
                Z = KIN1DZ (M)
                BATCH (M) = SCALE (M) * (X + Y + Z)
             END DO

         ELSE IF (SHELLB.EQ.0) THEN

             CALL    OED__KIN_INT1D_TO_A0
     +
     +                    ( SHELLA,
     +                      ATOMIC,
     +                      MIJ,
     +                      NXYZA,
     +                      KIN1DX,KIN1DY,KIN1DZ,
     +                      OVL1DX,OVL1DY,OVL1DZ,
     +                      EA,EB,
     +                      SCALE,
     +
     +                                BATCH )
     +
     +
         ELSE
             CALL    OED__KIN_INT1D_TO_AB
     +
     +                    ( SHELLA,SHELLB,SHELLP,
     +                      ATOMIC,
     +                      MIJ,
     +                      NXYZA,NXYZB,
     +                      KIN1DX,KIN1DY,KIN1DZ,
     +                      OVL1DX,OVL1DY,OVL1DZ,
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
