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
         SUBROUTINE  OED__KIN_DERV_INT1D_TO_00
     +
     +                    ( NEXP,
     +                      KIN1DX,KIN1DY,KIN1DZ,
     +                      OVL1DX,OVL1DY,OVL1DZ,
     +                      DIFFX,DIFFY,DIFFZ,
     +                      TEMP1,TEMP2,
     +                      SCALE,
     +
     +                                BATCH )
     +
C------------------------------------------------------------------------
C  OPERATION   : OED__KIN_DERV_INT1D_TO_00
C  MODULE      : ONE ELECTRON INTEGRALS DIRECT
C  MODULE-ID   : OED
C  SUBROUTINES : none
C  DESCRIPTION : This routine assembles the set of batches of cartesian
C                kinetic derivative integrals [0|0], adding up all
C                contributions from all the kinetic and overlap 1D
C                derivative integrals:
C
C                        X = KIN (x) * OVL (y) * OVL (z)
C                        Y = OVL (x) * KIN (y) * OVL (z)
C                        Z = OVL (x) * OVL (y) * KIN (z)
C
C                       [0|0] = [x=0,y=0,z=0] = X + Y + Z
C
C
C                  Input:
C
C                    NEXP        =  current # of exponent pairs
C                    KIN1Dx      =  all current kinetic 1D derivative
C                                   integrals for each cartesian
C                                   component x = X,Y,Z
C                    OVL1Dx      =  all current overlap 1D derivative
C                                   integrals for each cartesian
C                                   component x = X,Y,Z
C                    DIFFp       =  is true, if differentiation was
C                                   performed on the p = x,y,z
C                                   coordinates
C                    TEMP1(2)    =  scratch arrays holding intermediate
C                                   1DX kinetic and overlap integral
C                                   times scaling factor products
C                    SCALE       =  the NEXP scaling factors
C
C
C                  Output:
C
C                    BATCH       =  batch of primitive cartesian [0|0]
C                                   kinetic derivative integrals
C                                   corresponding to all current
C                                   exponent pairs
C
C
C  AUTHOR      : Norbert Flocke
C------------------------------------------------------------------------
C
C             ...include files and declare variables.
C
C
         IMPLICIT    NONE

         LOGICAL     DIFFX,DIFFY,DIFFZ

         INTEGER     N
         INTEGER     NEXP

         DOUBLE PRECISION  X,Y,Z

         DOUBLE PRECISION  SCALE (1:NEXP)
         DOUBLE PRECISION  TEMP1 (1:NEXP)
         DOUBLE PRECISION  TEMP2 (1:NEXP)

         DOUBLE PRECISION  BATCH (1:NEXP)

         DOUBLE PRECISION  KIN1DX (1:NEXP)
         DOUBLE PRECISION  KIN1DY (1:NEXP)
         DOUBLE PRECISION  KIN1DZ (1:NEXP)

         DOUBLE PRECISION  OVL1DX (1:NEXP,0:1,0:1)
         DOUBLE PRECISION  OVL1DY (1:NEXP,0:1,0:1)
         DOUBLE PRECISION  OVL1DZ (1:NEXP,0:1,0:1)
C
C
C------------------------------------------------------------------------
C
C
C             ...skip multiplications by 1 comming from original
C                (not derivated) 1D pure s-type overlap integrals.
C
C
         IF (DIFFX) THEN
             DO N = 1,NEXP
                TEMP1 (N) = SCALE (N) * KIN1DX (N)
                TEMP2 (N) = SCALE (N) * OVL1DX (N,0,0)
             END DO
         ELSE
             DO N = 1,NEXP
                TEMP1 (N) = SCALE (N) * KIN1DX (N)
                TEMP2 (N) = SCALE (N)
             END DO
         END IF

         IF (DIFFY .AND. DIFFZ) THEN

             DO N = 1,NEXP
                X = TEMP1 (N) * OVL1DY (N,0,0)
     +                        * OVL1DZ (N,0,0)
                Y = TEMP2 (N) * KIN1DY (N)
     +                        * OVL1DZ (N,0,0)
                Z = TEMP2 (N) * OVL1DY (N,0,0)
     +                        * KIN1DZ (N)
                BATCH (N) = X + Y + Z
             END DO

         ELSE IF (DIFFY) THEN

             DO N = 1,NEXP
                X = TEMP1 (N) * OVL1DY (N,0,0)
                Y = TEMP2 (N) * KIN1DY (N)
                Z = TEMP2 (N) * OVL1DY (N,0,0)
     +                        * KIN1DZ (N)
                BATCH (N) = X + Y + Z
             END DO

         ELSE IF (DIFFZ) THEN

             DO N = 1,NEXP
                X = TEMP1 (N) * OVL1DZ (N,0,0)
                Y = TEMP2 (N) * KIN1DY (N)
     +                        * OVL1DZ (N,0,0)
                Z = TEMP2 (N) * KIN1DZ (N)
                BATCH (N) = X + Y + Z
             END DO

         ELSE

             DO N = 1,NEXP
                X = TEMP1 (N)
                Y = TEMP2 (N) * KIN1DY (N)
                Z = TEMP2 (N) * KIN1DZ (N)
                BATCH (N) = X + Y + Z
             END DO

         END IF
C
C
C             ...ready!
C
C
         RETURN
         END
