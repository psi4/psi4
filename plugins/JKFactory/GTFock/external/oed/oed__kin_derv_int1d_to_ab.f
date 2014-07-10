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
         SUBROUTINE  OED__KIN_DERV_INT1D_TO_AB
     +
     +                    ( SHELLA,SHELLB,
     +                      NEXP,
     +                      NXYZA,NXYZB,
     +                      KIN1DX,KIN1DY,KIN1DZ,
     +                      OVL1DX,OVL1DY,OVL1DZ,
     +                      DIFFX,DIFFY,DIFFZ,
     +                      TEMP1,TEMP2,
     +                      SCALE,
     +
     +                                BATCH )
     +
C------------------------------------------------------------------------
C  OPERATION   : OED__KIN_DERV_INT1D_TO_AB
C  MODULE      : ONE ELECTRON INTEGRALS DIRECT
C  MODULE-ID   : OED
C  SUBROUTINES : none
C  DESCRIPTION : This routine assembles the set of batches of cartesian
C                kinetic derivative integrals [A|B], adding up all
C                contributions from all the kinetic and overlap 1D
C                derivative integrals:
C
C                      X = KIN (XA,XB) * OVL (YA,YB) * OVL (ZA,ZB)
C                      Y = OVL (XA,XB) * KIN (YA,YB) * OVL (ZA,ZB)
C                      Z = OVL (XA,XB) * OVL (YA,YB) * KIN (ZA,ZB)
C
C                        [A|B] = [XA,YA,ZA|XB,YB,ZB] = X + Y + Z
C
C                The routine uses the reduced Rys multiplication scheme
C                as suggested in R.Lindh, U.Ryu and B.Liu, JCP 95, 5889.
C                This scheme reuses intermediate products between the
C                kinetic and overlap 1DX derivative integrals with the
C                scaling factors.
C
C
C                  Input:
C
C                    SHELLx      =  shell types for individual csh
C                                   x=A and B
C                    NEXP        =  current # of exponent pairs
C                    NXYZx       =  # of cartesian monomials for the
C                                   shells x=A,B
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
C                                   1DX kinetic and overlap derivative
C                                   integral times scaling factor
C                                   products
C                    SCALE       =  the NEXP scaling factors
C
C
C                  Output:
C
C                    BATCH       =  batch of primitive cartesian [A|B]
C                                   kinetic derivative integrals
C                                   corresponding to all current
C                                   exponent pairs
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

         LOGICAL     DIFFX,DIFFY,DIFFZ
         LOGICAL     YNOT1,ZNOT1

         INTEGER     I,J,N
         INTEGER     NEXP
         INTEGER     NXYZA,NXYZB
         INTEGER     SHELLA,SHELLB
         INTEGER     XA,YA,ZA,XB,YB,ZB
         INTEGER     XAP,XBP
         INTEGER     YAMAX,YBMAX

         DOUBLE PRECISION  X,Y,Z
         DOUBLE PRECISION  ZERO

         DOUBLE PRECISION  SCALE (1:NEXP)
         DOUBLE PRECISION  TEMP1 (1:NEXP)
         DOUBLE PRECISION  TEMP2 (1:NEXP)

         DOUBLE PRECISION  BATCH (1:NEXP,1:NXYZA,1:NXYZB)

         DOUBLE PRECISION  KIN1DX (1:NEXP,0:SHELLA  ,0:SHELLB  )
         DOUBLE PRECISION  KIN1DY (1:NEXP,0:SHELLA  ,0:SHELLB  )
         DOUBLE PRECISION  KIN1DZ (1:NEXP,0:SHELLA  ,0:SHELLB  )
         DOUBLE PRECISION  OVL1DX (1:NEXP,0:SHELLA+1,0:SHELLB+1)
         DOUBLE PRECISION  OVL1DY (1:NEXP,0:SHELLA+1,0:SHELLB+1)
         DOUBLE PRECISION  OVL1DZ (1:NEXP,0:SHELLA+1,0:SHELLB+1)

         PARAMETER  (ZERO  =  0.D0)
C
C
C------------------------------------------------------------------------
C
C
C             ...outer loop is over x contributions, inner loop over
C                y,z contributions. Skip multiplications by 1 comming
C                from original (not derivated) 1D pure s-type overlap
C                integrals.
C
C
         XBP = 0
         DO 100 XB = SHELLB,0,-1
            YBMAX = SHELLB - XB
            XAP = 0
            DO 110 XA = SHELLA,0,-1
               YAMAX = SHELLA - XA

               IF (DIFFX .OR. XA+XB.GT.0) THEN
                   DO N = 1,NEXP
                      TEMP1 (N) = SCALE (N) * KIN1DX (N,XA,XB)
                      TEMP2 (N) = SCALE (N) * OVL1DX (N,XA,XB)
                   END DO
               ELSE
                   DO N = 1,NEXP
                      TEMP1 (N) = SCALE (N) * KIN1DX (N,XA,XB)
                      TEMP2 (N) = SCALE (N)
                   END DO
               END IF

               J = XBP
               DO 120 YB = YBMAX,0,-1
                  J = J + 1
                  ZB = YBMAX - YB
                  I = XAP
                  DO 130 YA = YAMAX,0,-1
                     I = I + 1
                     ZA = YAMAX - YA

                     YNOT1 = DIFFY .OR. YA+YB.GT.0
                     ZNOT1 = DIFFZ .OR. ZA+ZB.GT.0

                     IF (YNOT1 .AND. ZNOT1) THEN

                         DO N = 1,NEXP
                            X = TEMP1 (N) * OVL1DY (N,YA,YB)
     +                                    * OVL1DZ (N,ZA,ZB)
                            Y = TEMP2 (N) * KIN1DY (N,YA,YB)
     +                                    * OVL1DZ (N,ZA,ZB)
                            Z = TEMP2 (N) * OVL1DY (N,YA,YB)
     +                                    * KIN1DZ (N,ZA,ZB)
                            BATCH (N,I,J) = X + Y + Z
                         END DO

                     ELSE IF (YNOT1) THEN

                         DO N = 1,NEXP
                            X = TEMP1 (N) * OVL1DY (N,YA,YB)
                            Y = TEMP2 (N) * KIN1DY (N,YA,YB)
                            Z = TEMP2 (N) * OVL1DY (N,YA,YB)
     +                                    * KIN1DZ (N,ZA,ZB)
                            BATCH (N,I,J) = X + Y + Z
                         END DO

                     ELSE IF (ZNOT1) THEN

                         DO N = 1,NEXP
                            X = TEMP1 (N) * OVL1DZ (N,ZA,ZB)
                            Y = TEMP2 (N) * KIN1DY (N,YA,YB)
     +                                    * OVL1DZ (N,ZA,ZB)
                            Z = TEMP2 (N) * KIN1DZ (N,ZA,ZB)
                            BATCH (N,I,J) = X + Y + Z
                         END DO

                     ELSE

                         DO N = 1,NEXP
                            X = TEMP1 (N)
                            Y = TEMP2 (N) * KIN1DY (N,YA,YB)
                            Z = TEMP2 (N) * KIN1DZ (N,ZA,ZB)
                            BATCH (N,I,J) = X + Y + Z
                         END DO

                     END IF

  130             CONTINUE
  120          CONTINUE

               XAP = I
  110       CONTINUE
            XBP = J
  100    CONTINUE
C
C
C             ...ready!
C
C
         RETURN
         END
