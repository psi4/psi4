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
         SUBROUTINE  OED__OVL_DERV_INT1D_TO_A0
     +
     +                    ( SHELLA,
     +                      NEXP,
     +                      NXYZA,
     +                      INT1DX,INT1DY,INT1DZ,
     +                      DIFFX,DIFFY,DIFFZ,
     +                      TEMP1,TEMP2,
     +                      SCALE,
     +
     +                                BATCH )
     +
C------------------------------------------------------------------------
C  OPERATION   : OED__OVL_DERV_INT1D_TO_A0
C  MODULE      : ONE ELECTRON INTEGRALS DIRECT
C  MODULE-ID   : OED
C  SUBROUTINES : none
C  DESCRIPTION : This routine assembles the set of batches of cartesian
C                derivative overlap integrals [A|0] or [0|A].
C
C                Simplified version of the general AB routine to reduce
C                loop overheads for those cases where there is at least
C                one s-shell. For comments and details see the general
C                AB routine.
C
C
C                  Input:
C
C                    SHELLA      =  shell types for contraction shell A
C                    NEXP        =  current # of exponent pairs
C                    NXYZA       =  # of cartesian monomials for shell A
C                    INT1Dx      =  all current 1D A0/0A derivative
C                                   nuclear attraction integrals for
C                                   each cartesian component (x = X,Y,Z)
C                    DIFFx       =  is true, if differentiation was
C                                   performed along the x = X,Y,Z
C                                   direction
C                    TEMP1(2)    =  scratch arrays holding intermediate
C                                   1D A0/0A derivative integral
C                                   products
C                    SCALE       =  the NEXP scaling factors
C
C
C                  Output:
C
C                    BATCH       =  batch of primitive cartesian
C                                   [A|0] or [0|A] derivative overlap
C                                   integrals corresponding to all
C                                   current exponent pairs
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

         INTEGER     I,N
         INTEGER     NEXP
         INTEGER     NXYZA
         INTEGER     SHELLA
         INTEGER     XA,YA,ZA
         INTEGER     YAMAX

         DOUBLE PRECISION  SCALE (1:NEXP)
         DOUBLE PRECISION  TEMP1 (1:NEXP)
         DOUBLE PRECISION  TEMP2 (1:NEXP)

         DOUBLE PRECISION  BATCH (1:NEXP,1:NXYZA)

         DOUBLE PRECISION  INT1DX (1:NEXP,0:SHELLA)
         DOUBLE PRECISION  INT1DY (1:NEXP,0:SHELLA)
         DOUBLE PRECISION  INT1DZ (1:NEXP,0:SHELLA)
C
C
C------------------------------------------------------------------------
C
C
C             ...outer loop over x-contribution. Skip the multiplication
C                of x-contributions, if no x-coordinate derivative was
C                formed and we have a 0-contribution, as then the 1DX
C                integrals are equal to 1.
C
C
         I = 0
         DO 100 XA = SHELLA,0,-1
            YAMAX = SHELLA - XA

            IF (.NOT.DIFFX .AND. XA.EQ.0) THEN
                DO N = 1,NEXP
                   TEMP1 (N) = SCALE (N)
                END DO
            ELSE
                DO N = 1,NEXP
                   TEMP1 (N) = SCALE (N) * INT1DX (N,XA)
                END DO
            END IF
C
C
C             ...inner loop over y-contribution. Skip the multiplication
C                of y-contributions, if no y-coordinate derivative
C                was formed and we have a 0-contribution, as then the
C                1DY integrals are equal to 1.
C
C
            DO 110 YA = YAMAX,0,-1
               I = I + 1
               ZA = YAMAX - YA

               IF (.NOT.DIFFY .AND. YA.EQ.0) THEN
                   DO N = 1,NEXP
                      TEMP2 (N) = TEMP1 (N)
                   END DO
               ELSE
                   DO N = 1,NEXP
                      TEMP2 (N) = TEMP1 (N) * INT1DY (N,YA)
                   END DO
               END IF
C
C
C             ...skip multiplication of z-contributions, if we
C                have a 0-type and no derivations were performed
C                on the z-coordinate, as then the 1DZ integrals
C                are equal to 1. All info concerning all three
C                x,y and z contributions have been collected for all
C                exponent pairs. Put the 1D X,Y,Z integral products
C                in the appropriate places of the [A|0] or [0|A] batch.
C
C
               IF (.NOT.DIFFZ .AND. ZA.EQ.0) THEN
                   DO N = 1,NEXP
                      BATCH (N,I) = TEMP2 (N)
                   END DO
               ELSE
                   DO N = 1,NEXP
                      BATCH (N,I) = TEMP2 (N) * INT1DZ (N,ZA)
                   END DO
               END IF
C
C
C             ...next z- and y-contributions.
C
C
  110       CONTINUE
C
C
C             ...next x-contribution.
C
C
  100    CONTINUE
C
C
C             ...ready!
C
C
         RETURN
         END
