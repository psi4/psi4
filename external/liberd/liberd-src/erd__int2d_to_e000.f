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
         SUBROUTINE  ERD__INT2D_TO_E000
     +
     +                    ( SHELLA,SHELLP,
     +                      NGQP,NEXQ,NGQEXQ,
     +                      NXYZET,NXYZP,
     +                      INT2DX,INT2DY,INT2DZ,
     +                      TEMP1,TEMP2,
     +                      SCALE,
     +
     +                                BATCH )
     +
C------------------------------------------------------------------------
C  OPERATION   : ERD__INT2D_TO_E000
C  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT
C  MODULE-ID   : ERD
C  SUBROUTINES : none
C  DESCRIPTION : This routine assembles the set of batches of cartesian
C                eris:
C
C                        [E0|00] or [00|E0] , E = A to P
C
C                adding up all the contributions from all the
C                respective 2D integrals:
C
C                                   P0 or 0P
C
C                Simplified version of the general E0F0 routine to
C                reduce loop overheads for those cases where there is
C                at least one s-shell on the bra or ket side. For
C                comments and details see the general E0F0 routine.
C
C
C                  Input:
C
C                    SHELLx      =  shell types for individual csh
C                                   x=A and csh sum P=A+B
C                    NGQP        =  # of gaussian quadrature points
C                                   (roots)
C                    NEXQ        =  current # of exponent quadruplets
C                    NGQEXQ      =  product of # of gaussian quadrature
C                                   points times exponent quadruplets
C                    NXYZET      =  sum of # of cartesian monomials
C                                   for all shells in the range
C                                   E = A,...,P=A+B
C                    NXYZP       =  # of cartesian monomials for the
C                                   P=A+B shell
C                    INT2Dx      =  all current 2D P0/0P integrals for
C                                   each cartesian component
C                                   (x = X,Y,Z)
C                    TEMP1(2)    =  scratch arrays holding intermediate
C                                   2D P0/0P integral products
C                    SCALE       =  the NGQEXQ scaling factors
C
C
C                  Output:
C
C                    BATCH       =  batch of primitive cartesian
C                                   [E0|00] integrals corresponding
C                                   to all current exponent quadruplets
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

         INTEGER     I,K,M,N
         INTEGER     NGQP,NEXQ,NGQEXQ
         INTEGER     NXYZE,NXYZET,NXYZP
         INTEGER     SE,SEEND
         INTEGER     SHELLA,SHELLP
         INTEGER     XE,YE,ZE
         INTEGER     XEMAX
         INTEGER     XEP,XYEP
         INTEGER     XYE
         INTEGER     YEEND

         DOUBLE PRECISION  SUM
         DOUBLE PRECISION  ZERO

         DOUBLE PRECISION  SCALE (1:NGQEXQ)
         DOUBLE PRECISION  TEMP1 (1:NGQEXQ)
         DOUBLE PRECISION  TEMP2 (1:NGQEXQ)

         DOUBLE PRECISION  BATCH (1:NEXQ,1:NXYZET)

         DOUBLE PRECISION  INT2DX (1:NGQEXQ,0:SHELLP)
         DOUBLE PRECISION  INT2DY (1:NGQEXQ,0:SHELLP)
         DOUBLE PRECISION  INT2DZ (1:NGQEXQ,0:SHELLP)

         PARAMETER  (ZERO  =  0.D0)
C
C
C------------------------------------------------------------------------
C
C
C             ...jump according to number of roots.
C
C
         GOTO  (1,2,3,4,5,6,7,8,9,10)  MIN (NGQP,10)
C
C
C                       ********************
C                       *  # of roots = 1  *
C                       ********************
C
C
    1    XEP = NXYZET + 3
         DO 100 XE = 0,SHELLP
            XEP = XEP + XE - 2
            XEMAX = XE * SHELLP
            YEEND = SHELLP - XE

            DO M = 1,NEXQ
               TEMP1 (M) = SCALE (M) * INT2DX (M,XE)
            END DO

            XYEP = XEP - XEMAX
            DO 110 YE = 0,YEEND
               XYE = XE + YE
               XYEP = XYEP - 1
               SEEND = MAX (SHELLA,XYE)

               IF (YE.EQ.0) THEN
                   DO M = 1,NEXQ
                      TEMP2 (M) = TEMP1 (M)
                   END DO
               ELSE
                   DO M = 1,NEXQ
                      TEMP2 (M) = TEMP1 (M) * INT2DY (M,YE)
                   END DO
               END IF

               I = XYEP
               NXYZE = NXYZP
               DO 120 SE = SHELLP,SEEND,-1
                  ZE = SE - XYE

                  IF (ZE.EQ.0) THEN
                      DO M = 1,NEXQ
                         BATCH (M,I) = TEMP2 (M)
                      END DO
                  ELSE
                      DO M = 1,NEXQ
                         BATCH (M,I) = TEMP2 (M) * INT2DZ (M,ZE)
                      END DO
                  END IF

                  I = I - NXYZE + XE
                  NXYZE = NXYZE - SE - 1
  120          CONTINUE

  110       CONTINUE
  100    CONTINUE

         RETURN
C
C
C                       ********************
C                       *  # of roots = 2  *
C                       ********************
C
C
    2    XEP = NXYZET + 3
         DO 200 XE = 0,SHELLP
            XEP = XEP + XE - 2
            XEMAX = XE * SHELLP
            YEEND = SHELLP - XE

            DO M = 1,NGQEXQ
               TEMP1 (M) = SCALE (M) * INT2DX (M,XE)
            END DO

            XYEP = XEP - XEMAX
            DO 210 YE = 0,YEEND
               XYE = XE + YE
               XYEP = XYEP - 1
               SEEND = MAX (SHELLA,XYE)

               IF (YE.EQ.0) THEN
                   DO N = 1,NGQEXQ
                      TEMP2 (N) = TEMP1 (N)
                   END DO
               ELSE
                   DO N = 1,NGQEXQ
                      TEMP2 (N) = TEMP1 (N) * INT2DY (N,YE)
                   END DO
               END IF

               I = XYEP
               NXYZE = NXYZP
               DO 220 SE = SHELLP,SEEND,-1
                  ZE = SE - XYE

                  IF (ZE.EQ.0) THEN
                      K = 1
                      DO M = 1,NEXQ
                         BATCH (M,I) =   TEMP2 (K)
     +                                 + TEMP2 (K+1)
                         K = K + 2
                      END DO
                  ELSE
                      K = 1
                      DO M = 1,NEXQ
                         BATCH (M,I) =   TEMP2 (K)
     +                                 * INT2DZ (K,ZE)
     +                                 + TEMP2 (K+1)
     +                                 * INT2DZ (K+1,ZE)
                         K = K + 2
                      END DO
                  END IF

                  I = I - NXYZE + XE
                  NXYZE = NXYZE - SE - 1
  220          CONTINUE

  210       CONTINUE
  200    CONTINUE

         RETURN
C
C
C                       ********************
C                       *  # of roots = 3  *
C                       ********************
C
C
    3    XEP = NXYZET + 3
         DO 300 XE = 0,SHELLP
            XEP = XEP + XE - 2
            XEMAX = XE * SHELLP
            YEEND = SHELLP - XE

            DO M = 1,NGQEXQ
               TEMP1 (M) = SCALE (M) * INT2DX (M,XE)
            END DO

            XYEP = XEP - XEMAX
            DO 310 YE = 0,YEEND
               XYE = XE + YE
               XYEP = XYEP - 1
               SEEND = MAX (SHELLA,XYE)

               IF (YE.EQ.0) THEN
                   DO N = 1,NGQEXQ
                      TEMP2 (N) = TEMP1 (N)
                   END DO
               ELSE
                   DO N = 1,NGQEXQ
                      TEMP2 (N) = TEMP1 (N) * INT2DY (N,YE)
                   END DO
               END IF

               I = XYEP
               NXYZE = NXYZP
               DO 320 SE = SHELLP,SEEND,-1
                  ZE = SE - XYE

                  IF (ZE.EQ.0) THEN
                      K = 1
                      DO M = 1,NEXQ
                         BATCH (M,I) =   TEMP2 (K)
     +                                 + TEMP2 (K+1)
     +                                 + TEMP2 (K+2)
                         K = K + 3
                      END DO
                  ELSE
                      K = 1
                      DO M = 1,NEXQ
                         BATCH (M,I) =   TEMP2 (K)
     +                                 * INT2DZ (K,ZE)
     +                                 + TEMP2 (K+1)
     +                                 * INT2DZ (K+1,ZE)
     +                                 + TEMP2 (K+2)
     +                                 * INT2DZ (K+2,ZE)
                         K = K + 3
                      END DO
                  END IF

                  I = I - NXYZE + XE
                  NXYZE = NXYZE - SE - 1
  320          CONTINUE

  310       CONTINUE
  300    CONTINUE

         RETURN
C
C
C                       ********************
C                       *  # of roots = 4  *
C                       ********************
C
C
    4    XEP = NXYZET + 3
         DO 400 XE = 0,SHELLP
            XEP = XEP + XE - 2
            XEMAX = XE * SHELLP
            YEEND = SHELLP - XE

            DO M = 1,NGQEXQ
               TEMP1 (M) = SCALE (M) * INT2DX (M,XE)
            END DO

            XYEP = XEP - XEMAX
            DO 410 YE = 0,YEEND
               XYE = XE + YE
               XYEP = XYEP - 1
               SEEND = MAX (SHELLA,XYE)

               IF (YE.EQ.0) THEN
                   DO N = 1,NGQEXQ
                      TEMP2 (N) = TEMP1 (N)
                   END DO
               ELSE
                   DO N = 1,NGQEXQ
                      TEMP2 (N) = TEMP1 (N) * INT2DY (N,YE)
                   END DO
               END IF

               I = XYEP
               NXYZE = NXYZP
               DO 420 SE = SHELLP,SEEND,-1
                  ZE = SE - XYE

                  IF (ZE.EQ.0) THEN
                      K = 1
                      DO M = 1,NEXQ
                         BATCH (M,I) =   TEMP2 (K)
     +                                 + TEMP2 (K+1)
     +                                 + TEMP2 (K+2)
     +                                 + TEMP2 (K+3)
                         K = K + 4
                      END DO
                  ELSE
                      K = 1
                      DO M = 1,NEXQ
                         BATCH (M,I) =   TEMP2 (K)
     +                                 * INT2DZ (K,ZE)
     +                                 + TEMP2 (K+1)
     +                                 * INT2DZ (K+1,ZE)
     +                                 + TEMP2 (K+2)
     +                                 * INT2DZ (K+2,ZE)
     +                                 + TEMP2 (K+3)
     +                                 * INT2DZ (K+3,ZE)
                         K = K + 4
                      END DO
                  END IF

                  I = I - NXYZE + XE
                  NXYZE = NXYZE - SE - 1
  420          CONTINUE

  410       CONTINUE
  400    CONTINUE

         RETURN
C
C
C                       ********************
C                       *  # of roots = 5  *
C                       ********************
C
C
    5    XEP = NXYZET + 3
         DO 500 XE = 0,SHELLP
            XEP = XEP + XE - 2
            XEMAX = XE * SHELLP
            YEEND = SHELLP - XE

            DO M = 1,NGQEXQ
               TEMP1 (M) = SCALE (M) * INT2DX (M,XE)
            END DO

            XYEP = XEP - XEMAX
            DO 510 YE = 0,YEEND
               XYE = XE + YE
               XYEP = XYEP - 1
               SEEND = MAX (SHELLA,XYE)

               IF (YE.EQ.0) THEN
                   DO N = 1,NGQEXQ
                      TEMP2 (N) = TEMP1 (N)
                   END DO
               ELSE
                   DO N = 1,NGQEXQ
                      TEMP2 (N) = TEMP1 (N) * INT2DY (N,YE)
                   END DO
               END IF

               I = XYEP
               NXYZE = NXYZP
               DO 520 SE = SHELLP,SEEND,-1
                  ZE = SE - XYE

                  IF (ZE.EQ.0) THEN
                      K = 1
                      DO M = 1,NEXQ
                         BATCH (M,I) =   TEMP2 (K)
     +                                 + TEMP2 (K+1)
     +                                 + TEMP2 (K+2)
     +                                 + TEMP2 (K+3)
     +                                 + TEMP2 (K+4)
                         K = K + 5
                      END DO
                  ELSE
                      K = 1
                      DO M = 1,NEXQ
                         BATCH (M,I) =   TEMP2 (K)
     +                                 * INT2DZ (K,ZE)
     +                                 + TEMP2 (K+1)
     +                                 * INT2DZ (K+1,ZE)
     +                                 + TEMP2 (K+2)
     +                                 * INT2DZ (K+2,ZE)
     +                                 + TEMP2 (K+3)
     +                                 * INT2DZ (K+3,ZE)
     +                                 + TEMP2 (K+4)
     +                                 * INT2DZ (K+4,ZE)
                         K = K + 5
                      END DO
                  END IF

                  I = I - NXYZE + XE
                  NXYZE = NXYZE - SE - 1
  520          CONTINUE

  510       CONTINUE
  500    CONTINUE

         RETURN
C
C
C                       ********************
C                       *  # of roots = 6  *
C                       ********************
C
C
    6    XEP = NXYZET + 3
         DO 600 XE = 0,SHELLP
            XEP = XEP + XE - 2
            XEMAX = XE * SHELLP
            YEEND = SHELLP - XE

            DO M = 1,NGQEXQ
               TEMP1 (M) = SCALE (M) * INT2DX (M,XE)
            END DO

            XYEP = XEP - XEMAX
            DO 610 YE = 0,YEEND
               XYE = XE + YE
               XYEP = XYEP - 1
               SEEND = MAX (SHELLA,XYE)

               IF (YE.EQ.0) THEN
                   DO N = 1,NGQEXQ
                      TEMP2 (N) = TEMP1 (N)
                   END DO
               ELSE
                   DO N = 1,NGQEXQ
                      TEMP2 (N) = TEMP1 (N) * INT2DY (N,YE)
                   END DO
               END IF

               I = XYEP
               NXYZE = NXYZP
               DO 620 SE = SHELLP,SEEND,-1
                  ZE = SE - XYE

                  IF (ZE.EQ.0) THEN
                      K = 1
                      DO M = 1,NEXQ
                         BATCH (M,I) =   TEMP2 (K)
     +                                 + TEMP2 (K+1)
     +                                 + TEMP2 (K+2)
     +                                 + TEMP2 (K+3)
     +                                 + TEMP2 (K+4)
     +                                 + TEMP2 (K+5)
                         K = K + 6
                      END DO
                  ELSE
                      K = 1
                      DO M = 1,NEXQ
                         BATCH (M,I) =   TEMP2 (K)
     +                                 * INT2DZ (K,ZE)
     +                                 + TEMP2 (K+1)
     +                                 * INT2DZ (K+1,ZE)
     +                                 + TEMP2 (K+2)
     +                                 * INT2DZ (K+2,ZE)
     +                                 + TEMP2 (K+3)
     +                                 * INT2DZ (K+3,ZE)
     +                                 + TEMP2 (K+4)
     +                                 * INT2DZ (K+4,ZE)
     +                                 + TEMP2 (K+5)
     +                                 * INT2DZ (K+5,ZE)
                         K = K + 6
                      END DO
                  END IF

                  I = I - NXYZE + XE
                  NXYZE = NXYZE - SE - 1
  620          CONTINUE

  610       CONTINUE
  600    CONTINUE

         RETURN
C
C
C                       ********************
C                       *  # of roots = 7  *
C                       ********************
C
C
    7    XEP = NXYZET + 3
         DO 700 XE = 0,SHELLP
            XEP = XEP + XE - 2
            XEMAX = XE * SHELLP
            YEEND = SHELLP - XE

            DO M = 1,NGQEXQ
               TEMP1 (M) = SCALE (M) * INT2DX (M,XE)
            END DO

            XYEP = XEP - XEMAX
            DO 710 YE = 0,YEEND
               XYE = XE + YE
               XYEP = XYEP - 1
               SEEND = MAX (SHELLA,XYE)

               IF (YE.EQ.0) THEN
                   DO N = 1,NGQEXQ
                      TEMP2 (N) = TEMP1 (N)
                   END DO
               ELSE
                   DO N = 1,NGQEXQ
                      TEMP2 (N) = TEMP1 (N) * INT2DY (N,YE)
                   END DO
               END IF

               I = XYEP
               NXYZE = NXYZP
               DO 720 SE = SHELLP,SEEND,-1
                  ZE = SE - XYE

                  IF (ZE.EQ.0) THEN
                      K = 1
                      DO M = 1,NEXQ
                         BATCH (M,I) =   TEMP2 (K)
     +                                 + TEMP2 (K+1)
     +                                 + TEMP2 (K+2)
     +                                 + TEMP2 (K+3)
     +                                 + TEMP2 (K+4)
     +                                 + TEMP2 (K+5)
     +                                 + TEMP2 (K+6)
                         K = K + 7
                      END DO
                  ELSE
                      K = 1
                      DO M = 1,NEXQ
                         BATCH (M,I) =   TEMP2 (K)
     +                                 * INT2DZ (K,ZE)
     +                                 + TEMP2 (K+1)
     +                                 * INT2DZ (K+1,ZE)
     +                                 + TEMP2 (K+2)
     +                                 * INT2DZ (K+2,ZE)
     +                                 + TEMP2 (K+3)
     +                                 * INT2DZ (K+3,ZE)
     +                                 + TEMP2 (K+4)
     +                                 * INT2DZ (K+4,ZE)
     +                                 + TEMP2 (K+5)
     +                                 * INT2DZ (K+5,ZE)
     +                                 + TEMP2 (K+6)
     +                                 * INT2DZ (K+6,ZE)
                         K = K + 7
                      END DO
                  END IF

                  I = I - NXYZE + XE
                  NXYZE = NXYZE - SE - 1
  720          CONTINUE

  710       CONTINUE
  700    CONTINUE

         RETURN
C
C
C                       ********************
C                       *  # of roots = 8  *
C                       ********************
C
C
    8    XEP = NXYZET + 3
         DO 800 XE = 0,SHELLP
            XEP = XEP + XE - 2
            XEMAX = XE * SHELLP
            YEEND = SHELLP - XE

            DO M = 1,NGQEXQ
               TEMP1 (M) = SCALE (M) * INT2DX (M,XE)
            END DO

            XYEP = XEP - XEMAX
            DO 810 YE = 0,YEEND
               XYE = XE + YE
               XYEP = XYEP - 1
               SEEND = MAX (SHELLA,XYE)

               IF (YE.EQ.0) THEN
                   DO N = 1,NGQEXQ
                      TEMP2 (N) = TEMP1 (N)
                   END DO
               ELSE
                   DO N = 1,NGQEXQ
                      TEMP2 (N) = TEMP1 (N) * INT2DY (N,YE)
                   END DO
               END IF

               I = XYEP
               NXYZE = NXYZP
               DO 820 SE = SHELLP,SEEND,-1
                  ZE = SE - XYE

                  IF (ZE.EQ.0) THEN
                      K = 1
                      DO M = 1,NEXQ
                         BATCH (M,I) =   TEMP2 (K)
     +                                 + TEMP2 (K+1)
     +                                 + TEMP2 (K+2)
     +                                 + TEMP2 (K+3)
     +                                 + TEMP2 (K+4)
     +                                 + TEMP2 (K+5)
     +                                 + TEMP2 (K+6)
     +                                 + TEMP2 (K+7)
                         K = K + 8
                      END DO
                  ELSE
                      K = 1
                      DO M = 1,NEXQ
                         BATCH (M,I) =   TEMP2 (K)
     +                                 * INT2DZ (K,ZE)
     +                                 + TEMP2 (K+1)
     +                                 * INT2DZ (K+1,ZE)
     +                                 + TEMP2 (K+2)
     +                                 * INT2DZ (K+2,ZE)
     +                                 + TEMP2 (K+3)
     +                                 * INT2DZ (K+3,ZE)
     +                                 + TEMP2 (K+4)
     +                                 * INT2DZ (K+4,ZE)
     +                                 + TEMP2 (K+5)
     +                                 * INT2DZ (K+5,ZE)
     +                                 + TEMP2 (K+6)
     +                                 * INT2DZ (K+6,ZE)
     +                                 + TEMP2 (K+7)
     +                                 * INT2DZ (K+7,ZE)
                         K = K + 8
                      END DO
                  END IF

                  I = I - NXYZE + XE
                  NXYZE = NXYZE - SE - 1
  820          CONTINUE

  810       CONTINUE
  800    CONTINUE

         RETURN
C
C
C                       ********************
C                       *  # of roots = 9  *
C                       ********************
C
C
    9    XEP = NXYZET + 3
         DO 900 XE = 0,SHELLP
            XEP = XEP + XE - 2
            XEMAX = XE * SHELLP
            YEEND = SHELLP - XE

            DO M = 1,NGQEXQ
               TEMP1 (M) = SCALE (M) * INT2DX (M,XE)
            END DO

            XYEP = XEP - XEMAX
            DO 910 YE = 0,YEEND
               XYE = XE + YE
               XYEP = XYEP - 1
               SEEND = MAX (SHELLA,XYE)

               IF (YE.EQ.0) THEN
                   DO N = 1,NGQEXQ
                      TEMP2 (N) = TEMP1 (N)
                   END DO
               ELSE
                   DO N = 1,NGQEXQ
                      TEMP2 (N) = TEMP1 (N) * INT2DY (N,YE)
                   END DO
               END IF

               I = XYEP
               NXYZE = NXYZP
               DO 920 SE = SHELLP,SEEND,-1
                  ZE = SE - XYE

                  IF (ZE.EQ.0) THEN
                      K = 1
                      DO M = 1,NEXQ
                         BATCH (M,I) =   TEMP2 (K)
     +                                 + TEMP2 (K+1)
     +                                 + TEMP2 (K+2)
     +                                 + TEMP2 (K+3)
     +                                 + TEMP2 (K+4)
     +                                 + TEMP2 (K+5)
     +                                 + TEMP2 (K+6)
     +                                 + TEMP2 (K+7)
     +                                 + TEMP2 (K+8)
                         K = K + 9
                      END DO
                  ELSE
                      K = 1
                      DO M = 1,NEXQ
                         BATCH (M,I) =   TEMP2 (K)
     +                                 * INT2DZ (K,ZE)
     +                                 + TEMP2 (K+1)
     +                                 * INT2DZ (K+1,ZE)
     +                                 + TEMP2 (K+2)
     +                                 * INT2DZ (K+2,ZE)
     +                                 + TEMP2 (K+3)
     +                                 * INT2DZ (K+3,ZE)
     +                                 + TEMP2 (K+4)
     +                                 * INT2DZ (K+4,ZE)
     +                                 + TEMP2 (K+5)
     +                                 * INT2DZ (K+5,ZE)
     +                                 + TEMP2 (K+6)
     +                                 * INT2DZ (K+6,ZE)
     +                                 + TEMP2 (K+7)
     +                                 * INT2DZ (K+7,ZE)
     +                                 + TEMP2 (K+8)
     +                                 * INT2DZ (K+8,ZE)
                         K = K + 9
                      END DO
                  END IF

                  I = I - NXYZE + XE
                  NXYZE = NXYZE - SE - 1
  920          CONTINUE

  910       CONTINUE
  900    CONTINUE

         RETURN
C
C
C                       ********************
C                       *  # of roots > 9  *
C                       ********************
C
C             ...outer loops over x-contributions. No skipping of the
C                x-contributions of 0-type can be done here, since
C                the 2DX integrals carry the Rys weight!
C
C
   10    XEP = NXYZET + 3
         DO 1000 XE = 0,SHELLP
            XEP = XEP + XE - 2
            XEMAX = XE * SHELLP
            YEEND = SHELLP - XE

            DO M = 1,NGQEXQ
               TEMP1 (M) = SCALE (M) * INT2DX (M,XE)
            END DO
C
C
C             ...middle loops over y-contributions. Skip multiplication
C                of y-contributions, if we have a 0-type, as then the
C                2DY integral is equal to 1.
C
C
            XYEP = XEP - XEMAX
            DO 1010 YE = 0,YEEND
               XYE = XE + YE
               XYEP = XYEP - 1
               SEEND = MAX (SHELLA,XYE)

               IF (YE.EQ.0) THEN
                   DO N = 1,NGQEXQ
                      TEMP2 (N) = TEMP1 (N)
                   END DO
               ELSE
                   DO N = 1,NGQEXQ
                      TEMP2 (N) = TEMP1 (N) * INT2DY (N,YE)
                   END DO
               END IF
C
C
C             ...inner loops over E-pairs. Skip multiplication
C                of z-contributions, if we have a 0-type, as
C                then the 2DZ integral is equal to 1.
C
C
               I = XYEP
               NXYZE = NXYZP
               DO 1020 SE = SHELLP,SEEND,-1
                  ZE = SE - XYE
C
C
C             ...all info concerning all x-, y- and z-contributions
C                have been collected for all exponent quadruplets at
C                once. Sum up the 2D X,Y,Z integral products to the
C                appropriate place of the batch.
C
C
                  IF (ZE.EQ.0) THEN
                      K = 0
                      DO M = 1,NEXQ
                         SUM = ZERO
                         DO N = 1,NGQP
                            SUM = SUM + TEMP2 (K+N)
                         END DO
                         K = K + NGQP
                         BATCH (M,I) = SUM
                      END DO
                  ELSE
                      K = 0
                      DO M = 1,NEXQ
                         SUM = ZERO
                         DO N = 1,NGQP
                            SUM = SUM + TEMP2 (K+N) * INT2DZ (K+N,ZE)
                         END DO
                         K = K + NGQP
                         BATCH (M,I) = SUM
                      END DO
                  END IF
C
C
C             ...next z-contribution.
C
C
                  I = I - NXYZE + XE
                  NXYZE = NXYZE - SE - 1
 1020          CONTINUE
C
C
C             ...next y- and x-contribution.
C
C
 1010       CONTINUE
 1000    CONTINUE
C
C
C             ...ready!
C
C
         RETURN
         END
