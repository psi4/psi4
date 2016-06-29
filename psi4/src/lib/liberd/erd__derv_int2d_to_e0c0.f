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
         SUBROUTINE  ERD__DERV_INT2D_TO_E0C0
     +
     +                    ( SHELLA,SHELLP,SHELLC,
     +                      NGQP,NEXQ,NGQEXQ,
     +                      NXYZET,NXYZP,NXYZC,
     +                      INT2DX,INT2DY,INT2DZ,
     +                      DIFFY,DIFFZ,
     +                      TEMP1,TEMP2,
     +                      SCALE,
     +
     +                                BATCH )
     +
C------------------------------------------------------------------------
C  OPERATION   : ERD__DERV_INT2D_TO_E0C0
C  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT
C  MODULE-ID   : ERD
C  SUBROUTINES : none
C  DESCRIPTION : This routine assembles the set of batches of cartesian
C                derivative eris:
C
C                         [E0|C0] or [E0|0C] , E = A to P
C
C                adding up the contributions from all the respective
C                derivative 2D integrals:
C
C                                   PC0 or P0C
C
C                Simplified version of the general E0CD routine to
C                reduce loop overheads for those cases where there is
C                at least one s-shell in the CD-part. For comments
C                and details see the general E0CD routine.
C
C
C                  Input:
C
C                    SHELLx      =  shell types for individual csh
C                                   x=A,C and csh sum P=A+B
C                    NGQP        =  # of gaussian quadrature points
C                                   (roots)
C                    NEXQ        =  current # of exponent quadruplets
C                    NGQEXQ      =  product of # of gaussian quadrature
C                                   points times exponent quadruplets
C                    NXYZET      =  sum of # of cartesian monomials
C                                   for all shells in the range
C                                   E = A,...,P=A+B
C                    NXYZy       =  # of cartesian monomials for
C                                   y = P,C shells
C                    INT2Dx      =  all current 2D PC0/P0C derivative
C                                   integrals for each cartesian
C                                   component (x = X,Y,Z)
C                    DIFFx       =  is true, if differentiation was
C                                   performed along the x=Y,Z direction
C                    TEMP1(2)    =  scratch arrays holding intermediate
C                                   2D PC0/P0C derivative integral
C                                   products
C                    SCALE       =  the NGQEXQ scaling factors
C
C
C                  Output:
C
C                    BATCH       =  batch of primitive cartesian [E0|C0]
C                                   or [E0|0C]derivative integrals
C                                   corresponding to all current
C                                   exponent quadruplets
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

         LOGICAL     DIFFY,DIFFZ

         INTEGER     I,J,M,N,R
         INTEGER     NGQP,NEXQ,NGQEXQ
         INTEGER     NXYZE,NXYZET,NXYZP,NXYZC
         INTEGER     SE
         INTEGER     SEEND
         INTEGER     SHELLA,SHELLP,SHELLC
         INTEGER     XC,YC,ZC
         INTEGER     XCP
         INTEGER     XE,YE,ZE
         INTEGER     XEMAX
         INTEGER     XEP,XYEP
         INTEGER     XYE
         INTEGER     YCMAX
         INTEGER     YEEND

         DOUBLE PRECISION  SUM
         DOUBLE PRECISION  ZERO

         DOUBLE PRECISION  SCALE (1:NGQEXQ)
         DOUBLE PRECISION  TEMP1 (1:NGQEXQ)
         DOUBLE PRECISION  TEMP2 (1:NGQEXQ)

         DOUBLE PRECISION  BATCH (1:NEXQ,1:NXYZET,1:NXYZC)

         DOUBLE PRECISION  INT2DX (1:NGQEXQ,0:SHELLP,0:SHELLC)
         DOUBLE PRECISION  INT2DY (1:NGQEXQ,0:SHELLP,0:SHELLC)
         DOUBLE PRECISION  INT2DZ (1:NGQEXQ,0:SHELLP,0:SHELLC)

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
    1    XCP = 0
         DO 100 XC = SHELLC,0,-1
            YCMAX = SHELLC - XC

            XEP = NXYZET + 3
            DO 102 XE = 0,SHELLP
               XEP = XEP + XE - 2
               XEMAX = XE * SHELLP
               YEEND = SHELLP - XE

               DO M = 1,NEXQ
                  TEMP1 (M) = SCALE (M) * INT2DX (M,XE,XC)
               END DO

               J = XCP
               DO 110 YC = YCMAX,0,-1
                  J = J + 1
                  ZC = YCMAX - YC

                  XYEP = XEP - XEMAX
                  DO 112 YE = 0,YEEND
                     XYE = XE + YE
                     XYEP = XYEP - 1
                     SEEND = MAX0 (SHELLA,XYE)

                     IF (.NOT.DIFFY .AND. YE+YC.EQ.0) THEN
                         DO M = 1,NEXQ
                            TEMP2 (M) = TEMP1 (M)
                         END DO
                     ELSE
                         DO M = 1,NEXQ
                            TEMP2 (M) = TEMP1 (M) * INT2DY (M,YE,YC)
                         END DO
                     END IF

                     I = XYEP
                     NXYZE = NXYZP
                     DO 120 SE = SHELLP,SEEND,-1
                        ZE = SE - XYE

                        IF (.NOT.DIFFZ .AND. ZE+ZC.EQ.0) THEN
                            DO M = 1,NEXQ
                               BATCH (M,I,J) = TEMP2 (M)
                            END DO
                        ELSE
                            DO M = 1,NEXQ
                               BATCH (M,I,J) =   TEMP2 (M)
     +                                         * INT2DZ (M,ZE,ZC)
                            END DO
                        END IF

                        I = I - NXYZE + XE
                        NXYZE = NXYZE - SE - 1
  120                CONTINUE
  112             CONTINUE
  110          CONTINUE
  102       CONTINUE
            XCP = J
  100    CONTINUE

         RETURN
C
C
C                       ********************
C                       *  # of roots = 2  *
C                       ********************
C
C
    2    XCP = 0
         DO 200 XC = SHELLC,0,-1
            YCMAX = SHELLC - XC

            XEP = NXYZET + 3
            DO 202 XE = 0,SHELLP
               XEP = XEP + XE - 2
               XEMAX = XE * SHELLP
               YEEND = SHELLP - XE

               DO M = 1,NGQEXQ
                  TEMP1 (M) = SCALE (M) * INT2DX (M,XE,XC)
               END DO

               J = XCP
               DO 210 YC = YCMAX,0,-1
                  J = J + 1
                  ZC = YCMAX - YC

                  XYEP = XEP - XEMAX
                  DO 212 YE = 0,YEEND
                     XYE = XE + YE
                     XYEP = XYEP - 1
                     SEEND = MAX0 (SHELLA,XYE)

                     IF (.NOT.DIFFY .AND. YE+YC.EQ.0) THEN
                         DO N = 1,NGQEXQ
                            TEMP2 (N) = TEMP1 (N)
                         END DO
                     ELSE
                         DO N = 1,NGQEXQ
                            TEMP2 (N) = TEMP1 (N) * INT2DY (N,YE,YC)
                         END DO
                     END IF

                     I = XYEP
                     NXYZE = NXYZP
                     DO 220 SE = SHELLP,SEEND,-1
                        ZE = SE - XYE

                        IF (.NOT.DIFFZ .AND. ZE+ZC.EQ.0) THEN
                            R = 1
                            DO M = 1,NEXQ
                               BATCH (M,I,J) =   TEMP2 (R)
     +                                         + TEMP2 (R+1)
                               R = R + 2
                            END DO
                        ELSE
                            R = 1
                            DO M = 1,NEXQ
                               BATCH (M,I,J) =   TEMP2 (R)
     +                                         * INT2DZ (R,ZE,ZC)
     +                                         + TEMP2 (R+1)
     +                                         * INT2DZ (R+1,ZE,ZC)
                               R = R + 2
                            END DO
                        END IF

                        I = I - NXYZE + XE
                        NXYZE = NXYZE - SE - 1
  220                CONTINUE
  212             CONTINUE
  210          CONTINUE
  202       CONTINUE
            XCP = J
  200    CONTINUE

         RETURN
C
C
C                       ********************
C                       *  # of roots = 3  *
C                       ********************
C
C
    3    XCP = 0
         DO 300 XC = SHELLC,0,-1
            YCMAX = SHELLC - XC

            XEP = NXYZET + 3
            DO 302 XE = 0,SHELLP
               XEP = XEP + XE - 2
               XEMAX = XE * SHELLP
               YEEND = SHELLP - XE

               DO M = 1,NGQEXQ
                  TEMP1 (M) = SCALE (M) * INT2DX (M,XE,XC)
               END DO

               J = XCP
               DO 310 YC = YCMAX,0,-1
                  J = J + 1
                  ZC = YCMAX - YC

                  XYEP = XEP - XEMAX
                  DO 312 YE = 0,YEEND
                     XYE = XE + YE
                     XYEP = XYEP - 1
                     SEEND = MAX0 (SHELLA,XYE)

                     IF (.NOT.DIFFY .AND. YE+YC.EQ.0) THEN
                         DO N = 1,NGQEXQ
                            TEMP2 (N) = TEMP1 (N)
                         END DO
                     ELSE
                         DO N = 1,NGQEXQ
                            TEMP2 (N) = TEMP1 (N) * INT2DY (N,YE,YC)
                         END DO
                     END IF

                     I = XYEP
                     NXYZE = NXYZP
                     DO 320 SE = SHELLP,SEEND,-1
                        ZE = SE - XYE

                        IF (.NOT.DIFFZ .AND. ZE+ZC.EQ.0) THEN
                            R = 1
                            DO M = 1,NEXQ
                               BATCH (M,I,J) =   TEMP2 (R)
     +                                         + TEMP2 (R+1)
     +                                         + TEMP2 (R+2)
                               R = R + 3
                            END DO
                        ELSE
                            R = 1
                            DO M = 1,NEXQ
                               BATCH (M,I,J) =   TEMP2 (R)
     +                                         * INT2DZ (R,ZE,ZC)
     +                                         + TEMP2 (R+1)
     +                                         * INT2DZ (R+1,ZE,ZC)
     +                                         + TEMP2 (R+2)
     +                                         * INT2DZ (R+2,ZE,ZC)
                               R = R + 3
                            END DO
                        END IF

                        I = I - NXYZE + XE
                        NXYZE = NXYZE - SE - 1
  320                CONTINUE
  312             CONTINUE
  310          CONTINUE
  302       CONTINUE
            XCP = J
  300    CONTINUE

         RETURN
C
C
C                       ********************
C                       *  # of roots = 4  *
C                       ********************
C
C
    4    XCP = 0
         DO 400 XC = SHELLC,0,-1
            YCMAX = SHELLC - XC

            XEP = NXYZET + 3
            DO 402 XE = 0,SHELLP
               XEP = XEP + XE - 2
               XEMAX = XE * SHELLP
               YEEND = SHELLP - XE

               DO M = 1,NGQEXQ
                  TEMP1 (M) = SCALE (M) * INT2DX (M,XE,XC)
               END DO

               J = XCP
               DO 410 YC = YCMAX,0,-1
                  J = J + 1
                  ZC = YCMAX - YC

                  XYEP = XEP - XEMAX
                  DO 412 YE = 0,YEEND
                     XYE = XE + YE
                     XYEP = XYEP - 1
                     SEEND = MAX0 (SHELLA,XYE)

                     IF (.NOT.DIFFY .AND. YE+YC.EQ.0) THEN
                         DO N = 1,NGQEXQ
                            TEMP2 (N) = TEMP1 (N)
                         END DO
                     ELSE
                         DO N = 1,NGQEXQ
                            TEMP2 (N) = TEMP1 (N) * INT2DY (N,YE,YC)
                         END DO
                     END IF

                     I = XYEP
                     NXYZE = NXYZP
                     DO 420 SE = SHELLP,SEEND,-1
                        ZE = SE - XYE

                        IF (.NOT.DIFFZ .AND. ZE+ZC.EQ.0) THEN
                            R = 1
                            DO M = 1,NEXQ
                               BATCH (M,I,J) =   TEMP2 (R)
     +                                         + TEMP2 (R+1)
     +                                         + TEMP2 (R+2)
     +                                         + TEMP2 (R+3)
                               R = R + 4
                            END DO
                        ELSE
                            R = 1
                            DO M = 1,NEXQ
                               BATCH (M,I,J) =   TEMP2 (R)
     +                                         * INT2DZ (R,ZE,ZC)
     +                                         + TEMP2 (R+1)
     +                                         * INT2DZ (R+1,ZE,ZC)
     +                                         + TEMP2 (R+2)
     +                                         * INT2DZ (R+2,ZE,ZC)
     +                                         + TEMP2 (R+3)
     +                                         * INT2DZ (R+3,ZE,ZC)
                               R = R + 4
                            END DO
                        END IF

                        I = I - NXYZE + XE
                        NXYZE = NXYZE - SE - 1
  420                CONTINUE
  412             CONTINUE
  410          CONTINUE
  402       CONTINUE
            XCP = J
  400    CONTINUE

         RETURN
C
C
C                       ********************
C                       *  # of roots = 5  *
C                       ********************
C
C
    5    XCP = 0
         DO 500 XC = SHELLC,0,-1
            YCMAX = SHELLC - XC

            XEP = NXYZET + 3
            DO 502 XE = 0,SHELLP
               XEP = XEP + XE - 2
               XEMAX = XE * SHELLP
               YEEND = SHELLP - XE

               DO M = 1,NGQEXQ
                  TEMP1 (M) = SCALE (M) * INT2DX (M,XE,XC)
               END DO

               J = XCP
               DO 510 YC = YCMAX,0,-1
                  J = J + 1
                  ZC = YCMAX - YC

                  XYEP = XEP - XEMAX
                  DO 512 YE = 0,YEEND
                     XYE = XE + YE
                     XYEP = XYEP - 1
                     SEEND = MAX0 (SHELLA,XYE)

                     IF (.NOT.DIFFY .AND. YE+YC.EQ.0) THEN
                         DO N = 1,NGQEXQ
                            TEMP2 (N) = TEMP1 (N)
                         END DO
                     ELSE
                         DO N = 1,NGQEXQ
                            TEMP2 (N) = TEMP1 (N) * INT2DY (N,YE,YC)
                         END DO
                     END IF

                     I = XYEP
                     NXYZE = NXYZP
                     DO 520 SE = SHELLP,SEEND,-1
                        ZE = SE - XYE

                        IF (.NOT.DIFFZ .AND. ZE+ZC.EQ.0) THEN
                            R = 1
                            DO M = 1,NEXQ
                               BATCH (M,I,J) =   TEMP2 (R)
     +                                         + TEMP2 (R+1)
     +                                         + TEMP2 (R+2)
     +                                         + TEMP2 (R+3)
     +                                         + TEMP2 (R+4)
                               R = R + 5
                            END DO
                        ELSE
                            R = 1
                            DO M = 1,NEXQ
                               BATCH (M,I,J) =   TEMP2 (R)
     +                                         * INT2DZ (R,ZE,ZC)
     +                                         + TEMP2 (R+1)
     +                                         * INT2DZ (R+1,ZE,ZC)
     +                                         + TEMP2 (R+2)
     +                                         * INT2DZ (R+2,ZE,ZC)
     +                                         + TEMP2 (R+3)
     +                                         * INT2DZ (R+3,ZE,ZC)
     +                                         + TEMP2 (R+4)
     +                                         * INT2DZ (R+4,ZE,ZC)
                               R = R + 5
                            END DO
                        END IF

                        I = I - NXYZE + XE
                        NXYZE = NXYZE - SE - 1
  520                CONTINUE
  512             CONTINUE
  510          CONTINUE
  502       CONTINUE
            XCP = J
  500    CONTINUE

         RETURN
C
C
C                       ********************
C                       *  # of roots = 6  *
C                       ********************
C
C
    6    XCP = 0
         DO 600 XC = SHELLC,0,-1
            YCMAX = SHELLC - XC

            XEP = NXYZET + 3
            DO 602 XE = 0,SHELLP
               XEP = XEP + XE - 2
               XEMAX = XE * SHELLP
               YEEND = SHELLP - XE

               DO M = 1,NGQEXQ
                  TEMP1 (M) = SCALE (M) * INT2DX (M,XE,XC)
               END DO

               J = XCP
               DO 610 YC = YCMAX,0,-1
                  J = J + 1
                  ZC = YCMAX - YC

                  XYEP = XEP - XEMAX
                  DO 612 YE = 0,YEEND
                     XYE = XE + YE
                     XYEP = XYEP - 1
                     SEEND = MAX0 (SHELLA,XYE)

                     IF (.NOT.DIFFY .AND. YE+YC.EQ.0) THEN
                         DO N = 1,NGQEXQ
                            TEMP2 (N) = TEMP1 (N)
                         END DO
                     ELSE
                         DO N = 1,NGQEXQ
                            TEMP2 (N) = TEMP1 (N) * INT2DY (N,YE,YC)
                         END DO
                     END IF

                     I = XYEP
                     NXYZE = NXYZP
                     DO 620 SE = SHELLP,SEEND,-1
                        ZE = SE - XYE

                        IF (.NOT.DIFFZ .AND. ZE+ZC.EQ.0) THEN
                            R = 1
                            DO M = 1,NEXQ
                               BATCH (M,I,J) =   TEMP2 (R)
     +                                         + TEMP2 (R+1)
     +                                         + TEMP2 (R+2)
     +                                         + TEMP2 (R+3)
     +                                         + TEMP2 (R+4)
     +                                         + TEMP2 (R+5)
                               R = R + 6
                            END DO
                        ELSE
                            R = 1
                            DO M = 1,NEXQ
                               BATCH (M,I,J) =   TEMP2 (R)
     +                                         * INT2DZ (R,ZE,ZC)
     +                                         + TEMP2 (R+1)
     +                                         * INT2DZ (R+1,ZE,ZC)
     +                                         + TEMP2 (R+2)
     +                                         * INT2DZ (R+2,ZE,ZC)
     +                                         + TEMP2 (R+3)
     +                                         * INT2DZ (R+3,ZE,ZC)
     +                                         + TEMP2 (R+4)
     +                                         * INT2DZ (R+4,ZE,ZC)
     +                                         + TEMP2 (R+5)
     +                                         * INT2DZ (R+5,ZE,ZC)
                               R = R + 6
                            END DO
                        END IF

                        I = I - NXYZE + XE
                        NXYZE = NXYZE - SE - 1
  620                CONTINUE
  612             CONTINUE
  610          CONTINUE
  602       CONTINUE
            XCP = J
  600    CONTINUE

         RETURN
C
C
C                       ********************
C                       *  # of roots = 7  *
C                       ********************
C
C
    7    XCP = 0
         DO 700 XC = SHELLC,0,-1
            YCMAX = SHELLC - XC

            XEP = NXYZET + 3
            DO 702 XE = 0,SHELLP
               XEP = XEP + XE - 2
               XEMAX = XE * SHELLP
               YEEND = SHELLP - XE

               DO M = 1,NGQEXQ
                  TEMP1 (M) = SCALE (M) * INT2DX (M,XE,XC)
               END DO

               J = XCP
               DO 710 YC = YCMAX,0,-1
                  J = J + 1
                  ZC = YCMAX - YC

                  XYEP = XEP - XEMAX
                  DO 712 YE = 0,YEEND
                     XYE = XE + YE
                     XYEP = XYEP - 1
                     SEEND = MAX0 (SHELLA,XYE)

                     IF (.NOT.DIFFY .AND. YE+YC.EQ.0) THEN
                         DO N = 1,NGQEXQ
                            TEMP2 (N) = TEMP1 (N)
                         END DO
                     ELSE
                         DO N = 1,NGQEXQ
                            TEMP2 (N) = TEMP1 (N) * INT2DY (N,YE,YC)
                         END DO
                     END IF

                     I = XYEP
                     NXYZE = NXYZP
                     DO 720 SE = SHELLP,SEEND,-1
                        ZE = SE - XYE

                        IF (.NOT.DIFFZ .AND. ZE+ZC.EQ.0) THEN
                            R = 1
                            DO M = 1,NEXQ
                               BATCH (M,I,J) =   TEMP2 (R)
     +                                         + TEMP2 (R+1)
     +                                         + TEMP2 (R+2)
     +                                         + TEMP2 (R+3)
     +                                         + TEMP2 (R+4)
     +                                         + TEMP2 (R+5)
     +                                         + TEMP2 (R+6)
                               R = R + 7
                            END DO
                        ELSE
                            R = 1
                            DO M = 1,NEXQ
                               BATCH (M,I,J) =   TEMP2 (R)
     +                                         * INT2DZ (R,ZE,ZC)
     +                                         + TEMP2 (R+1)
     +                                         * INT2DZ (R+1,ZE,ZC)
     +                                         + TEMP2 (R+2)
     +                                         * INT2DZ (R+2,ZE,ZC)
     +                                         + TEMP2 (R+3)
     +                                         * INT2DZ (R+3,ZE,ZC)
     +                                         + TEMP2 (R+4)
     +                                         * INT2DZ (R+4,ZE,ZC)
     +                                         + TEMP2 (R+5)
     +                                         * INT2DZ (R+5,ZE,ZC)
     +                                         + TEMP2 (R+6)
     +                                         * INT2DZ (R+6,ZE,ZC)
                               R = R + 7
                            END DO
                        END IF

                        I = I - NXYZE + XE
                        NXYZE = NXYZE - SE - 1
  720                CONTINUE
  712             CONTINUE
  710          CONTINUE
  702       CONTINUE
            XCP = J
  700    CONTINUE

         RETURN
C
C
C                       ********************
C                       *  # of roots = 8  *
C                       ********************
C
C
    8    XCP = 0
         DO 800 XC = SHELLC,0,-1
            YCMAX = SHELLC - XC

            XEP = NXYZET + 3
            DO 802 XE = 0,SHELLP
               XEP = XEP + XE - 2
               XEMAX = XE * SHELLP
               YEEND = SHELLP - XE

               DO M = 1,NGQEXQ
                  TEMP1 (M) = SCALE (M) * INT2DX (M,XE,XC)
               END DO

               J = XCP
               DO 810 YC = YCMAX,0,-1
                  J = J + 1
                  ZC = YCMAX - YC

                  XYEP = XEP - XEMAX
                  DO 812 YE = 0,YEEND
                     XYE = XE + YE
                     XYEP = XYEP - 1
                     SEEND = MAX0 (SHELLA,XYE)

                     IF (.NOT.DIFFY .AND. YE+YC.EQ.0) THEN
                         DO N = 1,NGQEXQ
                            TEMP2 (N) = TEMP1 (N)
                         END DO
                     ELSE
                         DO N = 1,NGQEXQ
                            TEMP2 (N) = TEMP1 (N) * INT2DY (N,YE,YC)
                         END DO
                     END IF

                     I = XYEP
                     NXYZE = NXYZP
                     DO 820 SE = SHELLP,SEEND,-1
                        ZE = SE - XYE

                        IF (.NOT.DIFFZ .AND. ZE+ZC.EQ.0) THEN
                            R = 1
                            DO M = 1,NEXQ
                               BATCH (M,I,J) =   TEMP2 (R)
     +                                         + TEMP2 (R+1)
     +                                         + TEMP2 (R+2)
     +                                         + TEMP2 (R+3)
     +                                         + TEMP2 (R+4)
     +                                         + TEMP2 (R+5)
     +                                         + TEMP2 (R+6)
     +                                         + TEMP2 (R+7)
                               R = R + 8
                            END DO
                        ELSE
                            R = 1
                            DO M = 1,NEXQ
                               BATCH (M,I,J) =   TEMP2 (R)
     +                                         * INT2DZ (R,ZE,ZC)
     +                                         + TEMP2 (R+1)
     +                                         * INT2DZ (R+1,ZE,ZC)
     +                                         + TEMP2 (R+2)
     +                                         * INT2DZ (R+2,ZE,ZC)
     +                                         + TEMP2 (R+3)
     +                                         * INT2DZ (R+3,ZE,ZC)
     +                                         + TEMP2 (R+4)
     +                                         * INT2DZ (R+4,ZE,ZC)
     +                                         + TEMP2 (R+5)
     +                                         * INT2DZ (R+5,ZE,ZC)
     +                                         + TEMP2 (R+6)
     +                                         * INT2DZ (R+6,ZE,ZC)
     +                                         + TEMP2 (R+7)
     +                                         * INT2DZ (R+7,ZE,ZC)
                               R = R + 8
                            END DO
                        END IF

                        I = I - NXYZE + XE
                        NXYZE = NXYZE - SE - 1
  820                CONTINUE
  812             CONTINUE
  810          CONTINUE
  802       CONTINUE
            XCP = J
  800    CONTINUE

         RETURN
C
C
C                       ********************
C                       *  # of roots = 9  *
C                       ********************
C
C
    9    XCP = 0
         DO 900 XC = SHELLC,0,-1
            YCMAX = SHELLC - XC

            XEP = NXYZET + 3
            DO 902 XE = 0,SHELLP
               XEP = XEP + XE - 2
               XEMAX = XE * SHELLP
               YEEND = SHELLP - XE

               DO M = 1,NGQEXQ
                  TEMP1 (M) = SCALE (M) * INT2DX (M,XE,XC)
               END DO

               J = XCP
               DO 910 YC = YCMAX,0,-1
                  J = J + 1
                  ZC = YCMAX - YC

                  XYEP = XEP - XEMAX
                  DO 912 YE = 0,YEEND
                     XYE = XE + YE
                     XYEP = XYEP - 1
                     SEEND = MAX0 (SHELLA,XYE)

                     IF (.NOT.DIFFY .AND. YE+YC.EQ.0) THEN
                         DO N = 1,NGQEXQ
                            TEMP2 (N) = TEMP1 (N)
                         END DO
                     ELSE
                         DO N = 1,NGQEXQ
                            TEMP2 (N) = TEMP1 (N) * INT2DY (N,YE,YC)
                         END DO
                     END IF

                     I = XYEP
                     NXYZE = NXYZP
                     DO 920 SE = SHELLP,SEEND,-1
                        ZE = SE - XYE

                        IF (.NOT.DIFFZ .AND. ZE+ZC.EQ.0) THEN
                            R = 1
                            DO M = 1,NEXQ
                               BATCH (M,I,J) =   TEMP2 (R)
     +                                         + TEMP2 (R+1)
     +                                         + TEMP2 (R+2)
     +                                         + TEMP2 (R+3)
     +                                         + TEMP2 (R+4)
     +                                         + TEMP2 (R+5)
     +                                         + TEMP2 (R+6)
     +                                         + TEMP2 (R+7)
     +                                         + TEMP2 (R+8)
                               R = R + 9
                            END DO
                        ELSE
                            R = 1
                            DO M = 1,NEXQ
                               BATCH (M,I,J) =   TEMP2 (R)
     +                                         * INT2DZ (R,ZE,ZC)
     +                                         + TEMP2 (R+1)
     +                                         * INT2DZ (R+1,ZE,ZC)
     +                                         + TEMP2 (R+2)
     +                                         * INT2DZ (R+2,ZE,ZC)
     +                                         + TEMP2 (R+3)
     +                                         * INT2DZ (R+3,ZE,ZC)
     +                                         + TEMP2 (R+4)
     +                                         * INT2DZ (R+4,ZE,ZC)
     +                                         + TEMP2 (R+5)
     +                                         * INT2DZ (R+5,ZE,ZC)
     +                                         + TEMP2 (R+6)
     +                                         * INT2DZ (R+6,ZE,ZC)
     +                                         + TEMP2 (R+7)
     +                                         * INT2DZ (R+7,ZE,ZC)
     +                                         + TEMP2 (R+8)
     +                                         * INT2DZ (R+8,ZE,ZC)
                               R = R + 9
                            END DO
                        END IF

                        I = I - NXYZE + XE
                        NXYZE = NXYZE - SE - 1
  920                CONTINUE
  912             CONTINUE
  910          CONTINUE
  902       CONTINUE
            XCP = J
  900    CONTINUE

         RETURN
C
C
C                       ********************
C                       *  # of roots > 9  *
C                       ********************
C
C             ...outer loops over x,x-pairs. No skipping of the
C                x,x-contribution of 0,0-type can be done here,
C                since the 2DX integrals carry the Rys weight!
C
C
   10    XCP = 0
         DO 1000 XC = SHELLC,0,-1
            YCMAX = SHELLC - XC

            XEP = NXYZET + 3
            DO 1002 XE = 0,SHELLP
               XEP = XEP + XE - 2
               XEMAX = XE * SHELLP
               YEEND = SHELLP - XE

               DO M = 1,NGQEXQ
                  TEMP1 (M) = SCALE (M) * INT2DX (M,XE,XC)
               END DO
C
C
C             ...middle loops over y,y-pairs. Skip multiplication
C                of y,y-contributions, if no y-coordinate derivative
C                was formed and we have a 0,0-pair, as then the 2DY
C                integrals are equal to 1.
C
C
               J = XCP
               DO 1010 YC = YCMAX,0,-1
                  J = J + 1
                  ZC = YCMAX - YC

                  XYEP = XEP - XEMAX
                  DO 1012 YE = 0,YEEND
                     XYE = XE + YE
                     XYEP = XYEP - 1
                     SEEND = MAX0 (SHELLA,XYE)

                     IF (.NOT.DIFFY .AND. YE+YC.EQ.0) THEN
                         DO N = 1,NGQEXQ
                            TEMP2 (N) = TEMP1 (N)
                         END DO
                     ELSE
                         DO N = 1,NGQEXQ
                            TEMP2 (N) = TEMP1 (N) * INT2DY (N,YE,YC)
                         END DO
                     END IF
C
C
C             ...inner loop over E shells. Skip multiplication of
C                z,z-contributions, if we have a 0,0-pair and
C                no derivations were performed on the z-coordinate,
C                as then the 2DZ integrals are equal to 1.
C                All info concerning all three x,x-, y,y- and
C                z,z-pairs have been collected for all exponent
C                quadruplets at once. Sum up the 2D X,Y,Z integral
C                products to the appropriate place of the batch.
C
C
                     I = XYEP
                     NXYZE = NXYZP
                     DO 1020 SE = SHELLP,SEEND,-1
                        ZE = SE - XYE

                        IF (.NOT.DIFFZ .AND. ZE+ZC.EQ.0) THEN
                            R = 0
                            DO M = 1,NEXQ
                               SUM = ZERO
                               DO N = 1,NGQP
                                  SUM = SUM + TEMP2 (R+N)
                               END DO
                               R = R + NGQP
                               BATCH (M,I,J) = SUM
                            END DO
                        ELSE
                            R = 0
                            DO M = 1,NEXQ
                               SUM = ZERO
                               DO N = 1,NGQP
                                  SUM = SUM + TEMP2 (R+N)
     +                                      * INT2DZ (R+N,ZE,ZC)
                               END DO
                               R = R + NGQP
                               BATCH (M,I,J) = SUM
                            END DO
                        END IF
C
C
C             ...next z,z-pair.
C
C
                        I = I - NXYZE + XE
                        NXYZE = NXYZE - SE - 1
 1020                CONTINUE
C
C
C             ...next y,y-pair and next x,x-pair.
C
C
 1012             CONTINUE
 1010          CONTINUE
 1002       CONTINUE
            XCP = J
 1000    CONTINUE
C
C
C             ...ready!
C
C
         RETURN
         END
