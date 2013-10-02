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
         SUBROUTINE  ERD__DERV_INT2D_TO_E0CD
     +
     +                    ( SHELLA,SHELLP,SHELLC,SHELLD,
     +                      NGQP,NEXQ,NGQEXQ,
     +                      NXYZET,NXYZP,NXYZC,NXYZD,
     +                      INT2DX,INT2DY,INT2DZ,
     +                      DIFFY,DIFFZ,
     +                      TEMP1,TEMP2,
     +                      SCALE,
     +
     +                                BATCH )
     +
C------------------------------------------------------------------------
C  OPERATION   : ERD__DERV_INT2D_TO_E0CD
C  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT
C  MODULE-ID   : ERD
C  SUBROUTINES : none
C  DESCRIPTION : This routine assembles the set of batches of cartesian
C                derivative eris [E0|CD] , E = A to P, adding up the
C                contributions from all the 2D PCD derivative integrals.
C
C                The routine uses the reduced Rys multiplication scheme
C                as suggested in R.Lindh, U.Ryu and B.Liu, JCP 95, 5889.
C                This scheme reuses intermediate products between
C                2DX and 2DY integrals, which can be achieved by having
C                the outer loops run over all possible x and y monomial
C                parts and the inner loops over all allowed E shell
C                combinations. The price to pay for such loop ordering
C                is the scattered addressing of locations within the
C                batch array, which has its rows and columns ordered
C                such that the E shells are increasing and within each
C                E shell the monomials are ordered such that x>y>z in
C                exponents.
C
C                An example follows:
C                -------------------
C
C                Let E = 0,2 and C = 0 and D = 1. Then we have the left
C                and right hand of the batch array ordered as follows:
C
C                          left xyz              right xyz
C
C                             000               000     100
C                             ---                       010
C                             100                       001
C                             010
C                             001
C                             ---
C                             200 -> -5
C                             110
C                             101 -> -3
C                             020
C                             011
C                             002 -> 0
C                        
C                The batch would thus have dimensions 10 x 1 x 3.
C                For the left side the reduced multiplication scheme
C                would have its outermost loop run over x=0,2, followed
C                by the next loop y=0,2-x. The innermost loop would then
C                run over the allowed shells E=E(max),max(E(min),x+y).
C                In this case all x,y-triples can be reused for all
C                appropriate shell combinations.
C
C                To find the address of a specific x,y,z,E combination
C                inside the batch array, we first note that the z-part
C                is dependent on the x,y-parts and is hence not needed.
C                Lets look at the E-part first. The E-part is evaluated
C                from its dimension formula (E+1)*(E+2)/2. Organizing
C                the inner E-loop to run from E(max) always, the
C                dimension for E(max) is passed as the argument NXYZP
C                and all lower E dimensions are calculated by the
C                formula relating dimensions between E and E+1:
C
C                           dim(E+1) = dim(E) + E + 2
C
C                In this way multiplications in the E-part can be
C                entirely avoided. The x,y-part is defined as the
C                part which has to be subtracted from dim(E) to
C                reach the xyz monomial position inside the E shell.
C                It can be divided into an x-part and a y-part. The
C                x-part is given by the fomula:
C
C                          x-part = - x*E + x(x-3)/2
C
C                and for the example above has been given for E=2
C                and marked with arrows ->. The last term of the x-part
C                involves 1 multiplication and division, however it
C                can be changed to:
C
C                                               x-1
C                          x-part = - x*E - x + sum i
C                                               i=0
C
C                and clever additions inside the x-loop avoid the use
C                of multiplications and divisions. The y-part is trivial
C                and is simply equal to -y. The overall conclusion is
C                thus that the location of a specific x,y,z,E quadruple
C                inside the batch comes at the cost of one x*E(max)
C                multiplication in the outermost x-loops, since the
C                other x*E ones can again be reached via stepwise
C                subtraction of x from x*E(max).
C
C                Due to the very computational intensive steps inside
C                the x,y,z,E loops, special sections of identical
C                x,y,z,E loop structures have been given for each
C                # of roots =< 9, thus saving considerable computing
C                time over the general case.
C
C                For comments on how the x,y,z,E loop structures are
C                coded please refer to the general root case.
C
C
C                  Input:
C
C                    SHELLx      =  shell types for individual csh
C                                   x=A,C,D and csh sum P=A+B
C                    NGQP        =  # of gaussian quadrature points
C                                   (roots)
C                    NEXQ        =  current # of exponent quadruplets
C                    NGQEXQ      =  product of # of gaussian quadrature
C                                   points times exponent quadruplets
C                    NXYZET      =  sum of # of cartesian monomials
C                                   for all shells in the range
C                                   E = A,...,P=A+B
C                    NXYZy       =  # of cartesian monomials for
C                                   y = P,C,D shells
C                    INT2Dx      =  all current 2D PCD derivative
C                                   integrals for each cartesian
C                                   component (x = X,Y,Z)
C                    DIFFx       =  is true, if differentiation was
C                                   performed along the x=Y,Z direction
C                    TEMP1(2)    =  scratch arrays holding intermediate
C                                   2D PCD derivative integral products
C                    SCALE       =  the NGQEXQ scaling factors
C
C
C                  Output:
C
C                    BATCH       =  batch of primitive cartesian
C                                   [E0|CD] derivative integrals
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

         INTEGER     I,J,K,M,N,R
         INTEGER     NGQP,NEXQ,NGQEXQ
         INTEGER     NXYZE,NXYZET,NXYZP,NXYZC,NXYZD
         INTEGER     SE
         INTEGER     SEEND
         INTEGER     SHELLA,SHELLP,SHELLC,SHELLD
         INTEGER     XC,YC,ZC,XD,YD,ZD
         INTEGER     XCP,XDP
         INTEGER     XE,YE,ZE
         INTEGER     XEMAX
         INTEGER     XEP,XYEP
         INTEGER     XYE
         INTEGER     YCMAX,YDMAX
         INTEGER     YEEND

         DOUBLE PRECISION  SUM
         DOUBLE PRECISION  ZERO

         DOUBLE PRECISION  SCALE (1:NGQEXQ)
         DOUBLE PRECISION  TEMP1 (1:NGQEXQ)
         DOUBLE PRECISION  TEMP2 (1:NGQEXQ)

         DOUBLE PRECISION  BATCH (1:NEXQ,1:NXYZET,1:NXYZC,1:NXYZD)

         DOUBLE PRECISION  INT2DX (1:NGQEXQ,0:SHELLP,0:SHELLC,0:SHELLD)
         DOUBLE PRECISION  INT2DY (1:NGQEXQ,0:SHELLP,0:SHELLC,0:SHELLD)
         DOUBLE PRECISION  INT2DZ (1:NGQEXQ,0:SHELLP,0:SHELLC,0:SHELLD)

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
    1    XDP = 0
         DO 100 XD = SHELLD,0,-1
            YDMAX = SHELLD - XD
            XCP = 0
            DO 102 XC = SHELLC,0,-1
               YCMAX = SHELLC - XC

               XEP = NXYZET + 3
               DO 104 XE = 0,SHELLP
                  XEP = XEP + XE - 2
                  XEMAX = XE * SHELLP
                  YEEND = SHELLP - XE

                  DO M = 1,NEXQ
                     TEMP1 (M) = SCALE (M) * INT2DX (M,XE,XC,XD)
                  END DO

                  K = XDP
                  DO 110 YD = YDMAX,0,-1
                     K = K + 1
                     ZD = YDMAX - YD
                     J = XCP
                     DO 112 YC = YCMAX,0,-1
                        J = J + 1
                        ZC = YCMAX - YC

                        XYEP = XEP - XEMAX
                        DO 114 YE = 0,YEEND
                           XYE = XE + YE
                           XYEP = XYEP - 1
                           SEEND = MAX0 (SHELLA,XYE)

                           IF (.NOT.DIFFY .AND. YE+YC+YD.EQ.0) THEN
                               DO M = 1,NEXQ
                                  TEMP2 (M) = TEMP1 (M)
                               END DO
                           ELSE
                               DO M = 1,NEXQ
                                  TEMP2 (M) = TEMP1 (M)
     +                                      * INT2DY (M,YE,YC,YD)
                               END DO
                           END IF

                           I = XYEP
                           NXYZE = NXYZP
                           DO 120 SE = SHELLP,SEEND,-1
                              ZE = SE - XYE

                              IF (.NOT.DIFFZ .AND. ZE+ZC+ZD.EQ.0) THEN
                                  DO M = 1,NEXQ
                                     BATCH (M,I,J,K) = TEMP2 (M)
                                  END DO
                              ELSE
                                  DO M = 1,NEXQ
                                     BATCH (M,I,J,K) =
     +                                            TEMP2 (M)
     +                                          * INT2DZ (M,ZE,ZC,ZD)
                                  END DO
                              END IF

                              I = I - NXYZE + XE
                              NXYZE = NXYZE - SE - 1
  120                      CONTINUE
  114                   CONTINUE
  112                CONTINUE
  110             CONTINUE
  104          CONTINUE
               XCP = J
  102       CONTINUE
            XDP = K
  100    CONTINUE

         RETURN
C
C
C                       ********************
C                       *  # of roots = 2  *
C                       ********************
C
C
    2    XDP = 0
         DO 200 XD = SHELLD,0,-1
            YDMAX = SHELLD - XD
            XCP = 0
            DO 202 XC = SHELLC,0,-1
               YCMAX = SHELLC - XC

               XEP = NXYZET + 3
               DO 204 XE = 0,SHELLP
                  XEP = XEP + XE - 2
                  XEMAX = XE * SHELLP
                  YEEND = SHELLP - XE

                  DO M = 1,NGQEXQ
                     TEMP1 (M) = SCALE (M) * INT2DX (M,XE,XC,XD)
                  END DO

                  K = XDP
                  DO 210 YD = YDMAX,0,-1
                     K = K + 1
                     ZD = YDMAX - YD
                     J = XCP
                     DO 212 YC = YCMAX,0,-1
                        J = J + 1
                        ZC = YCMAX - YC

                        XYEP = XEP - XEMAX
                        DO 214 YE = 0,YEEND
                           XYE = XE + YE
                           XYEP = XYEP - 1
                           SEEND = MAX0 (SHELLA,XYE)

                           IF (.NOT.DIFFY .AND. YE+YC+YD.EQ.0) THEN
                               DO N = 1,NGQEXQ
                                  TEMP2 (N) = TEMP1 (N)
                               END DO
                           ELSE
                               DO N = 1,NGQEXQ
                                  TEMP2 (N) = TEMP1 (N)
     +                                      * INT2DY (N,YE,YC,YD)
                               END DO
                           END IF

                           I = XYEP
                           NXYZE = NXYZP
                           DO 220 SE = SHELLP,SEEND,-1
                              ZE = SE - XYE

                              IF (.NOT.DIFFZ .AND. ZE+ZC+ZD.EQ.0) THEN
                                  R = 1
                                  DO M = 1,NEXQ
                                     BATCH (M,I,J,K) =   TEMP2 (R)
     +                                                 + TEMP2 (R+1)
                                     R = R + 2
                                  END DO
                              ELSE
                                  R = 1
                                  DO M = 1,NEXQ
                                     BATCH (M,I,J,K) =
     +                                            TEMP2 (R)
     +                                          * INT2DZ (R,ZE,ZC,ZD)
     +                                          + TEMP2 (R+1)
     +                                          * INT2DZ (R+1,ZE,ZC,ZD)
                                     R = R + 2
                                  END DO
                              END IF

                              I = I - NXYZE + XE
                              NXYZE = NXYZE - SE - 1
  220                      CONTINUE
  214                   CONTINUE
  212                CONTINUE
  210             CONTINUE
  204          CONTINUE
               XCP = J
  202       CONTINUE
            XDP = K
  200    CONTINUE

         RETURN
C
C
C                       ********************
C                       *  # of roots = 3  *
C                       ********************
C
C
    3    XDP = 0
         DO 300 XD = SHELLD,0,-1
            YDMAX = SHELLD - XD
            XCP = 0
            DO 302 XC = SHELLC,0,-1
               YCMAX = SHELLC - XC

               XEP = NXYZET + 3
               DO 304 XE = 0,SHELLP
                  XEP = XEP + XE - 2
                  XEMAX = XE * SHELLP
                  YEEND = SHELLP - XE

                  DO M = 1,NGQEXQ
                     TEMP1 (M) = SCALE (M) * INT2DX (M,XE,XC,XD)
                  END DO

                  K = XDP
                  DO 310 YD = YDMAX,0,-1
                     K = K + 1
                     ZD = YDMAX - YD
                     J = XCP
                     DO 312 YC = YCMAX,0,-1
                        J = J + 1
                        ZC = YCMAX - YC

                        XYEP = XEP - XEMAX
                        DO 314 YE = 0,YEEND
                           XYE = XE + YE
                           XYEP = XYEP - 1
                           SEEND = MAX0 (SHELLA,XYE)

                           IF (.NOT.DIFFY .AND. YE+YC+YD.EQ.0) THEN
                               DO N = 1,NGQEXQ
                                  TEMP2 (N) = TEMP1 (N)
                               END DO
                           ELSE
                               DO N = 1,NGQEXQ
                                  TEMP2 (N) = TEMP1 (N)
     +                                      * INT2DY (N,YE,YC,YD)
                               END DO
                           END IF

                           I = XYEP
                           NXYZE = NXYZP
                           DO 320 SE = SHELLP,SEEND,-1
                              ZE = SE - XYE

                              IF (.NOT.DIFFZ .AND. ZE+ZC+ZD.EQ.0) THEN
                                  R = 1
                                  DO M = 1,NEXQ
                                     BATCH (M,I,J,K) =   TEMP2 (R)
     +                                                 + TEMP2 (R+1)
     +                                                 + TEMP2 (R+2)
                                     R = R + 3
                                  END DO
                              ELSE
                                  R = 1
                                  DO M = 1,NEXQ
                                     BATCH (M,I,J,K) =
     +                                            TEMP2 (R)
     +                                          * INT2DZ (R,ZE,ZC,ZD)
     +                                          + TEMP2 (R+1)
     +                                          * INT2DZ (R+1,ZE,ZC,ZD)
     +                                          + TEMP2 (R+2)
     +                                          * INT2DZ (R+2,ZE,ZC,ZD)
                                     R = R + 3
                                  END DO
                              END IF

                              I = I - NXYZE + XE
                              NXYZE = NXYZE - SE - 1
  320                      CONTINUE
  314                   CONTINUE
  312                CONTINUE
  310             CONTINUE
  304          CONTINUE
               XCP = J
  302       CONTINUE
            XDP = K
  300    CONTINUE

         RETURN
C
C
C                       ********************
C                       *  # of roots = 4  *
C                       ********************
C
C
    4    XDP = 0
         DO 400 XD = SHELLD,0,-1
            YDMAX = SHELLD - XD
            XCP = 0
            DO 402 XC = SHELLC,0,-1
               YCMAX = SHELLC - XC

               XEP = NXYZET + 3
               DO 404 XE = 0,SHELLP
                  XEP = XEP + XE - 2
                  XEMAX = XE * SHELLP
                  YEEND = SHELLP - XE

                  DO M = 1,NGQEXQ
                     TEMP1 (M) = SCALE (M) * INT2DX (M,XE,XC,XD)
                  END DO

                  K = XDP
                  DO 410 YD = YDMAX,0,-1
                     K = K + 1
                     ZD = YDMAX - YD
                     J = XCP
                     DO 412 YC = YCMAX,0,-1
                        J = J + 1
                        ZC = YCMAX - YC

                        XYEP = XEP - XEMAX
                        DO 414 YE = 0,YEEND
                           XYE = XE + YE
                           XYEP = XYEP - 1
                           SEEND = MAX0 (SHELLA,XYE)

                           IF (.NOT.DIFFY .AND. YE+YC+YD.EQ.0) THEN
                               DO N = 1,NGQEXQ
                                  TEMP2 (N) = TEMP1 (N)
                               END DO
                           ELSE
                               DO N = 1,NGQEXQ
                                  TEMP2 (N) = TEMP1 (N)
     +                                      * INT2DY (N,YE,YC,YD)
                               END DO
                           END IF

                           I = XYEP
                           NXYZE = NXYZP
                           DO 420 SE = SHELLP,SEEND,-1
                              ZE = SE - XYE

                              IF (.NOT.DIFFZ .AND. ZE+ZC+ZD.EQ.0) THEN
                                  R = 1
                                  DO M = 1,NEXQ
                                     BATCH (M,I,J,K) =   TEMP2 (R)
     +                                                 + TEMP2 (R+1)
     +                                                 + TEMP2 (R+2)
     +                                                 + TEMP2 (R+3)
                                     R = R + 4
                                  END DO
                              ELSE
                                  R = 1
                                  DO M = 1,NEXQ
                                     BATCH (M,I,J,K) =
     +                                            TEMP2 (R)
     +                                          * INT2DZ (R,ZE,ZC,ZD)
     +                                          + TEMP2 (R+1)
     +                                          * INT2DZ (R+1,ZE,ZC,ZD)
     +                                          + TEMP2 (R+2)
     +                                          * INT2DZ (R+2,ZE,ZC,ZD)
     +                                          + TEMP2 (R+3)
     +                                          * INT2DZ (R+3,ZE,ZC,ZD)
                                     R = R + 4
                                  END DO
                              END IF

                              I = I - NXYZE + XE
                              NXYZE = NXYZE - SE - 1
  420                      CONTINUE
  414                   CONTINUE
  412                CONTINUE
  410             CONTINUE
  404          CONTINUE
               XCP = J
  402       CONTINUE
            XDP = K
  400    CONTINUE

         RETURN
C
C
C                       ********************
C                       *  # of roots = 5  *
C                       ********************
C
C
    5    XDP = 0
         DO 500 XD = SHELLD,0,-1
            YDMAX = SHELLD - XD
            XCP = 0
            DO 502 XC = SHELLC,0,-1
               YCMAX = SHELLC - XC

               XEP = NXYZET + 3
               DO 504 XE = 0,SHELLP
                  XEP = XEP + XE - 2
                  XEMAX = XE * SHELLP
                  YEEND = SHELLP - XE

                  DO M = 1,NGQEXQ
                     TEMP1 (M) = SCALE (M) * INT2DX (M,XE,XC,XD)
                  END DO

                  K = XDP
                  DO 510 YD = YDMAX,0,-1
                     K = K + 1
                     ZD = YDMAX - YD
                     J = XCP
                     DO 512 YC = YCMAX,0,-1
                        J = J + 1
                        ZC = YCMAX - YC

                        XYEP = XEP - XEMAX
                        DO 514 YE = 0,YEEND
                           XYE = XE + YE
                           XYEP = XYEP - 1
                           SEEND = MAX0 (SHELLA,XYE)

                           IF (.NOT.DIFFY .AND. YE+YC+YD.EQ.0) THEN
                               DO N = 1,NGQEXQ
                                  TEMP2 (N) = TEMP1 (N)
                               END DO
                           ELSE
                               DO N = 1,NGQEXQ
                                  TEMP2 (N) = TEMP1 (N)
     +                                      * INT2DY (N,YE,YC,YD)
                               END DO
                           END IF

                           I = XYEP
                           NXYZE = NXYZP
                           DO 520 SE = SHELLP,SEEND,-1
                              ZE = SE - XYE

                              IF (.NOT.DIFFZ .AND. ZE+ZC+ZD.EQ.0) THEN
                                  R = 1
                                  DO M = 1,NEXQ
                                     BATCH (M,I,J,K) =   TEMP2 (R)
     +                                                 + TEMP2 (R+1)
     +                                                 + TEMP2 (R+2)
     +                                                 + TEMP2 (R+3)
     +                                                 + TEMP2 (R+4)
                                     R = R + 5
                                  END DO
                              ELSE
                                  R = 1
                                  DO M = 1,NEXQ
                                     BATCH (M,I,J,K) =
     +                                            TEMP2 (R)
     +                                          * INT2DZ (R,ZE,ZC,ZD)
     +                                          + TEMP2 (R+1)
     +                                          * INT2DZ (R+1,ZE,ZC,ZD)
     +                                          + TEMP2 (R+2)
     +                                          * INT2DZ (R+2,ZE,ZC,ZD)
     +                                          + TEMP2 (R+3)
     +                                          * INT2DZ (R+3,ZE,ZC,ZD)
     +                                          + TEMP2 (R+4)
     +                                          * INT2DZ (R+4,ZE,ZC,ZD)
                                     R = R + 5
                                  END DO
                              END IF

                              I = I - NXYZE + XE
                              NXYZE = NXYZE - SE - 1
  520                      CONTINUE
  514                   CONTINUE
  512                CONTINUE
  510             CONTINUE
  504          CONTINUE
               XCP = J
  502       CONTINUE
            XDP = K
  500    CONTINUE

         RETURN
C
C
C                       ********************
C                       *  # of roots = 6  *
C                       ********************
C
C
    6    XDP = 0
         DO 600 XD = SHELLD,0,-1
            YDMAX = SHELLD - XD
            XCP = 0
            DO 602 XC = SHELLC,0,-1
               YCMAX = SHELLC - XC

               XEP = NXYZET + 3
               DO 604 XE = 0,SHELLP
                  XEP = XEP + XE - 2
                  XEMAX = XE * SHELLP
                  YEEND = SHELLP - XE

                  DO M = 1,NGQEXQ
                     TEMP1 (M) = SCALE (M) * INT2DX (M,XE,XC,XD)
                  END DO

                  K = XDP
                  DO 610 YD = YDMAX,0,-1
                     K = K + 1
                     ZD = YDMAX - YD
                     J = XCP
                     DO 612 YC = YCMAX,0,-1
                        J = J + 1
                        ZC = YCMAX - YC

                        XYEP = XEP - XEMAX
                        DO 614 YE = 0,YEEND
                           XYE = XE + YE
                           XYEP = XYEP - 1
                           SEEND = MAX0 (SHELLA,XYE)

                           IF (.NOT.DIFFY .AND. YE+YC+YD.EQ.0) THEN
                               DO N = 1,NGQEXQ
                                  TEMP2 (N) = TEMP1 (N)
                               END DO
                           ELSE
                               DO N = 1,NGQEXQ
                                  TEMP2 (N) = TEMP1 (N)
     +                                      * INT2DY (N,YE,YC,YD)
                               END DO
                           END IF

                           I = XYEP
                           NXYZE = NXYZP
                           DO 620 SE = SHELLP,SEEND,-1
                              ZE = SE - XYE

                              IF (.NOT.DIFFZ .AND. ZE+ZC+ZD.EQ.0) THEN
                                  R = 1
                                  DO M = 1,NEXQ
                                     BATCH (M,I,J,K) =   TEMP2 (R)
     +                                                 + TEMP2 (R+1)
     +                                                 + TEMP2 (R+2)
     +                                                 + TEMP2 (R+3)
     +                                                 + TEMP2 (R+4)
     +                                                 + TEMP2 (R+5)
                                     R = R + 6
                                  END DO
                              ELSE
                                  R = 1
                                  DO M = 1,NEXQ
                                     BATCH (M,I,J,K) =
     +                                            TEMP2 (R)
     +                                          * INT2DZ (R,ZE,ZC,ZD)
     +                                          + TEMP2 (R+1)
     +                                          * INT2DZ (R+1,ZE,ZC,ZD)
     +                                          + TEMP2 (R+2)
     +                                          * INT2DZ (R+2,ZE,ZC,ZD)
     +                                          + TEMP2 (R+3)
     +                                          * INT2DZ (R+3,ZE,ZC,ZD)
     +                                          + TEMP2 (R+4)
     +                                          * INT2DZ (R+4,ZE,ZC,ZD)
     +                                          + TEMP2 (R+5)
     +                                          * INT2DZ (R+5,ZE,ZC,ZD)
                                     R = R + 6
                                  END DO
                              END IF

                              I = I - NXYZE + XE
                              NXYZE = NXYZE - SE - 1
  620                      CONTINUE
  614                   CONTINUE
  612                CONTINUE
  610             CONTINUE
  604          CONTINUE
               XCP = J
  602       CONTINUE
            XDP = K
  600    CONTINUE

         RETURN
C
C
C                       ********************
C                       *  # of roots = 7  *
C                       ********************
C
C
    7    XDP = 0
         DO 700 XD = SHELLD,0,-1
            YDMAX = SHELLD - XD
            XCP = 0
            DO 702 XC = SHELLC,0,-1
               YCMAX = SHELLC - XC

               XEP = NXYZET + 3
               DO 704 XE = 0,SHELLP
                  XEP = XEP + XE - 2
                  XEMAX = XE * SHELLP
                  YEEND = SHELLP - XE

                  DO M = 1,NGQEXQ
                     TEMP1 (M) = SCALE (M) * INT2DX (M,XE,XC,XD)
                  END DO

                  K = XDP
                  DO 710 YD = YDMAX,0,-1
                     K = K + 1
                     ZD = YDMAX - YD
                     J = XCP
                     DO 712 YC = YCMAX,0,-1
                        J = J + 1
                        ZC = YCMAX - YC

                        XYEP = XEP - XEMAX
                        DO 714 YE = 0,YEEND
                           XYE = XE + YE
                           XYEP = XYEP - 1
                           SEEND = MAX0 (SHELLA,XYE)

                           IF (.NOT.DIFFY .AND. YE+YC+YD.EQ.0) THEN
                               DO N = 1,NGQEXQ
                                  TEMP2 (N) = TEMP1 (N)
                               END DO
                           ELSE
                               DO N = 1,NGQEXQ
                                  TEMP2 (N) = TEMP1 (N)
     +                                      * INT2DY (N,YE,YC,YD)
                               END DO
                           END IF

                           I = XYEP
                           NXYZE = NXYZP
                           DO 720 SE = SHELLP,SEEND,-1
                              ZE = SE - XYE

                              IF (.NOT.DIFFZ .AND. ZE+ZC+ZD.EQ.0) THEN
                                  R = 1
                                  DO M = 1,NEXQ
                                     BATCH (M,I,J,K) =   TEMP2 (R)
     +                                                 + TEMP2 (R+1)
     +                                                 + TEMP2 (R+2)
     +                                                 + TEMP2 (R+3)
     +                                                 + TEMP2 (R+4)
     +                                                 + TEMP2 (R+5)
     +                                                 + TEMP2 (R+6)
                                     R = R + 7
                                  END DO
                              ELSE
                                  R = 1
                                  DO M = 1,NEXQ
                                     BATCH (M,I,J,K) =
     +                                            TEMP2 (R)
     +                                          * INT2DZ (R,ZE,ZC,ZD)
     +                                          + TEMP2 (R+1)
     +                                          * INT2DZ (R+1,ZE,ZC,ZD)
     +                                          + TEMP2 (R+2)
     +                                          * INT2DZ (R+2,ZE,ZC,ZD)
     +                                          + TEMP2 (R+3)
     +                                          * INT2DZ (R+3,ZE,ZC,ZD)
     +                                          + TEMP2 (R+4)
     +                                          * INT2DZ (R+4,ZE,ZC,ZD)
     +                                          + TEMP2 (R+5)
     +                                          * INT2DZ (R+5,ZE,ZC,ZD)
     +                                          + TEMP2 (R+6)
     +                                          * INT2DZ (R+6,ZE,ZC,ZD)
                                     R = R + 7
                                  END DO
                              END IF

                              I = I - NXYZE + XE
                              NXYZE = NXYZE - SE - 1
  720                      CONTINUE
  714                   CONTINUE
  712                CONTINUE
  710             CONTINUE
  704          CONTINUE
               XCP = J
  702       CONTINUE
            XDP = K
  700    CONTINUE

         RETURN
C
C
C                       ********************
C                       *  # of roots = 8  *
C                       ********************
C
C
    8    XDP = 0
         DO 800 XD = SHELLD,0,-1
            YDMAX = SHELLD - XD
            XCP = 0
            DO 802 XC = SHELLC,0,-1
               YCMAX = SHELLC - XC

               XEP = NXYZET + 3
               DO 804 XE = 0,SHELLP
                  XEP = XEP + XE - 2
                  XEMAX = XE * SHELLP
                  YEEND = SHELLP - XE

                  DO M = 1,NGQEXQ
                     TEMP1 (M) = SCALE (M) * INT2DX (M,XE,XC,XD)
                  END DO

                  K = XDP
                  DO 810 YD = YDMAX,0,-1
                     K = K + 1
                     ZD = YDMAX - YD
                     J = XCP
                     DO 812 YC = YCMAX,0,-1
                        J = J + 1
                        ZC = YCMAX - YC

                        XYEP = XEP - XEMAX
                        DO 814 YE = 0,YEEND
                           XYE = XE + YE
                           XYEP = XYEP - 1
                           SEEND = MAX0 (SHELLA,XYE)

                           IF (.NOT.DIFFY .AND. YE+YC+YD.EQ.0) THEN
                               DO N = 1,NGQEXQ
                                  TEMP2 (N) = TEMP1 (N)
                               END DO
                           ELSE
                               DO N = 1,NGQEXQ
                                  TEMP2 (N) = TEMP1 (N)
     +                                      * INT2DY (N,YE,YC,YD)
                               END DO
                           END IF

                           I = XYEP
                           NXYZE = NXYZP
                           DO 820 SE = SHELLP,SEEND,-1
                              ZE = SE - XYE

                              IF (.NOT.DIFFZ .AND. ZE+ZC+ZD.EQ.0) THEN
                                  R = 1
                                  DO M = 1,NEXQ
                                     BATCH (M,I,J,K) =   TEMP2 (R)
     +                                                 + TEMP2 (R+1)
     +                                                 + TEMP2 (R+2)
     +                                                 + TEMP2 (R+3)
     +                                                 + TEMP2 (R+4)
     +                                                 + TEMP2 (R+5)
     +                                                 + TEMP2 (R+6)
     +                                                 + TEMP2 (R+7)
                                     R = R + 8
                                  END DO
                              ELSE
                                  R = 1
                                  DO M = 1,NEXQ
                                     BATCH (M,I,J,K) =
     +                                            TEMP2 (R)
     +                                          * INT2DZ (R,ZE,ZC,ZD)
     +                                          + TEMP2 (R+1)
     +                                          * INT2DZ (R+1,ZE,ZC,ZD)
     +                                          + TEMP2 (R+2)
     +                                          * INT2DZ (R+2,ZE,ZC,ZD)
     +                                          + TEMP2 (R+3)
     +                                          * INT2DZ (R+3,ZE,ZC,ZD)
     +                                          + TEMP2 (R+4)
     +                                          * INT2DZ (R+4,ZE,ZC,ZD)
     +                                          + TEMP2 (R+5)
     +                                          * INT2DZ (R+5,ZE,ZC,ZD)
     +                                          + TEMP2 (R+6)
     +                                          * INT2DZ (R+6,ZE,ZC,ZD)
     +                                          + TEMP2 (R+7)
     +                                          * INT2DZ (R+7,ZE,ZC,ZD)
                                     R = R + 8
                                  END DO
                              END IF

                              I = I - NXYZE + XE
                              NXYZE = NXYZE - SE - 1
  820                      CONTINUE
  814                   CONTINUE
  812                CONTINUE
  810             CONTINUE
  804          CONTINUE
               XCP = J
  802       CONTINUE
            XDP = K
  800    CONTINUE

         RETURN
C
C
C                       ********************
C                       *  # of roots = 9  *
C                       ********************
C
C
    9    XDP = 0
         DO 900 XD = SHELLD,0,-1
            YDMAX = SHELLD - XD
            XCP = 0
            DO 902 XC = SHELLC,0,-1
               YCMAX = SHELLC - XC

               XEP = NXYZET + 3
               DO 904 XE = 0,SHELLP
                  XEP = XEP + XE - 2
                  XEMAX = XE * SHELLP
                  YEEND = SHELLP - XE

                  DO M = 1,NGQEXQ
                     TEMP1 (M) = SCALE (M) * INT2DX (M,XE,XC,XD)
                  END DO

                  K = XDP
                  DO 910 YD = YDMAX,0,-1
                     K = K + 1
                     ZD = YDMAX - YD
                     J = XCP
                     DO 912 YC = YCMAX,0,-1
                        J = J + 1
                        ZC = YCMAX - YC

                        XYEP = XEP - XEMAX
                        DO 914 YE = 0,YEEND
                           XYE = XE + YE
                           XYEP = XYEP - 1
                           SEEND = MAX0 (SHELLA,XYE)

                           IF (.NOT.DIFFY .AND. YE+YC+YD.EQ.0) THEN
                               DO N = 1,NGQEXQ
                                  TEMP2 (N) = TEMP1 (N)
                               END DO
                           ELSE
                               DO N = 1,NGQEXQ
                                  TEMP2 (N) = TEMP1 (N)
     +                                      * INT2DY (N,YE,YC,YD)
                               END DO
                           END IF

                           I = XYEP
                           NXYZE = NXYZP
                           DO 920 SE = SHELLP,SEEND,-1
                              ZE = SE - XYE

                              IF (.NOT.DIFFZ .AND. ZE+ZC+ZD.EQ.0) THEN
                                  R = 1
                                  DO M = 1,NEXQ
                                     BATCH (M,I,J,K) =   TEMP2 (R)
     +                                                 + TEMP2 (R+1)
     +                                                 + TEMP2 (R+2)
     +                                                 + TEMP2 (R+3)
     +                                                 + TEMP2 (R+4)
     +                                                 + TEMP2 (R+5)
     +                                                 + TEMP2 (R+6)
     +                                                 + TEMP2 (R+7)
     +                                                 + TEMP2 (R+8)
                                     R = R + 9
                                  END DO
                              ELSE
                                  R = 1
                                  DO M = 1,NEXQ
                                     BATCH (M,I,J,K) =
     +                                            TEMP2 (R)
     +                                          * INT2DZ (R,ZE,ZC,ZD)
     +                                          + TEMP2 (R+1)
     +                                          * INT2DZ (R+1,ZE,ZC,ZD)
     +                                          + TEMP2 (R+2)
     +                                          * INT2DZ (R+2,ZE,ZC,ZD)
     +                                          + TEMP2 (R+3)
     +                                          * INT2DZ (R+3,ZE,ZC,ZD)
     +                                          + TEMP2 (R+4)
     +                                          * INT2DZ (R+4,ZE,ZC,ZD)
     +                                          + TEMP2 (R+5)
     +                                          * INT2DZ (R+5,ZE,ZC,ZD)
     +                                          + TEMP2 (R+6)
     +                                          * INT2DZ (R+6,ZE,ZC,ZD)
     +                                          + TEMP2 (R+7)
     +                                          * INT2DZ (R+7,ZE,ZC,ZD)
     +                                          + TEMP2 (R+8)
     +                                          * INT2DZ (R+8,ZE,ZC,ZD)
                                     R = R + 9
                                  END DO
                              END IF

                              I = I - NXYZE + XE
                              NXYZE = NXYZE - SE - 1
  920                      CONTINUE
  914                   CONTINUE
  912                CONTINUE
  910             CONTINUE
  904          CONTINUE
               XCP = J
  902       CONTINUE
            XDP = K
  900    CONTINUE

         RETURN
C
C
C                       ********************
C                       *  # of roots > 9  *
C                       ********************
C
C             ...outer loops over x,x,x-triples. No skipping of the
C                x,x,x-contribution of 0,0,0-type can be done here,
C                since the 2DX integrals carry the Rys weight!
C
C
   10    XDP = 0
         DO 1000 XD = SHELLD,0,-1
            YDMAX = SHELLD - XD
            XCP = 0
            DO 1002 XC = SHELLC,0,-1
               YCMAX = SHELLC - XC

               XEP = NXYZET + 3
               DO 1004 XE = 0,SHELLP
                  XEP = XEP + XE - 2
                  XEMAX = XE * SHELLP
                  YEEND = SHELLP - XE

                  DO M = 1,NGQEXQ
                     TEMP1 (M) = SCALE (M) * INT2DX (M,XE,XC,XD)
                  END DO
C
C
C             ...middle loops over y,y,y-triples. Skip multiplication
C                of y,y,y-contributions, if no y-coordinate derivative
C                was formed and we have a 0,0,0-triple, as then the 2DY
C                PCD integrals are equal to 1.
C
C
                  K = XDP
                  DO 1010 YD = YDMAX,0,-1
                     K = K + 1
                     ZD = YDMAX - YD
                     J = XCP
                     DO 1012 YC = YCMAX,0,-1
                        J = J + 1
                        ZC = YCMAX - YC

                        XYEP = XEP - XEMAX
                        DO 1014 YE = 0,YEEND
                           XYE = XE + YE
                           XYEP = XYEP - 1
                           SEEND = MAX0 (SHELLA,XYE)

                           IF (.NOT.DIFFY .AND. YE+YC+YD.EQ.0) THEN
                               DO N = 1,NGQEXQ
                                  TEMP2 (N) = TEMP1 (N)
                               END DO
                           ELSE
                               DO N = 1,NGQEXQ
                                  TEMP2 (N) = TEMP1 (N)
     +                                      * INT2DY (N,YE,YC,YD)
                               END DO
                           END IF
C
C
C             ...inner loop over E shells. Skip multiplication of
C                z,z,z-contributions, if we have a 0,0,0-triple and
C                no derivations were performed on the z-coordinate,
C                as then the 2DZ PCD integrals are equal to 1.
C                All info concerning all three x,x,x-, y,y,y- and
C                z,z,z-triples have been collected for all exponent
C                quadruplets at once. Sum up the 2D X,Y,Z integral
C                products to the appropriate place of the [E0|CD]
C                batch.
C
C
                           I = XYEP
                           NXYZE = NXYZP
                           DO 1020 SE = SHELLP,SEEND,-1
                              ZE = SE - XYE

                              IF (.NOT.DIFFZ .AND. ZE+ZC+ZD.EQ.0) THEN
                                  R = 0
                                  DO M = 1,NEXQ
                                     SUM = ZERO
                                     DO N = 1,NGQP
                                        SUM = SUM + TEMP2 (R+N)
                                     END DO
                                     R = R + NGQP
                                     BATCH (M,I,J,K) = SUM
                                  END DO
                              ELSE
                                  R = 0
                                  DO M = 1,NEXQ
                                     SUM = ZERO
                                     DO N = 1,NGQP
                                        SUM = SUM + TEMP2 (R+N)
     +                                           * INT2DZ (R+N,ZE,ZC,ZD)
                                     END DO
                                     R = R + NGQP
                                     BATCH (M,I,J,K) = SUM
                                  END DO
                              END IF
C
C
C             ...next z,z,z-triple.
C
C
                              I = I - NXYZE + XE
                              NXYZE = NXYZE - SE - 1
 1020                      CONTINUE
C
C
C             ...next y,y,y-triple and next x,x,x-triple.
C
C
 1014                   CONTINUE
 1012                CONTINUE
 1010             CONTINUE

 1004          CONTINUE
               XCP = J
 1002       CONTINUE
            XDP = K
 1000    CONTINUE
C
C
C             ...ready!
C
C
         RETURN
         END
