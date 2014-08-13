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
         SUBROUTINE  ERD__ATOM_INT2D_TO_E0F0
     +
     +                    ( SHELLA,SHELLP,SHELLC,SHELLQ,
     +                      NGQP,NEXQ,NGQEXQ,
     +                      NXYZET,NXYZFT,NXYZP,NXYZQ,
     +                      INT2DX,INT2DY,INT2DZ,
     +                      TEMP1,TEMP2,
     +                      SCALE,
     +
     +                                BATCH )
     +
C------------------------------------------------------------------------
C  OPERATION   : ERD__ATOM_INT2D_TO_E0F0
C  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT
C  MODULE-ID   : ERD
C  SUBROUTINES : none
C  DESCRIPTION : This routine assembles the set of atomic batches of
C                cartesian eris [E0|F0] , E = A to P, F = C to Q,
C                adding up all the contributions from all the atomic
C                2D PQ integrals.
C
C                The routine is a special routine derived from the
C                corresponding nonatomic one. Its structure is the
C                same except for the fact that for an atom individual
C                x,x-,y,y- and z,z-pair contributions are zero if
C                their monomial exponent sum is odd (vanishing integral
C                due to atomic D2h symmetry!). Hence checks of such
C                cases are placed inside the x,y,z loops.
C                
C                For comments on how the x,y,z,E loop structures are
C                coded please refer to the nonatomic routine.
C
C
C                  Input:
C
C                    SHELLx      =  shell types for individual csh
C                                   x=A,C and csh sums P=A+B,Q=C+D
C                    NGQP        =  # of gaussian quadrature points
C                                   (roots)
C                    NEXQ        =  current # of exponent quadruplets
C                    NGQEXQ      =  product of # of gaussian quadrature
C                                   points times exponent quadruplets
C                    NXYZE(F)T   =  sum of # of cartesian monomials
C                                   for all shells in the range
C                                   E = A,...,P=A+B and in the range
C                                   F = C,...,Q=C+D
C                    NXYZy       =  # of cartesian monomials for
C                                   y = P,Q shells
C                    INT2Dx      =  all current atomic 2D PQ integrals
C                                   for each cartesian component
C                                   (x = X,Y,Z)
C                    TEMP1(2)    =  scratch arrays holding intermediate
C                                   atomic 2D PQ integral products
C                    SCALE       =  the NGQEXQ scaling factors
C
C
C                  Output:
C
C                    BATCH       =  batch of atomic primitive cartesian
C                                   [E0|F0] integrals corresponding
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

         INTEGER     I,J,K,M,N
         INTEGER     NGQP,NEXQ,NGQEXQ
         INTEGER     NXYZE,NXYZF,NXYZET,NXYZFT,NXYZP,NXYZQ
         INTEGER     SE,SF
         INTEGER     SEEND,SFEND
         INTEGER     SHELLA,SHELLP,SHELLC,SHELLQ
         INTEGER     XE,YE,ZE
         INTEGER     XEMAX,XFMAX
         INTEGER     XEP,XFP,XYEP,XYFP
         INTEGER     XF,YF,ZF
         INTEGER     XYE,XYF
         INTEGER     YEEND,YFEND
         INTEGER     YTOT,ZTOT

         DOUBLE PRECISION  SUM
         DOUBLE PRECISION  ZERO

         DOUBLE PRECISION  SCALE (1:NGQEXQ)
         DOUBLE PRECISION  TEMP1 (1:NGQEXQ)
         DOUBLE PRECISION  TEMP2 (1:NGQEXQ)

         DOUBLE PRECISION  BATCH (1:NEXQ,1:NXYZET,1:NXYZFT)

         DOUBLE PRECISION  INT2DX (1:NGQEXQ,0:SHELLP,0:SHELLQ)
         DOUBLE PRECISION  INT2DY (1:NGQEXQ,0:SHELLP,0:SHELLQ)
         DOUBLE PRECISION  INT2DZ (1:NGQEXQ,0:SHELLP,0:SHELLQ)

         PARAMETER  (ZERO  =  0.D0)
C
C
C------------------------------------------------------------------------
C
C
C             ...set first all [E0|F0] batch elements equal to zero.
C                This is necessary, since not all elements of the
C                batch will be addressed inside the x,y,z loops.
C
C
         DO J = 1,NXYZFT
         DO I = 1,NXYZET
         DO M = 1,NEXQ
            BATCH (M,I,J) = ZERO
         END DO
         END DO
         END DO
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
    1    XFP = NXYZFT + 3
         DO 100 XF = 0,SHELLQ
            XFP = XFP + XF - 2
            XFMAX = XF * SHELLQ
            YFEND = SHELLQ - XF
            XEP = NXYZET + 3
            DO 102 XE = 0,SHELLP
               XEP = XEP + XE - 2
               XEMAX = XE * SHELLP
               YEEND = SHELLP - XE

               IF (MOD (XE+XF,2).EQ.1) GOTO 102

               DO M = 1,NEXQ
                  TEMP1 (M) = SCALE (M) * INT2DX (M,XE,XF)
               END DO

               XYFP = XFP - XFMAX
               DO 110 YF = 0,YFEND
                  XYF = XF + YF
                  XYFP = XYFP - 1
                  SFEND = MAX (SHELLC,XYF)
                  XYEP = XEP - XEMAX
                  DO 112 YE = 0,YEEND
                     XYE = XE + YE
                     XYEP = XYEP - 1
                     SEEND = MAX (SHELLA,XYE)

                     YTOT = YE + YF
                     IF (MOD (YTOT,2).EQ.1) GOTO 112

                     IF (YTOT.EQ.0) THEN
                         DO M = 1,NEXQ
                            TEMP2 (M) = TEMP1 (M)
                         END DO
                     ELSE
                         DO M = 1,NEXQ
                            TEMP2 (M) = TEMP1 (M) * INT2DY (M,YE,YF)
                         END DO
                     END IF

                     J = XYFP
                     NXYZF = NXYZQ
                     DO 120 SF = SHELLQ,SFEND,-1
                        ZF = SF - XYF
                        I = XYEP
                        NXYZE = NXYZP
                        DO 122 SE = SHELLP,SEEND,-1
                           ZE = SE - XYE
                           ZTOT = ZE + ZF

                           IF (MOD (ZTOT,2).EQ.0) THEN
                               IF (ZTOT.EQ.0) THEN
                                   DO M = 1,NEXQ
                                      BATCH (M,I,J) = TEMP2 (M)
                                   END DO
                               ELSE
                                   DO M = 1,NEXQ
                                      BATCH (M,I,J) = TEMP2 (M)
     +                                                * INT2DZ (M,ZE,ZF)
                                   END DO
                               END IF
                           END IF

                           I = I - NXYZE + XE
                           NXYZE = NXYZE - SE - 1
  122                   CONTINUE
                        J = J - NXYZF + XF
                        NXYZF = NXYZF - SF - 1
  120                CONTINUE
  112             CONTINUE
  110          CONTINUE
  102       CONTINUE
  100    CONTINUE

         RETURN
C
C
C                       ********************
C                       *  # of roots = 2  *
C                       ********************
C
C
    2    XFP = NXYZFT + 3
         DO 200 XF = 0,SHELLQ
            XFP = XFP + XF - 2
            XFMAX = XF * SHELLQ
            YFEND = SHELLQ - XF
            XEP = NXYZET + 3
            DO 202 XE = 0,SHELLP
               XEP = XEP + XE - 2
               XEMAX = XE * SHELLP
               YEEND = SHELLP - XE

               IF (MOD (XE+XF,2).EQ.1) GOTO 202

               DO M = 1,NGQEXQ
                  TEMP1 (M) = SCALE (M) * INT2DX (M,XE,XF)
               END DO

               XYFP = XFP - XFMAX
               DO 210 YF = 0,YFEND
                  XYF = XF + YF
                  XYFP = XYFP - 1
                  SFEND = MAX (SHELLC,XYF)
                  XYEP = XEP - XEMAX
                  DO 212 YE = 0,YEEND
                     XYE = XE + YE
                     XYEP = XYEP - 1
                     SEEND = MAX (SHELLA,XYE)

                     YTOT = YE + YF
                     IF (MOD (YTOT,2).EQ.1) GOTO 212

                     IF (YTOT.EQ.0) THEN
                         DO N = 1,NGQEXQ
                            TEMP2 (N) = TEMP1 (N)
                         END DO
                     ELSE
                         DO N = 1,NGQEXQ
                            TEMP2 (N) = TEMP1 (N) * INT2DY (N,YE,YF)
                         END DO
                     END IF

                     J = XYFP
                     NXYZF = NXYZQ
                     DO 220 SF = SHELLQ,SFEND,-1
                        ZF = SF - XYF
                        I = XYEP
                        NXYZE = NXYZP
                        DO 222 SE = SHELLP,SEEND,-1
                           ZE = SE - XYE
                           ZTOT = ZE + ZF

                           IF (MOD (ZTOT,2).EQ.0) THEN
                            IF (ZTOT.EQ.0) THEN
                               K = 1
                               DO M = 1,NEXQ
                                  BATCH (M,I,J) =   TEMP2 (K)
     +                                            + TEMP2 (K+1)
                                  K = K + 2
                               END DO
                            ELSE
                               K = 1
                               DO M = 1,NEXQ
                                  BATCH (M,I,J) =   TEMP2 (K)
     +                                            * INT2DZ (K,ZE,ZF)
     +                                            + TEMP2 (K+1)
     +                                            * INT2DZ (K+1,ZE,ZF)
                                  K = K + 2
                               END DO
                            END IF
                           END IF

                           I = I - NXYZE + XE
                           NXYZE = NXYZE - SE - 1
  222                   CONTINUE
                        J = J - NXYZF + XF
                        NXYZF = NXYZF - SF - 1
  220                CONTINUE
  212             CONTINUE
  210          CONTINUE
  202       CONTINUE
  200    CONTINUE

         RETURN
C
C
C                       ********************
C                       *  # of roots = 3  *
C                       ********************
C
C
    3    XFP = NXYZFT + 3
         DO 300 XF = 0,SHELLQ
            XFP = XFP + XF - 2
            XFMAX = XF * SHELLQ
            YFEND = SHELLQ - XF
            XEP = NXYZET + 3
            DO 302 XE = 0,SHELLP
               XEP = XEP + XE - 2
               XEMAX = XE * SHELLP
               YEEND = SHELLP - XE

               IF (MOD (XE+XF,2).EQ.1) GOTO 302

               DO M = 1,NGQEXQ
                  TEMP1 (M) = SCALE (M) * INT2DX (M,XE,XF)
               END DO

               XYFP = XFP - XFMAX
               DO 310 YF = 0,YFEND
                  XYF = XF + YF
                  XYFP = XYFP - 1
                  SFEND = MAX (SHELLC,XYF)
                  XYEP = XEP - XEMAX
                  DO 312 YE = 0,YEEND
                     XYE = XE + YE
                     XYEP = XYEP - 1
                     SEEND = MAX (SHELLA,XYE)

                     YTOT = YE + YF
                     IF (MOD (YTOT,2).EQ.1) GOTO 312

                     IF (YTOT.EQ.0) THEN
                         DO N = 1,NGQEXQ
                            TEMP2 (N) = TEMP1 (N)
                         END DO
                     ELSE
                         DO N = 1,NGQEXQ
                            TEMP2 (N) = TEMP1 (N) * INT2DY (N,YE,YF)
                         END DO
                     END IF

                     J = XYFP
                     NXYZF = NXYZQ
                     DO 320 SF = SHELLQ,SFEND,-1
                        ZF = SF - XYF
                        I = XYEP
                        NXYZE = NXYZP
                        DO 322 SE = SHELLP,SEEND,-1
                           ZE = SE - XYE
                           ZTOT = ZE + ZF

                           IF (MOD (ZTOT,2).EQ.0) THEN
                            IF (ZTOT.EQ.0) THEN
                               K = 1
                               DO M = 1,NEXQ
                                  BATCH (M,I,J) =   TEMP2 (K)
     +                                            + TEMP2 (K+1)
     +                                            + TEMP2 (K+2)
                                  K = K + 3
                               END DO
                            ELSE
                               K = 1
                               DO M = 1,NEXQ
                                  BATCH (M,I,J) =   TEMP2 (K)
     +                                            * INT2DZ (K,ZE,ZF)
     +                                            + TEMP2 (K+1)
     +                                            * INT2DZ (K+1,ZE,ZF)
     +                                            + TEMP2 (K+2)
     +                                            * INT2DZ (K+2,ZE,ZF)
                                  K = K + 3
                               END DO
                            END IF
                           END IF

                           I = I - NXYZE + XE
                           NXYZE = NXYZE - SE - 1
  322                   CONTINUE
                        J = J - NXYZF + XF
                        NXYZF = NXYZF - SF - 1
  320                CONTINUE
  312             CONTINUE
  310          CONTINUE
  302       CONTINUE
  300    CONTINUE

         RETURN
C
C
C                       ********************
C                       *  # of roots = 4  *
C                       ********************
C
C
    4    XFP = NXYZFT + 3
         DO 400 XF = 0,SHELLQ
            XFP = XFP + XF - 2
            XFMAX = XF * SHELLQ
            YFEND = SHELLQ - XF
            XEP = NXYZET + 3
            DO 402 XE = 0,SHELLP
               XEP = XEP + XE - 2
               XEMAX = XE * SHELLP
               YEEND = SHELLP - XE

               IF (MOD (XE+XF,2).EQ.1) GOTO 402

               DO M = 1,NGQEXQ
                  TEMP1 (M) = SCALE (M) * INT2DX (M,XE,XF)
               END DO

               XYFP = XFP - XFMAX
               DO 410 YF = 0,YFEND
                  XYF = XF + YF
                  XYFP = XYFP - 1
                  SFEND = MAX (SHELLC,XYF)
                  XYEP = XEP - XEMAX
                  DO 412 YE = 0,YEEND
                     XYE = XE + YE
                     XYEP = XYEP - 1
                     SEEND = MAX (SHELLA,XYE)

                     YTOT = YE + YF
                     IF (MOD (YTOT,2).EQ.1) GOTO 412

                     IF (YTOT.EQ.0) THEN
                         DO N = 1,NGQEXQ
                            TEMP2 (N) = TEMP1 (N)
                         END DO
                     ELSE
                         DO N = 1,NGQEXQ
                            TEMP2 (N) = TEMP1 (N) * INT2DY (N,YE,YF)
                         END DO
                     END IF

                     J = XYFP
                     NXYZF = NXYZQ
                     DO 420 SF = SHELLQ,SFEND,-1
                        ZF = SF - XYF
                        I = XYEP
                        NXYZE = NXYZP
                        DO 422 SE = SHELLP,SEEND,-1
                           ZE = SE - XYE
                           ZTOT = ZE + ZF

                           IF (MOD (ZTOT,2).EQ.0) THEN
                            IF (ZTOT.EQ.0) THEN
                               K = 1
                               DO M = 1,NEXQ
                                  BATCH (M,I,J) =   TEMP2 (K)
     +                                            + TEMP2 (K+1)
     +                                            + TEMP2 (K+2)
     +                                            + TEMP2 (K+3)
                                  K = K + 4
                               END DO
                            ELSE
                               K = 1
                               DO M = 1,NEXQ
                                  BATCH (M,I,J) =   TEMP2 (K)
     +                                            * INT2DZ (K,ZE,ZF)
     +                                            + TEMP2 (K+1)
     +                                            * INT2DZ (K+1,ZE,ZF)
     +                                            + TEMP2 (K+2)
     +                                            * INT2DZ (K+2,ZE,ZF)
     +                                            + TEMP2 (K+3)
     +                                            * INT2DZ (K+3,ZE,ZF)
                                  K = K + 4
                               END DO
                            END IF
                           END IF

                           I = I - NXYZE + XE
                           NXYZE = NXYZE - SE - 1
  422                   CONTINUE
                        J = J - NXYZF + XF
                        NXYZF = NXYZF - SF - 1
  420                CONTINUE
  412             CONTINUE
  410          CONTINUE
  402       CONTINUE
  400    CONTINUE

         RETURN
C
C
C                       ********************
C                       *  # of roots = 5  *
C                       ********************
C
C
    5    XFP = NXYZFT + 3
         DO 500 XF = 0,SHELLQ
            XFP = XFP + XF - 2
            XFMAX = XF * SHELLQ
            YFEND = SHELLQ - XF
            XEP = NXYZET + 3
            DO 502 XE = 0,SHELLP
               XEP = XEP + XE - 2
               XEMAX = XE * SHELLP
               YEEND = SHELLP - XE

               IF (MOD (XE+XF,2).EQ.1) GOTO 502

               DO M = 1,NGQEXQ
                  TEMP1 (M) = SCALE (M) * INT2DX (M,XE,XF)
               END DO

               XYFP = XFP - XFMAX
               DO 510 YF = 0,YFEND
                  XYF = XF + YF
                  XYFP = XYFP - 1
                  SFEND = MAX (SHELLC,XYF)
                  XYEP = XEP - XEMAX
                  DO 512 YE = 0,YEEND
                     XYE = XE + YE
                     XYEP = XYEP - 1
                     SEEND = MAX (SHELLA,XYE)

                     YTOT = YE + YF
                     IF (MOD (YTOT,2).EQ.1) GOTO 512

                     IF (YTOT.EQ.0) THEN
                         DO N = 1,NGQEXQ
                            TEMP2 (N) = TEMP1 (N)
                         END DO
                     ELSE
                         DO N = 1,NGQEXQ
                            TEMP2 (N) = TEMP1 (N) * INT2DY (N,YE,YF)
                         END DO
                     END IF

                     J = XYFP
                     NXYZF = NXYZQ
                     DO 520 SF = SHELLQ,SFEND,-1
                        ZF = SF - XYF
                        I = XYEP
                        NXYZE = NXYZP
                        DO 522 SE = SHELLP,SEEND,-1
                           ZE = SE - XYE
                           ZTOT = ZE + ZF

                           IF (MOD (ZTOT,2).EQ.0) THEN
                            IF (ZTOT.EQ.0) THEN
                               K = 1
                               DO M = 1,NEXQ
                                  BATCH (M,I,J) =   TEMP2 (K)
     +                                            + TEMP2 (K+1)
     +                                            + TEMP2 (K+2)
     +                                            + TEMP2 (K+3)
     +                                            + TEMP2 (K+4)
                                  K = K + 5
                               END DO
                            ELSE
                               K = 1
                               DO M = 1,NEXQ
                                  BATCH (M,I,J) =   TEMP2 (K)
     +                                            * INT2DZ (K,ZE,ZF)
     +                                            + TEMP2 (K+1)
     +                                            * INT2DZ (K+1,ZE,ZF)
     +                                            + TEMP2 (K+2)
     +                                            * INT2DZ (K+2,ZE,ZF)
     +                                            + TEMP2 (K+3)
     +                                            * INT2DZ (K+3,ZE,ZF)
     +                                            + TEMP2 (K+4)
     +                                            * INT2DZ (K+4,ZE,ZF)
                                  K = K + 5
                               END DO
                            END IF
                           END IF

                           I = I - NXYZE + XE
                           NXYZE = NXYZE - SE - 1
  522                   CONTINUE
                        J = J - NXYZF + XF
                        NXYZF = NXYZF - SF - 1
  520                CONTINUE
  512             CONTINUE
  510          CONTINUE
  502       CONTINUE
  500    CONTINUE

         RETURN
C
C
C                       ********************
C                       *  # of roots = 6  *
C                       ********************
C
C
    6    XFP = NXYZFT + 3
         DO 600 XF = 0,SHELLQ
            XFP = XFP + XF - 2
            XFMAX = XF * SHELLQ
            YFEND = SHELLQ - XF
            XEP = NXYZET + 3
            DO 602 XE = 0,SHELLP
               XEP = XEP + XE - 2
               XEMAX = XE * SHELLP
               YEEND = SHELLP - XE

               IF (MOD (XE+XF,2).EQ.1) GOTO 602

               DO M = 1,NGQEXQ
                  TEMP1 (M) = SCALE (M) * INT2DX (M,XE,XF)
               END DO

               XYFP = XFP - XFMAX
               DO 610 YF = 0,YFEND
                  XYF = XF + YF
                  XYFP = XYFP - 1
                  SFEND = MAX (SHELLC,XYF)
                  XYEP = XEP - XEMAX
                  DO 612 YE = 0,YEEND
                     XYE = XE + YE
                     XYEP = XYEP - 1
                     SEEND = MAX (SHELLA,XYE)

                     YTOT = YE + YF
                     IF (MOD (YTOT,2).EQ.1) GOTO 612

                     IF (YTOT.EQ.0) THEN
                         DO N = 1,NGQEXQ
                            TEMP2 (N) = TEMP1 (N)
                         END DO
                     ELSE
                         DO N = 1,NGQEXQ
                            TEMP2 (N) = TEMP1 (N) * INT2DY (N,YE,YF)
                         END DO
                     END IF

                     J = XYFP
                     NXYZF = NXYZQ
                     DO 620 SF = SHELLQ,SFEND,-1
                        ZF = SF - XYF
                        I = XYEP
                        NXYZE = NXYZP
                        DO 622 SE = SHELLP,SEEND,-1
                           ZE = SE - XYE
                           ZTOT = ZE + ZF

                           IF (MOD (ZTOT,2).EQ.0) THEN
                            IF (ZTOT.EQ.0) THEN
                               K = 1
                               DO M = 1,NEXQ
                                  BATCH (M,I,J) =   TEMP2 (K)
     +                                            + TEMP2 (K+1)
     +                                            + TEMP2 (K+2)
     +                                            + TEMP2 (K+3)
     +                                            + TEMP2 (K+4)
     +                                            + TEMP2 (K+5)
                                  K = K + 6
                               END DO
                            ELSE
                               K = 1
                               DO M = 1,NEXQ
                                  BATCH (M,I,J) =   TEMP2 (K)
     +                                            * INT2DZ (K,ZE,ZF)
     +                                            + TEMP2 (K+1)
     +                                            * INT2DZ (K+1,ZE,ZF)
     +                                            + TEMP2 (K+2)
     +                                            * INT2DZ (K+2,ZE,ZF)
     +                                            + TEMP2 (K+3)
     +                                            * INT2DZ (K+3,ZE,ZF)
     +                                            + TEMP2 (K+4)
     +                                            * INT2DZ (K+4,ZE,ZF)
     +                                            + TEMP2 (K+5)
     +                                            * INT2DZ (K+5,ZE,ZF)
                                  K = K + 6
                               END DO
                            END IF
                           END IF

                           I = I - NXYZE + XE
                           NXYZE = NXYZE - SE - 1
  622                   CONTINUE
                        J = J - NXYZF + XF
                        NXYZF = NXYZF - SF - 1
  620                CONTINUE
  612             CONTINUE
  610          CONTINUE
  602       CONTINUE
  600    CONTINUE

         RETURN
C
C
C                       ********************
C                       *  # of roots = 7  *
C                       ********************
C
C
    7    XFP = NXYZFT + 3
         DO 700 XF = 0,SHELLQ
            XFP = XFP + XF - 2
            XFMAX = XF * SHELLQ
            YFEND = SHELLQ - XF
            XEP = NXYZET + 3
            DO 702 XE = 0,SHELLP
               XEP = XEP + XE - 2
               XEMAX = XE * SHELLP
               YEEND = SHELLP - XE

               IF (MOD (XE+XF,2).EQ.1) GOTO 702

               DO M = 1,NGQEXQ
                  TEMP1 (M) = SCALE (M) * INT2DX (M,XE,XF)
               END DO

               XYFP = XFP - XFMAX
               DO 710 YF = 0,YFEND
                  XYF = XF + YF
                  XYFP = XYFP - 1
                  SFEND = MAX (SHELLC,XYF)
                  XYEP = XEP - XEMAX
                  DO 712 YE = 0,YEEND
                     XYE = XE + YE
                     XYEP = XYEP - 1
                     SEEND = MAX (SHELLA,XYE)

                     YTOT = YE + YF
                     IF (MOD (YTOT,2).EQ.1) GOTO 712

                     IF (YTOT.EQ.0) THEN
                         DO N = 1,NGQEXQ
                            TEMP2 (N) = TEMP1 (N)
                         END DO
                     ELSE
                         DO N = 1,NGQEXQ
                            TEMP2 (N) = TEMP1 (N) * INT2DY (N,YE,YF)
                         END DO
                     END IF

                     J = XYFP
                     NXYZF = NXYZQ
                     DO 720 SF = SHELLQ,SFEND,-1
                        ZF = SF - XYF
                        I = XYEP
                        NXYZE = NXYZP
                        DO 722 SE = SHELLP,SEEND,-1
                           ZE = SE - XYE
                           ZTOT = ZE + ZF

                           IF (MOD (ZTOT,2).EQ.0) THEN
                            IF (ZTOT.EQ.0) THEN
                               K = 1
                               DO M = 1,NEXQ
                                  BATCH (M,I,J) =   TEMP2 (K)
     +                                            + TEMP2 (K+1)
     +                                            + TEMP2 (K+2)
     +                                            + TEMP2 (K+3)
     +                                            + TEMP2 (K+4)
     +                                            + TEMP2 (K+5)
     +                                            + TEMP2 (K+6)
                                  K = K + 7
                               END DO
                            ELSE
                               K = 1
                               DO M = 1,NEXQ
                                  BATCH (M,I,J) =   TEMP2 (K)
     +                                            * INT2DZ (K,ZE,ZF)
     +                                            + TEMP2 (K+1)
     +                                            * INT2DZ (K+1,ZE,ZF)
     +                                            + TEMP2 (K+2)
     +                                            * INT2DZ (K+2,ZE,ZF)
     +                                            + TEMP2 (K+3)
     +                                            * INT2DZ (K+3,ZE,ZF)
     +                                            + TEMP2 (K+4)
     +                                            * INT2DZ (K+4,ZE,ZF)
     +                                            + TEMP2 (K+5)
     +                                            * INT2DZ (K+5,ZE,ZF)
     +                                            + TEMP2 (K+6)
     +                                            * INT2DZ (K+6,ZE,ZF)
                                  K = K + 7
                               END DO
                            END IF
                           END IF

                           I = I - NXYZE + XE
                           NXYZE = NXYZE - SE - 1
  722                   CONTINUE
                        J = J - NXYZF + XF
                        NXYZF = NXYZF - SF - 1
  720                CONTINUE
  712             CONTINUE
  710          CONTINUE
  702       CONTINUE
  700    CONTINUE

         RETURN
C
C
C                       ********************
C                       *  # of roots = 8  *
C                       ********************
C
C
    8    XFP = NXYZFT + 3
         DO 800 XF = 0,SHELLQ
            XFP = XFP + XF - 2
            XFMAX = XF * SHELLQ
            YFEND = SHELLQ - XF
            XEP = NXYZET + 3
            DO 802 XE = 0,SHELLP
               XEP = XEP + XE - 2
               XEMAX = XE * SHELLP
               YEEND = SHELLP - XE

               IF (MOD (XE+XF,2).EQ.1) GOTO 802

               DO M = 1,NGQEXQ
                  TEMP1 (M) = SCALE (M) * INT2DX (M,XE,XF)
               END DO

               XYFP = XFP - XFMAX
               DO 810 YF = 0,YFEND
                  XYF = XF + YF
                  XYFP = XYFP - 1
                  SFEND = MAX (SHELLC,XYF)
                  XYEP = XEP - XEMAX
                  DO 812 YE = 0,YEEND
                     XYE = XE + YE
                     XYEP = XYEP - 1
                     SEEND = MAX (SHELLA,XYE)

                     YTOT = YE + YF
                     IF (MOD (YTOT,2).EQ.1) GOTO 812

                     IF (YTOT.EQ.0) THEN
                         DO N = 1,NGQEXQ
                            TEMP2 (N) = TEMP1 (N)
                         END DO
                     ELSE
                         DO N = 1,NGQEXQ
                            TEMP2 (N) = TEMP1 (N) * INT2DY (N,YE,YF)
                         END DO
                     END IF

                     J = XYFP
                     NXYZF = NXYZQ
                     DO 820 SF = SHELLQ,SFEND,-1
                        ZF = SF - XYF
                        I = XYEP
                        NXYZE = NXYZP
                        DO 822 SE = SHELLP,SEEND,-1
                           ZE = SE - XYE
                           ZTOT = ZE + ZF

                           IF (MOD (ZTOT,2).EQ.0) THEN
                            IF (ZTOT.EQ.0) THEN
                               K = 1
                               DO M = 1,NEXQ
                                  BATCH (M,I,J) =   TEMP2 (K)
     +                                            + TEMP2 (K+1)
     +                                            + TEMP2 (K+2)
     +                                            + TEMP2 (K+3)
     +                                            + TEMP2 (K+4)
     +                                            + TEMP2 (K+5)
     +                                            + TEMP2 (K+6)
     +                                            + TEMP2 (K+7)
                                  K = K + 8
                               END DO
                            ELSE
                               K = 1
                               DO M = 1,NEXQ
                                  BATCH (M,I,J) =   TEMP2 (K)
     +                                            * INT2DZ (K,ZE,ZF)
     +                                            + TEMP2 (K+1)
     +                                            * INT2DZ (K+1,ZE,ZF)
     +                                            + TEMP2 (K+2)
     +                                            * INT2DZ (K+2,ZE,ZF)
     +                                            + TEMP2 (K+3)
     +                                            * INT2DZ (K+3,ZE,ZF)
     +                                            + TEMP2 (K+4)
     +                                            * INT2DZ (K+4,ZE,ZF)
     +                                            + TEMP2 (K+5)
     +                                            * INT2DZ (K+5,ZE,ZF)
     +                                            + TEMP2 (K+6)
     +                                            * INT2DZ (K+6,ZE,ZF)
     +                                            + TEMP2 (K+7)
     +                                            * INT2DZ (K+7,ZE,ZF)
                                  K = K + 8
                               END DO
                            END IF
                           END IF

                           I = I - NXYZE + XE
                           NXYZE = NXYZE - SE - 1
  822                   CONTINUE
                        J = J - NXYZF + XF
                        NXYZF = NXYZF - SF - 1
  820                CONTINUE
  812             CONTINUE
  810          CONTINUE
  802       CONTINUE
  800    CONTINUE

         RETURN
C
C
C                       ********************
C                       *  # of roots = 9  *
C                       ********************
C
C
    9    XFP = NXYZFT + 3
         DO 900 XF = 0,SHELLQ
            XFP = XFP + XF - 2
            XFMAX = XF * SHELLQ
            YFEND = SHELLQ - XF
            XEP = NXYZET + 3
            DO 902 XE = 0,SHELLP
               XEP = XEP + XE - 2
               XEMAX = XE * SHELLP
               YEEND = SHELLP - XE

               IF (MOD (XE+XF,2).EQ.1) GOTO 902

               DO M = 1,NGQEXQ
                  TEMP1 (M) = SCALE (M) * INT2DX (M,XE,XF)
               END DO

               XYFP = XFP - XFMAX
               DO 910 YF = 0,YFEND
                  XYF = XF + YF
                  XYFP = XYFP - 1
                  SFEND = MAX (SHELLC,XYF)
                  XYEP = XEP - XEMAX
                  DO 912 YE = 0,YEEND
                     XYE = XE + YE
                     XYEP = XYEP - 1
                     SEEND = MAX (SHELLA,XYE)

                     YTOT = YE + YF
                     IF (MOD (YTOT,2).EQ.1) GOTO 912

                     IF (YTOT.EQ.0) THEN
                         DO N = 1,NGQEXQ
                            TEMP2 (N) = TEMP1 (N)
                         END DO
                     ELSE
                         DO N = 1,NGQEXQ
                            TEMP2 (N) = TEMP1 (N) * INT2DY (N,YE,YF)
                         END DO
                     END IF

                     J = XYFP
                     NXYZF = NXYZQ
                     DO 920 SF = SHELLQ,SFEND,-1
                        ZF = SF - XYF
                        I = XYEP
                        NXYZE = NXYZP
                        DO 922 SE = SHELLP,SEEND,-1
                           ZE = SE - XYE
                           ZTOT = ZE + ZF

                           IF (MOD (ZTOT,2).EQ.0) THEN
                            IF (ZTOT.EQ.0) THEN
                               K = 1
                               DO M = 1,NEXQ
                                  BATCH (M,I,J) =   TEMP2 (K)
     +                                            + TEMP2 (K+1)
     +                                            + TEMP2 (K+2)
     +                                            + TEMP2 (K+3)
     +                                            + TEMP2 (K+4)
     +                                            + TEMP2 (K+5)
     +                                            + TEMP2 (K+6)
     +                                            + TEMP2 (K+7)
     +                                            + TEMP2 (K+8)
                                  K = K + 9
                               END DO
                            ELSE
                               K = 1
                               DO M = 1,NEXQ
                                  BATCH (M,I,J) =   TEMP2 (K)
     +                                            * INT2DZ (K,ZE,ZF)
     +                                            + TEMP2 (K+1)
     +                                            * INT2DZ (K+1,ZE,ZF)
     +                                            + TEMP2 (K+2)
     +                                            * INT2DZ (K+2,ZE,ZF)
     +                                            + TEMP2 (K+3)
     +                                            * INT2DZ (K+3,ZE,ZF)
     +                                            + TEMP2 (K+4)
     +                                            * INT2DZ (K+4,ZE,ZF)
     +                                            + TEMP2 (K+5)
     +                                            * INT2DZ (K+5,ZE,ZF)
     +                                            + TEMP2 (K+6)
     +                                            * INT2DZ (K+6,ZE,ZF)
     +                                            + TEMP2 (K+7)
     +                                            * INT2DZ (K+7,ZE,ZF)
     +                                            + TEMP2 (K+8)
     +                                            * INT2DZ (K+8,ZE,ZF)
                                  K = K + 9
                               END DO
                            END IF
                           END IF

                           I = I - NXYZE + XE
                           NXYZE = NXYZE - SE - 1
  922                   CONTINUE
                        J = J - NXYZF + XF
                        NXYZF = NXYZF - SF - 1
  920                CONTINUE
  912             CONTINUE
  910          CONTINUE
  902       CONTINUE
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
   10    XFP = NXYZFT + 3
         DO 1000 XF = 0,SHELLQ
            XFP = XFP + XF - 2
            XFMAX = XF * SHELLQ
            YFEND = SHELLQ - XF
            XEP = NXYZET + 3
            DO 1002 XE = 0,SHELLP
               XEP = XEP + XE - 2
               XEMAX = XE * SHELLP
               YEEND = SHELLP - XE

               IF (MOD (XE+XF,2).EQ.1) GOTO 1002

               DO M = 1,NGQEXQ
                  TEMP1 (M) = SCALE (M) * INT2DX (M,XE,XF)
               END DO
C
C
C             ...middle loops over y,y-pairs. Skip multiplication
C                of y,y-contributions, if we have a 0,0-pair, as
C                then the 2DY integral is equal to 1.
C
C
               XYFP = XFP - XFMAX
               DO 1010 YF = 0,YFEND
                  XYF = XF + YF
                  XYFP = XYFP - 1
                  SFEND = MAX (SHELLC,XYF)
                  XYEP = XEP - XEMAX
                  DO 1012 YE = 0,YEEND
                     XYE = XE + YE
                     XYEP = XYEP - 1
                     SEEND = MAX (SHELLA,XYE)

                     YTOT = YE + YF
                     IF (MOD (YTOT,2).EQ.1) GOTO 1012

                     IF (YTOT.EQ.0) THEN
                         DO N = 1,NGQEXQ
                            TEMP2 (N) = TEMP1 (N)
                         END DO
                     ELSE
                         DO N = 1,NGQEXQ
                            TEMP2 (N) = TEMP1 (N) * INT2DY (N,YE,YF)
                         END DO
                     END IF
C
C
C             ...inner loops over E,F-pairs. Skip multiplication
C                of z,z-contributions, if we have a 0,0-pair, as
C                then the 2DZ integral is equal to 1.
C
C
                     J = XYFP
                     NXYZF = NXYZQ
                     DO 1020 SF = SHELLQ,SFEND,-1
                        ZF = SF - XYF
                        I = XYEP
                        NXYZE = NXYZP
                        DO 1022 SE = SHELLP,SEEND,-1
                           ZE = SE - XYE
                           ZTOT = ZE + ZF
C
C
C             ...all info concerning all three x,x-, y,y- and z,z-pairs
C                have been collected for all exponent quadruplets at
C                once. Sum up the 2D X,Y,Z integral products to the
C                appropriate place of the [E0|F0] batch.
C
C
                           IF (MOD (ZTOT,2).EQ.0) THEN
                            IF (ZTOT.EQ.0) THEN
                               K = 0
                               DO M = 1,NEXQ
                                  SUM = ZERO
                                  DO N = 1,NGQP
                                     SUM = SUM + TEMP2 (K+N)
                                  END DO
                                  K = K + NGQP
                                  BATCH (M,I,J) = SUM
                               END DO
                            ELSE
                               K = 0
                               DO M = 1,NEXQ
                                  SUM = ZERO
                                  DO N = 1,NGQP
                                     SUM = SUM + TEMP2 (K+N)
     +                                         * INT2DZ (K+N,ZE,ZF)
                                  END DO
                                  K = K + NGQP
                                  BATCH (M,I,J) = SUM
                               END DO
                            END IF
                           END IF
C
C
C             ...next z,z-pair.
C
C
                           I = I - NXYZE + XE
                           NXYZE = NXYZE - SE - 1
 1022                   CONTINUE
                        J = J - NXYZF + XF
                        NXYZF = NXYZF - SF - 1
 1020                CONTINUE
C
C
C             ...next y,y-pair and next x,x-pair.
C
C
 1012             CONTINUE
 1010          CONTINUE
 1002       CONTINUE
 1000    CONTINUE
C
C
C             ...ready!
C
C
         RETURN
         END
