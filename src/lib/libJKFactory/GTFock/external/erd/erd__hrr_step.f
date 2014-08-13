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
         SUBROUTINE  ERD__HRR_STEP
     +
     +                    ( NAB,NABO,MROWIN,MROWOUT,
     +                      NXYZX,NXYZP,NXYZA,NXYZB,NXYZAO,
     +                      SHELLX,SHELLP,SHELLB,
     +                      ABX,ABY,ABZ,
     +                      CPAIR,
     +                      NROWIN,ROWIN,
     +                      WIN,
     +
     +                                NROWOUT,ROWOUT,
     +                                WOUT )
     +
C------------------------------------------------------------------------
C  OPERATION   : ERD__HRR_STEP
C  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT
C  MODULE-ID   : ERD
C  SUBROUTINES : none
C  DESCRIPTION : This operation performs a single horizontal recurrence
C                relation (HRR) step on the input transformation matrix
C                WIN of formal dimension NXYZET x NAB, where NXYZET is
C                the dimension of the starting shell combination (e0|.
C                The columns of the matrix WIN are combined such that:
C
C                   (a(b+1i)| = ((a+1i)b| + (Ai-Bi)(ab| ; i=x,y,z
C
C                Note, that the shell symbol a stands for a range of
C                shells, starting from shell x and ending at shell p-1
C                for (a(b+1i)| and (ab| and ending at shell p for
C                ((a+1i)b|.
C
C                Since the input transformation matrix WIN contains
C                already all the HRR info from previous steps, the
C                output transformation matrix will have the complete
C                information for performing a HRR of the following
C                kind:
C
C                           (e0| --> (a(b+1i)|
C
C
C                Strategy used to perform the HRR step:
C                -------------------------------------
C
C                The HRR is split into two parts: part I) deals with
C                the (Ai-Bi)(ab| part and part II) with the ((a+1i)b|
C                part.
C
C                Part I):  In this case we have the following scheme,
C                          with a-shell range a = x to p-1:
C
C                           W = WIN (xa,ya,za,xb,yb,zb) ->
C
C                               + ABX*W to WOUT (xa,ya,za,xb+1,yb,zb)
C                               + ABY*W to WOUT (xa,ya,za,xb,yb+1,zb)
C                               + ABZ*W to WOUT (xa,ya,za,xb,yb,zb+1)
C
C                Part II): Here we have the following scheme, with
C                          a-shell range a = x+1 to p:
C
C                           W = WIN (xa,ya,za,xb,yb,zb) ->
C
C                               + W to WOUT (xa-1,ya,za,xb+1,yb,zb)
C                               + W to WOUT (xa,ya-1,za,xb,yb+1,zb)
C                               + W to WOUT (xa,ya,za-1,xb,yb,zb+1)
C
C
C                How to get from the b- to the (b+1)-shell monomials:
C                ---------------------------------------------------
C
C                To perform parts I) and II) of the HRR we need an
C                algorithm which creates a unique set of (b+1)-shell
C                monomials from those of the b-shell. The following
C                strategy is adopted:
C
C                     1) if x-exponent in b-shell monomial is > 0,
C                        add +1 to x-exponent.
C
C                     2) if x-exponent in b-shell monomial is = 0
C                        and y-exponent is > 0, add +1 to x-exponent
C                        and y-exponent.
C
C                     3) if x-exponent and y-exponent in b-shell
C                        monomial is = 0, add +1 to all exponents.
C
C                Example for b=2 --> b+1=3 case:
C
C                            200 --> 300
C                            110 --> 210
C                            101 --> 201
C                            020 --> 120,030
C                            011 --> 111,021
C                            002 --> 102,012,003
C
C
C                  Input:
C
C                    NAB          =  total # of monomials of the
C                                    input ((a+1i)b| part, i.e. #
C                                    of columns of input transformation
C                                    matrix
C                    NABO         =  total # of monomials of the
C                                    output (a(b+1i)| part, i.e. #
C                                    of columns of output transformation
C                                    matrix
C                    MROWIN       =  maximum # of nonzero row elements
C                                    per column in input transformation
C                                    matrix
C                    MROWOUT      =  maximum # of nonzero row elements
C                                    per column in output transformation
C                                    matrix
C                    NXYZX        =  total # of monomials for shell x
C                    NXYZP        =  total # of monomials for shell p
C                    NXYZA        =  sum of total # of monomials for
C                                    shell range a = x to p
C                    NXYZB        =  total # of monomials for shell b
C                    NXYZAO       =  sum of total # of monomials for
C                                    shell range a = x to p-1
C                    SHELLx       =  shell type for shells x = X,P,B
C                    ABm          =  the m=x,y,z-coordinate differences
C                                    between centers A and B
C                    CPAIR        =  integer array that will hold the
C                                    pair of column indices of the input
C                                    transformation matrix that have
C                                    to be combined to form each output
C                                    transformation matrix column
C                    NROWIN (J)   =  # of nonzero row elements of J-th
C                                    input transformation matrix column
C                    ROWIN (I,J)  =  nonzero row indices of J-th input
C                                    transformation matrix column
C                    WIN (I,J)    =  nonzero elements of the input
C                                    transformation matrix corresponding
C                                    to nonzero row index I and column J
C                  Output:
C
C                    NROWOUT (J)  =  # of nonzero row elements of J-th
C                                    output transformation matrix column
C                    ROWOUT (I,J) =  nonzero row indices of J-th output
C                                    transformation matrix column
C                    WOUT (I,J)   =  nonzero elements of the output
C                                    transformation matrix corresponding
C                                    to nonzero row index I and column J
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

         INTEGER     B,I,J,M,N
         INTEGER     C1,C2
         INTEGER     I1,I2
         INTEGER     MROWIN,MROWOUT
         INTEGER     NAB,NABO
         INTEGER     NOUT
         INTEGER     NROW1,NROW2
         INTEGER     NXBGT0
         INTEGER     NXYZAO
         INTEGER     NXYZA,NXYZB,NXYZX,NXYZP
         INTEGER     OFFYA,OFFYBO
         INTEGER     RDIFF
         INTEGER     ROW1,ROW2
         INTEGER     SHELLA,SHELLX,SHELLP,SHELLB
         INTEGER     X,Y,Z,XO,YO,ZO
         INTEGER     XA,YA

         INTEGER     NROWIN  (1:NAB)
         INTEGER     NROWOUT (1:NABO)

         INTEGER     CPAIR  (1:NABO,1:2)
         INTEGER     ROWIN  (1:MROWIN,1:NAB)
         INTEGER     ROWOUT (1:MROWOUT,1:NABO)

         DOUBLE PRECISION  ABX,ABY,ABZ
         DOUBLE PRECISION  COEFF
         DOUBLE PRECISION  ZERO

         DOUBLE PRECISION  WIN  (1:MROWIN,1:NAB)
         DOUBLE PRECISION  WOUT (1:MROWOUT,1:NABO)

         PARAMETER  (ZERO  =  0.D0)
C
C
C------------------------------------------------------------------------
C
C
C             ...determine some useful constants.
C
C
         NXBGT0 = NXYZB - SHELLB - 1
         OFFYBO = NXYZAO * (SHELLB + 1)
C
C
C             ...determine the column pairs to be 'added' to define
C                the output tranformation matrix.
C
C                      ------- Part I) section ---------
C
C                Outer loop over b- to (b+1)-shell monomials with
C                xb-exponent > 0. Inner loop over a-part contractions
C                for part I) of HRR for a-shell range a = x to p-1.
C
C
         X = - NXYZA
         XO = - NXYZAO
         DO 100 B = 1,NXBGT0
            X = X + NXYZA
            XO = XO + NXYZAO
            IF (DABS(ABX).GT.ZERO) THEN
                DO 102 J = 1,NXYZAO
                   CPAIR (XO+J,1) = X + J
  102           CONTINUE
            END IF
  100    CONTINUE
C
C
C             ...outer loop over b- to (b+1)-shell monomials with
C                xb-exponent = 0 and yb-exponent > 0. Inner loop over
C                a-part contractions for part I) of HRR for a-shell
C                range a = x to p-1.
C
C
         YO = XO + OFFYBO
         DO 110 B = 1,SHELLB
            X = X + NXYZA
            XO = XO + NXYZAO
            YO = YO + NXYZAO
            IF (DABS(ABX).GT.ZERO) THEN
                DO 112 J = 1,NXYZAO
                   CPAIR (XO+J,1) = X + J
  112           CONTINUE
            END IF
            IF (DABS(ABY).GT.ZERO) THEN
                DO 114 J = 1,NXYZAO
                   CPAIR (YO+J,1) = X + J
  114           CONTINUE
            END IF
  110    CONTINUE
C
C
C             ...last b- to (b+1)-shell monomial with xb-exponent = 0
C                and yb-exponent = 0. Inner loop over a-part
C                contractions for part I) of HRR for a-shell range
C                a = x to p-1.
C
C
         X = X + NXYZA
         XO = XO + NXYZAO
         YO = YO + NXYZAO
         ZO = YO + NXYZAO
         IF (DABS(ABX).GT.ZERO) THEN
             DO 120 J = 1,NXYZAO
                CPAIR (XO+J,1) = X + J
  120        CONTINUE
         END IF
         IF (DABS(ABY).GT.ZERO) THEN
             DO 122 J = 1,NXYZAO
                CPAIR (YO+J,1) = X + J
  122        CONTINUE
         END IF
         IF (DABS(ABZ).GT.ZERO) THEN
             DO 124 J = 1,NXYZAO
                CPAIR (ZO+J,1) = X + J
  124        CONTINUE
         END IF
C
C
C                      ------- Part II) section ---------
C
C             ...outer loop over b- to (b+1)-shell monomials with
C                xb-exponent > 0. Inner loop over a-part contractions
C                for part II) of HRR for a-shell range a = x+1 to p.
C
C
         X = 0
         XO = 0
         DO 200 B = 1,NXBGT0
            X = X + NXYZX
            DO 202 SHELLA = SHELLX,SHELLP-1
               DO 204 XA = SHELLA,0,-1
               DO 204 YA = SHELLA-XA,0,-1
                  X = X + 1
                  XO = XO + 1
                  CPAIR (XO,2) = X
  204          CONTINUE
               X = X + SHELLA + 2
  202       CONTINUE
  200    CONTINUE
C
C
C             ...outer loop over b- to (b+1)-shell monomials with
C                xb-exponent = 0 and yb-exponent > 0. Inner loop over
C                a-part contractions for part II) of HRR for a-shell
C                range a = x+1 to p.
C
C

         YO = XO + OFFYBO
         DO 210 B = 1,SHELLB
            X = X + NXYZX
            DO 212 SHELLA = SHELLX,SHELLP-1
               OFFYA = 0
               DO 214 XA = SHELLA,0,-1
                  OFFYA = OFFYA + 1
                  DO 216 YA = OFFYA-1,0,-1
                     X = X + 1
                     Y = X + OFFYA
                     XO = XO + 1
                     YO = YO + 1
                     CPAIR (XO,2) = X
                     CPAIR (YO,2) = Y
  216             CONTINUE
  214          CONTINUE
               X = Y + 1
  212       CONTINUE
  210    CONTINUE
C
C
C             ...last b- to (b+1)-shell monomial with xb-exponent = 0
C                and yb-exponent = 0. Inner loop over a-part
C                contractions for part II) of HRR for a-shell range
C                a = x+1 to p.
C
C
         X = X + NXYZX
         ZO = YO + NXYZAO
         DO 220 SHELLA = SHELLX,SHELLP-1
            OFFYA = 0
            DO 222 XA = SHELLA,0,-1
               OFFYA = OFFYA + 1
               DO 224 YA = OFFYA-1,0,-1
                  X = X + 1
                  Y = X + OFFYA
                  Z = Y + 1
                  XO = XO + 1
                  YO = YO + 1
                  ZO = ZO + 1
                  CPAIR (XO,2) = X
                  CPAIR (YO,2) = Y
                  CPAIR (ZO,2) = Z
  224          CONTINUE
  222       CONTINUE
            X = Z
  220    CONTINUE
C
C
C             ...the column pairs are ready. Construct the new
C                transformation matrix.
C
C
         DO 300 N = 1,NABO

            IF (N.LE.XO) THEN
                COEFF = ABX
            ELSE IF (N.LE.YO) THEN
                COEFF = ABY
            ELSE
                COEFF = ABZ
            END IF

            IF (DABS(COEFF).GT.ZERO) THEN

                C1 = CPAIR (N,1)
                C2 = CPAIR (N,2)
                NROW1 = NROWIN (C1)
                NROW2 = NROWIN (C2)

                M = MIN (NROW1,NROW2)
                I1 = 1
                I2 = 1
                NOUT = 0

                DO 310 I = 2,M
                   ROW1 = ROWIN (I1,C1)
                   ROW2 = ROWIN (I2,C2)
                   RDIFF = ROW1 - ROW2
                   IF (RDIFF.EQ.0) THEN
                       NOUT = NOUT + 1
                       WOUT (NOUT,N) = COEFF * WIN (I1,C1) + WIN (I2,C2)
                       ROWOUT (NOUT,N) = ROW1
                       I1 = I1 + 1
                       I2 = I2 + 1
                   ELSE IF (RDIFF.LT.0) THEN
                       NOUT = NOUT + 1
                       WOUT (NOUT,N) = COEFF * WIN (I1,C1)
                       ROWOUT (NOUT,N) = ROW1
                       I1 = I1 + 1
                   ELSE
                       NOUT = NOUT + 1
                       WOUT (NOUT,N) = WIN (I2,C2)
                       ROWOUT (NOUT,N) = ROW2
                       I2 = I2 + 1
                   END IF
  310           CONTINUE

 3000           ROW1 = ROWIN (I1,C1)
                ROW2 = ROWIN (I2,C2)
                RDIFF = ROW1 - ROW2
                IF (RDIFF.EQ.0) THEN
                    NOUT = NOUT + 1
                    WOUT (NOUT,N) = COEFF * WIN (I1,C1) + WIN (I2,C2)
                    ROWOUT (NOUT,N) = ROW1
                    I1 = I1 + 1
                    I2 = I2 + 1
                ELSE IF (RDIFF.LT.0) THEN
                    NOUT = NOUT + 1
                    WOUT (NOUT,N) = COEFF * WIN (I1,C1)
                    ROWOUT (NOUT,N) = ROW1
                    I1 = I1 + 1
                ELSE
                    NOUT = NOUT + 1
                    WOUT (NOUT,N) = WIN (I2,C2)
                    ROWOUT (NOUT,N) = ROW2
                    I2 = I2 + 1
                END IF

                IF (I1.GT.NROW1) THEN
                    DO 320 I = I2,NROW2
                       ROW2 = ROWIN (I,C2)
                       NOUT = NOUT + 1
                       WOUT (NOUT,N) = WIN (I,C2)
                       ROWOUT (NOUT,N) = ROW2
  320               CONTINUE
                ELSE IF (I2.GT.NROW2) THEN
                    DO 330 I = I1,NROW1
                       ROW1 = ROWIN (I,C1)
                       NOUT = NOUT + 1
                       WOUT (NOUT,N) = COEFF * WIN (I,C1)
                       ROWOUT (NOUT,N) = ROW1
  330               CONTINUE
                ELSE
                    GOTO 3000
                END IF

                NROWOUT (N) = NOUT

            ELSE

                C2 = CPAIR (N,2)
                NROW2 = NROWIN (C2)

                DO 340 I2 = 1,NROW2
                   WOUT (I2,N) = WIN (I2,C2)
                   ROWOUT (I2,N) = ROWIN (I2,C2)
  340           CONTINUE

                NROWOUT (N) = NROW2

            END IF

  300    CONTINUE
C
C
C             ...ready!
C
C
         RETURN
         END
