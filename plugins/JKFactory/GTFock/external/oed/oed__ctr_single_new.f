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
         SUBROUTINE  OED__CTR_SINGLE_NEW
     +
     +                    ( N,
     +                      MI,
     +                      NP,NC,
     +                      CC,
     +                      CCBEG,CCEND,
     +                      PRIM,
     +                      PCTR,PLOC,
     +                      X,
     +
     +                             Y )
     +
C------------------------------------------------------------------------
C  OPERATION   : OED__CTR_SINGLE_NEW
C  MODULE      : ONE ELECTRON INTEGRALS DIRECT
C  MODULE-ID   : OED
C  SUBROUTINES : none
C  DESCRIPTION : This operation performs a single contraction step
C                on the incomming integrals over primitives in one
C                shot (nonblocked) over invariant indices n and
C                generates new transformed integrals:
C
C                       y (n,r) = sum  cc (r,i) * x (n,i)
C                                  i
C
C                where cc is the array containing the contraction
C                coefficients. The sum is over the i primitives which
C                are transmitted in array PRIM and may constitute only
C                a subset of the full range of these primitives.
C
C                        --- SEGMENTED CONTRACTIONS ---
C
C                Segmented contractions are those defined to be
C                within a certain consecutive i-range of primitives.
C                The segmented limits for the contraction are sitting
C                in CCBEG (lowest limit) and CCEND (highest limit)
C                and they determine which of the actual i's from the
C                PRIM list have to be considered for each contraction
C                index.
C
C                The code is written in such form that the i-order of
C                the primitive integrals is immaterial.
C
C                The code also allows efficient contractions in case
C                there is only one contraction coefficient present
C                in a certain contraction and its value is equal to 1.
C                In such cases we can save lots of multiplications by
C                1 and since these cases are quite common for certain
C                types of basis functions it is worth including an
C                extra IF inside the contraction loop to gain speed.
C
C
C                  Input:
C
C                    N            =  # of invariant indices
C                    MI           =  # of primitives to be transformed
C                    NP           =  total # of i primitives
C                    NC           =  # of contractions for the i -> R
C                                    primitives
C                    CC           =  full set (including zeros) of
C                                    contraction coefficients for the
C                                    R contractions
C                    CCBEG        =  lowest nonzero primitive i index
C                                    for R contractions
C                    CCEND        =  highest nonzero primitive i index
C                                    for R contractions
C                    PRIM         =  primitive i indices to be
C                                    transformed
C                    Pxxx         =  intermediate storage arrays for
C                                    primitive labels to bundle
C                                    contraction steps in do loops
C                                    (xxx = CTR,LOC)
C                    X            =  input array containing the
C                                    primitive integrals
C
C                  Output:
C
C                    Y            =  contains the new transformed
C                                    integrals
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

         INTEGER     I,M,N,P,R
         INTEGER     MI,NI
         INTEGER     NC,NP
         INTEGER     NIBASE,NILEFT,NIREST,NISTEP
         INTEGER     P1,P2,P3,P4,P5,P6,P7,P8
         INTEGER     PMIN,PMAX

         INTEGER     CCBEG (1:NC)
         INTEGER     CCEND (1:NC)

         INTEGER     PRIM  (1:MI)
         INTEGER     PCTR  (1:NP)
         INTEGER     PLOC  (1:NP)

         DOUBLE PRECISION  C1,C2,C3,C4,C5,C6,C7,C8
         DOUBLE PRECISION  ZERO,ONE

         DOUBLE PRECISION  CC (1:NP,1:NC)

         DOUBLE PRECISION  X (1:N,1:MI)
         DOUBLE PRECISION  Y (1:N,1:NC)

         PARAMETER  (ZERO  =  0.D0)
         PARAMETER  (ONE   =  1.D0)
C
C
C------------------------------------------------------------------------
C
C
C             ...perform the single index contraction.
C
C
         DO 10 R = 1,NC
            PMIN = CCBEG (R)
            PMAX = CCEND (R)

            IF (PMIN.EQ.PMAX .AND. CC (PMIN,R).EQ.ONE) THEN
C
C
C             ...only one contraction coefficient equal to 1.0 is
C                present.
C
C
                DO M = 1,N
                   Y (M,R) = X (M,PMIN)
                END DO
            ELSE
C
C
C             ...the general contraction case with many coefficients.
C                Perform contractions for all primitives simultaneously.
C
C
                NI = 0
                DO 20 P = 1,MI
                   I = PRIM (P)
                   IF (I.GE.PMIN .AND. I.LE.PMAX) THEN
                       NI = NI + 1
                       PCTR (NI) = I
                       PLOC (NI) = P
                   END IF
   20           CONTINUE

                IF (NI.EQ.0) THEN
                    DO M = 1,N
                       Y (M,R) = ZERO
                    END DO
                ELSE IF (NI.EQ.1) THEN
                    P1 = PLOC (1)
                    C1 = CC (PCTR (1),R)
                    DO M = 1,N
                       Y (M,R) = C1 * X (M,P1)
                    END DO
                ELSE IF (NI.EQ.2) THEN
                    P1 = PLOC (1)
                    P2 = PLOC (2)
                    C1 = CC (PCTR (1),R)
                    C2 = CC (PCTR (2),R)
                    DO M = 1,N
                       Y (M,R) =   C1 * X (M,P1)
     +                           + C2 * X (M,P2)
                    END DO
                ELSE IF (NI.EQ.3) THEN
                    P1 = PLOC (1)
                    P2 = PLOC (2)
                    P3 = PLOC (3)
                    C1 = CC (PCTR (1),R)
                    C2 = CC (PCTR (2),R)
                    C3 = CC (PCTR (3),R)
                    DO M = 1,N
                       Y (M,R) =   C1 * X (M,P1)
     +                           + C2 * X (M,P2)
     +                           + C3 * X (M,P3)
                    END DO
                ELSE IF (NI.EQ.4) THEN
                    P1 = PLOC (1)
                    P2 = PLOC (2)
                    P3 = PLOC (3)
                    P4 = PLOC (4)
                    C1 = CC (PCTR (1),R)
                    C2 = CC (PCTR (2),R)
                    C3 = CC (PCTR (3),R)
                    C4 = CC (PCTR (4),R)
                    DO M = 1,N
                       Y (M,R) =   C1 * X (M,P1)
     +                           + C2 * X (M,P2)
     +                           + C3 * X (M,P3)
     +                           + C4 * X (M,P4)
                    END DO
                ELSE IF (NI.EQ.5) THEN
                    P1 = PLOC (1)
                    P2 = PLOC (2)
                    P3 = PLOC (3)
                    P4 = PLOC (4)
                    P5 = PLOC (5)
                    C1 = CC (PCTR (1),R)
                    C2 = CC (PCTR (2),R)
                    C3 = CC (PCTR (3),R)
                    C4 = CC (PCTR (4),R)
                    C5 = CC (PCTR (5),R)
                    DO M = 1,N
                       Y (M,R) =   C1 * X (M,P1)
     +                           + C2 * X (M,P2)
     +                           + C3 * X (M,P3)
     +                           + C4 * X (M,P4)
     +                           + C5 * X (M,P5)
                    END DO
                ELSE IF (NI.EQ.6) THEN
                    P1 = PLOC (1)
                    P2 = PLOC (2)
                    P3 = PLOC (3)
                    P4 = PLOC (4)
                    P5 = PLOC (5)
                    P6 = PLOC (6)
                    C1 = CC (PCTR (1),R)
                    C2 = CC (PCTR (2),R)
                    C3 = CC (PCTR (3),R)
                    C4 = CC (PCTR (4),R)
                    C5 = CC (PCTR (5),R)
                    C6 = CC (PCTR (6),R)
                    DO M = 1,N
                       Y (M,R) =   C1 * X (M,P1)
     +                           + C2 * X (M,P2)
     +                           + C3 * X (M,P3)
     +                           + C4 * X (M,P4)
     +                           + C5 * X (M,P5)
     +                           + C6 * X (M,P6)
                    END DO
                ELSE IF (NI.EQ.7) THEN
                    P1 = PLOC (1)
                    P2 = PLOC (2)
                    P3 = PLOC (3)
                    P4 = PLOC (4)
                    P5 = PLOC (5)
                    P6 = PLOC (6)
                    P7 = PLOC (7)
                    C1 = CC (PCTR (1),R)
                    C2 = CC (PCTR (2),R)
                    C3 = CC (PCTR (3),R)
                    C4 = CC (PCTR (4),R)
                    C5 = CC (PCTR (5),R)
                    C6 = CC (PCTR (6),R)
                    C7 = CC (PCTR (7),R)
                    DO M = 1,N
                       Y (M,R) =   C1 * X (M,P1)
     +                           + C2 * X (M,P2)
     +                           + C3 * X (M,P3)
     +                           + C4 * X (M,P4)
     +                           + C5 * X (M,P5)
     +                           + C6 * X (M,P6)
     +                           + C7 * X (M,P7)
                    END DO
                ELSE IF (NI.EQ.8) THEN
                    P1 = PLOC (1)
                    P2 = PLOC (2)
                    P3 = PLOC (3)
                    P4 = PLOC (4)
                    P5 = PLOC (5)
                    P6 = PLOC (6)
                    P7 = PLOC (7)
                    P8 = PLOC (8)
                    C1 = CC (PCTR (1),R)
                    C2 = CC (PCTR (2),R)
                    C3 = CC (PCTR (3),R)
                    C4 = CC (PCTR (4),R)
                    C5 = CC (PCTR (5),R)
                    C6 = CC (PCTR (6),R)
                    C7 = CC (PCTR (7),R)
                    C8 = CC (PCTR (8),R)
                    DO M = 1,N
                       Y (M,R) =   C1 * X (M,P1)
     +                           + C2 * X (M,P2)
     +                           + C3 * X (M,P3)
     +                           + C4 * X (M,P4)
     +                           + C5 * X (M,P5)
     +                           + C6 * X (M,P6)
     +                           + C7 * X (M,P7)
     +                           + C8 * X (M,P8)
                    END DO
                ELSE IF (NI.GT.8) THEN
                    P1 = PLOC (1)
                    P2 = PLOC (2)
                    P3 = PLOC (3)
                    P4 = PLOC (4)
                    P5 = PLOC (5)
                    P6 = PLOC (6)
                    P7 = PLOC (7)
                    P8 = PLOC (8)
                    C1 = CC (PCTR (1),R)
                    C2 = CC (PCTR (2),R)
                    C3 = CC (PCTR (3),R)
                    C4 = CC (PCTR (4),R)
                    C5 = CC (PCTR (5),R)
                    C6 = CC (PCTR (6),R)
                    C7 = CC (PCTR (7),R)
                    C8 = CC (PCTR (8),R)
                    DO M = 1,N
                       Y (M,R) =   C1 * X (M,P1)
     +                           + C2 * X (M,P2)
     +                           + C3 * X (M,P3)
     +                           + C4 * X (M,P4)
     +                           + C5 * X (M,P5)
     +                           + C6 * X (M,P6)
     +                           + C7 * X (M,P7)
     +                           + C8 * X (M,P8)
                    END DO

                    NIBASE = 8
                    NILEFT = NI - 8
                    NISTEP = NILEFT / 8
                    NIREST = MOD (NILEFT,8)

                    DO I = 1,NISTEP
                       P1 = PLOC (NIBASE+1)
                       P2 = PLOC (NIBASE+2)
                       P3 = PLOC (NIBASE+3)
                       P4 = PLOC (NIBASE+4)
                       P5 = PLOC (NIBASE+5)
                       P6 = PLOC (NIBASE+6)
                       P7 = PLOC (NIBASE+7)
                       P8 = PLOC (NIBASE+8)
                       C1 = CC (PCTR (NIBASE+1),R)
                       C2 = CC (PCTR (NIBASE+2),R)
                       C3 = CC (PCTR (NIBASE+3),R)
                       C4 = CC (PCTR (NIBASE+4),R)
                       C5 = CC (PCTR (NIBASE+5),R)
                       C6 = CC (PCTR (NIBASE+6),R)
                       C7 = CC (PCTR (NIBASE+7),R)
                       C8 = CC (PCTR (NIBASE+8),R)
                       DO M = 1,N
                          Y (M,R) = Y (M,R) + C1 * X (M,P1)
     +                                      + C2 * X (M,P2)
     +                                      + C3 * X (M,P3)
     +                                      + C4 * X (M,P4)
     +                                      + C5 * X (M,P5)
     +                                      + C6 * X (M,P6)
     +                                      + C7 * X (M,P7)
     +                                      + C8 * X (M,P8)
                       END DO
                       NIBASE = NIBASE + 8
                    END DO

                    DO I = 1,NIREST
                       P1 = PLOC (NIBASE+I)
                       C1 = CC (PCTR (NIBASE+I),R)
                       DO M = 1,N
                          Y (M,R) = Y (M,R) + C1 * X (M,P1)
                       END DO
                    END DO

                END IF

            END IF

   10    CONTINUE
C
C
C             ...ready!
C
C
         RETURN
         END
