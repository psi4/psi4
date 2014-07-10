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
         SUBROUTINE  OED__NAI_1D_CENDERV_INTEGRALS
     +
     +                    ( MIJ,MGQPIJ,
     +                      NGQP,
     +                      SHELLA,SHELLB,
     +                      EXP2AB,
     +                      P,C,
     +                      ORDER,
     +                      RTS,
     +                      TEMP,
     +                      INT1D,
     +
     +                               OUT1D )
     +
C------------------------------------------------------------------------
C  OPERATION   : OED__NAI_1D_CENDERV_INTEGRALS
C  MODULE      : ONE ELECTRON INTEGRALS DIRECT
C  MODULE-ID   : OED
C  SUBROUTINES : none
C  DESCRIPTION : This operation performs a single derivation step on
C                a nuclear attraction center C for input 1D nuclear
C                attraction integrals:
C
C            I''(n,a,b,C) =   r * [  a*I'(n,a-1,b,C) + b*I'(n,a,b-1,C) ]
C                         + 2pr * [(P-C)*I'(n,a,b,C) + (1-m)*I(n,a,b,C)]
C
C                and returns the result in a separate array. The
C                explanation of the symbols in the above equation is
C                as follows:
C
C                           m  =  current differentiation order.
C                    I'',I',I  =  m-th,(m-1)-th and (m-2)-th  
C                                 differentiated 1D nuclear
C                                 attraction integral. Note, that
C                                 if m = 1, i.e. when the first order
C                                 differentiated 1D integrals are
C                                 evaluated, the (m-1) term does
C                                 not exist. It only comes into play
C                                 for 2nd and higher derivatives.
C                           r  =  the set of all quadrature roots
C                           p  =  the set of all exponent sums between
C                                 shells A and B
C                         P,C  =  gaussian product center coordinate
C                                 and nuclear attraction center
C                                 coordinate
C                                 
C
C                The derivatives of the 1D integrals are calculated for
C                all roots, all nuclear centers and the present set of
C                exponent pairs.
C
C
C                  Input:
C
C                    MIJ         =  current # of ij primitive index
C                                   pairs corresponding to the
C                                   contracted shell pairs A,B
C                    MGQPIJ      =  # of roots times # of ij primitive
C                                   index pairs
C                    NGQP        =  # of gaussian quadrature points
C                                   (roots)
C                    SHELLx      =  shell types for centers x = A and B
C                    EXP2AB      =  the MIJ distinct exponent sums
C                                   between both centers x = A and B.
C                    P           =  the MIJ coordinates for the
C                                   gaussian product centers P=A+B
C                    C           =  value of the coordinate for the
C                                   nuclear attraction center
C                    ORDER       =  current order of differentiation
C                    RTS         =  all current MGQPIJ quadrature
C                                   roots
C                    TEMP        =  scratch array that will hold all
C                                   necessary 2pr*(P-C) values, where
C                                   r are all the MGQPIJ roots
C                    INT1D       =  all (ORDER-1)-th order
C                                   differentiated input 1D nuclear
C                                   attraction integrals
C                    OUT1D       =  if ORDER >= 2, this array contains
C                                   initially all (ORDER-2)-th order
C                                   differentiated input 1D nuclear
C                                   attraction integrals
C
C
C                  Output:
C
C                    OUT1D       =  all ORDER-th order differentiated
C                                   1D nuclear attraction integrals
C
C
C  AUTHOR      : Norbert Flocke
C------------------------------------------------------------------------
C
C             ...include files and declare variables.
C
C
         IMPLICIT  NONE

         INTEGER   A,B,K,M,N
         INTEGER   AM,AP,BM,BP
         INTEGER   IJ,NG
         INTEGER   MIJ,MGQPIJ
         INTEGER   NGQP
         INTEGER   ORDER
         INTEGER   SHELLA,SHELLB

         DOUBLE PRECISION  C,X,Y
         DOUBLE PRECISION  F1,F2
         DOUBLE PRECISION  ONE

         DOUBLE PRECISION  EXP2AB (1:MIJ)
         DOUBLE PRECISION  P      (1:MIJ)
         DOUBLE PRECISION  RTS    (1:MGQPIJ)
         DOUBLE PRECISION  TEMP   (1:MGQPIJ)

         DOUBLE PRECISION  INT1D (1:MGQPIJ,0:SHELLA,0:SHELLB)
         DOUBLE PRECISION  OUT1D (1:MGQPIJ,0:SHELLA,0:SHELLB)

         PARAMETER  (ONE    =  1.D0)
C
C
C------------------------------------------------------------------------
C
C
C             ...proceed.
C
C
         IF (ORDER.LT.2) THEN

             M = 0
             DO IJ = 1,MIJ
                X = (P (IJ) - C) * EXP2AB (IJ)
                DO NG = 1,NGQP
                   TEMP (M+NG) = RTS (M+NG) * X
                END DO
                M = M + NGQP
             END DO

             DO B = 0,SHELLB
             DO A = 0,SHELLA
             DO M = 1,MGQPIJ
                OUT1D (M,A,B) = TEMP (M) * INT1D (M,A,B)
             END DO
             END DO
             END DO

         ELSE

             M = 0
             X = ONE - DFLOAT (ORDER)
             DO IJ = 1,MIJ
                Y = X * EXP2AB (IJ)
                DO NG = 1,NGQP
                   TEMP (M+NG) = RTS (M+NG) * Y
                END DO
                M = M + NGQP
             END DO

             DO B = 0,SHELLB
             DO A = 0,SHELLA
             DO M = 1,MGQPIJ
                OUT1D (M,A,B) = TEMP (M) * OUT1D (M,A,B)
             END DO
             END DO
             END DO

             M = 0
             DO IJ = 1,MIJ
                X = (P (IJ) - C) * EXP2AB (IJ)
                DO NG = 1,NGQP
                   TEMP (M+NG) = RTS (M+NG) * X
                END DO
                M = M + NGQP
             END DO

             DO B = 0,SHELLB
             DO A = 0,SHELLA
             DO M = 1,MGQPIJ
                OUT1D (M,A,B) = OUT1D (M,A,B) + TEMP (M) * INT1D (M,A,B)
             END DO
             END DO
             END DO

         END IF

         IF (SHELLA.GT.0) THEN
             F1 = ONE
             DO A = 1,SHELLA
                AM = A - 1
                DO M = 1,MGQPIJ
                   OUT1D (M,A,0) = OUT1D (M,A,0)
     +                           + F1 * RTS (M) * INT1D (M,AM,0)
                END DO
                F1 = F1 + ONE
             END DO
         END IF

         IF (SHELLB.GT.0) THEN
             F2 = ONE
             DO B = 1,SHELLB
                BM = B - 1
                DO M = 1,MGQPIJ
                   OUT1D (M,0,B) = OUT1D (M,0,B)
     +                           + F2 * RTS (M) * INT1D (M,0,BM)
                END DO

                IF (SHELLA.GT.0) THEN
                    F1 = ONE
                    DO A = 1,SHELLA
                       AM = A - 1
                       DO M = 1,MGQPIJ
                          OUT1D (M,A,B) = OUT1D (M,A,B)
     +                                  + F1 * RTS (M) * INT1D (M,AM,B )
     +                                  + F2 * RTS (M) * INT1D (M,A ,BM)
                       END DO
                       F1 = F1 + ONE
                    END DO
                END IF

                F2 = F2 + ONE
             END DO
         END IF
C
C
C             ...ready!
C
C
         RETURN
         END
