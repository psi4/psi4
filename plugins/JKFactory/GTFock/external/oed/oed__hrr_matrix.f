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
         SUBROUTINE  OED__HRR_MATRIX
     +
     +                    ( NROTHRR,NCOLHRR,
     +                      NXYZET,NXYZA,NXYZP,
     +                      SHELLA,SHELLB,SHELLP,
     +                      NABCOOR,
     +                      ABX,ABY,ABZ,
     +                      WORK,
     +
     +                               IN1,IN2,
     +                               NROWOUT,
     +                               NROW,
     +                               ROW,
     +                               T )
     +
C------------------------------------------------------------------------
C  OPERATION   : OED__HRR_MATRIX
C  MODULE      : ONE ELECTRON INTEGRALS DIRECT
C  MODULE-ID   : OED
C  SUBROUTINES : OED__HRR_STEP
C  DESCRIPTION : This operation constructs the whole HRR transformation
C                matrix T, which will contain the information for
C                the complete transformation in xyz-basis:
C
C                             (e0| --> (ab|  ;  e = a+b ; a>=b
C
C                The matrix is constructed stepwise, starting from
C                a unit matrix and operating on its columns with
C                elementary HRR steps.
C
C
C                  Input:
C
C                    NROTHRR   = maximum # of elements of T and ROW
C                                matrix expected during construction
C                                of final T and ROW
C                    NCOLHRR   = maximum # of columns of T and ROW
C                                matrix expected during construction
C                                of final T and ROW
C                    NXYZET    = monomial dimension of the (e0| part.
C                    NXYZA     = monomial dimension for shell a
C                    NXYZP     = monomial dimension for shell p=a+b
C                    SHELLx    = shell types for x=a,b,p
C                    NABCOOR   = # of nonzero coordinate differences
C                                between nuclear centers A and B
C                    ABx       = the coordinate differences between
C                                nuclear centers A and B for x=X,Y,Z
C                    WORK      = scratch space used for assembly of
C                                final NROW vector and T and ROW
C                                matrices
C
C                  Output:
C
C                    IN1,IN2   = starting index position of final
C                                NROW vector (IN1) and final T and
C                                ROW matrices (IN2) in the big NROW
C                                T and ROW arrays (which are 2x the
C                                maximum size)
C                    NROWOUT   = maximum # of nonzero row labels
C                                of matrix T and ROW
C                    NROW      = vector containing # of nonzero
C                                entries in columns of T and ROW matrix
C                    ROW       = the nonzero row labels of matrix T
C                    T         = the HRR transformation matrix
C
C
C                The xyz-basis for the a- and b-parts in columns
C                of matrix T will be ordered such that a preceeds b.
C
C
C  AUTHOR      : Norbert Flocke
C------------------------------------------------------------------------
C
C             ...include files and declare variables.
C
C
         IMPLICIT    NONE

         INTEGER     BASE1,BASE2
         INTEGER     I,J,L,M,N
         INTEGER     IN1,IN2,OUT1,OUT2
         INTEGER     NABCOOR
         INTEGER     NCOLHRR,NROTHRR
         INTEGER     NGH,NGHO
         INTEGER     NROWIN,NROWOUT
         INTEGER     NXYZGO,NXYZHO
         INTEGER     NXYZA,NXYZP,NXYZG,NXYZH,NXYZI,NXYZET
         INTEGER     SHELLA,SHELLB,SHELLG,SHELLH,SHELLP
         INTEGER     TLEAP

         INTEGER     ADD  (0:2)
         INTEGER     NROW (1:NCOLHRR+NCOLHRR)
         INTEGER     WORK (1:NCOLHRR+NCOLHRR)

         INTEGER     ROW  (1:NROTHRR+NROTHRR)

         DOUBLE PRECISION  ABX,ABY,ABZ
         DOUBLE PRECISION  ZERO,ONE

         DOUBLE PRECISION  T (1:NROTHRR+NROTHRR)

         PARAMETER  (ZERO  =  0.D0)
         PARAMETER  (ONE   =  1.D0)

         DATA  ADD  /0,0,1/
C
C
C------------------------------------------------------------------------
C
C
C             ...accumulate T.
C
C
         IN1 = 1
         IN2 = 1

         OUT1 = IN1 + NCOLHRR
         OUT2 = IN2 + NROTHRR
C
C
C             ...form initial 'unit' T.
C
C
         DO 10 I = 1,NXYZET
            NROW (I) = 1
            ROW (I) = I
            T (I) = ONE
   10    CONTINUE
C
C
C             ...build up the HRR transformation matrix + data.
C
C
         NGH = NXYZET
         NXYZG = NXYZET
         NXYZH = 1
         NXYZGO = NXYZET
         NXYZHO = 1
         NXYZI = NXYZP
         SHELLG = SHELLP
         NROWIN = 1
         NROWOUT = 1

         DO 20 SHELLH = 1,SHELLB

            NXYZGO = NXYZGO - NXYZI
            NXYZHO = NXYZHO + SHELLH + 1
            NGHO   = NXYZGO * NXYZHO

            IF (NABCOOR.EQ.3) THEN
                M = 1 + SHELLH/3
                NROWOUT = NROWOUT + M*(M + ADD (MOD(SHELLH,3)))
            ELSE IF (NABCOOR.EQ.2) THEN
                NROWOUT = NROWOUT + SHELLH/2 + 1
            ELSE IF (NABCOOR.EQ.1) THEN
                NROWOUT = NROWOUT + 1
            END IF

            CALL  OED__HRR_STEP
     +
     +                 ( NGH,NGHO,NROWIN,NROWOUT,
     +                   NXYZA,NXYZI,NXYZG,NXYZH,NXYZGO,
     +                   SHELLA,SHELLG,SHELLH-1,
     +                   ABX,ABY,ABZ,
     +                   WORK,
     +                   NROW (IN1), ROW (IN2),
     +                   T (IN2),
     +
     +                           NROW (OUT1), ROW (OUT2),
     +                           T (OUT2) )
     +
     +
            NXYZH  = NXYZHO

            IF (SHELLH.NE.SHELLB) THEN
                NGH = NGHO
                NXYZG  = NXYZGO
                NXYZI  = NXYZI - SHELLG - 1
                SHELLG = SHELLG - 1
                NROWIN = NROWOUT
            END IF

            I = IN1
            IN1 = OUT1
            OUT1 = I
            I = IN2
            IN2 = OUT2
            OUT2 = I

   20    CONTINUE
C
C
C             ...the resulting T matrix of dimension NROWOUT x (NXYZA
C                x NXYZH) has the following structure: the columns
C                are such that they come in NXYZA identical copies
C                NXYZH times. Hence we can condense the T matrix into
C                only NROWOUT x NXYZH distinct elements. The same
C                applies to the array NROW of size NXYZA x NXYZH, which
C                also allows for condensation into only NXYZH elements.
C
C
         L = 1
         M = 0
         N = 0

         BASE1 = IN1 - 1
         BASE2 = IN2 - 1

         TLEAP = NROWOUT * NXYZA

         DO 30 J = 1,NXYZH
            NROW (BASE1+J) = NROW (BASE1+L)
            DO 40 I = 1,NROWOUT
               T (BASE2+M+I) = T (BASE2+N+I)
   40       CONTINUE
            L = L + NXYZA
            M = M + NROWOUT
            N = N + TLEAP
   30    CONTINUE
C
C
C             ...ready!
C
C
         RETURN
         END
