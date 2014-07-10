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
         SUBROUTINE  OED__CTR_RS_EXPAND
     +
     +                    ( N,NRS,
     +                      NR,NS,
     +                      X,
     +
     +                             Y )
     +
C------------------------------------------------------------------------
C  OPERATION   : OED__CTR_RS_EXPAND
C  MODULE      : ONE ELECTRON INTEGRALS DIRECT
C  MODULE-ID   : OED
C  SUBROUTINES : none
C  DESCRIPTION : This routine expands the rs contraction indices:
C
C                          x (n,r>=s) --> y (n,rs)
C
C
C                  Input:
C
C                      N          =  total # of invariant indices
C                      NRS        =  total # of contraction index pairs
C                                    for contraction shells rs before
C                                    expansion
C                      NR(S)      =  # of contractions for contraction
C                                    shells r(s)
C                      X          =  original integral batch with
C                                    contraction index restriction
C                                    r>=s
C
C                  Output:
C
C                      Y          =  new integral batch with complete
C                                    set of contraction indices rs
C
C
C  AUTHOR      : Norbert Flocke
C------------------------------------------------------------------------
C
C
C
C             ...include files and declare variables.
C
C
         IMPLICIT    NONE

         INTEGER     M,N,R,S
         INTEGER     NR,NS
         INTEGER     NRS
         INTEGER     RS

         DOUBLE PRECISION  X (1:N,1:NRS)
         DOUBLE PRECISION  Y (1:N,1:NR,1:NS)
C
C
C------------------------------------------------------------------------
C
C
C             ...do the expansion.
C
C
         RS = NRS + 1

         DO 100 S = NS,1,-1

            DO 200 R = NR,S+1,-1

               RS = RS - 1

               DO M = 1,N
                  Y (M,R,S) = X (M,RS)
               END DO

               DO M = 1,N
                  Y (M,S,R) = X (M,RS)
               END DO

  200       CONTINUE

            RS = RS - 1
            DO M = 1,N
               Y (M,S,S) = X (M,RS)
            END DO

  100    CONTINUE
C
C
C             ...ready!
C
C
         RETURN
         END
