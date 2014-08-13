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
         SUBROUTINE  OED__PRINT_BATCH
     +
     +                    ( N1,N2,N3,N4,
     +                      UNITID,
     +                      BATCH )
     +
C------------------------------------------------------------------------
C  OPERATION   : OED__PRINT_BATCH
C  MODULE      : ONE ELECTRON INTEGRALS DIRECT
C  MODULE-ID   : OED
C  DESCRIPTION : This routine prints a batch of integrals for testing
C                purposes. Up to 4-dimensional batches are accepted.
C
C
C                  Input:
C
C                    Nx           =  individual array dimensions for
C                                    all four indices x = 1,2,3,4
C                    UNITID       =  the unit identification # to which
C                                    the batch will be printed
C                    BATCH        =  the batch of integrals
C
C
C  AUTHOR      : Norbert Flocke
C------------------------------------------------------------------------
C
C             ...include files and declare variables.
C
C
         INTEGER     N1,N2,N3,N4
         INTEGER     I,J,K,L
         INTEGER     UNITID

         DOUBLE PRECISION  BATCH  (1:N1,1:N2,1:N3,1:N4)
C
C
C------------------------------------------------------------------------
C
C
C             ...print the batch to the unit specified.
C
C
         DO 10 L = 1,N4
         DO 10 K = 1,N3
         DO 10 J = 1,N2
         DO 10 I = 1,N1
            WRITE (UNITID,9999) ' [',I,J,'|',K,L,'] = ',BATCH (I,J,K,L)
   10    CONTINUE
 9999    FORMAT (A2,2I3,A1,2I3,A4,F20.10)
C
C
C             ...ready!
C
C
         RETURN
         END
