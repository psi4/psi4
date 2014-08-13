C  Copyright (c) 1997-1999, 2003 Massachusetts Institute of Technology
C 
C  This program is free software; you can redistribute it and/or modify
C  it under the terms of the GNU General Public License as published by
C  the Free Software Foundation; either version 2 of the License, or
C  (at your option) any later version.

C  This program is distributed in the hope that it will be useful,
C  but WITHOUT ANY WARRANTY; without even the implied warranty of
C  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C  GNU General Public License for more details.

C  You should have received a copy of the GNU General Public License
C  along with this program; if not, write to the Free Software
C  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
C  USA
         SUBROUTINE  OED__PRINT_2IND_BATCH
     +
     +                    ( N1,N2,N3,
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
         INTEGER     N1,N2,N3
         INTEGER     I,J,K,L
         INTEGER     UNITID

         DOUBLE PRECISION  BATCH  (1:N1,1:N2,1:N3)
C
C
C------------------------------------------------------------------------
C
C
C             ...print the batch to the unit specified.
C
C
         DO 10 N = 1,N1
         DO 10 I = 1,N2
         DO 10 J = 1, N3
            WRITE (UNITID,9999) N,' [',I,'|',J,'] = ',BATCH (N,I,J)
   10    CONTINUE
 9999    FORMAT (I3,A2,I3,A1,I3,A4,F20.10)
C
C
C             ...ready!
C
C
         RETURN
         END
