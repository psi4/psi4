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
         SUBROUTINE  OED__XYZ_INIT_00_INTS
     +
     +                    ( SHELLP,SHELLB,
     +                      NEXP,
     +                      XP,YP,ZP,
     +                      PINVHF,
     +                      MOMENTX,MOMENTY,MOMENTZ,
     +                      XMINUS,YMINUS,ZMINUS,
     +                      SCRMTX,SCRMTY,SCRMTZ,
     +
     +                               INT1DX,
     +                               INT1DY,
     +                               INT1DZ )
     +
C------------------------------------------------------------------------
C  OPERATION   : OED__XYZ_INIT_00_INTS 
C  MODULE      : ONE ELECTRON INTEGRALS DIRECT
C  MODULE-ID   : OED                 
C  SUBROUTINES : none                                 p q r               
C  DESCRIPTION : This operation computes all the (0| X Y Z |0) moment
C                integrals
C
C
C                  Input:
C
C                    SHELLx       =  the shell type for csh x = A,B,P
C                    NEXP         =  # of exponent pairs
C                    xP           =  current NEXP coordinate x=X,Y,Z
C                                    Gaussian Product centers
C                    PINVHF       =  current NEXP values of 1/(2*P),
C                                    where P are the exponent sums
C                                    for contraction shells A and B
C                    MOMENTx      =  tells which power of X,Y, or Z
C                                    integrals to compute
C                    SCRMTx       =  holds all previous integrals 
C                                    needed for the recursion relations
C                    xMINUS       =  holds previous (0|mu-2|0) integrals
C
C
C                  Output:
C
C                    INT1Dx       =  all 1D moment integrals for each
C                                    cartesian component x = X,Y,Z
C                                    for the current exponent pairs
C
C  
C  AUTHOR      : Norbert Flocke
C                  - Wrote original OED package
C
C  MODIFIER    : Thomas Watson Jr.                   p  q  r
C                  - Modified OED package to handle X, Y, Z integrals
C------------------------------------------------------------------------
C                      
C             
C             ...include files and declare variables.
C
C        
         IMPLICIT   NONE

         INTEGER     MOM,N,I,J
         INTEGER     SHELLP,SHELLB
         INTEGER     NEXP
         INTEGER     MOMENTX,MOMENTY,MOMENTZ

         DOUBLE PRECISION  F,P1
         DOUBLE PRECISION  ZERO,ONE

         DOUBLE PRECISION  PINVHF (1:NEXP)

         DOUBLE PRECISION  XP (1:NEXP)
         DOUBLE PRECISION  YP (1:NEXP)
         DOUBLE PRECISION  ZP (1:NEXP)

         DOUBLE PRECISION  XMINUS (1:NEXP)
         DOUBLE PRECISION  YMINUS (1:NEXP)
         DOUBLE PRECISION  ZMINUS (1:NEXP)

         DOUBLE PRECISION  SCRMTX (1:NEXP,0:SHELLP,0:SHELLB)
         DOUBLE PRECISION  SCRMTY (1:NEXP,0:SHELLP,0:SHELLB) 
         DOUBLE PRECISION  SCRMTZ (1:NEXP,0:SHELLP,0:SHELLB) 

         DOUBLE PRECISION  INT1DX (1:NEXP,0:SHELLP,0:SHELLB) 
         DOUBLE PRECISION  INT1DY (1:NEXP,0:SHELLP,0:SHELLB) 
         DOUBLE PRECISION  INT1DZ (1:NEXP,0:SHELLP,0:SHELLB) 

         PARAMETER  (ZERO = 0.D0)
         PARAMETER  (ONE  = 1.D0)
C
C
C------------------------------------------------------------------------
C
C
C             ...initialize 1D moment integrals.
C
C
         DO 200 N = 1,NEXP
            INT1DX (N,0,0) = ONE
            INT1DY (N,0,0) = ONE
            INT1DZ (N,0,0) = ONE
            SCRMTX (N,0,0) = ZERO
            SCRMTY (N,0,0) = ZERO
            SCRMTZ (N,0,0) = ZERO
  200    CONTINUE

         IF (SHELLP .EQ. 0) THEN

             IF (MOMENTX .GE. 1) THEN
                 F = ZERO
                 DO 205 MOM = 1,MOMENTX
                    DO N = 1,NEXP
                       P1 = F * PINVHF (N)
                       XMINUS (N)     = SCRMTX (N,0,0)
                       SCRMTX (N,0,0) = INT1DX (N,0,0)
                       INT1DX (N,0,0) = XP (N) * SCRMTX (N,0,0)
     +                                   +  P1 * XMINUS (N)
                    END DO
                 F = F + ONE
  205            CONTINUE
             END IF

             IF (MOMENTY .GE. 1) THEN
                 F = ZERO
                 DO 215 MOM = 1,MOMENTY
                    DO N = 1,NEXP
                       P1 = F * PINVHF (N)
                       YMINUS (N)     = SCRMTY (N,0,0)
                       SCRMTY (N,0,0) = INT1DY (N,0,0)
                       INT1DY (N,0,0) = YP (N) * SCRMTY (N,0,0)
     +                                   +  P1 * YMINUS (N)
                    END DO
                 F = F + ONE
  215            CONTINUE
             END IF

             IF (MOMENTZ .GE. 1) THEN
                 F = ZERO
                 DO 225 MOM = 1,MOMENTZ
                    DO N = 1,NEXP
                       P1 = F * PINVHF (N)
                       ZMINUS (N)     = SCRMTZ (N,0,0)
                       SCRMTZ (N,0,0) = INT1DZ (N,0,0)
                       INT1DZ (N,0,0) = ZP (N) * SCRMTZ (N,0,0)
     +                                   +  P1 * ZMINUS (N)
                    END DO
                 F = F + ONE
  225            CONTINUE
             END IF
C
C
C             ...handle the cases where SHELLP > 0
C
C                For moment integrals we need the set
C                of (0|XYZ|0) integrals as well as 
C                the set (0|XYZ-1|0) of integrals.
C                The latter go into SCRMTx.
C
C
         ELSE

             IF (MOMENTX .GT. 0) THEN
                 DO 300 N = 1,NEXP
                    INT1DX (N,0,0) = XP (N)
                    SCRMTX (N,0,0) = ONE
  300            CONTINUE
             END IF

             IF (MOMENTY .GT. 0) THEN
                 DO 310 N = 1,NEXP
                    INT1DY (N,0,0) = YP (N)
                    SCRMTY (N,0,0) = ONE
  310            CONTINUE  
             END IF

             IF (MOMENTZ .GT. 0) THEN
                 DO 320 N = 1,NEXP
                    INT1DZ (N,0,0) = ZP (N)
                    SCRMTZ (N,0,0) = ONE
  320            CONTINUE  
             END IF

         END IF
C
C             ...ready!
C
C
         RETURN
         END
