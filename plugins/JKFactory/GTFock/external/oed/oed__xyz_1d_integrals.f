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
         SUBROUTINE  OED__XYZ_1D_INTEGRALS
     +
     +                    ( SHELLP,SHELLB,
     +                      ATOMIC,
     +                      NEXP,
     +                      XA,YA,ZA,XB,YB,ZB,
     +                      PAX,PAY,PAZ,
     +                      XP,YP,ZP,
     +                      PINVHF,
     +                      ABX,ABY,ABZ,
     +                      MOMENTX,MOMENTY,MOMENTZ,
     +                      XMINUS,YMINUS,ZMINUS,
     +                      SCRMTX,SCRMTY,SCRMTZ,
     +
     +                               INT1DX,
     +                               INT1DY,
     +                               INT1DZ )
     +
C------------------------------------------------------------------------
C  OPERATION   : OED__XYZ_1D_INTEGRALS
C  MODULE      : ONE ELECTRON INTEGRALS DIRECT
C  MODULE-ID   : OED
C  SUBROUTINES : OED__XYZ_INIT_00_INTS
C                OED__XYZ_INIT_OVRLP
C                OED__XYZ_MOM_INTEGRALS
C
C  DESCRIPTION : This operation calculates a full table of 1-dimensional
C                moment integrals of x,y and z type using the Rys
C                vertical recurrence scheme RI.
C
C                The vertical recurrence scheme RI is due to Rys, Dupuis
C                and King, J. Comp. Chem. 4, p.154-157 (1983).
C
C
C                  Input:
C
C                    SHELLx       =  the shell types for shells
C                                    x=P=A+B and x=B
C                    ATOMIC       =  indicates, if purely atomic 1D
C                                    integrals will be evaluated
C                    NEXP         =  # of exponent pairs
C                    PAx          =  current NEXP coordinate x=X,Y,Z
C                                    differences P-A between centers
C                                    P and A
C                    PINVHF       =  current NEXP values of 1/(2*P),
C                                    where P are the exponent sums
C                                    for contraction shells A and B
C                    ABx          =  coordinate x=X,Y,Z differences
C                                    A-B between centers A and B
C                    Xx,Yx,Zx     =  X,Y,Z coordinates of centers
C                                    x = A,B,P
C                    MOMENTx      =  tells which power of X,Y, or Z
C                                    integrals to compute
C                    SCRMTX       =  holds all previous integrals 
C                                    needed for the recursion relations
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
C  MODIFIED    : Thomas Watson Jr.                   p  q  r
C                  - Modified OED package to handle X, Y, Z integrals
C------------------------------------------------------------------------
C
C
C             ...include files and declare variables.
C
C
         IMPLICIT   NONE

         LOGICAL    ATOMIC

         INTEGER    I,J,N
         INTEGER    I1,I2,J1
         INTEGER    IEND,JEND
         INTEGER    NEXP,MOM
         INTEGER    SHELLP,SHELLB
         INTEGER    MOMENTX,MOMENTY,MOMENTZ

         DOUBLE PRECISION  XA,YA,ZA,XB,YB,ZB
         DOUBLE PRECISION  ABM,ABX,ABY,ABZ
         DOUBLE PRECISION  F
         DOUBLE PRECISION  ZERO,ONE,TWO
         DOUBLE PRECISION  P1

         DOUBLE PRECISION  XP     (1:NEXP)
         DOUBLE PRECISION  YP     (1:NEXP)
         DOUBLE PRECISION  ZP     (1:NEXP)
         DOUBLE PRECISION  PAX    (1:NEXP)
         DOUBLE PRECISION  PAY    (1:NEXP)
         DOUBLE PRECISION  PAZ    (1:NEXP)
         DOUBLE PRECISION  XMINUS (1:NEXP)
         DOUBLE PRECISION  YMINUS (1:NEXP)
         DOUBLE PRECISION  ZMINUS (1:NEXP)
         DOUBLE PRECISION  PINVHF (1:NEXP)

         DOUBLE PRECISION  INT1DX (1:NEXP,0:SHELLP,0:SHELLB)
         DOUBLE PRECISION  INT1DY (1:NEXP,0:SHELLP,0:SHELLB)
         DOUBLE PRECISION  INT1DZ (1:NEXP,0:SHELLP,0:SHELLB)

         DOUBLE PRECISION  SCRMTX (1:NEXP,0:SHELLP,0:SHELLB)
         DOUBLE PRECISION  SCRMTY (1:NEXP,0:SHELLP,0:SHELLB)
         DOUBLE PRECISION  SCRMTZ (1:NEXP,0:SHELLP,0:SHELLB)

         PARAMETER  (ZERO  =  0.D0)
         PARAMETER  (ONE   =  1.D0)
         PARAMETER  (TWO   =  2.D0)
C
C
C------------------------------------------------------------------------
C
C
C             ...initialize 1D overlap integrals.
C
C
         CALL  OED__XYZ_INIT_00_INTS
     +
     +            ( SHELLP,SHELLB,
     +              NEXP,
     +              XP,YP,ZP,
     +              PINVHF,
     +              MOMENTX,MOMENTY,MOMENTZ,
     +              XMINUS,YMINUS,ZMINUS,
     +              SCRMTX,SCRMTY,SCRMTZ,
     +
     +                       INT1DX,
     +                       INT1DY,
     +                       INT1DZ )
     +

         IF (SHELLP .EQ. 0) RETURN
C
C
C             ...Watson
C                 Zero out the following because there are
C                 not enough stringent consistency checks 
C                 to determine if a batch of zero integrals
C                 is to be expected
C
C
         DO N = 1,NEXP
         DO I = 1,SHELLP
            INT1DX (N,I,0) = ZERO
            INT1DY (N,I,0) = ZERO
            INT1DZ (N,I,0) = ZERO
            SCRMTX (N,I,0) = ZERO
            SCRMTY (N,I,0) = ZERO 
            SCRMTZ (N,I,0) = ZERO 
         END DO
         END DO
C
C
C             ...handle the X direction.
C
C
         IF (MOMENTX .EQ. 0) THEN
             CALL  OED__XYZ_INIT_OVRLP
     +
     +              ( SHELLP,0,
     +                ATOMIC,
     +                NEXP,
     +                PAX,
     +                XP,
     +                PINVHF,
     +                ABX,
     +
     +                         INT1DX )
         ELSE

             CALL  OED__XYZ_INIT_OVRLP
     +
     +               ( SHELLP,0,
     +                 ATOMIC,
     +                 NEXP,
     +                 PAX,
     +                 XP,
     +                 PINVHF,
     +                 ABX,
     +
     +                          SCRMTX )

             F = ONE
             DO 100 MOM = 1,MOMENTX

                 CALL  OED__XYZ_MOM_INTEGRALS
     +
     +                  ( SHELLP,0,
     +                    MOM,
     +                    ATOMIC,
     +                    NEXP,
     +                    PAX,
     +                    PINVHF,
     +                    ABX,XA,XP,
     +                    SCRMTX,
     +
     +                             INT1DX )
     +

                 IF (MOM .NE. MOMENTX) THEN
                     DO N = 1,NEXP
                        XMINUS (N) = SCRMTX (N,0,0)
                        DO I = 0,SHELLP
                           SCRMTX (N,I,0) = INT1DX (N,I,0)
                        END DO
                        P1 = F * PINVHF (N)
                        INT1DX (N,0,0) = XP (N) * SCRMTX (N,0,0)
     +                                     + P1 * XMINUS (N)
                     END DO
                     F = F + ONE
                 END IF

  100        CONTINUE

         END IF
C
C
C             ...handle the Y direction.
C
C
         IF (MOMENTY .EQ. 0) THEN
             CALL  OED__XYZ_INIT_OVRLP
     +
     +              ( SHELLP,0,
     +                ATOMIC,
     +                NEXP,
     +                PAY,
     +                YP,
     +                PINVHF,
     +                ABY,
     +
     +                         INT1DY )
         ELSE

             CALL  OED__XYZ_INIT_OVRLP
     +
     +               ( SHELLP,0,
     +                 ATOMIC,
     +                 NEXP,
     +                 PAY,
     +                 YP,
     +                 PINVHF,
     +                 ABY,
     +
     +                          SCRMTY )

             F = ONE
             DO 110 MOM = 1,MOMENTY

                 CALL  OED__XYZ_MOM_INTEGRALS
     +
     +                  ( SHELLP,0,
     +                    MOM,
     +                    ATOMIC,
     +                    NEXP,
     +                    PAY,
     +                    PINVHF,
     +                    ABY,YA,YP,
     +                    SCRMTY,
     +
     +                             INT1DY )
     +

                 IF (MOM .NE. MOMENTY) THEN
                     DO N = 1,NEXP
                        YMINUS (N) = SCRMTY (N,0,0)
                        DO I = 0,SHELLP
                           SCRMTY (N,I,0) = INT1DY (N,I,0)
                        END DO
                        P1 = F * PINVHF (N)
                        INT1DY (N,0,0) = YP (N) * SCRMTY (N,0,0)
     +                                     + P1 * YMINUS (N)
                     END DO
                     F = F + ONE
                 END IF

  110        CONTINUE

         END IF
C
C
C             ...handle the Z direction.
C
C
         IF (MOMENTZ .EQ. 0) THEN
             CALL  OED__XYZ_INIT_OVRLP
     +
     +              ( SHELLP,0,
     +                ATOMIC,
     +                NEXP,
     +                PAZ,
     +                ZP,
     +                PINVHF,
     +                ABZ,
     +
     +                         INT1DZ )
         ELSE

             CALL  OED__XYZ_INIT_OVRLP
     +
     +               ( SHELLP,0,
     +                 ATOMIC,
     +                 NEXP,
     +                 PAZ,
     +                 ZP,
     +                 PINVHF,
     +                 ABZ,
     +
     +                          SCRMTZ )

             F = ONE
             DO 120 MOM = 1,MOMENTZ

                 CALL  OED__XYZ_MOM_INTEGRALS
     +
     +                  ( SHELLP,0,
     +                    MOM,
     +                    ATOMIC,
     +                    NEXP,
     +                    PAZ,
     +                    PINVHF,
     +                    ABZ,ZA,ZP,
     +                    SCRMTZ,
     +
     +                             INT1DZ )
     +

                 IF (MOM .NE. MOMENTZ) THEN
                     DO N = 1,NEXP
                        ZMINUS (N) = SCRMTZ (N,0,0)
                        DO I = 0,SHELLP
                           SCRMTZ (N,I,0) = INT1DZ (N,I,0)
                        END DO
                        P1 = F * PINVHF (N)
                        INT1DZ (N,0,0) = ZP (N) * SCRMTZ (N,0,0)
     +                                     + P1 * ZMINUS (N)
                     END DO
                     F = F + ONE
                 END IF

  120        CONTINUE

         END IF
C
C
C             ...ready!
C
C
         RETURN
         END
