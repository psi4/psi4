C*** November 1999 version ***)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION AZ(30000000),TYPE(200),U(200,20),IA(200,6),IU(200,0:20) INT00030
      CHARACTER TYPE*5
      COMMON/MDIM/M1,M2,M3,M4,M5,M6,M7,M8,M9,M10,M11,M12,M13
 1    FORMAT(A10)                                                       INT00060
 2    FORMAT(2I5,25X,I5)                                                INT00070
 3    FORMAT(1X,'NO INPUT FOUND.')                                      INT00080
      MAXCOR=30000000                                                   INT00090
      IIN1=5                                                            INT00100
      IOUT=6                                                            INT00110
      CALL TSTART(IOUT)                                                 INT00130
      CALL LOCATE(IIN1,'# INTDER #',IERR)                               INT00140
      IF(IERR.NE.0) GO TO 1000                                          INT00150
      READ(IIN1,2) NA,NS,NDUM                                           INT00160
      NAD=NA+NDUM                                                       INT00170
      NC=3*NA                                                           INT00180
      NCD=3*NAD                                                         INT00190
      M1=1                                                              INT00200
      M2=M1+NCD                                                         INT00210
      M3=M2+NA                                                          INT00220
      M4=M3+NS                                                          INT00230
      M5=M4+NS                                                          INT00240
      M6=M5+NS*NC                                                       INT00250
      M7=M6+NC*NC                                                       INT00260
      M8=M7+NC*NC                                                       INT00270
      M9=M8+NC*NC                                                       INT00280
      M10=M9+NC*NC                                                      INT00290
      M11=M10+NC*NC                                                     INT00300
      M12=M11+3*NC
      M13=M12+3*NA
      NCORE=MAXCOR-M13+1                                                INT00320
      CALL INTDER(NA,NAD,NC,NS,
     $ AZ(M1),AZ(M2),AZ(M3),AZ(M4),AZ(M5),AZ(M6),
     $ AZ(M7),AZ(M8),AZ(M9),AZ(M10),AZ(M11),AZ(M12),AZ(M13),            INT00340
     $ TYPE,U,IA,IU,NCORE,IFLAG)                                        INT00350
      CALL TSTOP(3)                                                     INT00360
      CALL TSTOP(IOUT)                                                  INT00370
      STOP                                                              INT00380
 1000 WRITE(IOUT,3)                                                     INT00390
      STOP                                                              INT00400
      END                                                               INT00410
C     ////////////////////////////////////////////////////////////      INT14360
      SUBROUTINE INTDER(NA,NAD,NC,NS,
     $        XA,XMASS,S,SS,B,BS,A,UGF,XS,XT,XKI,XAR,
     $        Z,TYPE,U,IA,IU,NCORE,IFLAG)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER LABEL*10,TYPE*5,NUMTST*5,ASTAT*6,LREF*4
      LOGICAL FLAG,RFLAG
      INTEGER QQ
      DIMENSION IOPT(30),TYPE(NS),IA(NS,6),IU(NS,0:1),S(NS),SS(NS)
      DIMENSION Z(1),U(NS,1),IR(4),XR(4),XA(NAD,3),B(NS,NC)
      DIMENSION A(NC,NC),UGF(NC,1),BS(NC,NC),XS(NC,NC)
      DIMENSION XT(NC,NC),XMASS(1),DIP(3),XKI(NC,3),XAR(NA,3),CM(3)
CWAx
      DIMENSION JSPF(100),SPFREF(100)
      PARAMETER(RAD=57.29577951308232D0)
      PARAMETER(ZERO=0.0D0,ONE=1.0D0,TOLINV=1.0D-10)
      COMMON /IO/ IIN1,IOUT,IIN2,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25,
     $   I31,I32,I33,I35,I36,I37,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10,
     $   ISCR11,ISCR12,ISCR13,ISCR14
      COMMON /DIPOLE/DIP,QQ
      COMMON /PHYCON/ BOHR,DEBYE,HART,WAVE0,CINT
CWAx
      COMMON /SPF1/NSPF,JSPF
      COMMON /SPF2/SPFREF    
 1    FORMAT(A10)
 3    FORMAT(16I5)
 4    FORMAT(A5,5I5,A5)
 5    FORMAT(I5,4(I4,F14.10))
 6    FORMAT(20X,3F20.10)
 7    FORMAT(1X,'B MATRIX DETERMINANT = 0 .')
 8    FORMAT(1X,'AT EQUILIBRIUM ALL INTERNAL GRADIENTS = 0 .')
 9    FORMAT(1X,'NOT ENOUGH CORE'/1X,I10,' WORDS NEEDED AND ONLY',I10,
     $          ' WORDS AVAILABLE')
 10   FORMAT(3X,'PROGRAM OPTIONS'//
     $'   NA',5X,' NSYM',5X,'  NEQ',5X,' NINV',5X,'NTEST',5X,'NFREQ',
     $5X,' NVEC',5X,'NDISP'/
     $5X,'   NS',5X,' NDER',5X,' NPRT',5X,' NDUM',5X,'NGEOM',5X,'IRINT',
     $5X,'NSTOP',5X,'NMODE')
 11   FORMAT(//1X,'NUCLEAR CARTESIAN COORDINATES (BOHR)'/)
 12   FORMAT(1X,3F20.10)
 13   FORMAT(//,1X,'DEFINITION OF INTERNAL COORDINATES')
 14   FORMAT(/,1X,'SIMPLE INTERNAL COORDINATES'/)
 15   FORMAT(/,1X,'SYMMETRY INTERNAL COORDINATES'/)
 16   FORMAT('S(',I2,')=',2X,4(F13.10,' L(',I2,')',3X),
     $     5(/,8X,4(F13.10,' L(',I2,')',3X)))
 17   FORMAT('L(',I2,')=',2X,A5,4I5,A5)
 18   FORMAT(///,1X,'VALUES OF SIMPLE INTERNAL COORDINATES ',
     $      '(ANG. OR DEG.)'/)
 19   FORMAT(4(I4,F16.10))
 20   FORMAT(/,1X,'B MATRIX FOR (SYMMETRY) INTERNAL COORDINATES'/)
 21   FORMAT(6F12.6)
 22   FORMAT(//' THIS VERSION OF THE PROGRAM USES NUMERICAL '/
     $  ' BR(I,J,K,L) MATRICES TO COMPUTE CONTRIBUTIONS TO THE'/
     $  ' QUARTIC FORCE CONSTANTS ARISING FROM NONZERO GRADIENTS.')
 23   FORMAT(/' B MATRIX INVERSION UNSUCCESSFUL FOR GIVEN TOLERANCE.')
 24   FORMAT(/' B*A MATRIX'/)
 25   FORMAT(/' NA FROM INPUT IS INCONSISTENT WITH THE HEADER IN A'/
     $    ' CARTESIAN DERIVATIVE FILE.   IFLAG=',I5)
 26   FORMAT(F20.10)
 27   FORMAT(///,1X,'VALUES OF SYMMETRY INTERNAL COORDINATES ',
     $      '(ANG. OR RAD.)'/)
 28   FORMAT(3F20.10)
 29   FORMAT(/,1X,'EXECUTION STOPPED DUE TO MASS INPUT ERROR IN ',
     $      'MASSIN.')
 30   FORMAT(//,1X,'NORMAL MODE ANALYSIS IN INTERNAL COORDINATES')
 31   FORMAT(//,1X,'NORMAL MODE ANALYSIS IN CARTESIAN COORDINATES')
 32   FORMAT(/,1X,'INTDER2000'//
     $ ' A GENERAL PROGRAM DEVELOPED BY WESLEY D. ALLEN AND CO-WORKERS'/
     $ ' WHICH PERFORMS VARIOUS VIBRATIONAL ANALYSES AND HIGHER-ORDER'/
     $ ' NONLINEAR TRANSFORMATIONS AMONG FORCE FIELD REPRESENTATIONS'/) 
 33   FORMAT(//'NUCLEAR MASSES'/)
 34   FORMAT(4(I4,F12.6,4X))
 35   FORMAT(I5,F20.10)
 36   FORMAT(/,1X,'EXECUTION STOPPED DUE TO ERROR IN GFMAT.',
     $      '  IFLAG=',I5)
 37   FORMAT(/,1X,'EXECUTION STOPPED DUE TO ERROR IN NORMCO.',
     $      '  IFLAG=',I5)
 38   FORMAT(1X,  'REFERENCE VALUE IS ',F20.10,' A.')
 39   FORMAT(//,1X,   'NDER=0   NO TRANSFORMATION REQUESTED.')
 40   FORMAT(//,1X,'SCALED QUANTUM MECHANICAL FORCE CONSTANT ANALYSIS')
 41   FORMAT(/,1X,'EXECUTION STOPPED DUE TO ERROR IN SQMFC',
     $      '  IFLAG=',I5)
 42   FORMAT(2I5,5X,I5)
 43   FORMAT(/1X,'DUMMY ATOM VECTORS (BOHR)'/)
 44   FORMAT(/1X,'EXECUTION STOPPED, NOT ENOUGH CORE TO PERFORM',
     $       ' SQM ANALYSIS.'/)
 45   FORMAT(//,'EXECUTION STOPPED.  REFERENCE GEOMETRY NOT FOUND.')
 46   FORMAT(/1X,'REFERENCE GEOMETRY (BOHR)'/)
 47   FORMAT(/,' Bending S ', 3F20.10/)
      IFLAG=0
      IIN1=5
      IOUT=6
      IIN2=7
      ICHECK=3
      I11=11
      I12=12
      I15=15
      I16=16
      I17=17
      I18=18
      I20=20
      I21=21
      I24=24
      I25=25
CSA
      I31=31
      I32=32
      I33=33
      I35=35
      I36=36
      I37=37
CSA
      ISCR1=91
      ISCR2=92
      ISCR3=93
      ISCR4=94
      ISCR5=95
      ISCR6=96
      ISCR7=97
      ISCR8=98
      ISCR9=80
      ISCR10=81
      ISCR11=85
      ISCR12=86
      ISCR13=87
      ISCR14=88
CSA
      BOHR=0.529177249D0
      DEBYE=2.54176548D0
      HART=4.3597482D0
      WAVE0=1302.7910D0
      CINT=42.25472D0
CSA
      NCHUNK=1000
      ENERGY=ZERO
      CALL NOUNFL
      CALL RFILE(ISCR1)
      CALL RFILE(ISCR2)
      CALL RFILE(ISCR3)
      CALL RFILE(ISCR4)
      CALL RFILE(ISCR6)
      CALL RFILE(ISCR7)
      CALL RFILE(ISCR8)
      CALL RFILE(ISCR9)
      CALL RFILE(ISCR10)
      
      CALL OPENFF(I11,0)
      CALL OPENFF(I12,0)
      CALL OPENFF(I15,0)
      CALL OPENFF(I16,0)
      CALL OPENFF(I17,0)
      CALL OPENFF(I18,0)
      CALL OPENFF(I20,0)
      CALL OPENFF(I21,0)
      CALL OPENFF(I24,0)
      CALL OPENFF(I25,0)
      CALL OPENFF(ISCR5,0)
      CALL OPENFF(ICHECK,0)
      CALL OPENFF(IIN2,0)

      CALL LOCATE(IIN1,'# INTDER #',IERR)
COMMENT111111...................
 105  READ(IIN1,3) (IOPT(K),K=1,16)
      NA=IOPT(1)
      NS=IOPT(2)
      NSYM=IOPT(3)
      NDER=IOPT(4)
      NEQ=IOPT(5)
      NPRT=IOPT(6)
      NINV=IOPT(7)
      NDUM=IOPT(8)
      NTEST=IOPT(9)
      NGEOM=IOPT(10)
      NFREQ=IOPT(11)
      IRINT=IOPT(12)
      NVEC=IOPT(13)
      NSTOP=IOPT(14)
      NDISP=IOPT(15)
      NMODE=IOPT(16)
      NREF=0
      NCOM=0 
      IF(NVEC.NE.0) THEN
         IOPT(5)=1
         NEQ=1
      END IF
      NC=3*NA
      NAD=NA+NDUM
CWAx
      NSPF=0
CSA
      IF(NTEST.EQ.5) THEN
          ITST=2
          NTEST=0
      ELSE IF (ABS(NTEST).EQ.3) THEN
          NSTOP=1
          NINV=ISIGN(1,NTEST)
      ELSE
          ITST=1
      END IF
CSA
      IF (ABS(NINV).EQ.3) THEN
         NEQ=0
         NVEC=0
         GO TO 117
      END IF
CSA
      DO 100  J=1,NS
      IA(J,6)=0
CSA   IF(QTEST) IA(J,6)=NTEST
      IF(ABS(NTEST).LE.2.AND.NTEST.NE.0) IA(J,6)=NTEST
 100  CONTINUE
      DO 110  J=1,NS
      READ(IIN1,4) TYPE(J),(IA(J,K),K=1,5),NUMTST
      IF(TYPE(J).EQ.'STRE ') TYPE(J)=' STRE'
      IF(TYPE(J).EQ.'BEND ') TYPE(J)=' BEND'
      IF(TYPE(J).EQ.'LIN1 ') TYPE(J)=' LIN1'
      IF(TYPE(J).EQ.'LINX ') TYPE(J)=' LINX'
      IF(TYPE(J).EQ.'LINY ') TYPE(J)=' LINY'
      IF(TYPE(J).EQ.'TORS ') TYPE(J)=' TORS'
      IF(TYPE(J).EQ.' OUT '.OR.TYPE(J).EQ.'OUT  ') TYPE(J)='  OUT'
      IF(TYPE(J).EQ.' SPF '.OR.TYPE(J).EQ.'SPF  ') TYPE(J)='  SPF'
      IF(TYPE(J).EQ.'RCOM ') TYPE(J)=' RCOM'
      IF(NUMTST.EQ.'  ST '.OR.NUMTST.EQ.' ST  '.OR.NUMTST.EQ.'ST   ')
     $      NUMTST='   ST'
      IF(NTEST.NE.0.AND.NUMTST.EQ.'   ST') IA(J,6)=0
      IF(TYPE(J).EQ.'  SPF') THEN
            IF(IA(J,3).NE.0)  READ(IIN1,26) S(J)
CWAx
            NSPF=NSPF+1
            JSPF(NSPF)=J
            SPFREF(NSPF)=S(J)
      END IF
      IF(TYPE(J).EQ.' RCOM') NCOM=1
 110  CONTINUE
      IF(NSYM.EQ.0) GO TO 125
      DO 115  K=1,NSYM
 115  IU(K,0)=0
 116  READ(IIN1,5) L,(IR(K),XR(K),K=1,4)
      IF(L.EQ.0) GO TO 125
      IPOINT=IU(L,0)
      DO 120  K=1,4
      M=IR(K)
      IF(M.GT.0) THEN
         IPOINT=IPOINT+1
         IU(L,IPOINT)=M
         U(L,IPOINT)=XR(K)
      END IF
 120  CONTINUE
      IU(L,0)=IPOINT
      GO TO 116
 125  DO 121  K=1,NSYM
      M=IU(K,0)
      UNORM=ZERO
      DO 122  L=1,M
 122  UNORM=UNORM+U(K,L)*U(K,L)
      UNORM=ONE/DSQRT(UNORM)
      DO 123  L=1,M
 123  U(K,L)=U(K,L)*UNORM
 121  CONTINUE
CSA
 117  IF(NGEOM.EQ.0) THEN
CSA
 126      READ(I11,1,END=127) LABEL
          GO TO 126
 127      IK=2*NA+2
          DO 128  K=1,IK
 128      BACKSPACE I11
          READ(I11,35) NA11,ENERGY
          DO 130  J=1,NA
 130      READ(I11,6) (XA(J,K),K=1,3)
          IF(NDUM.GT.0) THEN
          DO 132  J=NA+1,NAD
 132      READ(IIN1,28) (XA(J,K),K=1,3)
          END IF
      ELSE IF (NGEOM.EQ.1) THEN
          DO 210  J=1,NAD
 210      READ(IIN1,28) (XA(J,K),K=1,3)
      END IF
C
      IF (NREF.NE.0) THEN
         READ(IIN1,FMT='(A4)') LREF 
         IF (LREF.NE.'REFG') THEN
            WRITE (IOUT,45)
            GO TO 750
         END IF
         DO 211 J=1,NA
 211     READ (IIN1,28) (XAR(J,K),K=1,3)
      END IF
CSA
      IF(ABS(NINV).EQ.3) GO TO 215
CWAx
      DO 214  J=1,NSPF
      K=JSPF(J)
      IF(IA(K,3).EQ.0) THEN
         CALL VECT1(NAD,IA(K,1),IA(K,2),XT(1,1),XA,S(K))
         S(K)=S(K)*BOHR
         SPFREF(J)=S(K)
      END IF
 214  CONTINUE
C
 215  CONTINUE
CSA
      IF((NFREQ.NE.0.AND.ABS(NFREQ).NE.5).OR.NVEC.GT.0.OR.NDISP.LT.0
     $    .OR.NCOM.GT.0.OR.NINV.LT.0) THEN
          CALL MASSIN (XMASS,NA,IFLAG)
          IF (IFLAG.NE.0) GO TO 1018
      END IF
      DO 212  I=1,NA
      IF (NVEC.GT.0.OR.NINV.LT.0) THEN
         XT(I,1) = XMASS(I)
      ELSE
         XT(I,1) = ONE
      END IF
      IF (NINV.EQ.3)    XMASS(I)= ONE
 212  CONTINUE
C
      IF ((ABS(NINV).EQ.3.OR.ABS(NTEST).EQ.3).AND.NREF.EQ.0) THEN
         CM(1)=ZERO
         CM(2)=ZERO
         CM(3)=ZERO
C
         TM=ZERO
         DO 240 J=1,NA
            TM=TM+XMASS(J)
         DO 240 K=1,3
 240        CM(K)=CM(K)+XA(J,K)*XMASS(J)
C
         DO 241 K=1,3
 241        CM(K)=CM(K)/TM
C
         DO 242 J=1,NA
         DO 242 K=1,3
 242        XAR(J,K)=XA(J,K)-CM(K)
C
      END IF
C
C     *********************************
      WRITE(IOUT,32)
      WRITE(IOUT,10)
      WRITE(IOUT,3)(IOPT(K),K=1,16)
      WRITE(IOUT,11)
      DO 133  J=1,NA
 133  WRITE(IOUT,12) (XA(J,K),K=1,3)
      IF(NDUM.GT.0) THEN
      WRITE(IOUT,43)
      DO 143  J=NA+1,NAD
 143  WRITE(IOUT,12) (XA(J,K),K=1,3)
      END IF
      DO 131  I=1,NAD
      DO 131  J=1,3
 131  XA(I,J)=XA(I,J)*BOHR
CSA
      IF (ABS(NINV).EQ.3.OR.ABS(NTEST).EQ.3) THEN
         WRITE(IOUT,46)
         DO 144 J=1,NA
            WRITE (IOUT,12) (XAR(J,K),K=1,3)
         DO 144 K=1,3
 144     XAR(J,K)=XAR(J,K)*BOHR
      END IF
CSA
      IF(ABS(NINV).EQ.3) GO TO 700
CSA
      WRITE(IOUT,13)
      WRITE(IOUT,14)
      DO 134  J=1,NS
      NUMTST='     '
      IF(NTEST.NE.0.AND.ABS(NTEST).LE.2.AND.IA(J,6).EQ.0) NUMTST='   ST'
      WRITE(IOUT,17) J,TYPE(J),(IA(J,K),K=1,4),NUMTST
      IF(TYPE(J).EQ.'  SPF') WRITE(IOUT,38) S(J)
 134  CONTINUE
      DO 136  J=1,NS
      DO 136  I=1,NC
 136  UGF(I,J)=ZERO
      IF(NSYM.EQ.0) GO TO 145
      WRITE(IOUT,15)
      DO 137  I=1,NSYM
      LL=IU(I,0)
      WRITE(IOUT,16) I,(U(I,L),IU(I,L),L=1,LL)
      DO 137  L=1,LL
      K=IU(I,L)
 137  UGF(I,K)=U(I,L)
 145  IF(LPRT(4,NPRT).LT.4) GO TO 150
      DO 141  J=1,NAD
 141  WRITE(ICHECK,142) (XA(J,K),K=1,3)
 142  FORMAT(3F10.7)
      DO 138  I=1,NSYM
 138  WRITE(ICHECK,139) (UGF(I,K),K=1,NS)
 139  FORMAT(8F10.7)
C     **********************************
COMMENT2222222222.......................
 150  CONTINUE
      IF(NDISP.NE.0) THEN
         N1=1
         N2=N1+3*NAD
         N3=N2+3*NA
         NT=N3+NS
         CALL DISP(NAD,NA,NC,NS,NDISP,IOPT,XA,XMASS,S,SS,B,BS,A,
     $     TYPE,U,IA,IU,Z(N1),Z(N2),Z(N3),UGF,XS,XT,IFLAG)
         GO TO 750
      END IF
      CALL MACHB(NAD,NC,NS,XA,XMASS,TYPE,IA,B,S)
      IF(NSYM.LE.0) THEN
         DO 154  I=1,NS
         SS(I)=S(I)
         DO 154  J=1,NC
 154     BS(I,J)=B(I,J)
      ELSE
         DO 156  J=1,NC
         DO 156  I=1,NSYM
           BS(I,J)=ZERO
           DO 158  L=1,IU(I,0)
           K=IU(I,L)
 158       BS(I,J)=BS(I,J)+U(I,L)*B(K,J)
 156     CONTINUE
         DO 160  I=1,NSYM
           SS(I)=ZERO
           DO 160  L=1,IU(I,0)
           K=IU(I,L)
 160       SS(I)=SS(I)+U(I,L)*S(K)
      END IF
      IF(NVEC.GT.0.OR.NINV.LT.0) THEN
           WRITE(IOUT,33)
           WRITE(IOUT,34) (I,XMASS(I),I=1,NA)
      END IF
      NSX=NS
      IF(NSYM.GT.0) NSX=NSYM
      NSY=NSX
      DO 162  I=1,NS
      UGF(I,1)=S(I)
      IF(TYPE(I).EQ.' BEND'.OR.TYPE(I).EQ.' LIN1'.OR.TYPE(I).EQ.
     $   ' TORS'.OR.TYPE(I).EQ.'  OUT') THEN
      UGF(I,1)=UGF(I,1)*RAD
      END IF
 162  CONTINUE
      WRITE(IOUT,18)
      WRITE(IOUT,19) (I,UGF(I,1),I=1,NS)
      IF(NSYM.GT.0) THEN
      WRITE(IOUT,27)
      WRITE(IOUT,19) (I,SS(I),I=1,NSYM)
      END IF
      IF(LPRT(1,NPRT).GE.2) THEN
         WRITE(IOUT,20)
         DO 170  I=1,NSX
         WRITE(IOUT,*) 'I= ',I
 170     WRITE(IOUT,21) (BS(I,J),J=1,NC)
      END IF
      IF(NFREQ.LT.0) GO TO 600
      IF(NINV.GT.0) THEN
           NSY=NC
CSA        IF(NTEST.EQ.0) GO TO 175
      END IF
      CALL BINVRT(NSX,NC,XT(1,1),BS,UGF,A,XS,XT(1,2),IFLAG,TOLINV,ITST)
      IF(IFLAG.NE.0) GO TO 1002
COMMENT44444444444........................
 175  N1=1
      IF(NDER.EQ.0) GO TO 1014
      IF(NDER.EQ.1.AND.NEQ.EQ.0) GO TO 1004
      IF(NDER.EQ.1) THEN
         N2=1
         N3=1
         N4=1
         N5=N1+NC
         GO TO 180
      END IF
      IF(NEQ.EQ.0) THEN
         N2=1
      ELSE
         N2=NC+1
      END IF
      IF(NDER.EQ.2) THEN
         N3=1
         N4=1
         N5=N2+NC**2
         GO TO 180
      END IF
      N3=N2+NC**2
      IF(NDER.EQ.3) THEN
         N4=1
         N5=N3+NC**3
         GO TO 180
      END IF
      IF(NDER.EQ.4) THEN
         N4=N3+NC**3
         N5=N4+NC**4
         GO TO 180
      END IF
 180  NMM=NCHUNK*3
      NT=N5+NMM
      MDER=NDER+NEQ
      IF(MDER.GE.4) THEN
           NNN=NC**3
           N6=N5+NNN
           NMM=MAX0(NMM,2*NNN)
           NT=N5+NMM
      END IF
      IF(NT.GT.NCORE) GO TO 1006
      IF(NINV.EQ.2) THEN
        CALL INPFKM(NC,NDER,NEQ,NSX,Z(N1),Z(N2),Z(N3),Z(N4))
      END IF
      IF(MDER.GE.3) THEN
         IF(ABS(NTEST).EQ.1) THEN
            CALL NUMX(NAD,NC,NS,NSX,XA,XMASS,TYPE,IA,XT(1,1),A,S,UGF,XS)
         END IF
         CALL MACHX(NAD,NC,NS,NSX,IOPT,XA,XMASS,TYPE,IA,A,S,U,IU,UGF,XS)
         IF(NTEST.EQ.1) THEN
            CALL SRTST1(NC,NS,NSX,NSYM,IA,XS,XT)
         END IF
      END IF
      IF (MDER.GE.4) THEN
         IF(ABS(NTEST).EQ.2) THEN
            CALL NUMY(NAD,NC,NS,NSX,XA,TYPE,IA,A,S,Z(N5),Z(N6),UGF)
         END IF
       CALL MACHY(NAD,NC,NS,NSX,IOPT,XA,TYPE,IA,A,S,U,IU,Z(N5),Z(N6))
         IF(NTEST.EQ.2) THEN
            CALL SRTST2(NC,NS,NSX,NSYM,IA,Z(N5),Z(N6))
         END IF
      END IF
      IF (MDER.GE.5) THEN
       CALL NUMZ(NAD,NC,NS,NSX,IOPT,XA,TYPE,IA,A,S,Z(N5),Z(N6))
       CALL MACHZ(NAD,NC,NS,NSX,IOPT,XA,TYPE,IA,A,S,U,IU,Z(N5),Z(N6))
      END IF
      IF(NSTOP.EQ.1) GO TO 700
COMMENT5555555555...................................
      NRUN=0
      IF(NVEC.GT.0) THEN
         CALL VCXKI(NA,NAD,NC,XMASS,XA,XKI)
      END IF
 200  NRUN=NRUN+1
 220  IF(NINV.LE.0) THEN
      CALL LINTR(NC,NSX,NC,IOPT,A,XS,Z(N1),Z(N2),Z(N3),Z(N4),Z(N5),
     $         NMM,NRUN,IFLAG)
      ELSE
      CALL LINTR(NC,NC,NSX,IOPT,BS,XS,Z(N1),Z(N2),Z(N3),Z(N4),Z(N5),
     $         NMM,NRUN,IFLAG)
      END IF
 230  IF(IFLAG.NE.0) GO TO 1008
      IF(NVEC.GT.0) THEN
      CALL VCDER1(NA,NC,NSX,NINV,XMASS,XKI,A,Z(N1),Z(N2),
     $        Z(N3),Z(N4),Z(N5),NRUN)
      END IF
      IF(NDER.LE.1) GO TO 500
      IF(NDER.EQ.2.AND.NEQ.EQ.0) GO TO 500
      IF(NEQ.NE.0) THEN
         CALL MACHF2(NC,NSX,NINV,Z(N1),Z(N2),UGF,Z(N5))
      END IF
      IF(NDER.LE.2) GO TO 500
C      ADD X-V(M,N) CONTRIBUTIONS TO 3RD AND PERHAPS 4TH DERIVATIVES.
      CALL XF2(NDER,NC,NSX,NINV,BS,Z(N2),Z(N3),Z(N4),Z(N5),UGF,XS,XT)
      IF(NDER.EQ.3.AND.NEQ.EQ.0) GO TO 500
C      FORM Y MATRICES AND WRITE TO DISK.
C      COMPLETE V(M,N,P) COMPUTATION.
      IF(NEQ.NE.0) THEN
           CALL MACHF3(NC,NSX,NINV,Z(N1),Z(N3),Z(N5),Z(N6))
      END IF
      IF(NDER.LE.3) GO TO 500
C      ADD X-V(M,N,P) CONTRIBUTIONS TO 4TH DERIVATIVES.
      CALL XF3(NC,NSX,NINV,BS,UGF,XS,XT,Z(N3),Z(N4),Z(N5))
C      ADD Y-V(M,N) CONTRIBUTIONS TO 4TH DERIVATIVES.
      CALL YF2(NC,NSX,NINV,BS,UGF,XS,Z(N2),Z(N4),Z(N5))
      IF(NDER.EQ.4.AND.NEQ.EQ.0) GO TO 500
      WRITE(IOUT,22)
C      FORM Z MATRIX ELEMENTS NUMERICALLY AND
C        COMPLETE THE V(M,N,P,Q) COMPUTATION.
      CALL MACHF4(NC,NSX,NINV,Z(N1),Z(N4),Z(N5),Z(N6))
 500  CALL FCOUT(NC,NS,NSY,NEQ,NDER,NINV,NVEC,NRUN,SS,ENERGY,
     $           Z(N1),Z(N2),Z(N3),Z(N4))
      IF(NVEC.EQ.1.AND.NRUN.LT.3) GO TO 200
 600  NFREQ=IABS(NFREQ)
      IF(NFREQ.NE.0.AND.NFREQ.NE.2.AND.NFREQ.NE.5.AND.NFREQ.NE.7) THEN
         WRITE(IOUT,30)
         CALL GFMAT(NA,NAD,NC,NSX,NFREQ,IRINT,NMODE,IA,
     $           XMASS,XA,S,BS,UGF,A,XS,XT,XS,IFLAG)
         IF(IFLAG.NE.0) GO TO 1010
      END IF
      IF (NFREQ.EQ.2.OR.NFREQ.EQ.3.OR.NFREQ.EQ.7.OR.NFREQ.EQ.8) THEN
         WRITE(IOUT,31)
         CALL NORMCO(NA,NAD,NC,NFREQ,IRINT,
     $            XMASS,XA,BS,UGF,A,XS,XT,XS,IFLAG)
         IF(IFLAG.NE.0) GO TO 1012
      END IF
      IF (NFREQ.EQ.5) THEN
         WRITE(IOUT,40)
         CALL LOCATE(IIN1,'# SQMFC ##',IERR)
         READ(IIN1,42) NSF,NISO,NH
         N1=1
         N2=N1+NA*NISO
         N3=N2+NSX*NISO
         N4=N3+NSX*NISO
         N5=N4+NSX*NSX*NISO
         N6=N5+NSX*NSX
         N7=N6+NSF
         N8=N7+NSF
         N9=N8+NSF
         N10=N9+NSF
         N11=N10+NSF*NSF
         N12=N11+NSF*NSF
         IF(ABS(NH).GE.2) THEN
           N13=N12+NSX*NSX*NSF
           NT=N13+NSX*NSF*NISO
         ELSE
           N13=N12
           NT=N12
         END IF
         IF(NT.GT.NCORE) GO TO 1020
         CALL SQMFC(NA,NAD,NC,NSX,NSF,NISO,NMODE,XA,BS,
     $   Z(N1),Z(N2),Z(N3),Z(N4),Z(N5),Z(N6),Z(N7),Z(N8),Z(N9),Z(N10),
     $   S,SS,XKI(1,1),XKI(1,2),XKI(1,3),A,XS,UGF,XT,UGF,Z(N11),
     $   Z(N12),Z(N13),IFLAG)
         IF(IFLAG.NE.0) GO TO 1016
      END IF
CSA
      GO TO 750
 700  CONTINUE
      IF (ABS(NTEST).EQ.3) THEN
         N2=NC*NC+1
         N3=N2+NC**3
         IF (N3.GT.NCORE) GO TO 1006
         CALL ORTHOG(NA,NC,NSX,IOPT,XAR,BS,Z(1),Z(N2),FLAG)
         IF (FLAG) GOTO 750
      END IF
      IF (NSTOP.EQ.1) GO TO 750
C
C IF ABS(NINV).EQ.3 CALCULATE THE PROJECTION 
C MATRIX TO SEPARATE THE ROTATION
C AND THE TRANSLATION FROM THE FORCE CONSTANT MATRICES CALCULATED AT
C A NON-STATIONARY REFERENCE GEOMETRY
C
C ALLOCATE MEMORY FOR PROJECTION ROUTINES
C NDER=2
C       DK1     6*NC            N1
C       P2      NC,NC           N2
C       F2P     NC,NC           N3
C NDER=3
C       DK2     3*NC*NC         N4
C       P3      NC*NC*NC        N5
C       F3P     NC*NC*NC        N6
C       SCRTCH  NC**3           N6A
C       SCRTCH  NC**2           N6B
C NDER=4
C       DK3      3*NC*NC*NC     N7
C       P4      NC*NC*NC*NC     N8
C       F4P     NC*NC*NC*NC     N9
C       SCRATCH NC**4           N10
      NMM=NCHUNK*3
      IF(NDER.EQ.2) THEN
        N1  = 1
        N2  = N1  + 6*NC
        N3  = N2  + NC*NC
        N4  = N3  + NC*NC
        N5  = N4  + NC*NC
        N6  = N5
        N6A = N6
        N6B = N6A
        N7  = N6B
        N8  = N7
        N9  = N8
        N10 = N9
        N11 = N10
      END IF
      IF(NDER.EQ.3) THEN
        N1  = 1
        N2  = N1  + 6*NC
        N3  = N2  + NC*NC
        N4  = N3  + NC*NC
        N5  = N4  + 3*NC*NC
        N6  = N5  + NC*NC*NC
        N6A = N6  + NC*NC*NC
        N6B = N6A + NC*NC*NC
        N7  = N6B + NC**2
        N8  = N7
        N9  = N8
        N10 = N9
        N11 = N10
      END IF
      IF(NDER.EQ.4) THEN
        N1  = 1
        N2  = N1  + 6*NC
        N3  = N2  + NC*NC
        N4  = N3  + NC*NC
        N5  = N4  + 3*NC*NC
        N6  = N5  + NC*NC*NC
        N6A = N6  + NC*NC*NC
        N6B = N6A + NC*NC*NC
        N7  = N6B + NC**2
        N8  = N7  + 3*NC**3
        N9  = N8  + NC**4
        N10 = N9  + NC**4
        N11 = N10 + NC**4
      END IF
CSA
      NT=N11+NMM
      IF(NT.GT.NCORE) GO TO 1006
CSA
CSA     CALL ROTIFLIN (NA,NC,XAR,XMASS,IOPT,
CSA  $                  Z(N1),Z(N3),Z(N6),Z(N9),RFLAG)
CSA     IF (NSTOP.EQ.4) GOTO 750
CSA
        RFLAG=.FALSE.
        CALL FORMDK0 (NA,NC,XAR,XMASS,RFLAG,IOPT,Z(N5),
     $                     Z(N6),Z(N8),Z(N9),
     $                     FLAG,Z(N1),Z(N4),Z(N7))
        IF (FLAG.OR.NSTOP.EQ.2) GOTO 750
        CALL FORMP (IOPT,XMASS,NA,NC,XAR,Z(N1),Z(N4),Z(N7),Z(N6),
     $              Z(N6A),Z(N6B),Z(N9),Z(N2),Z(N5),Z(N8))
        IF (NSTOP.EQ.3) GOTO 750
CSA     IF (ABS(NDISP).EQ.3) THEN
CSA     CALL PRJDISP(NC,IOPT,RFLAG,Z(N2),Z(N5),Z(N8),
CSA  $                  Z(N1), Z(N3), FLAG, Z(N4))
CSA     GOTO 750
CSA     END IF
        CALL LINTR (NC,NC,NC,IOPT,Z(N2),Z(N4),Z(N1),Z(N3),Z(N6),
     $              Z(N9),Z(N10),3000,1,IFLAG)
        CALL PROJV (NC,IOPT,Z(N2),Z(N5),Z(N8),Z(N1),Z(N4),Z(N6A),
     $              Z(N7),Z(N10),Z(N3),Z(N6),Z(N9))
CSA     CALL UNROTATE (NC,IOPT,Z(N3),Z(N6),Z(N9))
        CALL FCOUT (NC,NC,NC,0,NDER,3,0,1,Z(N1),ENERGY,Z(N1),
     $              Z(N3),Z(N6),Z(N9))
 720  CONTINUE
  750 CONTINUE
      ISTAT=4
      IF (NSTOP.EQ.1) ISTAT=3
      CALL RCLOSE(ISCR1, ISTAT)
      CALL RCLOSE(ISCR2, ISTAT)
      CALL RCLOSE(ISCR3, ISTAT)
      CALL RCLOSE(ISCR4, ISTAT)
      CALL RCLOSE(ISCR6, ISTAT)
      CALL RCLOSE(ISCR7, ISTAT)
      CALL RCLOSE(ISCR8, ISTAT)
      CALL RCLOSE(ISCR9, ISTAT)
      CALL RCLOSE(ISCR10,ISTAT)
      CALL RCLOSE(ISCR11,ISTAT)
      CALL RCLOSE(ISCR12,ISTAT)
      CALL RCLOSE(ISCR13,ISTAT)
      CALL RCLOSE(ISCR14,ISTAT)
CWA 4/28/99
      CALL NCLOSE(ICHECK)
      CALL NCLOSE(IIN2)
      CALL NCLOSE(I11)
      CALL NCLOSE(I12)
      CALL NCLOSE(I15)
      CALL NCLOSE(I16)
      CALL NCLOSE(I17)
      CALL NCLOSE(I18)
      CALL NCLOSE(I20)
      CALL NCLOSE(I21)
      CALL NCLOSE(I24)
      CALL NCLOSE(I25)
CSA
      ASTAT='DELETE'
      CLOSE(ISCR5,STATUS=ASTAT)
      IF (NSTOP.EQ.2) ASTAT='KEEP'
         CLOSE(I31,STATUS=ASTAT)
         CLOSE(I32,STATUS=ASTAT)
         CLOSE(I33,STATUS=ASTAT)
      ASTAT='DELETE'
      IF (NSTOP.EQ.3) ASTAT='KEEP'
         CLOSE(I35,STATUS=ASTAT)
         CLOSE(I36,STATUS=ASTAT)
         CLOSE(I37,STATUS=ASTAT)
      ASTAT='DELETE'
      IF (NSTOP.EQ.4) ASTAT='KEEP'
         CLOSE(ISCR11,STATUS=ASTAT)
         CLOSE(ISCR12,STATUS=ASTAT)
         CLOSE(ISCR13,STATUS=ASTAT)
         CLOSE(ISCR14,STATUS=ASTAT)
      RETURN
CSA
 1002 IF(IFLAG.EQ.1) THEN
         WRITE(IOUT,7)
      ELSE IF(IFLAG.EQ.2) THEN
         WRITE(IOUT,23)
         WRITE(IOUT,24)
         DO 1003  I=1,NSX
 1003    WRITE(IOUT,21) (UGF(I,J),J=1,NSX)
      END IF
      RETURN
 1004 WRITE(IOUT,8)
      RETURN
 1006 WRITE(IOUT,9)  NT,NCORE
      RETURN
 1008 WRITE(IOUT,25) IFLAG
      RETURN
 1010 WRITE(IOUT,36) IFLAG
      RETURN
 1012 WRITE(IOUT,37) IFLAG
      RETURN
 1014 WRITE(IOUT,39)
      RETURN
 1016 WRITE(IOUT,41) IFLAG
      RETURN
 1018 WRITE(IOUT,29)
      RETURN
 1020 WRITE(IOUT,44)
      RETURN
 1022 WRITE(IOUT,46)
      RETURN
      END
C     ///////////////////////////////////////////////////////////////   INT23800
C STRE    (OPT2) (NUMB)                                                 INT23810
      SUBROUTINE VECT1(NAD,K1,K2,V1,XA,W)                               INT23820
      IMPLICIT REAL*8 (A-H,O-Z)                                         INT23830
      DIMENSION XA(NAD,3),V1(3)                                         INT23840
      PARAMETER(ZERO=0.0D0,ONE=1.0D0)                                   INT23850
      W=ZERO                                                            INT23860
      DO 5  I=1,3                                                       INT23870
      V1(I)=XA(K1,I)-XA(K2,I)                                           INT23880
 5    W=W+V1(I)*V1(I)                                                   INT23890
      W=DSQRT(W)                                                        INT23900
      C=ONE/W                                                           INT23910
      DO 10  I=1,3                                                      INT23920
 10   V1(I)=V1(I)*C                                                     INT23930
      RETURN                                                            INT23940
      END                                                               INT23950
C     //////////////////////////////////////////////////////////////    INT23960
C BEND    (OPT2) (NUMB)                                                 INT23970
      SUBROUTINE VECT2(NAD,K1,K2,K3,V1,V2,V3,XA,W)                      INT23980
      IMPLICIT REAL*8 (A-H,O-Z)                                         INT23990
      DIMENSION XA(NAD,3),V1(3),V2(3),V3(3),E21(3),E23(3)               INT24000
      PARAMETER(ONE=1.0D0)                                              INT24010
      CALL VECT1(NAD,K1,K2,E21,XA,T12)                                  INT24020
      CALL VECT1(NAD,K3,K2,E23,XA,T32)                                  INT24030
      CALL SCAPRO(E21,E23,W)                                            INT24040
      CP=W                                                              INT24050
      SP=DSQRT(ONE-W*W)                                                 INT24060
      W=DACOS(W)                                                        INT24070
      C1=ONE/(T12*SP)                                                   INT24080
      C2=ONE/(T32*SP)                                                   INT24090
      DO 5  I=1,3                                                       INT24100
      V1(I)=(CP*E21(I)-E23(I))*C1                                       INT24110
      V3(I)=(CP*E23(I)-E21(I))*C2                                       INT24120
  5   V2(I)=-V1(I)-V3(I)                                                INT24130
      WRITE(6, 47) V3(1) V3(2) V3(3)
      RETURN                                                            INT24140
      END                                                               INT24150
C     //////////////////////////////////////////////////////////////    INT24160
C LIN1   (OPT2) (NUMB)                                                  INT24170
      SUBROUTINE VECT3(NAD,K1,K2,K3,K4,V1,V2,V3,XA,W)                   INT24180
      IMPLICIT REAL*8 (A-H,O-Z)                                         INT24190
      REAL*8 V1(3),V2(3),V3(3),XA(NAD,3),E21(3),E23(3)                  INT24200
      DIMENSION E2M(3),EA(3),V4(3),V5(3)                                INT24210
      PARAMETER(ONE=1.0D0)                                              INT24220
      CALL VECT1(NAD,K1,K2,E21,XA,T21)                                  INT24230
      CALL VECT1(NAD,K3,K2,E23,XA,T23)                                  INT24240
      DO 5  J=1,3                                                       INT24250
  5   EA(J)=XA(K4,J)                                                    INT24260
      CALL SCAPRO(EA,EA,D)                                              INT24270
      D=ONE/DSQRT(D)                                                    INT24280
      DO 10  I=1,3                                                      INT24290
 10   EA(I)=D*EA(I)                                                     INT24300
      CALL VECPRO(E23,E21,E2M)                                          INT24310
      CALL SCAPRO(EA,E2M,STHETA)                                        INT24320
      W=DASIN(STHETA)                                                   INT24330
      CTHETA=DCOS(W)                                                    INT24340
      TTHETA=STHETA/CTHETA                                              INT24350
      CALL VECPRO(EA,E23,V4)                                            INT24360
      CALL VECPRO(EA,E21,V5)                                            INT24370
      C1=ONE/(CTHETA*T21)                                               INT24380
      C2=TTHETA/T21                                                     INT24390
      C3=ONE/(CTHETA*T23)                                               INT24400
      C4=TTHETA/T23                                                     INT24410
      DO 15  I=1,3                                                      INT24420
      V1(I)=C1*V4(I)-C2*E21(I)                                          INT24430
      V3(I)=-(C3*V5(I)+C4*E23(I))                                       INT24440
 15   V2(I)=-(V1(I)+V3(I))                                              INT24450
      RETURN                                                            INT24460
      END                                                               INT24470
C     ////////////////////////////////////////////////////////////      INT24480
C LIN2                                                                  INT24490
      SUBROUTINE VECT4(NAD,K1,K2,K3,K4,V1,V2,V3,XA,W)                   INT24500
      IMPLICIT REAL*8 (A-H,O-Z)                                         INT24510
      REAL*8 V1(3),V2(3),V3(3),XA(NAD,3),E21(3),E23(3)                  INT24520
      DIMENSION E2M(3),EA(3),V4(3),V5(3)                                INT24530
      PARAMETER(ONE=1.0D0)                                              INT24540
      CALL VECT1(NAD,K1,K2,E21,XA,T21)                                  INT24550
      CALL VECT1(NAD,K3,K2,E23,XA,T23)                                  INT24560
      CALL VECT1(NAD,K4,K2,E2M,XA,W)                                    INT24570
      CALL SCAPRO(E21,E2M,W1)                                           INT24580
      W2=DSQRT(ONE-W1*W1)                                               INT24590
      DO 5  I=1,3                                                       INT24600
  5   EA(I)=(W1*E21(I)-E2M(I))/W2                                       INT24610
      CALL VECPRO(E23,E21,E2M)                                          INT24620
      CALL SCAPRO(EA,E2M,STHETA)                                        INT24630
      W=DASIN(STHETA)                                                   INT24640
      CTHETA=DCOS(W)                                                    INT24650
      TTHETA=STHETA/CTHETA                                              INT24660
      CALL VECPRO(EA,E23,V4)                                            INT24670
      CALL VECPRO(EA,E21,V5)                                            INT24680
      DO 15  I=1,3                                                      INT24690
      V1(I)=(V4(I)/CTHETA-E21(I)*TTHETA)/T21                            INT24700
      V3(I)=-(V5(I)/CTHETA+E23(I)*TTHETA)/T23                           INT24710
 15   V2(I)=-(V1(I)+V3(I))                                              INT24720
      RETURN                                                            INT24730
      END                                                               INT24740
C     ///////////////////////////////////////////////////////////////   INT24750
C OUT     (OPT2)  (NUMB)                                                INT24760
      SUBROUTINE VECT5(NAD,K1,K2,K3,K4,V1,V2,V3,V4,XA,W)                INT24770
      IMPLICIT REAL*8 (A-H,O-Z)                                         INT24780
      DIMENSION V1(3),V2(3),V3(3),V4(3),V5(3),V6(3),V7(3)               INT24790
      DIMENSION XA(NAD,3),E23(3),E21(3),E24(3),B3P(3),B4P(3)            INT24800
      PARAMETER(ZERO=0.0D0,ONE=1.0D0,
     $PI=3.14159265358979323846264D0)
      CALL VECT1(NAD,K1,K2,E21,XA,T21)                                  INT24820
      CALL VECT1(NAD,K3,K2,E23,XA,T23)                                  INT24830
      CALL VECT1(NAD,K4,K2,E24,XA,T24)                                  INT24840
      CALL VECPRO(E23,E24,V5)                                           INT24850
      CALL SCAPRO(E21,E23,W1)
      CALL SCAPRO(E21,E24,W2)
      CALL VECPRO(E24,E21,V6)                                           INT24860
      CALL VECPRO(E21,E23,V7)                                           INT24870
      CALL VECT2(NAD,K3,K2,K4,B3P,V2,B4P,XA,PHI)                        INT24880
      SPHI=DSIN(PHI)                                                    INT24890
      CALL SCAPRO(E21,V5,W)                                             INT24900
      W=DASIN(W/SPHI)                                                   INT24910
      IF(W1+W2.GT.0) W=DSIGN(PI,W)-W
      CG=DCOS(W)                                                        INT24920
      SG=DSIN(W)                                                        INT24930
      TG=SG/CG                                                          INT24940
      W1=CG*SPHI                                                        INT24950
      W2=ONE/(T21*W1)                                                   INT24960
      W3=TG/T21                                                         INT24970
      W4=ONE/(T23*W1)                                                   INT24980
      W5=T24*SG*W4                                                      INT24990
      W6=ONE/(T24*W1)                                                   INT25000
      W7=T23*SG*W6                                                      INT25010
      DO 10  I=1,3                                                      INT25020
      V1(I)=V5(I)*W2-E21(I)*W3                                          INT25030
      V3(I)=V6(I)*W4+B4P(I)*W5                                          INT25040
      V4(I)=V7(I)*W6+B3P(I)*W7                                          INT25050
 10   V2(I)=-V1(I)-V3(I)-V4(I)                                          INT25060
      RETURN                                                            INT25070
      END                                                               INT25080
C     //////////////////////////////////////////////////////////////    INT25090
C TORS    (OPT2) (NUMB)                                                 INT25100
      SUBROUTINE VECT6(NAD,K1,K2,K3,K4,V1,V2,V3,V4,XA,W)                INT25110
      IMPLICIT REAL*8 (A-H,O-Z)                                         INT25120
      DIMENSION XA(NAD,3),V1(3),V2(3),V3(3),V4(3),V5(3),V6(3)           INT25130
      DIMENSION E21(3),E32(3),E43(3)                                    INT25140
      PARAMETER(ZERO=0.0D0,ONE=1.0D0,
     $PI=3.14159265358979323846264D0)
      CALL VECT1(NAD,K1,K2,E21,XA,T21)                                  INT25160
      CALL VECT1(NAD,K2,K3,E32,XA,T32)                                  INT25170
      CALL VECT1(NAD,K3,K4,E43,XA,T43)                                  INT25180
      CALL VECPRO(E21,E32,V5)                                           INT25190
      CALL VECPRO(E43,E32,V6)                                           INT25200
      CALL SCAPRO(E21,E32,W2)                                           INT25210
      CALL SCAPRO(E32,E43,W3)                                           INT25220
      CP2=-W2                                                           INT25230
      CP3=-W3                                                           INT25240
      SP2=DSQRT(ONE-CP2*CP2)                                            INT25250
      SP3=DSQRT(ONE-CP3*CP3)                                            INT25260
      CALL SCAPRO(E21,V6,W2)                                            INT25270
      CALL SCAPRO(V5,V6,W3)                                             INT25280
      W3=-W3                                                            INT25290
      W=DASIN(W2/(SP2*SP3))                                             INT25300
      IF(W3.LT.ZERO) W=PI-W                                             INT25310
      W1=ONE/(T21*SP2*SP2)                                              INT25320
      W2=ONE/(T43*SP3*SP3)                                              INT25330
      DO 5  I=1,3                                                       INT25340
      V1(I)=-W1*V5(I)                                                   INT25350
 5    V4(I)=-W2*V6(I)                                                   INT25360
      W3=(T32-T21*CP2)*W1/T32                                           INT25370
      W4=CP3/(T32*SP3*SP3)                                              INT25380
      W5=(T32-T43*CP3)*W2/T32                                           INT25390
      W6=CP2/(T32*SP2*SP2)                                              INT25400
      DO 10  I=1,3                                                      INT25410
      V2(I)=W3*V5(I)+W4*V6(I)                                           INT25420
 10   V3(I)=W5*V6(I)+W6*V5(I)                                           INT25430
      RETURN                                                            INT25440
      END
C     ////////////////////////////////////////////////////////////      INT25460
C LINX    (OPT2) (NUMB)                                                 INT25100
      SUBROUTINE VECT8(NAD,K1,K2,K3,K4,V1,V2,V3,V4,XA,W)                INT25110
      IMPLICIT REAL*8 (A-H,O-Z)                                         INT25120
      DIMENSION XA(NAD,3),V1(3),V2(3),V3(3),V4(3)
      DIMENSION E1(3),E2(3),E3(3),E32(3),E34(3),H3(3,3),H11(3,3)
      DIMENSION H21(3,3),H31(3,3),H22(3,3),H32(3,3),H33(3,3)
      PARAMETER(ZERO=0.0D0,ONE=1.0D0,
     $          PI=3.14159265358979323846264D0)
      CALL VECT1(NAD,K2,K3,E32,XA,T32)                                  INT25170
      CALL VECT1(NAD,K4,K3,E34,XA,T34)                                  INT25180
      CALL VECT2(NAD,K1,K2,K3,E1,E2,E3,XA,T123)                         INT25160
      CALL SCAPRO(E34,E3,T)                                             INT25220
      W=-T32*T                                                          INT25230
      CALL HIJS1(NAD,K3,K4,XA,H3) 
      CALL HIJS2(NAD,K1,K2,K3,XA,H11,H21,H31,H22,H32,H33) 
      DO 5  I=1,3                                                       INT25340
      V1(I)=ZERO                                                        INT25350
      V4(I)=ZERO                                                        INT25350
      V2(I)=W*E32(I)/T32 
      DO 1 J=1,3
      V1(I)=V1(I)-T32*E34(J)*H31(J,I)
      V2(I)=V2(I)-T32*E34(J)*H32(J,I)
 1    V4(I)=V4(I)-T32*E3(J)*H3(J,I)
 5    V3(I)=-V1(I)-V2(I)-V4(I)                                          INT25360
      RETURN                                                            INT25440
      END                                                               INT25450
C     ////////////////////////////////////////////////////////////      INT25460
C LINY    (OPT2) (NUMB)                                                 INT25100
      SUBROUTINE VECT9(NAD,K1,K2,K3,K4,V1,V2,V3,V4,XA,W)                INT25110
      IMPLICIT REAL*8 (A-H,O-Z)                                         INT25120
      DIMENSION XA(NAD,3),V1(3),V2(3),V3(3),V4(3)
      DIMENSION E4(3),E3(3),E2(3),E1(3)
      PARAMETER(ZERO=0.0D0,ONE=1.0D0,
     $PI=3.14159265358979323846264D0)
      CALL VECT5(NAD,K4,K3,K2,K1,E4,E3,E2,E1,XA,TOUT)
      W=-DSIN(TOUT)
      COSY=DCOS(TOUT)                                                   INT25310
      DO 5  I=1,3                                                       INT25340
      V2(I)=-COSY*E2(I)
      V3(I)=-COSY*E3(I)
      V4(I)=-COSY*E4(I)   
 5    V1(I)=-COSY*E1(I)   
      RETURN                                                            INT25440
      END                                                               INT25450
C     ////////////////////////////////////////////////////////////
C RCOM
      SUBROUTINE VECT10(NAD,K1,K2,K3,K4,V1,XA,XMASS,XMA,XMB,W)
      IMPLICIT REAL*8 (A-H,O-Z)                                         INT23830
      DIMENSION XA(NAD,3),XMASS(1),V1(3),RA(3),RB(3)                    INT23840
      PARAMETER(ZERO=0.0D0,ONE=1.0D0)                                   INT23850
      W=ZERO
      XMA=ZERO
      XMB=ZERO
      DO 2  J=1,3
      RA(J)=ZERO
  2   RB(J)=ZERO
      DO 5  K = K1,K2
        XMA=XMA+XMASS(K)
        DO 5  J = 1,3
  5     RA(J)=RA(J)+XA(K,J)*XMASS(K)
      DO 10  K = K3,K4
        XMB=XMB+XMASS(K)
        DO 10  J = 1,3
 10     RB(J)=RB(J)+XA(K,J)*XMASS(K)
      DO 15  J = 1,3
        RA(J)=RA(J)/XMA
        RB(J)=RB(J)/XMB
 15     V1(J)=RA(J)-RB(J)
      CALL SCAPRO(V1,V1,W)
      W=DSQRT(W)
      DO 20  J = 1,3
 20   V1(J)=V1(J)/W
      RETURN
      END
C     ////////////////////////////////////////////////////////////      INT25460
C STRE    (OPT2)  (NUMB)                                                INT25470
      SUBROUTINE HIJS1(NAD,K1,K2,XA,H11)                                INT25480
      IMPLICIT REAL*8 (A-H,O-Z)                                         INT25490
      DIMENSION XA(NAD,3),V1(3),H11(3,3)                                INT25500
      PARAMETER(ONE=1.0D0)                                              INT25510
      CALL VECT1(NAD,K1,K2,V1,XA,T21)                                   INT25520
      DO 5  J=1,3                                                       INT25530
      DO 5  I=J,3                                                       INT25540
  5   H11(I,J)=-V1(I)*V1(J)                                             INT25550
      DO 10  I=1,3                                                      INT25560
  10  H11(I,I)=H11(I,I)+ONE                                             INT25570
      W=ONE/T21                                                         INT25580
      DO 15  J=1,3                                                      INT25590
      DO 15  I=J,3                                                      INT25600
  15  H11(I,J)=H11(I,J)*W                                               INT25610
      DO 20  J=1,2                                                      INT25620
      DO 20  I=J+1,3                                                    INT25630
  20  H11(J,I)=H11(I,J)                                                 INT25640
      RETURN                                                            INT25650
      END                                                               INT25660
C     //////////////////////////////////////////////////////////////    INT25670
C BEND    (OPT2)  (NUMB)                                                INT25680
      SUBROUTINE HIJS2(NAD,K1,K2,K3,XA,H11,H21,H31,H22,H32,H33)         INT25690
      IMPLICIT REAL*8 (A-H,O-Z)                                         INT25700
      DIMENSION XA(NAD,3),V1(3),V2(3),V3(3),E21(3),E23(3)               INT25710
      DIMENSION H11(3,3),H21(3,3),H31(3,3),H22(3,3),H32(3,3),H33(3,3)   INT25720
      DIMENSION H11A(3,3),H33A(3,3)                                     INT25730
      PARAMETER(ONE=1.0D0)                                              INT25740
         CALL VECT2(NAD,K1,K2,K3,V1,V2,V3,XA,PHI)                       INT25750
         CALL VECT1(NAD,K1,K2,E21,XA,T21)                               INT25760
         CALL VECT1(NAD,K3,K2,E23,XA,T23)                               INT25770
         CALL HIJS1(NAD,K1,K2,XA,H11A)                                  INT25780
         CALL HIJS1(NAD,K3,K2,XA,H33A)                                  INT25790
      SPHI=DSIN(PHI)                                                    INT25800
      CTPHI=DCOS(PHI)/SPHI                                              INT25810
      W1=CTPHI                                                          INT25820
      W2=ONE/T21                                                        INT25830
      W3=W1*W2                                                          INT25840
      W4=ONE/T23                                                        INT25850
      W5=W1*W4                                                          INT25860
         DO 5  J=1,3                                                    INT25870
         DO 5  I=J,3                                                    INT25880
         H11(I,J)=H11A(I,J)*W3-V1(I)*V1(J)*W1                           INT25890
     $                        -(E21(I)*V1(J)+V1(I)*E21(J))*W2           INT25900
         H33(I,J)=H33A(I,J)*W5-V3(I)*V3(J)*W1                           INT25910
     $                        -(E23(I)*V3(J)+V3(I)*E23(J))*W4           INT25920
  5      CONTINUE                                                       INT25930
         DO 10  J=1,2                                                   INT25940
         DO 10  I=J+1,3                                                 INT25950
         H11(J,I)=H11(I,J)                                              INT25960
  10     H33(J,I)=H33(I,J)                                              INT25970
         W3=ONE/(T21*SPHI)                                              INT25980
         DO 15  J=1,3                                                   INT25990
         W4=W2*E21(J)+W1*V1(J)                                          INT26000
         DO 15  I=1,3                                                   INT26010
         H31(I,J)=-H33A(I,J)*W3-V3(I)*W4                                INT26020
         H21(I,J)=-(H11(I,J)+H31(I,J))                                  INT26030
         H32(I,J)=-(H31(I,J)+H33(I,J))                                  INT26040
 15      CONTINUE                                                       INT26050
         DO 20  J=1,3                                                   INT26060
         DO 20  I=1,3                                                   INT26070
         H22(I,J)=-(H21(J,I)+H32(I,J))                                  INT26080
 20      CONTINUE                                                       INT26090
      RETURN                                                            INT26100
      END                                                               INT26110
C     //////////////////////////////////////////////////////////////    INT26120
C LIN1   (OPT2) (NUMB)                                                  INT26130
      SUBROUTINE HIJS3(NAD,K1,K2,K3,K4,XA,H11,H21,H31,H22,H32,H33)      INT26140
      IMPLICIT REAL*8 (A-H,O-Z)                                         INT26150
      DIMENSION XA(NAD,3),V1(3),V2(3),V3(3),E21(3),E23(3),EA(3)         INT26160
      DIMENSION H11(3,3),H21(3,3),H31(3,3),H22(3,3),H32(3,3),H33(3,3)   INT26170
      DIMENSION H11A(3,3),H33A(3,3),EM(3,3)                             INT26180
      PARAMETER(ZERO=0.0D0,ONE=1.0D0)                                   INT26190
      CALL VECT3(NAD,K1,K2,K3,K4,V1,V2,V3,XA,TH)                        INT26200
      CALL VECT1(NAD,K1,K2,E21,XA,T21)                                  INT26210
      CALL VECT1(NAD,K3,K2,E23,XA,T23)                                  INT26220
      CALL HIJS1(NAD,K1,K2,XA,H11A)                                     INT26230
      CALL HIJS1(NAD,K3,K2,XA,H33A)                                     INT26240
         DO 5  J=1,3                                                    INT26250
  5      EA(J)=XA(K4,J)                                                 INT26260
           CALL SCAPRO(EA,EA,D)                                         INT26270
           D=ONE/DSQRT(D)                                               INT26280
           DO 10  I=1,3                                                 INT26290
  10       EA(I)=D*EA(I)                                                INT26300
      TANTH=DTAN(TH)                                                    INT26310
      COSTH=DCOS(TH)                                                    INT26320
      EM(2,1)=EA(3)                                                     INT26330
      EM(3,1)=-EA(2)                                                    INT26340
      EM(3,2)=EA(1)                                                     INT26350
      EM(1,2)=-EM(2,1)                                                  INT26360
      EM(1,3)=-EM(3,1)                                                  INT26370
      EM(2,3)=-EM(3,2)                                                  INT26380
      DO 15  I=1,3                                                      INT26390
 15   EM(I,I)=ZERO                                                      INT26400
      DO 20  J=1,3                                                      INT26410
      DO 20  I=1,3                                                      INT26420
      H22(I,J)=ZERO                                                     INT26430
      DO 20  K=1,3                                                      INT26440
 20   H22(I,J)=H22(I,J)+EM(I,K)*H33A(K,J)                               INT26450
      W1=ONE/T21                                                        INT26460
      W2=ONE/T23                                                        INT26470
         DO 25  J=1,3                                                   INT26480
         DO 25  I=1,3                                                   INT26490
         H11(I,J)=(-H11A(I,J)*W1+V1(I)*V1(J))*TANTH-(E21(I)*V1(J)       INT26500
     $                  +V1(I)*E21(J))*W1                               INT26510
         H31(I,J)=(H22(J,I)/COSTH-V3(I)*E21(J))/T21+V3(I)*V1(J)*TANTH   INT26520
         H33(I,J)=(-H33A(I,J)*W2+V3(I)*V3(J))*TANTH-(E23(I)*V3(J)       INT26530
     $                  +V3(I)*E23(J))*W2                               INT26540
         H21(I,J)=-(H11(I,J)+H31(I,J))                                  INT26550
         H32(I,J)=-(H31(I,J)+H33(I,J))                                  INT26560
 25      CONTINUE                                                       INT26570
         DO 30  J=1,3                                                   INT26580
         DO 30  I=1,3                                                   INT26590
         H22(I,J)=-(H21(J,I)+H32(I,J))                                  INT26600
 30      CONTINUE                                                       INT26610
      RETURN                                                            INT26620
      END                                                               INT26630
C     //////////////////////////////////////////////////////////////    INT26640
C TORS   (OPT2)   (NUMB)                                                INT26650
      SUBROUTINE HIJS6(NAD,K1,K2,K3,K4,XA,H11,H21,H31,H41,              INT26660
     $               H22,H32,H42,H33,H43,H44)                           INT26670
      IMPLICIT REAL*8 (A-H,O-Z)                                         INT26680
      DIMENSION XA(NAD,3),V1(3),V2(3),V3(3),V4(3),E21(3),E23(3),E34(3)  INT26690
      DIMENSION BP21(3),BP22(3),BP23(3),BP32(3),BP34(3)                 INT26700
      DIMENSION H11(3,3),H21(3,3),H31(3,3),H41(3,3),H22(3,3),H32(3,3)   INT26710
      DIMENSION H42(3,3),H33(3,3),H43(3,3),H44(3,3)                     INT26720
      PARAMETER(ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0)                         INT26730
         CALL VECT6(NAD,K1,K2,K3,K4,V1,V2,V3,V4,XA,W)                   INT26740
         CALL VECT1(NAD,K1,K2,E21,XA,T21)                               INT26750
         CALL VECT1(NAD,K3,K2,E23,XA,T23)                               INT26760
         CALL VECT1(NAD,K4,K3,E34,XA,T34)                               INT26770
         CALL VECT2(NAD,K1,K2,K3,BP21,BP22,BP23,XA,P2)                  INT26780
         CALL VECT2(NAD,K2,K3,K4,BP32,V2,BP34,XA,P3)                    INT26790
      CALL MAT1(H11,E23)                                                INT26800
      CALL MAT1(H31,E21)                                                INT26810
      CALL MAT1(H44,E23)                                                INT26820
      CALL MAT1(H42,E34)                                                INT26830
      XX=DSIN(P2)                                                       INT26840
      XY=DSIN(P3)                                                       INT26850
      XX=T21*XX*XX                                                      INT26860
      XY=T34*XY*XY                                                      INT26870
      W1=ONE/(T21*XX)                                                   INT26880
      W2=ONE/(T23*XX)                                                   INT26890
      W3=ONE/(T34*XY)                                                   INT26900
      W4=ONE/(T23*XY)                                                   INT26910
         DO 5  J=1,3                                                    INT26920
         DO 5  I=1,3                                                    INT26930
         H11(I,J)=-H11(I,J)*W1                                          INT26940
         H31(I,J)=H31(I,J)*W2                                           INT26950
         H44(I,J)=H44(I,J)*W3                                           INT26960
  5      H42(I,J)=-H42(I,J)*W4                                          INT26970
      XX=DCOS(P2)/DSIN(P2)                                              INT26980
      XY=DCOS(P3)/DSIN(P3)                                              INT26990
         DO 10  I=1,3                                                   INT27000
         W1=TWO*(E21(I)/T21+BP21(I)*XX)                                 INT27010
         W2=(E23(I)/T23+TWO*BP23(I)*XX)                                 INT27020
         W3=TWO*(E34(I)/T34+BP34(I)*XY)                                 INT27030
         W4=(E23(I)/T23-TWO*BP32(I)*XY)                                 INT27040
         DO 10  J=1,3                                                   INT27050
         H11(I,J)=H11(I,J)-W1*V1(J)                                     INT27060
         H31(I,J)=H31(I,J)-W2*V1(J)                                     INT27070
         H44(I,J)=H44(I,J)-W3*V4(J)                                     INT27080
 10      H42(J,I)=H42(J,I)+W4*V4(J)                                     INT27090
         DO 15  J=1,3                                                   INT27100
         DO 15  I=1,3                                                   INT27110
         H41(I,J)=ZERO                                                  INT27120
         H21(I,J)=-(H11(I,J)+H31(I,J))                                  INT27130
 15      H43(I,J)=-(H44(I,J)+H42(I,J))                                  INT27140
      X1=T21/T23                                                        INT27150
      Y1=T34/T23                                                        INT27160
      X2=DCOS(P2)                                                       INT27170
      Y2=DSIN(P2)                                                       INT27180
      X3=DCOS(P3)                                                       INT27190
      Y3=DSIN(P3)                                                       INT27200
      C1=X1*X2-ONE                                                      INT27210
      C2=-X3*Y1                                                         INT27220
      C3=-X2/T23                                                        INT27230
      C4=-X1*Y2                                                         INT27240
      C5=X1*X2/T23                                                      INT27250
      C6=Y1*Y3                                                          INT27260
      C7=-Y1*X3/T23                                                     INT27270
         DO 20  I=1,3                                                   INT27280
         W1=C3*E21(I)+C4*BP22(I)+C5*E23(I)                              INT27290
         W2=C6*BP32(I)+C7*E23(I)                                        INT27300
         DO 20  J=1,3                                                   INT27310
 20      H22(I,J)=C1*H21(I,J)+C2*H42(J,I)+W1*V1(J)+W2*V4(J)             INT27320
         DO 25  J=1,3                                                   INT27330
         DO 25  I=1,3                                                   INT27340
 25      H32(I,J)=-(H21(J,I)+H22(I,J)+H42(I,J))                         INT27350
         DO 30  J=1,3                                                   INT27360
         DO 30  I=1,3                                                   INT27370
 30      H33(I,J)=-(H31(I,J)+H32(I,J)+H43(J,I))                         INT27380
      RETURN                                                            INT27390
      END                                                               INT27400
C     //////////////////////////////////////////////////////////////    INT27410
C SUBROUTINE HIJS7 WRITTEN BY DOUG GIBSON 4/89                          YY 02350
C OUT                                                                   YY 02360
      SUBROUTINE HIJS7(NAD,K1,K2,K3,K4,XA,H11,H21,H31,H41,              YY 02370
     $               H22,H32,H42,H33,H43,H44)                           YY 02380
      IMPLICIT REAL*8 (A-H,O-Z)                                         YY 02390
      DIMENSION XA(NAD,3),V1(3),V2(3),V3(3),V4(3),V5(3),V6(3)           YY 02400
      DIMENSION BP3(3),BP4(3),E21(3),E23(3),E24(3)                      YY 02410
      DIMENSION H11(3,3),H21(3,3),H31(3,3),H41(3,3),H22(3,3)            YY 02420
      DIMENSION H32(3,3),H42(3,3),H33(3,3),H43(3,3),H44(3,3)            YY 02430
      DIMENSION HP43(3,3),HP44(3,3),CP31(3,3),CP41(3,3),CP43(3,3)       YY 02440
      PARAMETER(ZERO=0.0D0,ONE=1.0D0)                                   YY 02450
         CALL VECT1(NAD,K1,K2,E21,XA,T21)                               YY 02460
         CALL VECT1(NAD,K3,K2,E23,XA,T23)                               YY 02470
         CALL VECT1(NAD,K4,K2,E24,XA,T24)                               YY 02480
         CALL VECT2(NAD,K3,K2,K4,BP3,V2,BP4,XA,PHI)                     YY 02490
         CALL VECT5(NAD,K1,K2,K3,K4,V1,V2,V3,V4,XA,GAMMA)               YY 02500
         CALL HIJS2(NAD,K3,K2,K4,XA,H11,H21,HP43,H22,H32,HP44)          YY 02510
      CALL VECPRO(E23,E24,V5)                                           YY 02520
      CALL VECPRO(E24,E21,V6)                                           YY 02530
      CALL MAT1(CP31,E24)                                               YY 02540
      CALL MAT1(CP41,E23)                                               YY 02550
      CALL MAT1(CP43,E21)                                               YY 02560
      SP=DSIN(PHI)                                                      YY 02570
      CP=DCOS(PHI)                                                      YY 02580
      TP=SP/CP                                                          YY 02590
      SG=DSIN(GAMMA)                                                    YY 02600
      CG=DCOS(GAMMA)                                                    YY 02610
      TG=SG/CG                                                          YY 02620
      C21=ONE/T21                                                       YY 02630
      C23=ONE/T23                                                       YY 02640
      C24=ONE/T24                                                       YY 02650
      CTP=ONE/TP                                                        YY 02660
      C11=TG*C21*C21                                                    YY 02670
      C312=C21/(CG*SP)                                                  YY 02680
      C311=C312*C23                                                     YY 02690
      C313=C312*CTP                                                     YY 02700
      C411=C312*C24                                                     YY 02710
      C3=C23/SP                                                         YY 02720
      C4=C24/SP                                                         YY 02730
      C331=T24*C3                                                       YY 02740
      C332=C331*TG                                                      YY 02750
      C441=T23*C4                                                       YY 02760
      C442=C441*TG                                                      YY 02770
      C431=C3*C24/CG                                                    YY 02780
      C432=TG                                                           YY 02790
      C434=TG*C3                                                        YY 02800
      C435=T24*C3                                                       YY 02810
      C436=C435*TG                                                      YY 02820
      DO 3 J=1,3                                                        YY 02830
      DO 3 I=J,3                                                        YY 02840
         H11(I,J)=V1(J)*(TG*V1(I)-E21(I)*C21)-V1(I)*E21(J)*C21          YY 02850
         H11(I,J)=H11(I,J)+E21(I)*E21(J)*C11                            YY 02860
         IF (I.EQ.J) H11(I,J)=H11(I,J)-C11                              YY 02870
         H11(J,I)=H11(I,J)                                              YY 02880
         H33(I,J)=V3(I)*BP4(J)*C331+HP43(J,I)*C332                      YY 02890
         H33(I,J)=H33(I,J)+V3(J)*(TG*V3(I)-E23(I)*C23-BP3(I)*CTP)       YY 02900
         H33(J,I)=H33(I,J)                                              YY 02910
         H44(I,J)=V4(I)*BP3(J)*C441+HP43(I,J)*C442                      YY 02920
         H44(I,J)=H44(I,J)+V4(J)*(TG*V4(I)-E24(I)*C24-BP4(I)*CTP)       YY 02930
 3       H44(J,I)=H44(I,J)                                              YY 02940
      DO 5  J=1,3                                                       YY 02950
      XJ=TG*V1(J)-E21(J)*C21                                            YY 02960
      DO 5  I=1,3                                                       YY 02970
         H31(I,J)=V3(I)*XJ-CP31(I,J)*C311                               YY 02980
         H31(I,J)=H31(I,J)-E23(I)*V5(J)*C311-BP3(I)*V5(J)*C313          YY 02990
         H41(I,J)=V4(I)*XJ+CP41(I,J)*C411                               YY 03000
         H41(I,J)=H41(I,J)-E24(I)*V5(J)*C411-BP4(I)*V5(J)*C313          YY 03010
         H21(I,J)=-(H11(I,J)+H31(I,J)+H41(I,J))                         YY 03020
         H43(I,J)=(CP43(J,I)-E24(I)*V6(J))*C431+V3(J)*V4(I)*C432        YY 03030
         H43(I,J)=H43(I,J)-V3(J)*BP4(I)*CTP+E24(I)*BP4(J)*C434          YY 03040
 5       H43(I,J)=H43(I,J)+V4(I)*BP4(J)*C435+HP44(I,J)*C436             YY 03050
      DO 10 I=1,3                                                       YY 03060
      DO 10 J=1,3                                                       YY 03070
         H32(I,J)=-(H31(I,J)+H33(I,J)+H43(J,I))                         YY 03080
         H42(I,J)=-(H41(I,J)+H43(I,J)+H44(I,J))                         YY 03090
 10      H22(I,J)=-(H21(J,I)+H32(I,J)+H42(I,J))                         YY 03100
      RETURN                                                            YY 03110
      END                                                               YY 03120
C     /////////////////////////////////////////////////////////////////
      SUBROUTINE HIJS8(NAD,K1,K2,K3,K4,XA,H11,H21,H31,H41,    
     $               H22,H32,H42,H33,H43,H44)                           INT26670
      IMPLICIT REAL*8 (A-H,O-Z)                                         INT26680
      DIMENSION XA(NAD,3),E2(3),E4(3),Q1(3),Q2(3),Q3(3)
      DIMENSION E22(3,3),E44(3,3)
      DIMENSION H11(3,3),H21(3,3),H31(3,3),H41(3,3),H22(3,3),H32(3,3)   INT26710
      DIMENSION H42(3,3),H33(3,3),H43(3,3),H44(3,3)                     INT26720
      DIMENSION Q11(3,3),Q21(3,3),Q31(3,3),Q33(3,3),Q22(3,3),Q32(3,3)   INT26710
      DIMENSION Q111(3,3,3),Q222(3,3,3),Q333(3,3,3),Q444(3,3,3)         INT26710
      DIMENSION Q112(3,3,3),Q223(3,3,3),Q331(3,3,3)                     INT26710
      DIMENSION Q221(3,3,3),Q332(3,3,3),Q113(3,3,3),Q123(3,3,3)         INT26710
      PARAMETER(ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0)                         INT26730
      CALL VECT1(NAD,K2,K3,E2,XA,T32)                                   INT27470
      CALL VECT1(NAD,K4,K3,E4,XA,T34)                      
      CALL VECT2(NAD,K1,K2,K3,Q1,Q2,Q3,XA,T123)
      CALL SCAPRO(E4,Q3,T)
      W=-T32*T
      CALL HIJS1(NAD,K4,K3,XA,E44)                                      INT27480
      CALL HIJS2(NAD,K1,K2,K3,XA,Q11,Q21,Q31,Q22,Q32,Q33) 
      CALL HIJS1(NAD,K2,K3,XA,E22)                                      INT27480
      CALL HIJKS1(NAD,K4,K3,XA,Q444)                                    INT27760
      DO 1 J=1,3
      DO 1 K=1,3
      H44(J,K)=ZERO
      DO 1 I=1,3
  1   H44(J,K)=H44(J,K)-T32*Q444(I,J,K)*Q3(I)
      CALL HIJKS2(NAD,K1,K2,K3,XA,Q111,Q112,Q113,Q123,Q221,             YY 03380
     $    Q222,Q223,Q331,Q332,Q333)                                     YY 03390
      DO 5 K=1,3
      DO 5 J=1,3
      H41(J,K)=ZERO
      H42(J,K)=ZERO
      H11(J,K)=ZERO
      H21(J,K)=ZERO
      H22(J,K)=W*E22(J,K)/T32
      DO 5 I=1,3
      H11(J,K)=H11(J,K)-T32*E4(I)*Q113(J,K,I) 
      H21(J,K)=H21(J,K)-E4(I)*( E2(J)*Q31(I,K)+T32*Q123(K,J,I) )
      H22(J,K)=H22(J,K)-E4(I)*( E2(J)*Q32(I,K)+E2(K)*Q32(I,J)+
     $      T32*Q223(J,K,I) )
      H41(J,K)=H41(J,K)-T32*E44(I,J)*Q31(I,K)   
  5   H42(J,K)=H42(J,K)-E44(I,J)*( T32*Q32(I,K)+E2(K)*Q3(I) )
      DO 3 J=1,3
      DO 3 K=1,3
      H31(J,K)=-H11(J,K)-H21(J,K)-H41(J,K)
      H32(J,K)=-H21(K,J)-H22(J,K)-H42(J,K)
C     H32(J,K)=ZERO
C     DO 2 I=1,3
C     H32(J,K)=H32(J,K)-E2(K)*( -Q3(I)*E44(I,J)+E4(I)*Q33(I,J) )   
C     H32(J,K)=H32(J,K)+E4(I)*(  Q3(I)*E22(K,J)+E2(J)*Q32(I,K) )   
C 2   H32(J,K)=H32(J,K)-T32*( E4(I)*Q332(I,J,K)-E44(J,I)*Q32(I,K) )  
  3   H43(J,K)=-H41(J,K)-H42(J,K)-H44(J,K)
      DO 4 J=1,3
      DO 4 K=1,3
  4   H33(J,K)=-H31(J,K)-H32(J,K)-H43(K,J)
      RETURN
      END                                                               INT27400
C     //////////////////////////////////////////////////////////////////
      SUBROUTINE HIJS9(NAD,K1,K2,K3,K4,XA,H11,H21,H31,H41,         
     $               H22,H32,H42,H33,H43,H44)                           INT26670
      IMPLICIT REAL*8 (A-H,O-Z)                                         INT26680
      DIMENSION XA(NAD,3),E1(3),E2(3),E3(3),E4(3)                       INT26690
      DIMENSION H11(3,3),H21(3,3),H31(3,3),H41(3,3),H22(3,3),H32(3,3)   INT26710
      DIMENSION H42(3,3),H33(3,3),H43(3,3),H44(3,3)                     INT26720
      DIMENSION Q11(3,3),Q12(3,3),Q13(3,3),Q14(3,3),Q22(3,3),Q23(3,3)   INT26710
      DIMENSION Q24(3,3),Q33(3,3),Q34(3,3),Q44(3,3)                     INT26720
      PARAMETER(ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0)                         INT26730
      CALL VECT5(NAD,K4,K3,K2,K1,E4,E3,E2,E1,XA,TOUT)                   INT25180
      W=-DSIN(TOUT)
      COSY=DCOS(TOUT)                                                   INT25310
      CALL HIJS7(NAD,K4,K3,K2,K1,XA,Q44,Q34,Q24,Q14,                    YY 00250
     $           Q33,Q23,Q13,Q22,Q12,Q11)                               YY 00260
      DO 1 K=1,3
      DO 1 J=1,3
      H22(J,K)=-W*E2(J)*E2(K)-COSY*Q22(K,J)        
      H32(J,K)=-W*E3(J)*E2(K)-COSY*Q23(K,J)        
      H42(J,K)=-W*E4(J)*E2(K)-COSY*Q24(K,J)        
      H33(J,K)=-W*E3(J)*E3(K)-COSY*Q33(K,J)        
      H43(J,K)=-W*E4(J)*E3(K)-COSY*Q34(K,J)        
      H44(J,K)=-W*E4(J)*E4(K)-COSY*Q44(K,J)                         
      H41(J,K)=-W*E4(J)*E1(K)-COSY*Q14(K,J)                         
      H31(J,K)=-W*E3(J)*E1(K)-COSY*Q13(K,J)                         
      H21(J,K)=-W*E2(J)*E1(K)-COSY*Q12(K,J)                         
  1   H11(J,K)=-W*E1(J)*E1(K)-COSY*Q11(K,J)                         
      RETURN
      END                                                               INT27400
C     /////////////////////////////////////////////////////////////////
C     RCOM
      SUBROUTINE HIJS10(NAD,K1,K2,K3,K4,XA,XMASS,XMA,XMB,H11)
      IMPLICIT REAL*8 (A-H,O-Z)                                         INT27440
      DIMENSION XA(NAD,3),XMASS(1),RA(3),RB(3),V1(3),H11(3,3)
      PARAMETER(ZERO=0.0D0,ONE=1.0D0)
      XMA=ZERO
      XMB=ZERO
      DO 2  J=1,3
      RA(J)=ZERO
 2    RB(J)=ZERO
      DO 5  K = K1,K2
      XMA=XMA+XMASS(K)
         DO 5  J = 1,3
 5       RA(J)=RA(J)+XMASS(K)*XA(K,J)
      DO 10  K = K3,K4
      XMB=XMB+XMASS(K)
         DO 10  J = 1,3
 10      RB(J)=RB(J)+XMASS(K)*XA(K,J)
      DO 15  J = 1,3
         RA(J)=RA(J)/XMA
         RB(J)=RB(J)/XMB
 15      V1(J)=RA(J)-RB(J)
      CALL SCAPRO(V1,V1,W)
      W=DSQRT(W)
      DO 20  I = 1,3
 20     V1(I)=V1(I)/W
      DO 25  J = 1,3
      DO 25  I = 1,3
 25     H11(I,J)=-V1(I)*V1(J)
      DO 30  I = 1,3
 30     H11(I,I)=H11(I,I)+ONE
      DO 35  J = 1,3
      DO 35  I = 1,3 
 35     H11(I,J)=H11(I,J)/W
      RETURN
      END
C     /////////////////////////////////////////////////////////////////
C STRE  (OPT2)   (NUMB)                                                 INT27420
      SUBROUTINE HIJKS1(NAD,K1,K2,XA,H111)                              INT27430
      IMPLICIT REAL*8 (A-H,O-Z)                                         INT27440
      DIMENSION XA(NAD,3),V1(3),H11(3,3),H111(3,3,3)                    INT27450
      PARAMETER(ONE=1.0D0)                                              INT27460
      CALL VECT1(NAD,K1,K2,V1,XA,T21)                                   INT27470
      CALL HIJS1(NAD,K1,K2,XA,H11)                                      INT27480
      W1=ONE/T21                                                        INT27490
      DO 5  K=1,3                                                       INT27500
      DO 5  J=K,3                                                       INT27510
      DO 5  I=J,3                                                       INT27520
  5   H111(I,J,K)=-(V1(I)*H11(K,J)+V1(J)*H11(K,I)+V1(K)*H11(J,I))*W1    INT27530
      CALL FILL3B(3,3,H111)                                             INT27540
      RETURN                                                            INT27550
      END                                                               INT27560
C     //////////////////////////////////////////////////////////////    INT27570
C BEND   (OPT2)  (NUMB)                                                 INT27580
      SUBROUTINE HIJKS2(NAD,K1,K2,K3,XA,H111,H112,H113,H123,H221,       INT27590
     $    H222,H223,H331,H332,H333)                                     INT27600
      IMPLICIT REAL*8 (A-H,O-Z)                                         INT27610
      DIMENSION XA(NAD,3),V1(3),V2(3),V3(3),E21(3),E23(3)               INT27620
      DIMENSION H11(3,3),H21(3,3),H31(3,3),H22(3,3),H32(3,3),H33(3,3)   INT27630
      DIMENSION H111(3,3,3),H112(3,3,3),H113(3,3,3),H123(3,3,3)         INT27640
      DIMENSION H221(3,3,3),H222(3,3,3),H223(3,3,3),H331(3,3,3)         INT27650
      DIMENSION H332(3,3,3),H333(3,3,3)                                 INT27660
      DIMENSION H111A(3,3,3),H333A(3,3,3),H11A(3,3),H33A(3,3)           INT27670
      PARAMETER(ONE=1.0D0)                                              INT27680
      CALL VECT2(NAD,K1,K2,K3,V1,V2,V3,XA,PHI)                          INT27690
      CALL VECT1(NAD,K1,K2,E21,XA,T21)                                  INT27700
      CALL VECT1(NAD,K3,K2,E23,XA,T23)                                  INT27710
      CALL HIJS1(NAD,K1,K2,XA,H11A)                                     INT27720
      CALL HIJS1(NAD,K3,K2,XA,H33A)                                     INT27730
      CALL HIJS2(NAD,K1,K2,K3,XA,H11,H21,H31,H22,H32,H33)               INT27740
      CALL HIJKS1(NAD,K1,K2,XA,H111A)                                   INT27750
      CALL HIJKS1(NAD,K3,K2,XA,H333A)                                   INT27760
      SPHI=DSIN(PHI)                                                    INT27770
      CTPHI=DCOS(PHI)/SPHI                                              INT27780
      W1=ONE/T21                                                        INT27790
      W2=ONE/T23                                                        INT27800
      W3=CTPHI*W1                                                       INT27810
      W4=CTPHI*W2                                                       INT27820
      DO 10  K=1,3                                                      INT27830
      W5=V1(K)*CTPHI+E21(K)*W1                                          INT27840
      W6=E21(K)*W3                                                      INT27850
      W7=V1(K)*W1                                                       INT27860
      W8=V3(K)*CTPHI+E23(K)*W2                                          INT27870
      W9=E23(K)*W4                                                      INT27880
      W10=V3(K)*W2                                                      INT27890
      DO 10  J=1,3                                                      INT27900
      DO 10  I=1,3                                                      INT27910
      H221(I,J,K)=W5*H11(I,J)+V1(I)*V1(J)*W6+H11A(I,J)*W7               INT27920
 10   H223(I,J,K)=W8*H33(I,J)+V3(I)*V3(J)*W9+H33A(I,J)*W10              INT27930
      DO 15  K=1,3                                                      INT27940
      DO 15  J=K,3                                                      INT27950
      DO 15  I=J,3                                                      INT27960
      H111(I,J,K)=-(H221(I,J,K)+H221(J,K,I)+H221(I,K,J))+               INT27970
     $    V1(I)*V1(J)*V1(K)+H111A(I,J,K)*W3                             INT27980
      H333(I,J,K)=-(H223(I,J,K)+H223(J,K,I)+H223(I,K,J))+               INT27990
     $    V3(I)*V3(J)*V3(K)+H333A(I,J,K)*W4                             INT28000
 15   CONTINUE                                                          INT28010
      CALL FILL3B(3,3,H111)                                             INT28020
      CALL FILL3B(3,3,H333)                                             INT28030
      DO 20  I=1,3                                                      INT28040
      W3=V1(I)*CTPHI+E21(I)*W1                                          INT28050
      W4=V3(I)*CTPHI+E23(I)*W2                                          INT28060
      DO 20  J=1,3                                                      INT28070
      DO 20  K=1,3                                                      INT28080
      H221(I,J,K)=W3*H31(K,J)                                           INT28090
 20   H223(I,J,K)=W4*H31(J,K)                                           INT28100
      W3=ONE/(SPHI*SPHI)                                                INT28110
      DO 25  K=1,3                                                      INT28120
      DO 25  J=1,3                                                      INT28130
      DO 25  I=1,3                                                      INT28140
      H113(I,J,K)=V3(K)*(V1(I)*V1(J)-H11A(I,J)*W1)*W3                   INT28150
     $        -H221(I,J,K)-H221(J,I,K)                                  INT28160
      H331(I,J,K)=V1(K)*(V3(I)*V3(J)-H33A(I,J)*W2)*W3                   INT28170
     $        -H223(I,J,K)-H223(J,I,K)                                  INT28180
 25   CONTINUE                                                          INT28190
      DO 30  K=1,3                                                      INT28200
      DO 30  J=1,3                                                      INT28210
      DO 30  I=1,3                                                      INT28220
      H123(I,J,K)=-(H331(J,K,I)+H113(I,J,K))                            INT28230
      H112(I,J,K)=-(H111(I,J,K)+H113(I,J,K))                            INT28240
 30   H332(I,J,K)=-(H333(I,J,K)+H331(I,J,K))                            INT28250
      DO 35  K=1,3                                                      INT28260
      DO 35  J=1,3                                                      INT28270
      DO 35  I=1,3                                                      INT28280
      H221(J,K,I)=-(H123(I,J,K)+H112(I,K,J))                            INT28290
 35   H223(J,K,I)=-(H332(I,J,K)+H123(J,K,I))                            INT28300
      DO 40  K=1,3                                                      INT28310
      DO 40  J=1,3                                                      INT28320
      DO 40  I=1,3                                                      INT28330
 40   H222(I,J,K)=-(H223(J,K,I)+H221(J,K,I))                            INT28340
      RETURN                                                            INT28350
      END                                                               INT28360
C     ///////////////////////////////////////////////////////////////   INT28370
C LIN1  (OPT1)                                                          INT28380
      SUBROUTINE HIJKS3(NAD,K1,K2,K3,K4,XA,                             INT28390
     $    H111,H112,H113,H123,H221,H222,H223,H331,H332,H333)            INT28400
      IMPLICIT REAL*8 (A-H,O-Z)                                         INT28410
      DIMENSION XA(NAD,3),V1(3),V2(3),V3(3),E21(3),E23(3),EA(3)         INT28420
      DIMENSION H11(3,3),H21(3,3),H31(3,3),H22(3,3),H32(3,3),H33(3,3)   INT28430
      DIMENSION H111(3,3,3),H112(3,3,3),H113(3,3,3),H123(3,3,3)         INT28440
      DIMENSION H221(3,3,3),H222(3,3,3),H223(3,3,3),H331(3,3,3)         INT28450
      DIMENSION H332(3,3,3),H333(3,3,3)                                 INT28460
      DIMENSION H111A(3,3,3),H333A(3,3,3),H11A(3,3),H33A(3,3)           INT28470
      PARAMETER(ONE=1.0D0)                                              INT28480
         CALL VECT3(NAD,K1,K2,K3,K4,V1,V2,V3,XA,TH)                     INT28490
         CALL VECT1(NAD,K1,K2,E21,XA,T21)                               INT28500
         CALL VECT1(NAD,K3,K2,E23,XA,T23)                               INT28510
         CALL HIJS1(NAD,K1,K2,XA,H11A)                                  INT28520
         CALL HIJS1(NAD,K3,K2,XA,H33A)                                  INT28530
         DO 4  J=1,3                                                    INT28540
  4      EA(J)=XA(K4,J)                                                 INT28550
           CALL SCAPRO(EA,EA,D)                                         INT28560
           D=ONE/DSQRT(D)                                               INT28570
           DO 6  I=1,3                                                  INT28580
  6        EA(I)=D*EA(I)                                                INT28590
         CALL HIJS3(NAD,K1,K2,K3,K4,XA,H11,H21,H31,H22,H32,H33)         INT28600
         CALL HIJKS1(NAD,K1,K2,XA,H111A)                                INT28610
         CALL HIJKS1(NAD,K3,K2,XA,H333A)                                INT28620
      TANTH=DTAN(TH)                                                    INT28630
      COSTH=DCOS(TH)                                                    INT28640
      W1=ONE/T21                                                        INT28650
      W2=ONE/T23                                                        INT28660
      W3=TANTH*W1                                                       INT28670
      W4=TANTH*W2                                                       INT28680
      DO 10  K=1,3                                                      INT28690
      DO 10  J=1,3                                                      INT28700
      DO 10  I=1,3                                                      INT28710
      H221(I,J,K)=H11(I,J)*(V1(K)*TANTH-E21(K)/T21)                     INT28720
      H221(I,J,K)=H221(I,J,K)+V1(K)*V1(J)*E21(I)*TANTH/T21              INT28730
      H221(I,J,K)=H221(I,J,K)-(H11A(I,J)*V1(K))/T21                     INT28740
      H223(I,J,K)=H33(I,J)*(V3(K)*TANTH-E23(K)/T23)                     INT28750
      H223(I,J,K)=H223(I,J,K)+V3(K)*V3(J)*E23(I)*TANTH/T23              INT28760
 10   H223(I,J,K)=H223(I,J,K)-(H33A(I,J)*V3(K))/T23                     INT28770
      DO 15  K=1,3                                                      INT28780
      DO 15  J=K,3                                                      INT28790
      DO 15  I=J,3                                                      INT28800
      H111(I,J,K)=(H221(I,J,K)+H221(J,K,I)+H221(K,I,J))+                INT28810
     $    V1(I)*V1(J)*V1(K)-H111A(I,J,K)*W3                             INT28820
      H333(I,J,K)=(H223(I,J,K)+H223(J,K,I)+H223(K,I,J))+                INT28830
     $    V3(I)*V3(J)*V3(K)-H333A(I,J,K)*W4                             INT28840
 15   CONTINUE                                                          INT28850
      CALL FILL3B(3,3,H111)                                             INT28860
      CALL FILL3B(3,3,H333)                                             INT28870
      DO 20  I=1,3                                                      INT28880
      W5=V1(I)*TANTH-E21(I)*W1                                          INT28890
      W6=V3(I)*TANTH-E23(I)*W2                                          INT28900
      DO 20  J=1,3                                                      INT28910
      DO 20  K=1,3                                                      INT28920
      H221(I,J,K)=W5*H31(K,J)                                           INT28930
 20   H223(I,J,K)=W6*H31(J,K)                                           INT28940
      W5=ONE/(COSTH*COSTH)                                              INT28950
      DO 25  K=1,3                                                      INT28960
      DO 25  J=1,3                                                      INT28970
      DO 25  I=1,3                                                      INT28980
      H113(I,J,K)=V3(K)*(V1(I)*V1(J)-H11A(I,J)*W1)*W5                   INT28990
     $        +H221(I,J,K)+H221(J,I,K)                                  INT29000
      H331(I,J,K)=V1(K)*(V3(I)*V3(J)-H33A(I,J)*W2)*W5                   INT29010
     $        +H223(I,J,K)+H223(J,I,K)                                  INT29020
 25   CONTINUE                                                          INT29030
      DO 30  K=1,3                                                      INT29040
      DO 30  J=1,3                                                      INT29050
      DO 30  I=1,3                                                      INT29060
      H123(I,J,K)=-(H331(J,K,I)+H113(I,J,K))                            INT29070
      H112(I,J,K)=-(H111(I,J,K)+H113(I,J,K))                            INT29080
 30   H332(I,J,K)=-(H333(I,J,K)+H331(I,J,K))                            INT29090
      DO 35  K=1,3                                                      INT29100
      DO 35  J=1,3                                                      INT29110
      DO 35  I=1,3                                                      INT29120
      H221(J,K,I)=-(H123(I,J,K)+H112(I,K,J))                            INT29130
 35   H223(J,K,I)=-(H332(I,J,K)+H123(J,K,I))                            INT29140
      DO 40  K=1,3                                                      INT29150
      DO 40  J=1,3                                                      INT29160
      DO 40  I=1,3                                                      INT29170
 40   H222(I,J,K)=-(H223(J,K,I)+H221(J,K,I))                            INT29180
      RETURN                                                            INT29190
      END                                                               INT29200
C     ///////////////////////////////////////////////////////////////   INT29210
C  TORS     (OPT2)                                                      INT29220
      SUBROUTINE HIJKS6(NAD,K1,K2,K3,K4,XA,H111,H112,H221,H222,H113,    INT29230
     $       H123,H223,H331,H332,H333,H411,H421,H422,H431,H432,H433,    INT29240
     $       H441,H442,H443,H444)                                       INT29250
      IMPLICIT REAL*8 (A-H,O-Z)                                         INT29260
      DIMENSION XA(NAD,3),V1(3),V2(3),V3(3),V4(3),E21(3),E23(3),E34(3)  INT29270
      DIMENSION BP21(3),BP22(3),BP23(3),BP32(3),BP33(3),BP34(3)         INT29280
      DIMENSION H11(3,3),H21(3,3),H31(3,3),H41(3,3),H22(3,3),H32(3,3)   INT29290
      DIMENSION H42(3,3),H33(3,3),H43(3,3),H44(3,3)                     INT29300
      DIMENSION H111(3,3,3),H112(3,3,3),H221(3,3,3),H222(3,3,3)         INT29310
      DIMENSION H113(3,3,3),H123(3,3,3),H223(3,3,3),H331(3,3,3)         INT29320
      DIMENSION H332(3,3,3),H333(3,3,3),H411(3,3,3),H421(3,3,3)         INT29330
      DIMENSION H422(3,3,3),H431(3,3,3),H432(3,3,3),H433(3,3,3)         INT29340
      DIMENSION H441(3,3,3),H442(3,3,3),H443(3,3,3),H444(3,3,3)         INT29350
      PARAMETER(ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0)                         INT29360
      CALL VECT6(NAD,K1,K2,K3,K4,V1,V2,V3,V4,XA,TAU)                    INT29370
      CALL VECT1(NAD,K1,K2,E21,XA,T21)                                  INT29380
      CALL VECT1(NAD,K3,K2,E23,XA,T23)                                  INT29390
      CALL VECT1(NAD,K4,K3,E34,XA,T34)                                  INT29400
      CALL VECT2(NAD,K1,K2,K3,BP21,BP22,BP23,XA,P2)                     INT29410
      CALL VECT2(NAD,K2,K3,K4,BP32,BP33,BP34,XA,P3)                     INT29420
      CALL MAT1(H32,E23)                                                INT29430
      CALL MAT1(H21,E21)                                                INT29440
      CALL MAT1(H43,E34)                                                INT29450
      C1=ONE/T21                                                        INT29460
      C2=ONE/T34                                                        INT29470
      C3=ONE/T23                                                        INT29480
      C4=DSIN(P2)                                                       INT29490
      C5=DCOS(P2)                                                       INT29500
      C6=C5/C4                                                          INT29510
      C7=DSIN(P3)                                                       INT29520
      C8=DCOS(P3)                                                       INT29530
      C9=C8/C7                                                          INT29540
      C10=ONE/C4**2                                                     INT29550
      C11=ONE/C7**2                                                     INT29560
      C12=C1*C1                                                         INT29570
      C13=C2*C2                                                         INT29580
      C14=C3*C3                                                         INT29590
      C15=T21*C3                                                        INT29600
      C16=T34*C3                                                        INT29610
      W1=TWO*C6                                                         INT29620
      W2=TWO*C9                                                         INT29630
      W3=TWO*C1                                                         INT29640
      W4=TWO*C2                                                         INT29650
      W5=C5*C3                                                          INT29660
      W6=C8*C3                                                          INT29670
      DO 10  K=1,3                                                      INT29680
      H411(1,1,K)=E21(K)*C1+BP21(K)*C6                                  INT29690
      H411(1,2,K)=E34(K)*C2+BP34(K)*C9                                  INT29700
      H411(1,3,K)=E23(K)*C3+BP23(K)*W1                                  INT29710
      H411(2,1,K)=-E23(K)*C3+BP32(K)*W2                                 INT29720
      H411(2,2,K)=E21(K)*W3+E23(K)*C3-BP22(K)*W1                        INT29730
      H411(2,3,K)=E34(K)*W4-E23(K)*C3-BP33(K)*W2                        INT29740
      H411(3,1,K)=E23(K)*W5+BP23(K)*C4                                  INT29750
 10   H411(3,2,K)=-E23(K)*W6+BP32(K)*C7                                 INT29760
      DO 15  K=1,3                                                      INT29770
        DO 16  M=1,3                                                    INT29780
 16     V2(M)=ZERO                                                      INT29790
      V2(K)=ONE                                                         INT29800
      CALL MAT1(H22,V2)                                                 INT29810
      DO 20  J=1,3                                                      INT29820
      DO 20  I=1,3                                                      INT29830
 20   H421(I,J,K)=H22(I,J)                                              INT29840
 15   CONTINUE                                                          INT29850
      W1=TWO*C10*C12                                                    INT29860
      W2=TWO*C11*C13                                                    INT29870
      DO 100 K=1,3                                                      INT29880
      W3=W1*H411(1,1,K)                                                 INT29890
      W4=W2*H411(1,2,K)                                                 INT29900
        DO 100  J=1,K                                                   INT29910
        DO 100  I=1,J                                                   INT29920
        H111(I,J,K)=W3*H32(I,J)                                         INT29930
 100    H444(I,J,K)=-W4*H32(I,J)                                        INT29940
      W1=C10*C12                                                        INT29950
      W2=C11*C13                                                        INT29960
      W3=W1*C3                                                          INT29970
      W4=W2*C3                                                          INT29980
      DO 105 K=1,3                                                      INT29990
      W5=W1*H411(1,3,K)                                                 INT30000
      W6=W2*H411(2,1,K)                                                 INT30010
        DO 105  J=1,3                                                   INT30020
        DO 105  I=1,3                                                   INT30030
        H113(I,J,K)=W5*H32(I,J)-W3*H421(I,J,K)                          INT30040
 105    H442(I,J,K)=-W6*H32(I,J)-W4*H421(I,J,K)                         INT30050
      W1=C1*C3*C10                                                      INT30060
      W2=C2*C3*C11                                                      INT30070
      DO 110 K=1,3                                                      INT30080
        W5=W1*H411(2,2,K)                                               INT30090
        W6=W2*H411(2,3,K)                                               INT30100
        DO 110  J=1,3                                                   INT30110
        DO 110  I=1,3                                                   INT30120
        H123(I,K,J)=-W5*H21(I,J)+W3*H421(I,J,K)                         INT30130
 110    H432(I,K,J)=-W6*H43(I,J)+W4*H421(I,J,K)                         INT30140
      CALL HIJS2(NAD,K1,K2,K3,XA,H11,H21,H31,H22,H32,H33)               INT30150
      CALL HIJS1(NAD,K1,K2,XA,H44)                                      INT30160
      CALL HIJS1(NAD,K2,K3,XA,H42)                                      INT30170
      W1=TWO*C1                                                         INT30180
      W2=TWO*C12                                                        INT30190
      DO 115  K=1,3                                                     INT30200
        DO 116  I=1,K                                                   INT30210
 116    H43(I,K)=TWO*(W1*H44(I,K)+C6*H11(I,K)-C10*BP21(I)*BP21(K))      INT30220
 115  H43(K,K)=H43(K,K)-W2                                              INT30230
      DO 120  K=1,3                                                     INT30240
      DO 120  J=1,K                                                     INT30250
      DO 120  I=1,J                                                     INT30260
 120  H111(I,J,K)=H111(I,J,K)-V1(J)*H43(I,K)                            INT30270
        W1=TWO*C6                                                       INT30280
        W2=TWO*C10                                                      INT30290
        W3=TWO*C3                                                       INT30300
        DO 125  K=1,3                                                   INT30310
        DO 125  I=1,3                                                   INT30320
 125    H43(I,K)=H31(K,I)*W1-BP21(I)*BP23(K)*W2                         INT30330
        DO 130  K=1,3                                                   INT30340
        DO 130  J=1,3                                                   INT30350
        DO 130  I=1,3                                                   INT30360
 130    H113(I,J,K)=H113(I,J,K)-V1(J)*H43(I,K)                          INT30370
        DO 135  K=1,3                                                   INT30380
        DO 136  J=1,3                                                   INT30390
 136    H43(J,K)=W3*H42(J,K)-W1*H32(K,J)+W2*BP22(J)*BP23(K)             INT30400
 135    H43(K,K)=H43(K,K)-C14                                           INT30410
        DO 140  K=1,3                                                   INT30420
        DO 140  J=1,3                                                   INT30430
        DO 140  I=1,3                                                   INT30440
 140    H123(I,J,K)=H123(I,J,K)+V1(I)*H43(J,K)                          INT30450
      W1=C4*C3                                                          INT30460
      W2=C4*C15                                                         INT30470
      W3=C5*C15                                                         INT30480
      W4=W3*C3                                                          INT30490
      W5=C3*C15                                                         INT30500
      DO 145  K=1,3                                                     INT30510
      DO 145  I=1,3                                                     INT30520
      W6=-E21(I)*BP23(K)*W1+H32(K,I)*W2+BP22(I)*BP23(K)*W3-H42(I,K)*W4  INT30530
 145  H43(I,K)=W6                                                       INT30540
      DO 150  K=1,3                                                     INT30550
      DO 150  I=1,3                                                     INT30560
      W6=H43(I,K)+E23(I)*W5*H411(3,1,K)                                 INT30570
      DO 150  J=1,3                                                     INT30580
 150  H223(I,J,K)=-V1(J)*W6                                             INT30590
      W1=C3*C4*C15                                                      INT30600
      W2=C5*C14                                                         INT30610
      DO 155  K=1,3                                                     INT30620
      DO 155  I=1,3                                                     INT30630
      W3=-E23(K)*BP22(I)*W1+E23(K)*W2*(C15*E23(I)-E21(I))               INT30640
 155  H43(I,K)=H43(I,K)+W3                                              INT30650
      DO 160  K=1,3                                                     INT30660
      DO 160  J=1,3                                                     INT30670
      DO 160  I=1,3                                                     INT30680
 160  H332(I,J,K)=V1(J)*H43(K,I)                                        INT30690
      CALL HIJS2(NAD,K2,K3,K4,XA,H22,H32,H42,H33,H43,H44)               INT30700
      CALL HIJS1(NAD,K4,K3,XA,H11)                                      INT30710
      CALL HIJS1(NAD,K3,K2,XA,H31)                                      INT30720
      W1=TWO*C2                                                         INT30730
      W2=TWO*C13                                                        INT30740
      DO 165  K=1,3                                                     INT30750
      DO 166  I=1,K                                                     INT30760
 166  H21(I,K)=TWO*(W1*H11(I,K)+C9*H44(I,K)-C11*BP34(I)*BP34(K))        INT30770
 165  H21(K,K)=H21(K,K)-W2                                              INT30780
      DO 170  K=1,3                                                     INT30790
      DO 170  J=1,K                                                     INT30800
      DO 170  I=1,J                                                     INT30810
 170  H444(I,J,K)=H444(I,J,K)-V4(J)*H21(I,K)                            INT30820
      W1=TWO*C9                                                         INT30830
      W2=TWO*C11                                                        INT30840
      W3=TWO*C3                                                         INT30850
      DO 175  K=1,3                                                     INT30860
      DO 175  I=1,3                                                     INT30870
 175  H21(I,K)=W1*H42(I,K)-W2*BP34(I)*BP32(K)                           INT30880
      DO 180  K=1,3                                                     INT30890
      DO 180  J=1,3                                                     INT30900
      DO 180  I=1,3                                                     INT30910
 180  H442(I,J,K)=H442(I,J,K)-V4(J)*H21(I,K)                            INT30920
      DO 185  K=1,3                                                     INT30930
      DO 186  J=1,3                                                     INT30940
 186  H21(J,K)=W3*H31(J,K)-W1*H32(J,K)+W2*BP33(J)*BP32(K)               INT30950
 185  H21(K,K)=H21(K,K)-C14                                             INT30960
      DO 190  K=1,3                                                     INT30970
      DO 190  J=1,3                                                     INT30980
      DO 190  I=1,3                                                     INT30990
 190  H432(I,J,K)=H432(I,J,K)+V4(I)*H21(J,K)                            INT31000
      W1=C7*C3                                                          INT31010
      W2=C7*C16                                                         INT31020
      W3=C8*C16                                                         INT31030
      W4=W3*C3                                                          INT31040
      W5=T34*C14                                                        INT31050
      DO 195  K=1,3                                                     INT31060
      DO 195  I=1,3                                                     INT31070
      W6=-E34(I)*BP32(K)*W1+H32(I,K)*W2+BP33(I)*BP32(K)*W3-H31(I,K)*W4  INT31080
 195  H21(I,K)=W6                                                       INT31090
      DO 200  K=1,3                                                     INT31100
      DO 200  I=1,3                                                     INT31110
      W6=H21(I,K)-E23(I)*W5*H411(3,2,K)                                 INT31120
      DO 200  J=1,3                                                     INT31130
 200  H332(I,J,K)=H332(I,J,K)-V4(J)*W6                                  INT31140
      W1=C3*C7*C16                                                      INT31150
      W2=C8*C14                                                         INT31160
      DO 205  K=1,3                                                     INT31170
      DO 205  I=1,3                                                     INT31180
      W3=E23(K)*BP33(I)*W1+E23(K)*W2*(E34(I)+C16*E23(I))                INT31190
 205  H21(I,K)=H21(I,K)+W3                                              INT31200
      DO 210  K=1,3                                                     INT31210
      DO 210  J=1,3                                                     INT31220
      DO 210  I=1,3                                                     INT31230
 210  H223(I,J,K)=H223(I,J,K)+V4(J)*H21(K,I)                            INT31240
      CALL HIJS6(NAD,K1,K2,K3,K4,XA,H11,H21,H31,H41,H22,H32,H42,        INT31250
     $            H33,H43,H44)                                          INT31260
      DO 220  K=1,3                                                     INT31270
      DO 220  J=1,K                                                     INT31280
      DO 220  I=1,J                                                     INT31290
      H111(I,J,K)=H111(I,J,K)-TWO*H11(J,K)*H411(1,1,I)                  INT31300
 220  H444(I,J,K)=H444(I,J,K)-TWO*H44(J,K)*H411(1,2,I)                  INT31310
      CALL FILL3A(3,3,H111)                                             INT31320
      CALL FILL3A(3,3,H444)                                             INT31330
      DO 225  I=1,3                                                     INT31340
      W1=TWO*H411(1,1,I)                                                INT31350
      W2=TWO*H411(1,2,I)                                                INT31360
      DO 225  J=1,3                                                     INT31370
      DO 225  K=1,3                                                     INT31380
      H113(I,J,K)=H113(I,J,K)-W1*H31(K,J)                               INT31390
      H442(I,J,K)=H442(I,J,K)-W2*H42(J,K)                               INT31400
      H123(I,J,K)=H123(I,J,K)-H21(J,I)*H411(1,3,K)                      INT31410
 225  H432(I,J,K)=H432(I,J,K)-H43(I,J)*H411(2,1,K)                      INT31420
      W4=C5*C15                                                         INT31430
      W1=W4-ONE                                                         INT31440
      W2=C8*C16                                                         INT31450
      W3=W2-ONE                                                         INT31460
      DO 230  K=1,3                                                     INT31470
      DO 230  J=1,3                                                     INT31480
      DO 230  I=1,3                                                     INT31490
      H223(I,J,K)=H223(I,J,K)+W1*H123(J,I,K)-W2*H432(J,K,I)             INT31500
 230  H332(I,J,K)=H332(I,J,K)+W3*H432(J,I,K)-W4*H123(J,K,I)             INT31510
      DO 240  K=1,3                                                     INT31520
      DO 240  J=1,3                                                     INT31530
      DO 240  I=1,3                                                     INT31540
      H223(I,J,K)=H223(I,J,K)-C15*H21(I,J)*H411(3,1,K)                  INT31550
 240  H332(I,J,K)=H332(I,J,K)-C16*H43(J,I)*H411(3,2,K)                  INT31560
      DO 250  K=1,3                                                     INT31570
      DO 250  J=1,3                                                     INT31580
      W1=C16*(H43(J,K)-C3*V4(J)*E23(K))                                 INT31590
      W2=C15*(H21(K,J)+C3*V1(J)*E23(K))                                 INT31600
      DO 250  I=1,3                                                     INT31610
      H223(I,J,K)=H223(I,J,K)+W1*H411(3,2,I)                            INT31620
 250  H332(I,J,K)=H332(I,J,K)+W2*H411(3,1,I)                            INT31630
      W1=C5*C3                                                          INT31640
      W2=C4*C15                                                         INT31650
      W3=C5*T21*C14                                                     INT31660
      W4=C8*C3                                                          INT31670
      W5=C7*C16                                                         INT31680
      W6=C8*T34*C14                                                     INT31690
      DO 260  K=1,3                                                     INT31700
      H411(1,1,K)=W5*BP33(K)+W6*E23(K)+W4*E34(K)                        INT31710
      H411(1,2,K)=W2*BP22(K)-W3*E23(K)+W1*E21(K)                        INT31720
      H411(1,3,K)=-W1*E21(K)-W2*BP22(K)+W3*E23(K)                       INT31730
 260  H411(2,1,K)=-W4*E34(K)-W5*BP33(K)-W6*E23(K)                       INT31740
      DO 265  K=1,3                                                     INT31750
      DO 265  J=1,3                                                     INT31760
      DO 265  I=1,3                                                     INT31770
      H223(I,J,K)=H223(I,J,K)+H42(J,I)*H411(1,1,K)                      INT31780
 265  H332(I,J,K)=H332(I,J,K)+H31(I,J)*H411(1,2,K)                      INT31790
      DO 270  K=1,3                                                     INT31800
      DO 270  J=1,3                                                     INT31810
      W1=(H31(K,J)-C3*V1(J)*E23(K))                                     INT31820
      W2=(H42(J,K)+C3*V4(J)*E23(K))                                     INT31830
      DO 270  I=1,3                                                     INT31840
      H223(I,J,K)=H223(I,J,K)+W1*H411(1,3,I)                            INT31850
 270  H332(I,J,K)=H332(I,J,K)+W2*H411(2,1,I)                            INT31860
      DO 280  K=1,3                                                     INT31870
      DO 280  J=1,3                                                     INT31880
      DO 280  I=1,3                                                     INT31890
      H411(I,J,K)=ZERO                                                  INT31900
      H421(I,J,K)=ZERO                                                  INT31910
      H431(I,J,K)=ZERO                                                  INT31920
      H441(I,J,K)=ZERO                                                  INT31930
      H443(J,K,I)=-(H444(I,J,K)+H442(J,K,I))                            INT31940
 280  H112(I,J,K)=-(H111(I,J,K)+H113(I,J,K))                            INT31950
      DO 285  K=1,3                                                     INT31960
      DO 285  J=1,3                                                     INT31970
      DO 285  I=1,3                                                     INT31980
      H433(K,I,J)=-(H443(I,K,J)+H432(K,J,I))                            INT31990
 285  H221(I,J,K)=-(H112(I,K,J)+H123(K,J,I))                            INT32000
      DO 290  K=1,3                                                     INT32010
      DO 290  J=1,3                                                     INT32020
      DO 290  I=1,3                                                     INT32030
      H422(K,I,J)=-(H432(K,I,J)+H442(I,K,J))                            INT32040
 290  H331(I,J,K)=-(H123(K,I,J)+H113(I,K,J))                            INT32050
      DO 300  K=1,3                                                     INT32060
      DO 300  J=1,K                                                     INT32070
      DO 300  I=1,J                                                     INT32080
      H222(I,J,K)=-(H221(I,J,K)+H223(I,J,K)+H422(K,I,J))                INT32090
 300  H333(I,J,K)=-(H331(I,J,K)+H332(I,J,K)+H433(K,I,J))                INT32100
      CALL FILL3A(3,3,H222)                                             INT32110
      CALL FILL3A(3,3,H333)                                             INT32120
      RETURN                                                            INT32130
      END                                                               INT32140
C     ///////////////////////////////////////////////////////////////   YY 01540
C SUBROUTINE HIJKS7 WRITTEN BY DOUG GIBSON 5/89                         YY 02350
C OUT
      SUBROUTINE HIJKS7(NAD,K1,K2,K3,K4,XA,H111,H112,H221,H222,H113,    YY 00020
     $       H123,H223,H331,H332,H333,H411,H421,H422,H431,H432,H433,    YY 00030
     $       H441,H442,H443,H444)                                       YY 00040
      IMPLICIT REAL*8 (A-H,O-Z)                                         YY 00050
      DIMENSION XA(NAD,3),V1(3),V2(3),V3(3),V4(3),E21(3),E23(3),E24(3)  YY 00060
      DIMENSION BP3(3),BP4(3),HP33(3,3),HP43(3,3),HP44(3,3)             YY 00070
      DIMENSION HT11(3,3),HT111(3,3,3),HP334(3,3,3),HP443(3,3,3)        YY 00080
      DIMENSION CP21(3,3),CP24(3,3),CP2124(3),PROD(3,3,3)               YY 00090
      DIMENSION H11(3,3),H21(3,3),H31(3,3),H41(3,3),H22(3,3),H32(3,3)   YY 00100
      DIMENSION H42(3,3),H33(3,3),H43(3,3),H44(3,3)                     YY 00110
      DIMENSION H111(3,3,3),H112(3,3,3),H221(3,3,3),H222(3,3,3)         YY 00120
      DIMENSION H113(3,3,3),H123(3,3,3),H223(3,3,3),H331(3,3,3)         YY 00130
      DIMENSION H332(3,3,3),H333(3,3,3),H411(3,3,3),H421(3,3,3)         YY 00140
      DIMENSION H422(3,3,3),H431(3,3,3),H432(3,3,3),H433(3,3,3)         YY 00150
      DIMENSION H441(3,3,3),H442(3,3,3),H443(3,3,3),H444(3,3,3)         YY 00160
      PARAMETER(ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0)                         YY 00170
      CALL VECT1(NAD,K1,K2,E21,XA,T21)                                  YY 00180
      CALL VECT1(NAD,K3,K2,E23,XA,T23)                                  YY 00190
      CALL VECT1(NAD,K4,K2,E24,XA,T24)                                  YY 00200
      CALL VECT2(NAD,K3,K2,K4,BP3,V2,BP4,XA,PHI)                        YY 00210
      CALL VECT5(NAD,K1,K2,K3,K4,V1,V2,V3,V4,XA,GAMMA)                  YY 00220
      CALL HIJS1(NAD,K1,K2,XA,HT11)                                     YY 00230
      CALL HIJS2(NAD,K3,K2,K4,XA,HP33,H32,HP43,H22,H42,HP44)            YY 00240
      CALL HIJS7(NAD,K1,K2,K3,K4,XA,H11,H21,H31,H41,                    YY 00250
     $           H22,H32,H42,H33,H43,H44)                               YY 00260
      CALL HIJKS1(NAD,K1,K2,XA,HT111)                                   YY 00270
      CALL HIJKS2(NAD,K3,K2,K4,XA,H111,H112,HP334,H123,H221,            YY 00280
     $            H222,H223,HP443,H332,H333)                            YY 00290
      CALL MAT1(CP21,E21)                                               YY 00300
      CALL MAT1(CP24,E24)                                               YY 00310
      CALL VECPRO(E21,E24,CP2124)                                       YY 00320
      CALL TRIPRO(PROD)                                                 YY 00330
      CG=DCOS(GAMMA)                                                    YY 00340
      SG=DSIN(GAMMA)                                                    YY 00350
      TG=SG/CG                                                          YY 00360
      CP=DCOS(PHI)                                                      YY 00370
      SP=DSIN(PHI)                                                      YY 00380
      TP=SP/CP                                                          YY 00390
      CTP=CP/SP                                                         YY 00400
      C1=ONE/T21                                                        YY 00410
      S2G=ONE/(CG*CG)                                                   YY 00420
      C3=ONE/T23                                                        YY 00430
      C4=ONE/T24                                                        YY 00440
      C2P=ONE/(SP*SP)                                                   YY 00450
      C1111=TG*C1                                                       YY 00460
      C1112=S2G*C1                                                      YY 00470
      C3331=T24*C3/SP                                                   YY 00480
      C3332=C3331*TG                                                    YY 00490
      C3333=C3331*S2G                                                   YY 00500
      C3335=C3*C3                                                       YY 00510
      C3334=TWO*C3335                                                   YY 00520
      C4411=T23*C4/SP                                                   YY 00530
      C4412=C4411*S2G                                                   YY 00540
      C431=C3*C4/(CG*SP)                                                YY 00550
      C4311=C431*C1                                                     YY 00560
      C4312=C431*TG                                                     YY 00570
      C4313=C3333*C4                                                    YY 00580
      C4431=C4411*TG                                                    YY 00590
      C4442=C4*C4                                                       YY 00600
      C4441=TWO*C4442                                                   YY 00610
      DO 10 I=1,3                                                       YY 00620
      DO 10 J=1,I                                                       YY 00630
      DO 10 K=1,J                                                       YY 00640
      H111(I,J,K)=S2G*V1(I)*V1(J)*V1(K)+TG*H11(I,J)*V1(K)               YY 00650
     $    +TG*H11(I,K)*V1(J)+C1111*(E21(I)*V1(J)*V1(K)-HT111(I,J,K))    YY 00660
     $    -C1112*HT11(J,K)*V1(I)-C1*(E21(I)*H11(J,K)+E21(J)*H11(I,K)    YY 00670
     $       +E21(K)*H11(I,J)+HT11(I,J)*V1(K)+HT11(I,K)*V1(J))          YY 00680
      H333(I,J,K)=C3331*(H33(I,J)*BP4(K)+HP43(K,I)*V3(J)                YY 00690
     $       -(V3(J)*BP4(K)+TG*HP43(K,J))*(C3*E23(I)+CTP*BP3(I)))       YY 00700
     $    +C3332*HP334(I,J,K)+C3333*HP43(K,J)*V3(I)                     YY 00710
     $    +H33(I,K)*(TG*V3(J)-CTP*BP3(J)-C3*E23(J))                     YY 00720
     $    +V3(K)*(V3(I)*V3(J)*S2G+H33(I,J)*TG+BP3(I)*BP3(J)*C2P         YY 00730
     $       -HP33(I,J)*CTP+E23(I)*E23(J)*C3334)                        YY 00740
      H444(I,J,K)=C4411*(H44(I,J)*BP3(K)+HP43(I,K)*V4(J)                YY 00750
     $       -(V4(J)*BP3(K)+HP43(J,K)*TG)*(E24(I)*C4+BP4(I)*CTP))       YY 00760
     $    +C4431*HP443(I,J,K)+C4412*HP43(J,K)*V4(I)                     YY 00770
     $    +H44(I,K)*(TG*V4(J)-CTP*BP4(J)-C4*E24(J))                     YY 00780
     $    +V4(K)*(V4(I)*V4(J)*S2G+H44(I,J)*TG+BP4(I)*BP4(J)*C2P         YY 00790
     $       -HP44(I,J)*CTP+E24(I)*E24(J)*C4441)                        YY 00800
      IF (I.EQ.J) THEN                                                  YY 00810
         H333(I,J,K)=H333(I,J,K)-V3(K)*C3335                            YY 00820
         H444(I,J,K)=H444(I,J,K)-V4(K)*C4442                            YY 00830
      END IF                                                            YY 00840
10    CONTINUE                                                          YY 00850
      CALL FILL3B(3,3,H111)                                             YY 00860
      CALL FILL3B(3,3,H333)                                             YY 00870
      CALL FILL3B(3,3,H444)                                             YY 00880
      DO 20 K=1,3                                                       YY 00890
      DO 20 J=1,K                                                       YY 00900
      DO 20 I=1,3                                                       YY 00910
      H113(J,K,I)=S2G*V3(I)*V1(J)*V1(K)+TG*H31(I,J)*V1(K)               YY 00920
     $    +TG*H31(I,K)*V1(J)-C1*(E21(J)*H31(I,K)+E21(K)*H31(I,J))       YY 00930
     $    -C1112*HT11(J,K)*V3(I)                                        YY 00940
      H113(K,J,I)=H113(J,K,I)                                           YY 00950
      H411(I,J,K)=S2G*V4(I)*V1(J)*V1(K)+TG*H41(I,J)*V1(K)               YY 00960
     $    +TG*H41(I,K)*V1(J)-C1*(E21(J)*H41(I,K)+E21(K)*H41(I,J))       YY 00970
     $    -C1112*HT11(J,K)*V4(I)                                        YY 00980
      H411(I,K,J)=H411(I,J,K)                                           YY 00990
      H433(I,J,K)=C3331*(H43(I,J)*BP4(K)+HP44(I,K)*V3(J)                YY 01000
     $       +(C4*E24(I)-CTP*BP4(I))*(V3(J)*BP4(K)+TG*HP43(K,J)))       YY 01010
     $    +C3332*HP443(I,K,J)+C3333*HP43(K,J)*V4(I)                     YY 01020
     $    +H43(I,K)*(TG*V3(J)-CTP*BP3(J)-C3*E23(J))+V3(K)*(TG*H43(I,J)  YY 01030
     $       +S2G*V4(I)*V3(J)-CTP*HP43(I,J)+C2P*BP4(I)*BP3(J))          YY 01040
20    H433(I,K,J)=H433(I,J,K)                                           YY 01050
      DO 30 I=1,3                                                       YY 01060
      DO 30 J=1,I                                                       YY 01070
      DO 30 K=1,3                                                       YY 01080
      H331(I,J,K)=C3331*H31(I,K)*BP4(J)+C3333*HP43(J,I)*V1(K)           YY 01090
     $    +(TG*V3(I)-CTP*BP3(I)-C3*E23(I))*H31(J,K)+TG*H31(I,K)*V3(J)   YY 01100
     $    +S2G*V3(I)*V3(J)*V1(K)                                        YY 01110
      H331(J,I,K)=H331(I,J,K)                                           YY 01120
      H441(I,J,K)=C4411*H41(I,K)*BP3(J)+C4412*HP43(I,J)*V1(K)           YY 01130
     $    +(TG*V4(I)-CTP*BP4(I)-C4*E24(I))*H41(J,K)+TG*H41(I,K)*V4(J)   YY 01140
     $    +S2G*V4(I)*V4(J)*V1(K)                                        YY 01150
      H441(J,I,K)=H441(I,J,K)                                           YY 01160
      H443(I,J,K)=C4411*(H43(J,K)*BP3(I)+HP33(I,K)*V4(J)                YY 01170
     $       +(C3*E23(K)-CTP*BP3(K))*(BP3(I)*V4(J)+TG*HP43(J,I)))       YY 01180
     $    +C4431*HP334(I,K,J)+C4412*HP43(J,I)*V3(K)                     YY 01190
     $    +H43(I,K)*(TG*V4(J)-CTP*BP4(J)-C4*E24(J))+V4(I)*(TG*H43(J,K)  YY 01200
     $       +S2G*V4(J)*V3(K)-CTP*HP43(J,K)+C2P*BP4(J)*BP3(K))          YY 01210
30    H443(J,I,K)=H443(I,J,K)                                           YY 01220
      DO 40 I=1,3                                                       YY 01230
      DO 40 J=1,3                                                       YY 01240
      DO 40 K=1,3                                                       YY 01250
      H431(I,J,K)=C4311*(PROD(K,J,I)+E21(K)*CP21(I,J)+E24(I)*CP24(J,K)  YY 01260
     $    -E24(I)*E21(K)*CP2124(J))+C4312*V1(K)*(E24(I)*CP2124(J)       YY 01270
     $       -CP21(I,J))+TG*H41(I,K)*V3(J)+S2G*V4(I)*V3(J)*V1(K)        YY 01280
     $    +H31(J,K)*(TG*V4(I)-CTP*BP4(I))+C4313*E24(I)*BP4(J)*V1(K)     YY 01290
     $    +C3331*H41(I,K)*BP4(J)+C3333*HP44(I,J)*V1(K)                  YY 01300
40    CONTINUE                                                          YY 01310
      DO 50 I=1,3                                                       YY 01320
      DO 50 J=1,3                                                       YY 01330
      DO 50 K=1,3                                                       YY 01340
      H112(I,J,K)=-(H111(I,J,K)+H113(I,J,K)+H411(K,I,J))                YY 01350
      H421(I,J,K)=-(H411(I,J,K)+H431(I,J,K)+H441(I,J,K))                YY 01360
      H123(I,J,K)=-(H113(I,J,K)+H331(J,K,I)+H431(J,K,I))                YY 01370
      H332(I,J,K)=-(H331(I,J,K)+H333(I,J,K)+H433(K,I,J))                YY 01380
      H432(I,J,K)=-(H431(I,J,K)+H433(I,J,K)+H443(I,K,J))                YY 01390
50    H442(I,J,K)=-(H441(I,J,K)+H443(I,J,K)+H444(I,J,K))                YY 01400
      DO 60 I=1,3                                                       YY 01410
      DO 60 J=1,3                                                       YY 01420
      DO 60 K=1,3                                                       YY 01430
      H221(I,J,K)=-(H112(I,K,J)+H123(K,J,I)+H421(I,J,K))                YY 01440
      H223(I,J,K)=-(H123(I,J,K)+H332(I,K,J)+H432(I,K,J))                YY 01450
60    H422(I,J,K)=-(H421(I,J,K)+H432(I,K,J)+H442(I,K,J))                YY 01460
      DO 70 I=1,3                                                       YY 01470
      DO 70 J=1,I                                                       YY 01480
      DO 70 K=1,J                                                       YY 01490
70    H222(I,J,K)=-(H221(I,J,K)+H223(I,J,K)+H422(K,I,J))                YY 01500
      CALL FILL3B(3,3,H222)                                             YY 01510
      RETURN                                                            YY 01520
      END                                                                      
C     /////////////////////////////////////////////////////////////////
      SUBROUTINE HIJKS8(NAD,K1,K2,K3,K4,XA,H111,H112,H221,H222,H113,
     $       H123,H223,H331,H332,H333,H411,H421,H422,H431,H432,H433,    INT29240
     $       H441,H442,H443,H444)                                       INT29250
      IMPLICIT REAL*8 (A-H,O-Z)                                         INT29260
      DIMENSION XA(NAD,3),Q1(3),Q2(3),Q3(3),Q4(3),QB(3),QC(3)
      DIMENSION Q11(3,3),Q21(3,3),Q31(3,3),Q41(3,3)                     INT29300
      DIMENSION Q22(3,3),Q32(3,3)                                       INT29290
      DIMENSION Q42(3,3),Q33(3,3),Q43(3,3),Q44(3,3)                     INT29300
      DIMENSION H111(3,3,3),H112(3,3,3),H221(3,3,3),H222(3,3,3)         INT29310
      DIMENSION H113(3,3,3),H123(3,3,3),H223(3,3,3),H331(3,3,3)         INT29320
      DIMENSION H332(3,3,3),H333(3,3,3),H411(3,3,3),H421(3,3,3)         INT29330
      DIMENSION H422(3,3,3),H431(3,3,3),H432(3,3,3),H433(3,3,3)         INT29340
      DIMENSION H441(3,3,3),H442(3,3,3),H443(3,3,3),H444(3,3,3)         INT29350
      DIMENSION Q221(3,3,3),Q222(3,3,3),Q111(3,3,3),Q112(3,3,3)         INT29330
      DIMENSION Q113(3,3,3),Q123(3,3,3),Q223(3,3,3),Q331(3,3,3)         INT29340
      DIMENSION Q332(3,3,3),Q444(3,3,3)          
      DIMENSION Q1111(3,3,3,3),Q1112(3,3,3,3),Q1113(3,3,3,3)
      DIMENSION Q1122(3,3,3,3),Q1123(3,3,3,3),Q1133(3,3,3,3)
      DIMENSION Q1222(3,3,3,3),Q1223(3,3,3,3),Q1233(3,3,3,3)
      DIMENSION Q1333(3,3,3,3),Q2222(3,3,3,3),Q2223(3,3,3,3)
      DIMENSION Q2233(3,3,3,3),Q2333(3,3,3,3),Q4444(3,3,3,3)
      PARAMETER(ZERO=0.0D0,SIX=6.0D0,TWO=2.0D0)                         INT29360
      CALL VECT1(NAD,K2,K3,QB,XA,R23)  
      CALL VECT1(NAD,K4,K3,QC,XA,R34)  
      CALL VECT2(NAD,K1,K2,K3,Q1,Q2,Q3,XA,PHI)
      CALL HIJS2(NAD,K1,K2,K3,XA,Q11,Q21,Q31,Q22,Q32,Q33)  
      CALL HIJKS1(NAD,K4,K3,XA,Q444)
      CALL H4TH1(NAD,K4,K3,XA,Q4444)
      DO 1 I=1,3
      DO 1 J=1,3
      DO 1 K=1,3
      H444(I,J,K)=0.0D0
      H441(I,J,K)=0.0D0
      H442(I,J,K)=0.0D0
      DO 1 L=1,3
      H444(I,J,K)=H444(I,J,K)-R23*Q4444(L,I,J,K)*Q3(L)
      H441(I,J,K)=H441(I,J,K)-R23*Q444(L,I,J)*Q31(L,K)
  1   H442(I,J,K)=H442(I,J,K)-R23*Q444(L,I,J)*Q32(L,K)
      CALL HIJS1(NAD,K4,K3,XA,Q44)
      CALL HIJKS2(NAD,K1,K2,K3,XA,Q111,Q112,Q113,Q123,Q221,             
     $    Q222,Q223,Q331,Q332,Q444)                                    
      CALL H4TH2(NAD,K1,K2,K3,XA,Q1111,Q1112,Q1113,Q1122,
     $    Q1123,Q1133,Q1222,Q1223,Q1233,Q1333,Q2222,Q2223,Q2233,
     $    Q2333,Q4444)
      DO 2 I=1,3
      DO 2 J=1,3
      DO 2 K=1,3
      H111(I,J,K)=0.0D0
      H112(I,J,K)=0.0D0
      H221(I,J,K)=0.0D0
      H222(I,J,K)=0.0D0
      H421(I,J,K)=0.0D0
      H422(I,J,K)=0.0D0
      H411(I,J,K)=0.0D0
      DO 2 L=1,3
      H111(I,J,K)=H111(I,J,K)-R23*Q1113(I,J,K,L)*QC(L)
      H112(I,J,K)=H112(I,J,K)-R23*Q1123(I,J,K,L)*QC(L)
      H221(I,J,K)=H221(I,J,K)-R23*Q1223(K,J,I,L)*QC(L)
      H222(I,J,K)=H222(I,J,K)-R23*Q2223(I,J,K,L)*QC(L)
      H421(I,J,K)=H421(I,J,K)-R23*Q123(K,J,L)*Q44(L,I)
      H422(I,J,K)=H422(I,J,K)-R23*Q223(J,K,L)*Q44(L,I)
  2   H411(I,J,K)=H411(I,J,K)-R23*Q113(J,K,L)*Q44(L,I)
      CALL HIJS8(NAD,K1,K2,K3,K4,XA,Q11,Q21,Q31,Q41,
     $               Q22,Q32,Q42,Q33,Q43,Q44)
      DO 5 I=1,3
      DO 5 J=1,3
      DO 5 K=1,3
      H442(I,J,K)=H442(I,J,K)+Q44(I,J)*QB(K)/R23
      H421(I,J,K)=H421(I,J,K)+Q41(I,K)*QB(J)/R23
      H112(I,J,K)=H112(I,J,K)+Q11(I,J)*QB(K)/R23
      H222(I,J,K)=H222(I,J,K)+Q22(J,K)*QB(I)/R23
      H222(I,J,K)=H222(I,J,K)+(Q22(I,K)*QB(J)+Q22(I,J)*QB(K))/R23
      H221(I,J,K)=H221(I,J,K)+(Q21(I,K)*QB(J)+Q21(J,K)*QB(I))/R23
  5   H422(I,J,K)=H422(I,J,K)+(Q42(I,K)*QB(J)+Q42(I,J)*QB(K))/R23
      CALL VECT8(NAD,K1,K2,K3,K4,Q1,Q2,Q3,Q4,XA,W)
      CALL HIJS1(NAD,K2,K3,XA,Q22)
      CALL HIJKS1(NAD,K2,K3,XA,Q222)
      DO 6 I=1,3
      DO 6 J=1,3
      DO 6 K=1,3
      H422(I,J,K)=H422(I,J,K)+Q22(J,K)*Q4(I)/R23
      H422(I,J,K)=H422(I,J,K)-2.0D0*QB(J)*QB(K)*Q4(I)/R23/R23
      H221(I,J,K)=H221(I,J,K)+Q22(I,J)*Q1(K)/R23
      H221(I,J,K)=H221(I,J,K)-2.0D0*QB(I)*QB(J)*Q1(K)/R23/R23
      H222(I,J,K)=H222(I,J,K)+(Q22(I,J)*Q2(K)+Q22(J,K)*Q2(I))/R23
      H222(I,J,K)=H222(I,J,K)+Q22(I,K)*Q2(J)/R23
      H222(I,J,K)=H222(I,J,K)+SIX*W*QB(I)*QB(J)*QB(K)/R23/R23/R23 -
     $       (QB(I)*QB(J)*Q2(K)+QB(J)*QB(K)*Q2(I)+QB(I)*QB(K)*Q2(J))*
     $       2.0D0/R23/R23 + W*Q222(I,J,K)/R23 -
     $ TWO*W*(Q22(I,J)*QB(K)+Q22(I,K)*QB(J)+Q22(J,K)*QB(I))/R23/R23
  6   CONTINUE
      DO 7 I=1,3
      DO 7 J=1,3
      DO 7 K=1,3
      H223(I,J,K)=-H222(I,J,K)-H221(I,J,K)-H422(K,I,J)
      H113(I,J,K)=-H112(I,J,K)-H111(I,J,K)-H411(K,I,J)
      H123(I,J,K)=-H112(I,K,J)-H221(K,J,I)-H421(K,J,I)
      H443(I,J,K)=-H442(I,J,K)-H441(I,J,K)-H444(I,J,K)
      H431(I,J,K)=-H421(I,J,K)-H411(I,J,K)-H441(I,J,K)
  7   H432(I,J,K)=-H422(I,J,K)-H421(I,K,J)-H442(I,J,K)
      DO 8 K=1,3
      DO 8 I=1,3
      DO 8 J=1,3
      H331(I,J,K)=-H431(I,J,K)-H123(K,I,J)-H113(I,K,J)
      H332(I,J,K)=-H432(I,J,K)-H223(I,K,J)-H123(I,K,J)
  8   H433(I,J,K)=-H431(I,J,K)-H432(I,J,K)-H443(K,I,J)
      DO 9 K=1,3
      DO 9 I=1,3
      DO 9 J=1,3
  9   H333(I,J,K)=-H433(I,J,K)-H331(J,K,I)-H332(J,K,I)
      RETURN
      END
C     ////////////////////////////////////////////////////////////////
      SUBROUTINE HIJKS9(NAD,K1,K2,K3,K4,XA,H111,H112,H221,H222,H113,
     $       H123,H223,H331,H332,H333,H411,H421,H422,H431,H432,H433,    INT29240
     $       H441,H442,H443,H444)                                       INT29250
      IMPLICIT REAL*8 (A-H,O-Z)                                         INT29260
      DIMENSION XA(NAD,3),E1(3),E2(3),E3(3),E4(3)
      DIMENSION Q11(3,3),Q12(3,3),Q13(3,3),Q14(3,3)                     INT29300
      DIMENSION Q22(3,3),Q23(3,3)                                       INT29290
      DIMENSION Q24(3,3),Q33(3,3),Q34(3,3),Q44(3,3)                     INT29300
      DIMENSION Q444(3,3,3),Q443(3,3,3),Q334(3,3,3),Q333(3,3,3)         INT29310
      DIMENSION Q442(3,3,3),Q432(3,3,3),Q332(3,3,3),Q224(3,3,3)         INT29320
      DIMENSION Q223(3,3,3),Q222(3,3,3),Q144(3,3,3),Q134(3,3,3)         INT29330
      DIMENSION Q133(3,3,3),Q124(3,3,3),Q123(3,3,3),Q122(3,3,3)         INT29340
      DIMENSION Q114(3,3,3),Q113(3,3,3),Q112(3,3,3),Q111(3,3,3)         INT29350
      DIMENSION H111(3,3,3),H112(3,3,3),H221(3,3,3),H222(3,3,3)         INT29310
      DIMENSION H113(3,3,3),H123(3,3,3),H223(3,3,3),H331(3,3,3)         INT29320
      DIMENSION H332(3,3,3),H333(3,3,3),H411(3,3,3),H421(3,3,3)         INT29330
      DIMENSION H422(3,3,3),H431(3,3,3),H432(3,3,3),H433(3,3,3)         INT29340
      DIMENSION H441(3,3,3),H442(3,3,3),H443(3,3,3),H444(3,3,3)         INT29350
      PARAMETER(ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0)                         INT26730
      CALL VECT5(NAD,K4,K3,K2,K1,E4,E3,E2,E1,XA,TOUT)                   INT25180
      W=-DSIN(TOUT)
      COSY=DCOS(TOUT)                                                   INT25310
      CALL HIJS7(NAD,K4,K3,K2,K1,XA,Q44,Q34,Q24,Q14,                    YY 00250
     $           Q33,Q23,Q13,Q22,Q12,Q11)                               YY 00260
      CALL HIJKS7(NAD,K4,K3,K2,K1,XA,Q444,Q443,Q334,Q333,Q442,
     $       Q432,Q332,Q224,Q223,Q222,Q144,Q134,Q133,Q124,Q123,Q122,    INT29240
     $       Q114,Q113,Q112,Q111)                                       INT29250
      DO 1 K=1,3
      DO 1 I=1,3
      DO 1 J=1,3
      H222(I,J,K)=COSY*E2(I)*E2(J)*E2(K)-COSY*Q222(I,J,K)-
     $   W*(E2(I)*Q22(J,K)+E2(J)*Q22(I,K)+E2(K)*Q22(I,J))
      H223(I,J,K)=COSY*E2(I)*E2(J)*E3(K)-COSY*Q223(I,J,K)-
     $   W*(E2(I)*Q23(J,K)+E2(J)*Q23(I,K)+E3(K)*Q22(I,J))
      H422(I,J,K)=COSY*E4(I)*E2(J)*E2(K)-COSY*Q224(J,K,I)-
     $   W*(E2(K)*Q24(J,I)+E2(J)*Q24(K,I)+E4(I)*Q22(J,K))
      H333(I,J,K)=COSY*E3(I)*E3(J)*E3(K)-COSY*Q333(I,J,K)-
     $   W*(E3(K)*Q33(J,I)+E3(J)*Q33(K,I)+E3(I)*Q33(J,K))
      H433(I,J,K)=COSY*E4(I)*E3(J)*E3(K)-COSY*Q334(J,K,I)-
     $   W*(E3(K)*Q34(J,I)+E3(J)*Q34(K,I)+E4(I)*Q33(J,K))
      H332(I,J,K)=COSY*E3(I)*E3(J)*E2(K)-COSY*Q332(I,J,K)-
     $   W*(E3(I)*Q23(K,J)+E3(J)*Q23(K,I)+E2(K)*Q33(I,J))
      H432(I,J,K)=COSY*E4(I)*E3(J)*E2(K)-COSY*Q432(I,J,K)-
     $   W*(E4(I)*Q23(K,J)+E3(J)*Q24(K,I)+E2(K)*Q34(J,I))
      H444(I,J,K)=COSY*E4(I)*E4(J)*E4(K)-COSY*Q444(I,J,K)-
     $   W*(E4(I)*Q44(K,J)+E4(J)*Q44(K,I)+E4(K)*Q44(I,J))
      H443(I,J,K)=COSY*E4(I)*E4(J)*E3(K)-COSY*Q443(I,J,K)-
     $   W*(E4(I)*Q34(K,J)+E4(J)*Q34(K,I)+E3(K)*Q44(I,J))
      H442(I,J,K)=COSY*E4(I)*E4(J)*E2(K)-COSY*Q442(I,J,K)-
     $   W*(E4(I)*Q24(K,J)+E4(J)*Q24(K,I)+E2(K)*Q44(I,J))
  1   CONTINUE
      DO 2 K=1,3
      DO 2 I=1,3
      DO 2 J=1,3
      H221(I,J,K)=-H222(I,J,K)-H223(I,J,K)-H422(K,I,J)
      H331(I,J,K)=-H332(I,J,K)-H333(I,J,K)-H433(K,I,J)
      H123(I,J,K)=-H332(I,K,J)-H223(I,J,K)-H432(I,K,J)
      H441(I,J,K)=-H442(I,J,K)-H443(I,J,K)-H444(I,J,K)
      H431(I,J,K)=-H432(I,J,K)-H433(I,K,J)-H443(I,K,J)
  2   H421(I,J,K)=-H422(I,J,K)-H432(I,K,J)-H442(I,K,J)
      DO 3 K=1,3
      DO 3 I=1,3
      DO 3 J=1,3
      H112(I,J,K)=-H421(I,K,J)-H123(J,K,I)-H221(I,K,J)
      H113(I,J,K)=-H431(I,K,J)-H331(I,K,J)-H123(J,I,K)
  3   H411(I,J,K)=-H441(I,J,K)-H431(I,J,K)-H421(I,J,K)
      DO 4 K=1,3
      DO 4 I=1,3
      DO 4 J=1,3
  4   H111(I,J,K)=-H411(K,I,J)-H113(I,J,K)-H112(I,J,K)
      RETURN
      END
C     /////////////////////////////////////////////////////////////////
      SUBROUTINE H4TH1(NAD,K1,K2,XA,H1111)                              INT27430
      IMPLICIT REAL*8 (A-H,O-Z)                                         INT27440
      DIMENSION XA(NAD,3),V1(3),H11(3,3),H111(3,3,3),H1111(3,3,3,3)   
      PARAMETER(ONE=1.0D0)                                              INT27460
      CALL VECT1(NAD,K1,K2,V1,XA,T21)                                   INT27470
      CALL HIJS1(NAD,K1,K2,XA,H11)                                      INT27480
      CALL HIJKS1(NAD,K1,K2,XA,H111)                                    INT27480
      DO 5  L=1,3                                                       INT27500
      DO 5  K=1,L                                                       INT27500
      DO 5  J=1,K                                                       INT27510
      DO 5  I=1,J                                                       INT27520
      F=H11(I,L)*H11(K,J)+H11(J,L)*H11(K,I)+H11(I,J)*H11(K,L)
      F=F+V1(I)*H111(J,K,L)+V1(J)*H111(I,K,L)+V1(K)*H111(I,J,L)
      F=F+V1(L)*H111(I,J,K)
  5   H1111(I,J,K,L)=-F/T21
      CALL FILL4A(3,3,H1111)  
      RETURN                                                            INT27550
      END
C     //////////////////////////////////////////////////////////////    INT27570
      SUBROUTINE H4TH2(NAD,K1,K2,K3,XA,H1111,H1112,H1113,H1122,
     $    H1123,H1133,H1222,H1223,H1233,H1333,H2222,H2223,H2233,
     $    H2333,H3333)
      IMPLICIT REAL*8 (A-H,O-Z)                                         INT27610
      DIMENSION XA(NAD,3),V1(3),V2(3),V3(3)
      DIMENSION H11(3,3),H21(3,3),H31(3,3),H22(3,3),H32(3,3),H33(3,3)   INT27630
      DIMENSION H111(3,3,3),H112(3,3,3),H113(3,3,3),H123(3,3,3)         INT27640
      DIMENSION H221(3,3,3),H222(3,3,3),H223(3,3,3),H331(3,3,3)         INT27650
      DIMENSION H332(3,3,3),H333(3,3,3)                                 INT27660
      DIMENSION H1111(3,3,3,3),H1112(3,3,3,3),H1113(3,3,3,3)
      DIMENSION H1122(3,3,3,3),H1123(3,3,3,3),H1133(3,3,3,3)
      DIMENSION H1222(3,3,3,3),H1223(3,3,3,3),H1233(3,3,3,3)
      DIMENSION H1333(3,3,3,3),H2222(3,3,3,3),H2223(3,3,3,3)
      DIMENSION H2233(3,3,3,3),H2333(3,3,3,3),H3333(3,3,3,3)
      DIMENSION Q11111(3,3,3,3,3),Q33333(3,3,3,3,3)           
      DIMENSION Q1111(3,3,3,3),Q3333(3,3,3,3)           
      PARAMETER(ONE=1.0D0)                                              INT27680
      CALL VECT2(NAD,K1,K2,K3,V1,V2,V3,XA,PHI)                          INT27690
      CALL HIJS2(NAD,K1,K2,K3,XA,H11,H21,H31,H22,H32,H33)               INT27740
      CALL HIJKS2(NAD,K1,K2,K3,XA,H111,H112,H113,H123,H221,    
     $    H222,H223,H331,H332,H333)
      CSCP=ONE/DSIN(PHI)            
      COTP=DCOS(PHI)*CSCP  
      DO 10  L=1,3                                                      INT27830
      DO 10  K=1,3                                                      INT27830
      DO 10  J=1,3                                                      INT27900
      DO 10  I=1,3                                                      INT27910
      H1111(I,J,K,L)=V1(I)*V1(J)*V1(K)*V1(L) - H111(J,K,L)*V1(I) -
     $  H111(I,J,K)*V1(L) - H111(I,J,L)*V1(K) - H111(I,K,L)*V1(J) -
     $  H11(I,J)*H11(K,L) - H11(I,L)*H11(K,J) - H11(I,K)*H11(J,L)
      H1111(I,J,K,L)=H1111(I,J,K,L)*COTP +
     $  H11(I,J)*V1(K)*V1(L) + H11(I,K)*V1(J)*V1(L) +
     $  H11(I,L)*V1(J)*V1(K) + H11(J,K)*V1(I)*V1(L) +
     $  H11(J,L)*V1(K)*V1(I) + H11(K,L)*V1(J)*V1(I) 
      H1113(I,J,K,L)=V1(I)*V1(J)*V1(K)*V3(L) - H113(J,K,L)*V1(I) -
     $  H111(I,J,K)*V3(L) - H113(I,J,L)*V1(K) - H113(I,K,L)*V1(J) -
     $  H11(I,J)*H31(L,K) - H31(L,I)*H11(K,J) - H11(I,K)*H31(L,J)
      H1113(I,J,K,L)=H1113(I,J,K,L)*COTP +
     $  H11(I,J)*V1(K)*V3(L) + H11(I,K)*V1(J)*V3(L) +
     $  H31(L,I)*V1(J)*V1(K) + H11(J,K)*V1(I)*V3(L) +
     $  H31(L,J)*V1(K)*V1(I) + H31(L,K)*V1(J)*V1(I) 
      H1133(I,J,K,L)=V1(I)*V1(J)*V3(K)*V3(L) - H331(L,K,J)*V1(I) -
     $  H113(I,J,K)*V3(L) - H113(I,J,L)*V3(K) - H331(L,K,I)*V1(J) -
     $  H11(I,J)*H33(L,K) - H31(L,I)*H31(K,J) - H31(K,I)*H31(L,J)
      H1133(I,J,K,L)=H1133(I,J,K,L)*COTP +
     $  H11(I,J)*V3(K)*V3(L) + H31(K,I)*V1(J)*V3(L) +
     $  H31(L,I)*V1(J)*V3(K) + H31(K,J)*V1(I)*V3(L) +
     $  H31(L,J)*V3(K)*V1(I) + H33(L,K)*V1(J)*V1(I) 
      H1333(I,J,K,L)=V1(I)*V3(J)*V3(K)*V3(L) - H333(L,K,J)*V1(I) -
     $  H331(K,J,I)*V3(L) - H331(L,J,I)*V3(K) - H331(L,K,I)*V3(J) -
     $  H31(J,I)*H33(L,K) - H31(L,I)*H33(K,J) - H31(K,I)*H33(L,J)
      H1333(I,J,K,L)=H1333(I,J,K,L)*COTP +
     $  H31(J,I)*V3(K)*V3(L) + H31(K,I)*V3(J)*V3(L) +
     $  H31(L,I)*V3(J)*V3(K) + H33(K,J)*V1(I)*V3(L) +
     $  H33(L,J)*V3(K)*V1(I) + H33(L,K)*V3(J)*V1(I) 
      H3333(I,J,K,L)=V3(I)*V3(J)*V3(K)*V3(L) - H333(L,K,J)*V3(I) -
     $  H333(K,J,I)*V3(L) - H333(L,J,I)*V3(K) - H333(L,K,I)*V3(J) -
     $  H33(J,I)*H33(L,K) - H33(L,I)*H33(K,J) - H33(K,I)*H33(L,J)
      H3333(I,J,K,L)=H3333(I,J,K,L)*COTP +
     $  H33(J,I)*V3(K)*V3(L) + H33(K,I)*V3(J)*V3(L) +
     $  H33(L,I)*V3(J)*V3(K) + H33(K,J)*V3(I)*V3(L) +
     $  H33(L,J)*V3(K)*V3(I) + H33(L,K)*V3(J)*V3(I) 
 10   CONTINUE
      CALL VECT1(NAD,K1,K2,V1,XA,R1)                                    INT27700
      CALL VECT1(NAD,K3,K2,V3,XA,R3)                                    INT27710
      CALL HIJS1(NAD,K1,K2,XA,H11)                                      INT27720
      CALL HIJS1(NAD,K3,K2,XA,H33)                                      INT27730
      CALL HIJKS1(NAD,K1,K2,XA,H111)                                    INT27750
      CALL HIJKS1(NAD,K3,K2,XA,H333)                                    INT27760
      CALL H4TH1(NAD,K1,K2,XA,Q1111)                                    INT27750
      CALL H4TH1(NAD,K3,K2,XA,Q3333)                                    INT27760
      CALL H5TH1(NAD,K1,K2,XA,Q11111)                                   INT27750
      CALL H5TH1(NAD,K3,K2,XA,Q33333)                                   INT27760
      DO 15  M=1,3                                                      INT27940
      DO 15  L=1,3                                                      INT27940
      DO 15  K=1,3                                                      INT27940
      DO 15  J=1,3                                                      INT27950
      DO 15  I=1,3                                                      INT27960
      H1111(I,J,K,L)=H1111(I,J,K,L)-CSCP*Q11111(M,I,J,K,L)*V3(M)
      H1113(I,J,K,L)=H1113(I,J,K,L)-CSCP*Q1111(M,I,J,K)*H33(L,M)
      H1133(I,J,K,L)=H1133(I,J,K,L)-CSCP*H111(M,I,J)*H333(K,L,M)
      H1333(I,J,K,L)=H1333(I,J,K,L)-CSCP*H11(M,I)*Q3333(J,K,L,M)
      H3333(I,J,K,L)=H3333(I,J,K,L)-CSCP*V1(M)*Q33333(I,J,K,L,M)
 15   CONTINUE                                                          INT28010
      DO 16  L=1,3                                                      INT27940
      DO 16  K=1,3                                                      INT27940
      DO 16  J=1,3                                                      INT27950
      DO 16  I=1,3                                                      INT27960
      H1112(I,J,K,L)=-H1111(I,J,K,L)-H1113(I,J,K,L)
      H1123(I,J,K,L)=-H1113(I,J,K,L)-H1133(I,J,K,L)
      H1233(I,J,K,L)=-H1133(I,J,K,L)-H1333(I,J,K,L)
 16   H2333(I,J,K,L)=-H1333(I,J,K,L)-H3333(I,J,K,L)
      DO 17  L=1,3                                                      INT27940
      DO 17  K=1,3                                                      INT27940
      DO 17  J=1,3                                                      INT27950
      DO 17  I=1,3                                                      INT27960
      H1122(I,J,K,L)=-H1112(I,J,K,L)-H1123(I,J,L,K)
      H1223(I,J,K,L)=-H1123(I,J,K,L)-H1233(I,K,J,L)
 17   H2233(I,J,K,L)=-H1233(I,J,K,L)-H2333(J,I,K,L)
      DO 18  L=1,3                                                      INT27940
      DO 18  K=1,3   
      DO 18  J=1,3                                                      INT27950
      DO 18  I=1,3                                                      INT27960
      H1222(I,J,K,L)=-H1122(I,J,K,L)-H1223(I,K,L,J)
 18   H2223(I,J,K,L)=-H1223(I,J,K,L)-H2233(J,K,I,L)
      DO 19  L=1,3                                                      INT27940
      DO 19  K=1,3   
      DO 19  J=1,3                                                      INT27950
      DO 19  I=1,3                                                      INT27960
 19   H2222(I,J,K,L)=-H1222(I,J,K,L)-H2223(J,K,L,I)
      RETURN                                                            INT28350
      END
C     /////////////////////////////////////////////////////////////////
      SUBROUTINE H5TH1(NAD,K1,K2,XA,H11111)                             INT27430
      IMPLICIT REAL*8 (A-H,O-Z)                                         INT27440
      DIMENSION XA(NAD,3),V1(3),H11(3,3),H111(3,3,3),H1111(3,3,3,3)     INT27450
      DIMENSION H11111(3,3,3,3,3)
      PARAMETER(ONE=1.0D0)                                              INT27460
      CALL VECT1(NAD,K1,K2,V1,XA,T21)                                   INT27470
      CALL HIJS1(NAD,K1,K2,XA,H11)                                      INT27480
      CALL HIJKS1(NAD,K1,K2,XA,H111)                                    INT27480
      CALL H4TH1(NAD,K1,K2,XA,H1111)                                    INT27480
      DO 5  M=1,3                                                       INT27500
      DO 5  L=1,3                                                       INT27500
      DO 5  K=1,L                                                       INT27500
      DO 5  J=1,K                                                       INT27510
      DO 5  I=1,J                                                       INT27520
      A=H11(I,L)*H111(K,J,M)+H11(J,L)*H111(K,I,M)+H11(I,J)*H111(K,L,M)
      B=H111(I,L,M)*H11(K,J)+H111(J,L,M)*H11(K,I)+H111(I,J,M)*H11(K,L)
      C=H11(I,M)*H111(K,J,L)+H11(J,M)*H111(K,I,L)
      D=H11(K,M)*H111(I,J,L)+H11(L,M)*H111(K,I,J)
      E=V1(I)*H1111(J,K,L,M)+V1(J)*H1111(I,K,L,M)+V1(K)*H1111(I,J,L,M)
      F=V1(L)*H1111(I,J,K,M)+V1(M)*H1111(I,J,K,L)
  5   H11111(I,J,K,L,M)=-(A+B+C+D+E+F)/T21
      DO 6 M=1,3
      CALL FILL4A(3,3,H11111(1,1,1,1,M))
  6   CONTINUE
      RETURN                                                            INT27550
      END
C     ///////////////////////////////////////////////////////////////
      SUBROUTINE FILL3A(NX,NY,F3)                                       INT32160
      IMPLICIT REAL*8 (A-H,O-Z)                                         INT32170
      INTEGER M,N,P                                                     INT32180
      DIMENSION F3(NX,NX,NX)                                            INT32190
      DO 5  P=1,NY                                                      INT32200
      DO 6  N=1,P-1                                                     INT32210
      DO 7  M=1,N-1                                                     INT32220
      F3(N,M,P)=F3(M,N,P)                                               INT32230
      F3(N,P,M)=F3(M,N,P)                                               INT32240
      F3(M,P,N)=F3(M,N,P)                                               INT32250
      F3(P,M,N)=F3(M,N,P)                                               INT32260
   7  F3(P,N,M)=F3(M,N,P)                                               INT32270
      F3(N,P,N)=F3(N,N,P)                                               INT32280
   6  F3(P,N,N)=F3(N,N,P)                                               INT32290
      DO 8  M=1,P-1                                                     INT32300
      F3(P,M,P)=F3(M,P,P)                                               INT32310
   8  F3(P,P,M)=F3(M,P,P)                                               INT32320
   5  CONTINUE                                                          INT32330
      RETURN                                                            INT32340
      END                                                               INT32350
C     ///////////////////////////////////////////////////////////////   INT32360
      SUBROUTINE FILL3B(NX,NY,F3)                                       INT32370
      IMPLICIT REAL*8 (A-H,O-Z)                                         INT32380
      INTEGER M,N,P                                                     INT32390
      DIMENSION F3(NX,NX,NX)                                            INT32400
      DO 5   M=1,NY                                                     INT32410
      DO 10  N=1,M-1                                                    INT32420
      DO 15  P=1,N-1                                                    INT32430
      F3(N,M,P)=F3(M,N,P)                                               INT32440
      F3(N,P,M)=F3(M,N,P)                                               INT32450
      F3(M,P,N)=F3(M,N,P)                                               INT32460
      F3(P,M,N)=F3(M,N,P)                                               INT32470
  15  F3(P,N,M)=F3(M,N,P)                                               INT32480
      F3(N,M,N)=F3(M,N,N)                                               INT32490
  10  F3(N,N,M)=F3(M,N,N)                                               INT32500
      DO 20  P=1,M-1                                                    INT32510
      F3(M,P,M)=F3(M,M,P)                                               INT32520
  20  F3(P,M,M)=F3(M,M,P)                                               INT32530
   5  CONTINUE                                                          INT32540
      RETURN                                                            INT32550
      END                                                               INT32560
C     ///////////////////////////////////////////////////////////////   INT32570
      SUBROUTINE FILL4A(NX,NY,F4)                                       INT32580
      IMPLICIT REAL*8 (A-H,O-Z)                                         INT32590
      INTEGER M,N,P,Q                                                   INT32600
      DIMENSION F4(NX,NX,NX,NX)                                         INT32610
      DO 5  Q=1,NY                                                      INT32620
      DO 5  P=1,Q                                                       INT32630
      DO 5  N=1,P                                                       INT32640
      DO 5  M=1,N                                                       INT32650
      F4(N,M,P,Q)=F4(M,N,P,Q)                                           INT32660
      F4(N,P,M,Q)=F4(M,N,P,Q)                                           INT32670
      F4(N,P,Q,M)=F4(M,N,P,Q)                                           INT32680
      F4(M,P,N,Q)=F4(M,N,P,Q)                                           INT32690
      F4(P,M,N,Q)=F4(M,N,P,Q)                                           INT32700
      F4(P,N,M,Q)=F4(M,N,P,Q)                                           INT32710
      F4(P,N,Q,M)=F4(M,N,P,Q)                                           INT32720
      F4(M,P,Q,N)=F4(M,N,P,Q)                                           INT32730
      F4(P,M,Q,N)=F4(M,N,P,Q)                                           INT32740
      F4(P,Q,M,N)=F4(M,N,P,Q)                                           INT32750
      F4(P,Q,N,M)=F4(M,N,P,Q)                                           INT32760
      F4(M,N,Q,P)=F4(M,N,P,Q)                                           INT32770
      F4(N,M,Q,P)=F4(M,N,P,Q)                                           INT32780
      F4(N,Q,M,P)=F4(M,N,P,Q)                                           INT32790
      F4(N,Q,P,M)=F4(M,N,P,Q)                                           INT32800
      F4(M,Q,N,P)=F4(M,N,P,Q)                                           INT32810
      F4(Q,M,N,P)=F4(M,N,P,Q)                                           INT32820
      F4(Q,N,M,P)=F4(M,N,P,Q)                                           INT32830
      F4(Q,N,P,M)=F4(M,N,P,Q)                                           INT32840
      F4(M,Q,P,N)=F4(M,N,P,Q)                                           INT32850
      F4(Q,M,P,N)=F4(M,N,P,Q)                                           INT32860
      F4(Q,P,M,N)=F4(M,N,P,Q)                                           INT32870
      F4(Q,P,N,M)=F4(M,N,P,Q)                                           INT32880
   5  CONTINUE                                                          INT32890
      RETURN                                                            INT32900
      END                                                               INT32910
C///////////////////////////////////////////////////////////////////
	subroutine fill4b (NC, A)
	IMPLICIT REAL*8 (A-H,O-Z)
	dimension A(NC,NC,NC,NC)

	do 1 i = 1, NC
	do 1 j = 1, i
	do 1 k = 1, j
	do 1 l = 1, k
	A(i,j,l,k) = A(i,j,k,l)
        A(i,k,j,l) = A(i,j,k,l)
        A(i,k,l,j) = A(i,j,k,l)
        A(i,l,j,k) = A(i,j,k,l)
        A(i,l,k,j) = A(i,j,k,l)
        A(j,i,k,l) = A(i,j,k,l)
        A(j,i,l,k) = A(i,j,k,l)
        A(j,k,i,l) = A(i,j,k,l)
        A(j,k,l,i) = A(i,j,k,l)
        A(j,l,i,k) = A(i,j,k,l)
        A(j,l,k,i) = A(i,j,k,l)
        A(k,i,j,l) = A(i,j,k,l)
        A(k,i,l,j) = A(i,j,k,l)
        A(k,j,i,l) = A(i,j,k,l)
        A(k,j,l,i) = A(i,j,k,l)
        A(k,l,i,j) = A(i,j,k,l)
        A(k,l,j,i) = A(i,j,k,l)
        A(l,i,j,k) = A(i,j,k,l)
        A(l,i,k,j) = A(i,j,k,l)
        A(l,j,i,k) = A(i,j,k,l)
        A(l,j,k,i) = A(i,j,k,l)
        A(l,k,i,j) = A(i,j,k,l)
 1      A(l,k,j,i) = A(i,j,k,l)
	return
	END
C     ////////////////////////////////////////////////////////////      INT19670
      SUBROUTINE VCXKI(NA,NAD,NC,XMASS,XA,XKI)                          INT19680
      IMPLICIT REAL*8 (A-H,O-Z)                                         INT19690
      INTEGER QQ                                                        INT19700
      DIMENSION XMASS(NA),XA(NAD,3),DIP(3),XKI(NC,3)                    INT19710
      DIMENSION W(3,3),XX(3,6),B(3),RCM(3)                              INT19720
      PARAMETER(ZERO=0.0D0,ONE=1.0D0)                                   INT19730
      COMMON /IO/ IIN1,IOUT,IIN2,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25,
     $   I31,I32,I33,I35,I36,I37,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10,
     $   ISCR11,ISCR12,ISCR13,ISCR14
      COMMON /DIPOLE/DIP,QQ                                             INT19780
      RMT=ZERO                                                          INT19790
      DO 80  K=1,NA                                                     INT19800
  80  RMT=RMT+XMASS(K)                                                  INT19810
      DO 90  I=1,3                                                      INT19820
      RCM(I)=ZERO                                                       INT19830
      DO 95  K=1,NA                                                     INT19840
  95  RCM(I)=RCM(I)+XA(K,I)*XMASS(K)                                    INT19850
  90  RCM(I)=RCM(I)/RMT                                                 INT19860
      DO 100  I=1,3                                                     INT19870
      DO 100  J=1,3                                                     INT19880
      XX(I,J+3)=ZERO                                                    INT19890
 100  W(I,J)=-RMT*RCM(I)*RCM(J)                                         INT19900
      DO 110  I=1,3                                                     INT19910
      DO 110  J=I,3                                                     INT19920
      DO 110  K=1,NA                                                    INT19930
 110  W(I,J)=W(I,J)+XA(K,I)*XA(K,J)*XMASS(K)                            INT19940
      XX(1,1)=W(1,3)                                                    INT19950
      XX(2,1)=W(2,3)                                                    INT19960
      XX(3,1)=-(W(1,1)+W(2,2))                                          INT19970
      XX(1,2)=-W(1,2)                                                   INT19980
      XX(2,2)=W(1,1)+W(3,3)                                             INT19990
      XX(3,2)=-W(2,3)                                                   INT20000
      XX(1,3)=-(W(2,2)+W(3,3))                                          INT20010
      XX(2,3)=W(1,2)                                                    INT20020
      XX(3,3)=W(1,3)                                                    INT20030
      XX(1,4)=ONE                                                       INT20040
      XX(2,5)=ONE                                                       INT20050
      XX(3,6)=ONE                                                       INT20060
      CALL FLIN(XX,3,3,3,DD)                                            INT20070
      DO 1000  IX=1,NC                                                  INT20080
      II=(IX-1)/3+1                                                     INT20090
      IG=IX-3*(II-1)                                                    INT20100
      IF(IG.EQ.1) THEN                                                  INT20110
        B(1)=ZERO                                                       INT20120
        B(2)=XMASS(II)*(XA(II,3)-RCM(3))                                INT20130
        B(3)=-XMASS(II)*(XA(II,2)-RCM(2))                               INT20140
      ELSE IF(IG.EQ.2) THEN                                             INT20150
        B(1)=-XMASS(II)*(XA(II,3)-RCM(3))                               INT20160
        B(2)=ZERO                                                       INT20170
        B(3)=XMASS(II)*(XA(II,1)-RCM(1))                                INT20180
      ELSE IF(IG.EQ.3) THEN                                             INT20190
        B(1)=XMASS(II)*(XA(II,2)-RCM(2))                                INT20200
        B(2)=-XMASS(II)*(XA(II,1)-RCM(1))                               INT20210
        B(3)=ZERO                                                       INT20220
      END IF                                                            INT20230
      DO 1010 I=1,3                                                     INT20240
      XKI(IX,I)=ZERO                                                    INT20250
      DO 1010 J=1,3                                                     INT20260
 1010 XKI(IX,I)=XKI(IX,I)+XX(I,J+3)*B(J)                                INT20270
 1000 CONTINUE                                                          INT20280
      END                                                               INT20290
C     ////////////////////////////////////////////////////////////      INT20300
      SUBROUTINE VCDER1(NA,NC,NSX,NINV,XMASS,XKI,A,F1,F2,F3,F4,V,NRUN)  INT20310
      IMPLICIT REAL*8 (A-H,O-Z)                                         INT20320
      INTEGER QQ                                                        INT20330
      DIMENSION DIP(3),XKI(NC,3),A(NC,NC),V(NC)                         INT20340
      DIMENSION F1(NC),F2(NC,NC),F3(NC,NC,NC),F4(NC,NC,NC,NC)           INT20350
      DIMENSION AA(3,3),VK(3),XMASS(NA)                                 INT20360
      PARAMETER(ZERO=0.0D0)
      COMMON /IO/ IIN1,IOUT,IIN2,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25,
     $   I31,I32,I33,I35,I36,I37,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10,
     $   ISCR11,ISCR12,ISCR13,ISCR14
      COMMON /DIPOLE/DIP,QQ                                             INT20420
      COMMON /PHYCON/ BOHR,DEBYE,HART,WAVE0,CINT
      XMT=ZERO                                                          INT20430
      DO 10  I=1,NA                                                     INT20440
 10   XMT=XMT+XMASS(I)                                                  INT20450
      DO 15  I=1,NC                                                     INT20460
 15   V(I)=ZERO                                                         INT20470
      IF(NINV.NE.0) THEN                                                INT20480
      DO 100  I=1,NC                                                    INT20490
      K=(I-1)/3+1                                                       INT20500
 100  V(I)=DEBYE*QQ*XMASS(K)/XMT                                        INT20510
      ELSE                                                              INT20520
      DO 105  I=1,NSX                                                   INT20530
         XX=ZERO                                                        INT20540
         DO 110  J=1,NC                                                 INT20550
         K=(J-1)/3+1                                                    INT20560
 110     XX=XX+XMASS(K)*A(J,I)                                          INT20570
 105  V(I)=-DEBYE*QQ*XX/XMT                                             INT20580
      END IF                                                            INT20590
      IF(NINV.NE.0) THEN                                                INT20600
      DO 120  I=1,NC                                                    INT20610
         DO 125  KK=1,3                                                 INT20620
 125     VK(KK)=XKI(I,KK)                                               INT20630
         CALL MAT2(AA,VK)                                               INT20640
         DO 130 KK=1,3                                                  INT20650
 130     V(I)=V(I)+AA(NRUN,KK)*DIP(KK)                                  INT20660
 120  CONTINUE                                                          INT20670
      ELSE                                                              INT20680
      DO 140  I=1,NSX                                                   INT20690
         DO 145  KK=1,3                                                 INT20700
         VK(KK)=ZERO                                                    INT20710
         DO 145  J=1,NC                                                 INT20720
 145     VK(KK)=VK(KK)+A(J,I)*XKI(J,KK)                                 INT20730
         CALL MAT2(AA,VK)                                               INT20740
         DO 150 KK=1,3                                                  INT20750
 150     V(I)=V(I)-AA(NRUN,KK)*DIP(KK)                                  INT20760
 140  CONTINUE                                                          INT20770
      END IF                                                            INT20780
      DO 200 I=1,NC                                                     INT20790
 200  F1(I)=F1(I)+V(I)                                                  INT20800
      END                                                               INT20810
C     ////////////////////////////////////////////////////////////      INT14360
C  NOTE THE OVERLAP IN MEMORY OF THE XR AND SRX ARRAYS.                 DIS00010
      SUBROUTINE DISP(NAD,NA,NC,NS,NDISP,IOPT,                          DIS00010
     $        XA,XMASS,S,SS,B,BS,A,TYPE,U,IA,IU,                        DIS00020
     $        XD,XE,SD,XR,SRX,XT,IFLAG)                                 DIS00030
      IMPLICIT REAL*8 (A-H,O-Z)                                         DIS00040
      CHARACTER TYPE*5,LABEL*4                                          DIS00050
      DIMENSION XECK(6)                                                 DIS00060
      DIMENSION IOPT(30),TYPE(NS),IA(NS,6),IU(NS,0:1),S(NS),SS(NS)      DIS00070
      DIMENSION U(NS,1),XA(NAD,3),XMASS(1)                              DIS00080
      DIMENSION B(NS,NC),A(NC,NC),BS(NC,NC)                             DIS00090
      DIMENSION SRX(NC,NC),XR(NC,1),XT(NC,NC)                           DIS00100
      DIMENSION SD(NS),XD(NAD,3),XE(NA,3)                               DIS00110
CWAx
      DIMENSION JSPF(100),SPFREF(100)
      PARAMETER(RAD=57.29577951308232D0)
      PARAMETER(ZERO=0.0D0,ONE=1.0D0,HALF=0.50D0,SIX=6.0D0)             DIS00130
      PARAMETER(TOLINV=1.0D-10,TOLDISP=1.0D-14,MXITER=20)               DIS00140
      COMMON /IO/ IIN1,IOUT,IIN2,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25,
     $   I31,I32,I33,I35,I36,I37,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10,
     $   ISCR11,ISCR12,ISCR13,ISCR14
      COMMON /PHYCON/ BOHR,DEBYE,HART,WAVE0,CINT
CWAx
      COMMON /SPF1/NSPF,JSPF
      COMMON /SPF2/SPFREF 
  1   FORMAT(//'NUCLEAR MASSES'/)                                       DIS00200
  2   FORMAT(4(I4,F12.6,4X))                                            DIS00210
  3   FORMAT(A4,I4)                                                     DIS00220
  4   FORMAT(/,1X,'EXECUTION STOPPED DUE TO INPUT ERROR IN DISP.')      DIS00230
  5   FORMAT(///,1X,'VALUES OF SIMPLE INTERNAL COORDINATES ',           DIS00240
     $      '(ANG. OR DEG.) FOR REFERENCE GEOMETRY'/)                   DIS00250
  6   FORMAT(4(I4,F16.10))                                              DIS00260
  7   FORMAT(///,1X,'VALUES OF SYMMETRY INTERNAL COORDINATES ',         DIS00270
     $      '(ANG. OR RAD.) FOR REFERENCE GEOMETRY'/)                   DIS00280
  8   FORMAT(I5,F20.10)                                                 DIS00290
  9   FORMAT(//,1X,'INTERNAL (SYMMETRY) COORDINATE FINAL VALUES'/)      DIS00300
 10   FORMAT(I5,F20.10)                                                 DIS00310
 11   FORMAT(/,1X,' ITER = ',I5,'   MAX INTERNAL DEVIATION = ',E10.3)   DIS00320
 12   FORMAT(1X,'B MATRIX DETERMINANT = 0 .')                           DIS00330
 13   FORMAT(/' B MATRIX INVERSION UNSUCCESSFUL FOR GIVEN TOLERANCE: ', DIS00340
     $      E10.3)                                                      DIS00350
 14   FORMAT(/' NEW CARTESIAN GEOMETRY FOUND IN ',I5,' ITERATIONS.')    DIS00360
 15   FORMAT(/' CARTESIAN GEOMETRY NOT FOUND IN ',I5,' ITERATIONS.')    DIS00370
 16   FORMAT(/' NEW CARTESIAN GEOMETRY (BOHR)'/)                        DIS00380
 17   FORMAT(3F20.10)                                                   DIS00390
 20   FORMAT(//' INTERNAL DISPLACEMENT SUBROUTINE'/)                    DIS00400
 21   FORMAT(/' INTERNAL DISPLACEMENTS'/)                               DIS00410
 22   FORMAT(I5,F20.10)                                                 DIS00420
 23   FORMAT(15X,'  XDISP=',E10.3)                                      DIS00430
 25   FORMAT(/' DISPLACEMENT',I5)                                       DIS00440
 26   FORMAT('# GEOMUP #################')                              DIS00450
 27   FORMAT(3F20.10)                                                   DIS00460
 28   FORMAT(/'ECKART CONDITION VECTORS'/)                              DIS00470
 29   FORMAT(3F20.10)                                                   DIS00480
C                                                                       DIS00490
      IFLAG=0                                                           DIS00500
      NSYM=IOPT(3)                                                      DIS00510
      NSX=NS                                                            DIS00520
      IF(NSYM.GT.0) NSX=NSYM                                            DIS00530
      WRITE(IOUT,20)                                                    DIS00540
      IF(NDISP.LT.0) THEN                                               DIS00550
           WRITE(IOUT,1)                                                DIS00560
           WRITE(IOUT,2) (I,XMASS(I),I=1,NA)                            DIS00570
      END IF                                                            DIS00580
      READ(IIN1,3) LABEL,MDISP                                          DIS00590
      IF(LABEL.NE.'DISP') THEN                                          DIS00600
           WRITE(IOUT,4)                                                DIS00610
           IFLAG=1                                                      DIS00620
           RETURN                                                       DIS00630
      END IF                                                            DIS00640
C                                                                       DIS00650
      DO 1000 ID=1,MDISP                                                DIS00660
      ITER=0                                                            DIS00670
CWAx
      DO 102 I=1,NSPF
 102  S(JSPF(I))=SPFREF(I)
C
      CALL MACHB(NAD,NC,NS,XA,XMASS,TYPE,IA,B,S)                        DIS00680
      IF(NSYM.LE.0) THEN                                                DIS00690
         DO 100  I=1,NS                                                 DIS00700
 100     SD(I)=S(I)                                                     DIS00710
      ELSE                                                              DIS00720
         DO 110  I=1,NSYM                                               DIS00730
         SD(I)=ZERO                                                     DIS00740
         DO 110  L=1,IU(I,0)                                            DIS00750
         K=IU(I,L)                                                      DIS00760
 110     SD(I)=SD(I)+U(I,L)*S(K)                                        DIS00770
      END IF                                                            DIS00780
      DO 115 J=1,3                                                      DIS00790
      DO 115 I=1,NAD                                                    DIS00800
 115  XD(I,J)=XA(I,J)                                                   DIS00810
      IF(ID.NE.1) GO TO 150                                             DIS00820
CPRINT                                                                  DIS00830
      DO 120  I=1,NS                                                    DIS00840
      XR(I,1)=S(I)                                                      DIS00850
      IF(TYPE(I).EQ.' BEND'.OR.TYPE(I).EQ.' LIN1'.OR.TYPE(I).EQ.        DIS00860
     $   ' TORS'.OR.TYPE(I).EQ.'  OUT') THEN                            DIS00870
      XR(I,1)=XR(I,1)*RAD                                               DIS00880
      END IF                                                            DIS00890
 120  CONTINUE                                                          DIS00900
      WRITE(IOUT,5)                                                     DIS00910
      WRITE(IOUT,6) (I,XR(I,1),I=1,NS)                                  DIS00920
      IF(NSYM.GT.0) THEN                                                DIS00930
        WRITE(IOUT,7)                                                   DIS00940
        WRITE(IOUT,6) (I,SD(I),I=1,NSYM)                                DIS00950
      END IF                                                            DIS00960
C                                                                       DIS00970
 150  CONTINUE                                                          DIS00980
      WRITE(IOUT,25) ID                                                 DIS00990
      WRITE(IOUT,21)                                                    DIS01000
 130  READ(IIN1,8) IC,X                                                 DIS01010
      IF(IC.EQ.0) GO TO 135                                             DIS01020
      WRITE(IOUT,22) IC,X                                               DIS01030
      SD(IC)=SD(IC)+X                                                   DIS01040
      GO TO 130                                                         DIS01050
 135  WRITE(IOUT,9)                                                     DIS01060
      WRITE(IOUT,10) (I,SD(I),I=1,NSX)                                  DIS01070
C                                                                       DIS01080
 155  ITER=ITER+1                                                       DIS01090
CWAx
      DO 156 I=1,NSPF
 156  S(JSPF(I))=SPFREF(I)   
C
      CALL MACHB(NAD,NC,NS,XD,XMASS,TYPE,IA,B,S)                        DIS01100
      IF(NSYM.LE.0) THEN                                                DIS01110
         DO 160  I=1,NS                                                 DIS01120
         SS(I)=S(I)                                                     DIS01130
         DO 160  J=1,NC                                                 DIS01140
 160     BS(I,J)=B(I,J)                                                 DIS01150
      ELSE                                                              DIS01160
         DO 162  J=1,NC                                                 DIS01170
         DO 162  I=1,NSYM                                               DIS01180
           BS(I,J)=ZERO                                                 DIS01190
           DO 164  L=1,IU(I,0)                                          DIS01200
           K=IU(I,L)                                                    DIS01210
 164       BS(I,J)=BS(I,J)+U(I,L)*B(K,J)                                DIS01220
 162     CONTINUE                                                       DIS01230
         DO 166  I=1,NSYM                                               DIS01240
           SS(I)=ZERO                                                   DIS01250
           DO 166  L=1,IU(I,0)                                          DIS01260
           K=IU(I,L)                                                    DIS01270
 166       SS(I)=SS(I)+U(I,L)*S(K)                                      DIS01280
      END IF                                                            DIS01290
C                                                                       DIS01300
      ICONV=1                                                           DIS01310
      SDISP=ZERO                                                        DIS01320
      DO 170  I=1,NSX                                                   DIS01330
      X=DABS(SD(I)-SS(I))                                               DIS01340
 170  IF(X.GT.SDISP) SDISP=X                                            DIS01350
      IF(SDISP.LE.TOLDISP) ICONV=0                                      DIS01360
      WRITE(IOUT,11) ITER,SDISP                                         DIS01370
      IF(ICONV.EQ.0) GO TO 500                                          DIS01380
C                                                                       DIS01390
      IF(NDISP.LT.0) THEN                                               DIS01400
      DO 180 I=1,NA                                                     DIS01410
 180  XT(I,1)=XMASS(I)                                                  DIS01420
      ELSE                                                              DIS01430
      DO 182 I=1,NA                                                     DIS01440
 182  XT(I,1)=ONE                                                       DIS01450
      END IF                                                            DIS01460
      CALL BINVRT(NSX,NC,XT(1,1),BS,XR,A,SRX,
     $     XT(1,2),IFLAG,TOLINV,1)
      IF(IFLAG.EQ.1) THEN
         WRITE(IOUT,12)                                                 DIS01490
         RETURN                                                         DIS01500
      ELSE IF (IFLAG.EQ.2) THEN                                         DIS01510
         WRITE(IOUT,13) TOLINV                                          DIS01520
         RETURN                                                         DIS01530
      END IF                                                            DIS01540
C                                                                       DIS01550
      II=0                                                              DIS01560
      DO 200 I=1,NA                                                     DIS01570
      DO 200 J=1,3                                                      DIS01580
      II=II+1                                                           DIS01590
           XE(I,J)=ZERO                                                 DIS01600
           DO 205 K=1,NSX                                               DIS01610
 205       XE(I,J)=XE(I,J)+A(II,K)*(SD(K)-SS(K))                        DIS01620
 200  CONTINUE                                                          DIS01630
      IF(ABS(NDISP).EQ.1) GO TO 400                                     DIS01640
C                                                                       DIS01650
       CALL MACHX(NAD,NC,NS,NSX,IOPT,XD,XMASS,TYPE,IA,A,S,U,IU,XR,SRX)
       DO 300 K=1,NSX
       CALL XIN(NC,NC,SRX,-K,ISCR1)                                     DIS01680
          XR(K,1)=ZERO                                                  DIS01690
          DO 305 I=1,NC                                                 DIS01700
          II=(I-1)/3+1                                                  DIS01710
          JI=I-3*(II-1)                                                 DIS01720
          DO 305 J=1,NC                                                 DIS01730
          IJ=(J-1)/3+1                                                  DIS01740
          JJ=J-3*(IJ-1)                                                 DIS01750
 305      XR(K,1)=XR(K,1)+SRX(I,J)*XE(IJ,JJ)*XE(II,JI)                  DIS01760
 300      XR(K,1)=SD(K)-SS(K)-HALF*XR(K,1)                              DIS01770
C                                                                       DIS01780
 350   II=0                                                             DIS01790
       XDISP=ZERO                                                       DIS01800
       DO 360 I=1,NA                                                    DIS01810
       DO 360 J=1,3                                                     DIS01820
       II=II+1                                                          DIS01830
           X=ZERO                                                       DIS01840
           DO 365 K=1,NSX                                               DIS01850
 365       X=X+A(II,K)*XR(K,1)                                          DIS01860
           W=DABS(XE(I,J)-X)                                            DIS01870
           IF(W.GT.XDISP) XDISP=W                                       DIS01880
           XE(I,J)=X                                                    DIS01890
 360   CONTINUE                                                         DIS01900
       WRITE(IOUT,23) XDISP                                             DIS01910
C                                                                       DIS01920
 400  DO 405 J=1,3                                                      DIS01930
      DO 405 I=1,NA                                                     DIS01940
 405   XD(I,J)=XD(I,J)+XE(I,J)                                          DIS01950
      IF(ITER.LT.MXITER) GO TO 155                                      DIS01960
C   ROTATE COORDINATES, PRINT NEW CARTESIANS, WRITE OUT GEOM TO INPUT.  DIS01970
 500  IF(ICONV.EQ.0) THEN                                               DIS01980
         WRITE(IOUT,14) ITER                                            DIS01990
      ELSE                                                              DIS02000
         WRITE(IOUT,15) MXITER                                          DIS02010
         RETURN                                                         DIS02020
      END IF                                                            DIS02030
C                                                                       DIS02040
      IF(NDISP.GT.0) GO TO 510                                          DIS02050
      CALL ROTC(NA,NAD,XMASS,XA,XD,IFLAG)                               DIS02060
C                                                                       DIS02070
 510  DO 525 I=1,6                                                      DIS02080
 525  XECK(I)=ZERO                                                      DIS02090
      DO 530 I=1,NA                                                     DIS02100
      DO 535 J=1,3                                                      DIS02110
 535  XECK(J)=XECK(J)+XMASS(I)*(XD(I,J)-XA(I,J))                        DIS02120
      XECK(4)=XECK(4)+XMASS(I)*(XD(I,2)*XA(I,3)-XD(I,3)*XA(I,2))        DIS02130
      XECK(5)=XECK(5)+XMASS(I)*(XD(I,3)*XA(I,1)-XD(I,1)*XA(I,3))        DIS02140
 530  XECK(6)=XECK(6)+XMASS(I)*(XD(I,1)*XA(I,2)-XD(I,2)*XA(I,1))        DIS02150
      WRITE(IOUT,28)                                                    DIS02160
      WRITE(IOUT,29) (XECK(I),I=1,6)                                    DIS02170
C                                                                       DIS02180
      DO 605  J=1,3                                                     DIS02190
      DO 605  I=1,NA                                                    DIS02200
 605  XD(I,J)=XD(I,J)/BOHR                                              DIS02210
      WRITE(IOUT,16)                                                    DIS02220
      DO 610 I=1,NA                                                     DIS02230
 610  WRITE(IOUT,17) (XD(I,J),J=1,3)                                    DIS02240
      WRITE(IIN2,26)                                                    DIS02250
      DO 620 I=1,NA                                                     DIS02260
 620  WRITE(IIN2,27) (XD(I,J),J=1,3)                                    DIS02270
C                                                                       DIS02280
 1000  CONTINUE                                                         DIS02290
       RETURN                                                           DIS02300
       END                                                              DIS02310
C     //////////////////////////////////////////////////////////////////DIS02320
      SUBROUTINE ROTC(NA,NAD,XMASS,XA,XD,IFLAG)                         DIS02330
      IMPLICIT REAL*8 (A-H,O-Z)                                         DIS02340
      DIMENSION XMASS(NA),XA(NAD,3),XD(NAD,3)                           DIS02350
      DIMENSION W(3,3),XX(3,6),B(3),RK(3),RCM(3)                        DIS02360
      PARAMETER(ZERO=0.0D0,ONE=1.0D0)                                   DIS02370
      PARAMETER(MAXIT=40,TOLDISP=1.0D-12)                               DIS02380
      COMMON /IO/ IIN1,IOUT,IIN2,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25,
     $   I31,I32,I33,I35,I36,I37,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10,
     $   ISCR11,ISCR12,ISCR13,ISCR14
 1    FORMAT(/'ROTATION OF COORDINATES TO SATISFY ECKART CONDITIONS'/)  DIS02430
 2    FORMAT('CONVERGENCE NOT REACHED IN ROTC')                         DIS02440
 3    FORMAT('XDISP= ',E10.3,' TOLDISP= ',E10.3)                        DIS02450
 4    FORMAT('K VECTOR'/3E12.5)                                         DIS02460
 5    FORMAT('CENTER OF MASS'/3F20.10)                                  DIS02470
      ITER=0                                                            DIS02480
      IFLAG=0                                                           DIS02490
      RMT=ZERO                                                          DIS02500
      DO 15  I=1,NA                                                     DIS02510
 15   RMT=RMT+XMASS(I)                                                  DIS02520
      DO 25  J=1,3                                                      DIS02530
      RCM(J)=ZERO                                                       DIS02540
      DO 20  I=1,NA                                                     DIS02550
 20   RCM(J)=RCM(J)+XA(I,J)*XMASS(I)                                    DIS02560
 25   RCM(J)=RCM(J)/RMT                                                 DIS02570
      DO 30  I=1,NA                                                     DIS02580
      DO 30  J=1,3                                                      DIS02590
      XA(I,J)=XA(I,J)-RCM(J)                                            DIS02600
 30   XD(I,J)=XD(I,J)-RCM(J)                                            DIS02610
      WRITE(IOUT,1)                                                     DIS02620
      WRITE(IOUT,5) (RCM(I),I=1,3)                                      DIS02630
C                                                                       DIS02640
  50  ITER=ITER+1                                                       DIS02650
      DO 100  I=1,3                                                     DIS02660
      B(I)=ZERO                                                         DIS02670
      DO 100  J=1,3                                                     DIS02680
      XX(I,J+3)=ZERO                                                    DIS02690
 100  W(I,J)=ZERO                                                       DIS02700
      DO 110  I=1,3                                                     DIS02710
      DO 110  J=I,3                                                     DIS02720
      DO 110  K=1,NA                                                    DIS02730
 110  W(I,J)=W(I,J)+XA(K,I)*XD(K,J)*XMASS(K)                            DIS02740
      XX(1,1)=W(3,1)                                                    DIS02750
      XX(2,1)=W(3,2)                                                    DIS02760
      XX(3,1)=-(W(1,1)+W(2,2))                                          DIS02770
      XX(1,2)=-W(2,1)                                                   DIS02780
      XX(2,2)=W(1,1)+W(3,3)                                             DIS02790
      XX(3,2)=-W(2,3)                                                   DIS02800
      XX(1,3)=-(W(2,2)+W(3,3))                                          DIS02810
      XX(2,3)=W(1,2)                                                    DIS02820
      XX(3,3)=W(1,3)                                                    DIS02830
      XX(1,4)=ONE                                                       DIS02840
      XX(2,5)=ONE                                                       DIS02850
      XX(3,6)=ONE                                                       DIS02860
      CALL FLIN(XX,3,3,3,DD)                                            DIS02870
      DO 120  I=1,NA                                                    DIS02880
      B(1)=B(1)+XMASS(I)*(XD(I,2)*XA(I,3)-XD(I,3)*XA(I,2))              DIS02890
      B(2)=B(2)+XMASS(I)*(XD(I,3)*XA(I,1)-XD(I,1)*XA(I,3))              DIS02900
 120  B(3)=B(3)+XMASS(I)*(XD(I,1)*XA(I,2)-XD(I,2)*XA(I,1))              DIS02910
      DO 130  I=1,3                                                     DIS02920
      RK(I)=ZERO                                                        DIS02930
      DO 130  J=1,3                                                     DIS02940
 130  RK(I)=RK(I)+XX(I,J+3)*B(J)                                        DIS02950
      CALL EXPMAT(RK,W)                                                 DIS02960
      XDISP=ZERO                                                        DIS02970
      DO 200  I=1,NA                                                    DIS02980
      DO 205  J=1,3                                                     DIS02990
      XY=ZERO                                                           DIS03000
      DO 210  K=1,3                                                     DIS03010
  210 XY=XY+W(J,K)*XD(I,K)                                              DIS03020
      XZ=DABS(XY-XD(I,J))                                               DIS03030
      IF(XZ.GT.XDISP) XDISP=XZ                                          DIS03040
  205 B(J)=XY                                                           DIS03050
      DO 215  J=1,3                                                     DIS03060
  215 XD(I,J)=B(J)                                                      DIS03070
  200 CONTINUE                                                          DIS03080
      IF(XDISP.GT.TOLDISP.AND.ITER.LT.MAXIT) GO TO 50                   DIS03090
      DO 300 I=1,NA                                                     DIS03100
      DO 300 J=1,3                                                      DIS03110
      XA(I,J)=XA(I,J)+RCM(J)                                            DIS03120
 300  XD(I,J)=XD(I,J)+RCM(J)                                            DIS03130
      IF(XDISP.GT.TOLDISP) THEN                                         DIS03140
         WRITE(IOUT,2)                                                  DIS03150
         WRITE(IOUT,3) XDISP,TOLDISP                                    DIS03160
         WRITE(IOUT,4) (RK(I), I=1,3)                                   DIS03170
         IFLAG=1                                                        DIS03180
         RETURN                                                         DIS03190
      END IF                                                            DIS03200
      RETURN                                                            DIS03210
      END                                                               DIS03220
C////////////////////////////////////////////////////////////////////   INT35890
      SUBROUTINE INPFKM(NC,NDER,NEQ,NS,F1,F2,F3,F4)                     INT35900
      IMPLICIT REAL*8 (A-H,O-Z)                                         INT35910
      INTEGER P,Q,R,S                                                   INT35920
      DIMENSION F1(NC),F2(NC,NC),F3(NC,NC,NC),F4(NC,NC,NC,NC)           INT35930
      PARAMETER(ZERO=0.0D0)                                             INT35940
      COMMON /IO/ IIN1,IOUT,IIN2,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25,
     $   I31,I32,I33,I35,I36,I37,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10,
     $   ISCR11,ISCR12,ISCR13,ISCR14
 1    FORMAT(5I5)                                                       INT35990
 2    FORMAT(3F20.10)                                                   INT36000
 3    FORMAT(20X,F20.10)                                                INT36010
 5    FORMAT(20X,F20.10)                                                INT36020
      NA=NC/3                                                           INT36030
      IF(NEQ.EQ.0) GO TO 110                                            INT36040
           DO 108  I=1,NS                                               INT36050
 108       F1(I)=ZERO                                                   INT36060
      IF(NDER.LE.1) GO TO 190                                           INT36070
 110  DO 112  I=1,NS                                                    INT36080
      DO 112  J=1,NS                                                    INT36090
 112  F2(I,J)=ZERO                                                      INT36100
      IF(NDER.LE.2) GO TO 190                                           INT36110
      DO 114  I=1,NS                                                    INT36120
      DO 114  J=1,I                                                     INT36130
      DO 114  K=1,J                                                     INT36140
 114  F3(I,J,K)=ZERO                                                    INT36150
      IF(NDER.LE.3) GO TO 190                                           INT36160
      DO 116  I=1,NS                                                    INT36170
      DO 116  J=1,I                                                     INT36180
      DO 116  K=1,J                                                     INT36190
      DO 116  L=1,K                                                     INT36200
 116  F4(I,J,K,L)=ZERO                                                  INT36210
C                                                                       INT36220
 190  IF(NEQ.EQ.0) GO TO 200                                            INT36230
 192  READ(IIN1,1) M                                                    INT36240
      IF(M.EQ.0) GO TO 200                                              INT36250
      BACKSPACE IIN1                                                    INT36260
      READ(IIN1,5) F1(M)                                                INT36270
      GO TO 192                                                         INT36280
 200  READ(IIN1,1) M,N                                                  INT36290
      IF(M.EQ.0) GO TO 202                                              INT36300
      BACKSPACE IIN1                                                    INT36310
      READ(IIN1,5) F2(M,N)                                              INT36320
      GO TO 200                                                         INT36330
 202  IF(NDER.LE.2) GO TO 212                                           INT36340
 204  READ(IIN1,1) M,N,P                                                INT36350
      IF(M.EQ.0) GO TO 206                                              INT36360
      BACKSPACE IIN1                                                    INT36370
      READ(IIN1,5) F3(M,N,P)                                            INT36380
      GO TO 204                                                         INT36390
 206  IF(NDER.LE.3) GO TO 212                                           INT36400
 208  READ(IIN1,1) M,N,P,Q                                              INT36410
      IF(M.EQ.0) GO TO 212                                              INT36420
      BACKSPACE IIN1                                                    INT36430
      READ(IIN1,5) F4(M,N,P,Q)                                          INT36440
      GO TO 208                                                         INT36450
 212  CONTINUE                                                          INT36460
C                                                                       INT36470
      IF(NEQ.NE.0) THEN                                                 INT36480
      WRITE(I12,1) NA                                                   INT36490
      WRITE(I12,3) (F1(I),I=1,NS)                                       INT36500
      REWIND I12                                                        INT36510
      END IF                                                            INT36520
      DO 214  I=1,NS                                                    INT36530
      DO 214  J=1,I                                                     INT36540
 214  F2(J,I)=F2(I,J)                                                   INT36550
      WRITE (I16,1) NA,NA*6                                             INT36560
      WRITE (I16,2) ((F2(I,J),J=1,NS),I=1,NS)                           INT36570
      REWIND I16                                                        INT36580
      IF(NDER.LE.2) RETURN                                              INT36590
      WRITE (I21,1) NA,(NS*(NS+1)*(NS+2))/6                             INT36600
      WRITE (I21,2) (((F3(I,J,K),K=1,J),J=1,I),I=1,NS)                  INT36610
      REWIND I21                                                        INT36620
      IF(NDER.LE.3) RETURN                                              INT36630
      WRITE (I25,1) NA,(NS*(NS+1)*(NS+2)*(NS+3))/24                     INT36640
      WRITE (I25,2) ((((F4(I,J,K,L),L=1,K),K=1,J),J=1,I),I=1,NS)        INT36650
      REWIND I25                                                        INT36660
      RETURN                                                            INT36670
      END                                                               INT36680
C     ///////////////////////////////////////////////////////////       INT36690
      SUBROUTINE FCOUT(NC,NS,NSY,NEQ,NDER,NINV,NVEC,NRUN,SS,ENERGY,     INT36700
     $                 F1,F2,F3,F4)                                     INT36710
      IMPLICIT REAL*8 (A-H,O-Z)                                         INT36720
      INTEGER P,Q,QQ                                                    INT36730
      DIMENSION SS(NS)                                                  INT36740
      DIMENSION F1(NC),F2(NC,NC),F3(NC,NC,NC),F4(NC,NC,NC,NC)           INT36750
      DIMENSION DIP(3)                                                  INT36760
CSA   PARAMETER(CF1=0.121376050883D0,CF2=0.064229421739D0)              INT36770
CSA   PARAMETER(CF3=0.033988736554D0,CF4=0.017986059681D0)              INT36780
      PARAMETER(ONE=1.0D0)                                              INT36790
      PARAMETER(F1CUT=1.0D-8,F2CUT=1.0D-6,F3CUT=1.0D-5,F4CUT=1.0D-5)    INT36800
      COMMON /IO/ IIN1,IOUT,IIN2,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25,
     $   I31,I32,I33,I35,I36,I37,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10,
     $   ISCR11,ISCR12,ISCR13,ISCR14
      COMMON /PHYCON/ BOHR,DEBYE,HART,WAVE0,CINT
      COMMON /DIPOLE/ DIP,QQ
 1    FORMAT(/' GRADIENTS')                                             INT36860
 2    FORMAT(6F12.8)                                                    INT36870
 3    FORMAT(/' QUADRATIC FORCE CONSTANTS'/)                            INT36880
 4    FORMAT(6F12.6)                                                    INT36890
 5    FORMAT(/' CUBIC FORCE CONSTANTS'/)                                INT36900
 6    FORMAT('M=',I5)                                                   INT36910
 7    FORMAT(/' QUARTIC FORCE CONSTANTS'/)                              INT36920
 8    FORMAT('M=',I5,'  N=',I5)                                         INT36930
 9    FORMAT(//' FINAL INTERNAL (SYMMETRY) COORDINATE FORCE CONSTANTS') INT36940
 10   FORMAT(7F10.6)                                                    INT36950
 11   FORMAT(8F10.7)                                                    INT36960
 12   FORMAT(//' FINAL CARTESIAN COORDINATE FORCE CONSTANTS')           INT36970
 13   FORMAT(3F20.10)                                                   INT36980
 14   FORMAT(5I5)                                                       INT36990
 15   FORMAT(20X,3F20.10)                                               INT37000
 16   FORMAT(3F20.10)                                                   INT37010
 17   FORMAT(I5,12X,F20.10)                                             INT37020
 18   FORMAT(//' PRINTING OF CUBIC FORCE CONSTANTS',                    INT37030
     $         ' TO OUTPUT FILE SUPPRESSED.'/' LPRT(1,NPRT)=',I3)       INT37040
 19   FORMAT(//' PRINTING OF QUARTIC FORCE CONSTANTS',                  INT37050
     $         ' TO OUTPUT FILE SUPPRESSED.'/' LPRT(1,NPRT)=',I3)       INT37060
 20   FORMAT(I5,15X,F20.10)                                             INT37070
 21   FORMAT(2I5,10X,F20.10)                                            INT37080
 22   FORMAT(3I5,5X,F20.10)                                             INT37090
 23   FORMAT(4I5,F20.10)                                                INT37100
 24   FORMAT(2I5,3F20.10)                                               INT37110
      CF1=BOHR/HART
      CF2=CF1*BOHR
      CF3=CF2*BOHR
      CF4=CF3*BOHR
      NA=NC/3                                                           INT37120
      NNINV=NINV
CSA
      IF(ABS(NINV).NE.3) GO TO 55
         NNINV=0
         WRITE(IOUT,12)                                                 INT37210
         IOUT1=I12                                                      INT37150
         IOUT2=I16                                                      INT37160
         IOUT3=I21                                                      INT37170
         IOUT4=I25                                                      INT37180
      GO TO 60
CSA
 55   IF(NNINV.LE.0) THEN                                               INT37130
         WRITE(IOUT,9)                                                  INT37140
         IOUT1=I12                                                      INT37150
         IOUT2=I16                                                      INT37160
         IOUT3=I21                                                      INT37170
         IOUT4=I25                                                      INT37180
         IOUT5=I18                                                      INT37190
      ELSE                                                              INT37200
         WRITE(IOUT,12)                                                 INT37210
         IOUT1=I11                                                      INT37220
         IOUT2=I15                                                      INT37230
         IOUT3=I20                                                      INT37240
         IOUT4=I24                                                      INT37250
         IOUT5=I17                                                      INT37260
      END IF                                                            INT37270
  60  IF(NEQ.NE.0) THEN                                                 INT37280
         IF(NNINV.LE.0) THEN                                            INT37290
             XX=ONE                                                     INT37300
             IF(NVEC.EQ.1) XX=CF1                                       INT37310
         ELSE                                                           INT37320
             XX=CF1                                                     INT37330
             IF(NVEC.EQ.1) XX=ONE                                       INT37340
         END IF                                                         INT37350
         DO 200 M=1,NSY                                                 INT37360
 200     F1(M)=F1(M)*XX                                                 INT37370
      END IF                                                            INT37380
      IF(NNINV.LE.0) GO TO 220                                          INT37390
      IF(NDER.GE.2) THEN                                                INT37400
      DO 205  N=1,NSY                                                   INT37410
      DO 205  M=1,NSY                                                   INT37420
 205  F2(M,N)=F2(M,N)*CF2                                               INT37430
      END IF                                                            INT37440
      IF(NDER.GE.3) THEN                                                INT37450
      DO 210  P=1,NSY                                                   INT37460
      DO 210  N=1,NSY                                                   INT37470
      DO 210  M=1,NSY                                                   INT37480
 210  F3(M,N,P)=F3(M,N,P)*CF3                                           INT37490
      END IF                                                            INT37500
      IF(NDER.GE.4) THEN                                                INT37510
      DO 215  Q=1,NSY                                                   INT37520
      DO 215  P=1,NSY                                                   INT37530
      DO 215  N=1,NSY                                                   INT37540
      DO 215  M=1,NSY                                                   INT37550
 215  F4(M,N,P,Q)=F4(M,N,P,Q)*CF4                                       INT37560
      END IF                                                            INT37570
 220  CONTINUE                                                          INT37580
      IF(NEQ.NE.0) THEN                                                 INT37590
           WRITE(IOUT,1)                                                INT37600
           WRITE(IOUT,2) (F1(M),M=1,NSY)                                INT37610
           IF(NVEC.EQ.1) THEN                                           INT37620
                IF(NRUN.LE.1) WRITE(IOUT5,24) NA,QQ,(DIP(KK),KK=1,3)    INT37630
                WRITE(IOUT5,16) (F1(M),M=1,NSY)                         INT37640
                GO TO 225                                               INT37650
           END IF                                                       INT37660
           IF(NNINV.GT.0) THEN                                          INT37670
                WRITE(IOUT1,15) (F1(M),M=1,NSY)                         INT37680
           ELSE                                                         INT37690
                WRITE(IOUT1,17) NA,ENERGY                               INT37700
                DO 230  M=1,NSY                                         INT37710
 230            WRITE(IOUT1,13) SS(M),F1(M)                             INT37720
           END IF                                                       INT37730
      END IF                                                            INT37740
 225  IF(NDER.GE.2) THEN                                                INT37750
           WRITE(IOUT,3)                                                INT37760
           DO 100  M=1,NSY                                              INT37770
 100       WRITE(IOUT,4) (F2(M,N),N=1,NSY)                              INT37780
           WRITE (IOUT2,14) NA,NA*6                                     INT37790
           WRITE (IOUT2,13) ((F2(M,N),N=1,NSY),M=1,NSY)                 INT37800
           IF(LPRT(4,NPRT).GE.4)                                        INT37810
     $             WRITE(ICHECK,10)((F2(M,N),M=N,NSY),N=1,NSY)          INT37820
           IF(LPRT(4,NPRT).GE.3) THEN                                   INT37830
           DO 102  M=1,NSY                                              INT37840
 102       WRITE(ICHECK,11)(F2(M,N),N=1,M)                              INT37850
           END IF                                                       INT37860
      END IF                                                            INT37870
      IF(NDER.GE.3) THEN                                                INT37880
           LLL=LPRT(1,NPRT)
           IF(LLL.GE.1) THEN                                            INT37890
             WRITE(IOUT,5)                                              INT37900
             DO 105  M=1,NSY                                            INT37910
             WRITE(IOUT,6) M                                            INT37920
             DO 105  N=1,NSY                                            INT37930
 105         WRITE(IOUT,4) (F3(M,N,P),P=1,NSY)                          INT37940
           ELSE                                                         INT37950
             WRITE(IOUT,18) LPRT(1,NPRT)                                INT37960
           END IF                                                       INT37970
           WRITE (IOUT3,14) NA,(NSY*(NSY+1)*(NSY+2))/6                  INT37980
           WRITE (IOUT3,13) (((F3(M,N,P),P=1,N),N=1,M),M=1,NSY)         INT37990
      END IF                                                            INT38000
      IF(NDER.GE.4) THEN                                                INT38010
           IF(LPRT(1,NPRT).GE.1) THEN                                   INT38020
             WRITE(IOUT,7)                                              INT38030
             DO 110  M=1,NSY                                            INT38040
             DO 110  N=1,NSY                                            INT38050
             WRITE(IOUT,8) M,N                                          INT38060
             DO 110  P=1,NSY                                            INT38070
 110         WRITE(IOUT,4) (F4(M,N,P,Q),Q=1,NSY)                        INT38080
           ELSE                                                         INT38090
             WRITE(IOUT,19) LPRT(1,NPRT)                                INT38100
           END IF                                                       INT38110
      WRITE (IOUT4,14) NA,(NSY*(NSY+1)*(NSY+2)*(NSY+3))/24              INT38120
      WRITE (IOUT4,13) ((((F4(M,N,P,Q),Q=1,P),P=1,N),N=1,M),M=1,NSY)    INT38130
      END IF                                                            INT38140
      IF(LPRT(4,NPRT).NE.2) RETURN                                      INT38150
      NZ=0                                                              INT38160
      IF(NEQ.NE.0) THEN                                                 INT38170
         DO 300  M=1,NSY                                                INT38180
       IF(DABS(F1(M)).GT.F1CUT) WRITE(ICHECK,20) M,F1(M)                INT38190
 300     CONTINUE                                                       INT38200
         WRITE(ICHECK,20) NZ                                            INT38210
      END IF                                                            INT38220
      IF(NDER.GE.2) THEN                                                INT38230
         DO 310  M=1,NSY                                                INT38240
         DO 310  N=1,M                                                  INT38250
       IF(DABS(F2(M,N)).GT.F2CUT) WRITE(ICHECK,21) M,N,F2(M,N)          INT38260
 310     CONTINUE                                                       INT38270
         WRITE(ICHECK,20) NZ                                            INT38280
      END IF                                                            INT38290
      IF(NDER.GE.3) THEN                                                INT38300
         DO 320  M=1,NSY                                                INT38310
         DO 320  N=1,M                                                  INT38320
         DO 320  P=1,N                                                  INT38330
       IF(DABS(F3(M,N,P)).GT.F3CUT) WRITE(ICHECK,22) M,N,P,F3(M,N,P)    INT38340
 320     CONTINUE                                                       INT38350
         WRITE(ICHECK,20) NZ                                            INT38360
      END IF                                                            INT38370
      IF(NDER.GE.4) THEN                                                INT38380
         DO 330  M=1,NSY                                                INT38390
         DO 330  N=1,M                                                  INT38400
         DO 330  P=1,N                                                  INT38410
         DO 330  Q=1,P                                                  INT38420
       IF(DABS(F4(M,N,P,Q)).GT.F4CUT)                                   INT38430
     $     WRITE(ICHECK,23) M,N,P,Q,F4(M,N,P,Q)                         INT38440
 330     CONTINUE                                                       INT38450
         WRITE(ICHECK,20) NZ                                            INT38460
      END IF                                                            INT38470
      RETURN                                                            INT38480
      END                                                               INT38490
C     //////////////////////////////////////////////////////////////    INT34200
      SUBROUTINE XOUT(NC,NS,X,R,ISCR)                                   INT34210
      REAL*8 X(NC,NC)                                                   INT34220
      INTEGER R                                                         INT34230
      COMMON /IO/ IIN1,IOUT,IIN2,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25,
     $   I31,I32,I33,I35,I36,I37,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10,
     $   ISCR11,ISCR12,ISCR13,ISCR14
    1 FORMAT(1X,'SR MATRIX FOR COORDINATE ',I4,' WRITTEN ON FILE ',I4,  INT34280
     $    ' AT LOCATION ',I8)                                           INT34290
    2 FORMAT(1X,'X  MATRIX FOR COORDINATE ',I4,' WRITTEN ON FILE ',I4,  INT34300
     $    ' AT LOCATION ',I8)                                           INT34310
      NLEN=NC**2                                                        INT34320
      NNN=2*NLEN                                                        INT34330
      MMM=(ABS(R)-1)*NNN                                                INT34340
      IF(R.GT.0) MMM=MMM+NLEN                                           INT34350
      LLL=INTOWP(MMM)+1                                                 INT34360
      CALL WWRITW(ISCR,X,INTOWP(NLEN),LLL,KKK)                          INT34370
      IF(LPRT(4,NPRT).GE.1) RETURN                                      INT34380
      IF(R.LT.0) THEN                                                   INT34390
      WRITE(ICHECK,1) ABS(R),ISCR,LLL                                   INT34400
      ELSE                                                              INT34410
      WRITE(ICHECK,2) ABS(R),ISCR,LLL                                   INT34420
      END IF                                                            INT34430
      RETURN                                                            INT34440
      END                                                               INT34450
C     ////////////////////////////////////////////////////////////////  INT34460
      SUBROUTINE XIN(NC,NS,X,R,ISCR)                                    INT34470
      REAL*8 X(NC,NC)                                                   INT34480
      INTEGER R                                                         INT34490
      COMMON /IO/ IIN1,IOUT,IIN2,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25,
     $   I31,I32,I33,I35,I36,I37,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10,
     $   ISCR11,ISCR12,ISCR13,ISCR14
    1 FORMAT(1X,'SR MATRIX FOR COORDINATE ',I4,' READ FROM  FILE ',I4,  INT34540
     $    ' AT LOCATION ',I8)                                           INT34550
    2 FORMAT(1X,'X  MATRIX FOR COORDINATE ',I4,' READ FROM  FILE ',I4,  INT34560
     $    ' AT LOCATION ',I8)                                           INT34570
      NLEN=NC**2                                                        INT34580
      NNN=2*NLEN                                                        INT34590
      MMM=(ABS(R)-1)*NNN                                                INT34600
      IF(R.GT.0) MMM=MMM+NLEN                                           INT34610
      LLL=INTOWP(MMM)+1                                                 INT34620
      CALL WREADW(ISCR,X,INTOWP(NLEN),LLL,KKK)                          INT34630
      IF(LPRT(4,NPRT).GE.1) RETURN                                      INT34640
      IF(R.LT.0) THEN                                                   INT34650
      WRITE(ICHECK,1) ABS(R),ISCR,LLL                                   INT34660
      ELSE                                                              INT34670
      WRITE(ICHECK,2) ABS(R),ISCR,LLL                                   INT34680
      END IF                                                            INT34690
      RETURN                                                            INT34700
      END                                                               INT34710
C     ////////////////////////////////////////////////////////////////  INT34720
      SUBROUTINE YOUT(NC,NS,Y,R,ISCR)                                   INT34730
      REAL*8 Y(NC,NC,NC)                                                INT34740
      INTEGER R                                                         INT34750
      COMMON /IO/ IIN1,IOUT,IIN2,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25,
     $   I31,I32,I33,I35,I36,I37,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10,
     $   ISCR11,ISCR12,ISCR13,ISCR14
    1 FORMAT(1X,'SR MATRIX FOR COORDINATE ',I4,' WRITTEN ON FILE ',I4,  INT34800
     $    ' AT LOCATION ',I8)                                           INT34810
    2 FORMAT(1X,'Y  MATRIX FOR COORDINATE ',I4,' WRITTEN ON FILE ',I4,  INT34820
     $    ' AT LOCATION ',I8)                                           INT34830
      NLEN=NC**3                                                        INT34840
      NNN=2*NLEN                                                        INT34850
      MMM=(ABS(R)-1)*NNN                                                INT34860
      IF(R.GT.0) MMM=MMM+NLEN                                           INT34870
      LLL=INTOWP(MMM)+1                                                 INT34880
      CALL WWRITW(ISCR,Y,INTOWP(NLEN),LLL,KKK)                          INT34890
      IF(LPRT(4,NPRT).GE.1) RETURN                                      INT34900
      IF(R.LT.0) THEN                                                   INT34910
      WRITE(ICHECK,1) ABS(R),ISCR,LLL                                   INT34920
      ELSE                                                              INT34930
      WRITE(ICHECK,2) ABS(R),ISCR,LLL                                   INT34940
      END IF                                                            INT34950
      RETURN                                                            INT34960
      END                                                               INT34970
C     ////////////////////////////////////////////////////////////////  INT34720
      SUBROUTINE YIN(NC,NS,Y,R,ISCR)                                    INT34990
      REAL*8 Y(NC,NC,NC)                                                INT35000
      INTEGER R                                                         INT35010
      COMMON /IO/ IIN1,IOUT,IIN2,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25,
     $   I31,I32,I33,I35,I36,I37,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10,
     $   ISCR11,ISCR12,ISCR13,ISCR14
    1 FORMAT(1X,'SR MATRIX FOR COORDINATE ',I4,' READ FROM  FILE ',I4,  INT35060
     $    ' AT LOCATION ',I8)                                           INT35070
    2 FORMAT(1X,'Y  MATRIX FOR COORDINATE ',I4,' READ FROM  FILE ',I4,  INT35080
     $    ' AT LOCATION ',I8)                                           INT35090
      NLEN=NC**3                                                        INT35100
      NNN=2*NLEN                                                        INT35110
      MMM=(ABS(R)-1)*NNN                                                INT35120
      IF(R.GT.0) MMM=MMM+NLEN                                           INT35130
      LLL=INTOWP(MMM)+1                                                 INT35140
      CALL WREADW(ISCR,Y,INTOWP(NLEN),LLL,KKK)                          INT35150
      IF(LPRT(4,NPRT).GE.1) RETURN                                      INT35160
      IF(R.LT.0) THEN                                                   INT35170
      WRITE(ICHECK,1) ABS(R),ISCR,LLL                                   INT35180
      ELSE                                                              INT35190
      WRITE(ICHECK,2) ABS(R),ISCR,LLL                                   INT35200
      END IF                                                            INT35210
      RETURN                                                            INT35220
      END                                                               INT35230
C///////////////////////////////////////////////////////////////////////INT35240
      SUBROUTINE YIN2(NC,NS,Y,R,NB,ISCR)                                YY 02300
      REAL*8 Y(NC,NC,NC)                                                YY 02310
      INTEGER R,NB,LR                                                   YY 02320
      COMMON /IO/ IIN1,IOUT,IIN2,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25,
     $   I31,I32,I33,I35,I36,I37,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10,
     $   ISCR11,ISCR12,ISCR13,ISCR14
    1 FORMAT(1X,'YR BLOCK',I4,' FOR R=',I4,' READ ON FILE ',I4,         YY 02370
     $    ' AT LOCATION ',I8)                                           YY 02380
    2 FORMAT(1X,'R = 0 ENCOUNTERED IN SUBROUTINE YIN2.')                YY 02390
      IF(R.LT.0) THEN                                                   YY 02400
         LR=-2*(R+1)*NC+NB                                              YY 02410
      ELSE IF(R.GT.0) THEN                                              YY 02420
         LR=(2*R-1)*NC+NB                                               YY 02430
      ELSE                                                              YY 02440
         WRITE(ICHECK,2)                                                YY 02450
         STOP                                                           YY 02460
      END IF                                                            YY 02470
      NLEN=NC**3                                                        YY 02480
      MMM=(LR-1)*NLEN                                                   YY 02490
      LLL=INTOWP(MMM)+1                                                 YY 02500
      CALL WREADW(ISCR,Y,INTOWP(NLEN),LLL,KKK)                          YY 02510
      IF(LPRT(4,NPRT).GE.1) RETURN                                      YY 02520
      WRITE(ICHECK,1) NB,R,ISCR,LLL                                     YY 02530
      RETURN                                                            YY 02540
      END                                                               YY 02550
C     ///////////////////////////////////////////////////////////////// YY 02560
      SUBROUTINE YOUT2(NC,NS,Y,R,NB,ISCR)                               YY 02570
      REAL*8 Y(NC,NC,NC)                                                YY 02580
      INTEGER R,NB,LR                                                   YY 02590
      COMMON /IO/ IIN1,IOUT,IIN2,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25,
     $   I31,I32,I33,I35,I36,I37,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10,
     $   ISCR11,ISCR12,ISCR13,ISCR14
    1 FORMAT(1X,'YR BLOCK',I4,' FOR R=',I4,' WRITTEN ON FILE ',I4,      YY 02640
     $    ' AT LOCATION ',I8)                                           YY 02650
    2 FORMAT(1X,'R = 0 ENCOUNTERED IN SUBROUTINE YOUT2.')               YY 02660
      IF(R.LT.0) THEN                                                   YY 02670
         LR=-2*(R+1)*NC+NB                                              YY 02680
      ELSE IF(R.GT.0) THEN                                              YY 02690
         LR=(2*R-1)*NC+NB                                               YY 02700
      ELSE                                                              YY 02710
         WRITE(ICHECK,2)                                                YY 02720
         STOP                                                           YY 02730
      END IF                                                            YY 02740
      NLEN=NC**3                                                        YY 02750
      MMM=(LR-1)*NLEN                                                   YY 02760
      LLL=INTOWP(MMM)+1                                                 YY 02770
      CALL WWRITW(ISCR,Y,INTOWP(NLEN),LLL,KKK)                          YY 02780
      IF(LPRT(4,NPRT).GE.1) RETURN                                      YY 02790
      WRITE(ICHECK,1) NB,R,ISCR,LLL                                     YY 02800
      RETURN                                                            YY 02810
      END                                                               YY 02820
C     ///////////////////////////////////////////////////////////////// YY 02830
      SUBROUTINE ZOUT(NC,NS,Z,R,ISCR)                                   INT35250
      REAL*8 Z(NC,NC,NC,NC)                                             INT35260
      INTEGER R                                                         INT35270
      COMMON /IO/ IIN1,IOUT,IIN2,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25,
     $   I31,I32,I33,I35,I36,I37,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10,
     $   ISCR11,ISCR12,ISCR13,ISCR14
    1 FORMAT(1X,'SR MATRIX FOR COORDINATE ',I4,' WRITTEN ON FILE ',I4,  INT35320
     $    ' AT LOCATION ',I8)                                           INT35330
    2 FORMAT(1X,'Z  MATRIX FOR COORDINATE ',I4,' WRITTEN ON FILE ',I4,  INT35340
     $    ' AT LOCATION ',I8)                                           INT35350
      NLEN=NC**4                                                        INT35360
      NNN=2*NLEN                                                        INT35370
      MMM=(ABS(R)-1)*NNN                                                INT35380
      IF(R.GT.0) MMM=MMM+NLEN                                           INT35390
      LLL=INTOWP(MMM)+1                                                 INT35400
      CALL WWRITW(ISCR,Z,INTOWP(NLEN),LLL,KKK)                          INT35410
      IF(LPRT(4,NPRT).GE.1) RETURN                                      INT35420
      IF(R.LT.0) THEN                                                   INT35430
      WRITE(ICHECK,1) ABS(R),ISCR,LLL                                   INT35440
      ELSE                                                              INT35450
      WRITE(ICHECK,2) ABS(R),ISCR,LLL                                   INT35460
      END IF                                                            INT35470
      RETURN                                                            INT35480
      END                                                               INT35490
C///////////////////////////////////////////////////////////////////////INT35500
      SUBROUTINE ZIN(NC,NS,Z,R,ISCR)                                    INT35510
      REAL*8 Z(NC,NC,NC,NC)                                             INT35520
      INTEGER R                                                         INT35530
      COMMON /IO/ IIN1,IOUT,IIN2,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25,
     $   I31,I32,I33,I35,I36,I37,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10,
     $   ISCR11,ISCR12,ISCR13,ISCR14
    1 FORMAT(1X,'SR MATRIX FOR COORDINATE ',I4,' READ FROM FILE ',I4,   INT35580
     $    ' AT LOCATION ',I8)                                           INT35590
    2 FORMAT(1X,'Z  MATRIX FOR COORDINATE ',I4,' READ FROM FILE ',I4,   INT35600
     $    ' AT LOCATION ',I8)                                           INT35610
      NLEN=NC**4                                                        INT35620
      NNN=2*NLEN                                                        INT35630
      MMM=(ABS(R)-1)*NNN                                                INT35640
      IF(R.GT.0) MMM=MMM+NLEN                                           INT35650
      LLL=INTOWP(MMM)+1                                                 INT35660
      CALL WREADW(ISCR,Z,INTOWP(NLEN),LLL,KKK)                          INT35670
      IF(LPRT(4,NPRT).GE.1) RETURN                                      INT35680
      IF(R.LT.0) THEN                                                   INT35690
      WRITE(ICHECK,1) ABS(R),ISCR,LLL                                   INT35700
      ELSE                                                              INT35710
      WRITE(ICHECK,2) ABS(R),ISCR,LLL                                   INT35720
      END IF                                                            INT35730
      RETURN                                                            INT35740
      END                                                               INT35750
C////////////////////////////////////////////////////////////////////   INT35760
      INTEGER FUNCTION LPRT(K,NPRT)                                     INT35770
      DIMENSION N(4)                                                    INT35780
      NR=NPRT                                                           INT35790
      N(4)=NR/1000                                                      INT35800
      NR=NR-N(4)*1000                                                   INT35810
      N(3)=NR/100                                                       INT35820
      NR=NR-N(3)*100                                                    INT35830
      N(2)=NR/10                                                        INT35840
      NR=NR-N(2)*10                                                     INT35850
      N(1)=NR                                                           INT35860
      LPRT=N(K)                                                         INT35870
      END                                                               INT35880
C////////////////////////////////////////////////////////////////////// INT51150
      SUBROUTINE TABLE1(NN,N,A)                                         INT51050
      REAL*8 A(NN,NN)                                                   INT51060
      COMMON /IO/ IIN1,IOUT,IIN2,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25,
     $   I31,I32,I33,I35,I36,I37,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10,
     $   ISCR11,ISCR12,ISCR13,ISCR14
 1    FORMAT(6(I5,I3,F10.6))                                            INT51110
      DO 10 I=1,N                                                       INT51120
 10   WRITE(IOUT,1) (I,J,A(I,J),J=1,N)                                  INT51130
      END                                                               INT51140
C////////////////////////////////////////////////////////////////////// INT51150
      SUBROUTINE TABLE2(NN,N,M,W,A)                                     INT51160
      REAL*8 A(NN,NN),W(NN)                                             INT51170
      CHARACTER LINE*10                                                 INT51180
      COMMON /IO/ IIN1,IOUT,IIN2,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25,
     $   I31,I32,I33,I35,I36,I37,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10,
     $   ISCR11,ISCR12,ISCR13,ISCR14
 1    FORMAT(1X,11A10)                                                  INT51230
 2    FORMAT(11X,10(I6,4X))                                             INT51240
 3    FORMAT(1X,'COORDINATE',10F10.3)                                   INT51250
 4    FORMAT(I7,4X,10F10.6)                                             INT51260
 5    FORMAT(//)                                                        INT51270
      LINE='----------'                                                 INT51280
      N2=0                                                              INT51290
 100  N1=N2+1                                                           INT51300
      N2=N2+10                                                          INT51310
      IF(N2.GT.M) N2=M                                                  INT51320
      WRITE(IOUT,1) (LINE,I=N1,N2+1)                                    INT51330
      WRITE(IOUT,2) (I,I=N1,N2)                                         INT51340
      WRITE(IOUT,3) (W(I),I=N1,N2)                                      INT51350
      WRITE(IOUT,1) (LINE,I=N1,N2+1)                                    INT51360
      DO 105 I=1,N                                                      INT51370
 105  WRITE(IOUT,4) I,(A(I,J),J=N1,N2)                                  INT51380
      WRITE(IOUT,1) (LINE,I=N1,N2+1)                                    INT51390
      WRITE(IOUT,5)                                                     INT51400
      IF(N2.LT.M) GO TO 100                                             INT51410
      RETURN                                                            INT51420
      END                                                               INT51430
C////////////////////////////////////////////////////////////////////// INT51440
      SUBROUTINE TABLE3(NN,N,W,A)                                       INT51450
      REAL*8 A(NN,NN),W(NN)                                             INT51460
      CHARACTER LINE*10                                                 INT51470
      COMMON /IO/ IIN1,IOUT,IIN2,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25,
     $   I31,I32,I33,I35,I36,I37,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10,
     $   ISCR11,ISCR12,ISCR13,ISCR14
 1    FORMAT(1X,11A10)                                                  INT51520
 2    FORMAT(11X,10(I6,4X))                                             INT51530
 3    FORMAT(1X,'COORDINATE',10F10.3)                                   INT51540
 4    FORMAT(I7,4X,10F10.5)                                             INT51550
 5    FORMAT(//)                                                        INT51560
      LINE='----------'                                                 INT51570
      N2=0                                                              INT51580
 100  N1=N2+1                                                           INT51590
      N2=N2+10                                                          INT51600
      IF(N2.GT.N) N2=N                                                  INT51610
      WRITE(IOUT,1) (LINE,I=N1,N2+1)                                    INT51620
      WRITE(IOUT,2) (I,I=N1,N2)                                         INT51630
      WRITE(IOUT,3) (W(I),I=N1,N2)                                      INT51640
      WRITE(IOUT,1) (LINE,I=N1,N2+1)                                    INT51650
      DO 105 I=1,N                                                      INT51660
 105  WRITE(IOUT,4) I,(A(I,J),J=N1,N2)                                  INT51670
      WRITE(IOUT,1) (LINE,I=N1,N2+1)                                    INT51680
      WRITE(IOUT,5)                                                     INT51690
      IF(N2.LT.N) GO TO 100                                             INT51700
      RETURN                                                            INT51710
      END                                                               INT51720
C////////////////////////////////////////////////////////////////////// INT51150
      SUBROUTINE TABLE4(NN,N,M,W,A)                                     INT51160
      REAL*8 A(NN,NN),W(NN)                                             INT51170
      CHARACTER LINE*10                                                 INT51180
      COMMON /IO/ IIN1,IOUT,IIN2,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25,
     $   I31,I32,I33,I35,I36,I37,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10,
     $   ISCR11,ISCR12,ISCR13,ISCR14
 1    FORMAT(1X,11A10)                                                  INT51230
 2    FORMAT(11X,10(I6,4X))                                             INT51240
 3    FORMAT(1X,'COORDINATE',10E10.3)                                   INT51250
 4    FORMAT(I7,4X,10F10.5)                                             INT51260
 5    FORMAT(//)                                                        INT51270
      LINE='----------'                                                 INT51280
      N2=0                                                              INT51290
 100  N1=N2+1                                                           INT51300
      N2=N2+10                                                          INT51310
      IF(N2.GT.M) N2=M                                                  INT51320
      WRITE(IOUT,1) (LINE,I=N1,N2+1)                                    INT51330
      WRITE(IOUT,2) (I,I=N1,N2)                                         INT51340
      WRITE(IOUT,3) (W(I),I=N1,N2)                                      INT51350
      WRITE(IOUT,1) (LINE,I=N1,N2+1)                                    INT51360
      DO 105 I=1,N                                                      INT51370
 105  WRITE(IOUT,4) I,(A(I,J),J=N1,N2)                                  INT51380
      WRITE(IOUT,1) (LINE,I=N1,N2+1)                                    INT51390
      WRITE(IOUT,5)                                                     INT51400
      IF(N2.LT.M) GO TO 100                                             INT51410
      RETURN                                                            INT51420
      END                                                               INT51430
C////////////////////////////////////////////////////////////////////// INT51730
      SUBROUTINE TABLE5(NN,N,M,W,A)                                     INT51160
      REAL*8 A(NN,NN),W(NN)                                             INT51170
      CHARACTER LINE*10                                                 INT51180
      COMMON /IO/ IIN1,IOUT,IIN2,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25,
     $   I31,I32,I33,I35,I36,I37,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10,
     $   ISCR11,ISCR12,ISCR13,ISCR14
 1    FORMAT(1X,11A10)                                                  INT51230
 2    FORMAT(11X,10(I6,4X))                                             INT51240
 3    FORMAT(1X,'COORDINATE',10F10.6)                                   INT51250
 4    FORMAT(I7,4X,10F10.6)                                             INT51260
 5    FORMAT(//)                                                        INT51270
      LINE='----------'                                                 INT51280
      N2=0                                                              INT51290
 100  N1=N2+1                                                           INT51300
      N2=N2+10                                                          INT51310
      IF(N2.GT.M) N2=M                                                  INT51320
      WRITE(IOUT,1) (LINE,I=N1,N2+1)                                    INT51330
      WRITE(IOUT,2) (I,I=N1,N2)                                         INT51340
      WRITE(IOUT,3) (W(I),I=N1,N2)                                      INT51350
      WRITE(IOUT,1) (LINE,I=N1,N2+1)                                    INT51360
      DO 105 I=1,N                                                      INT51370
 105  WRITE(IOUT,4) I,(A(I,J),J=N1,N2)                                  INT51380
      WRITE(IOUT,1) (LINE,I=N1,N2+1)                                    INT51390
      WRITE(IOUT,5)                                                     INT51400
      IF(N2.LT.M) GO TO 100                                             INT51410
      RETURN                                                            INT51420
      END                                                               INT51430
C////////////////////////////////////////////////////////////////////// INT51150
      SUBROUTINE MASSIN (XMASS,NA,IFLAG)                                XX 04420
      IMPLICIT REAL*8 (A-H,O-Z)                                         XX 04430
      CHARACTER A*12,B*72                                               XX 04440
      LOGICAL P,Q                                                       XX 04450
      DIMENSION XMASS(NA)                                               XX 04460
      PARAMETER(ZERO=0.0D0)                                             XX 04470
      PARAMETER(W1=1.007825D0, W2=2.014D0, W3A=3.01605D0)               XX 04480
      PARAMETER(W3B=3.01603D0, W4=4.0026D0, W5=5.0125D0,W6=6.01512D0)   XX 04490
      PARAMETER(W7=7.01600D0, W9=9.01218D0,W10=10.0129D0)               XX 04500
      PARAMETER(W11=11.00931D0, W12=12.0000D0,W13=13.00335D0)           XX 04510
      PARAMETER(W14=14.00307D0, W15=15.00011D0,W16=15.99491D0)          XX 04520
      PARAMETER(W17=16.999131D0,W18=17.999160D0)                        XX 04530
      PARAMETER(W19=18.99840D0,W20=19.992435D0,W21=20.993843D0)         XX 04540
      PARAMETER(W22=21.991383D0,W23=22.989767D0,W24=23.985042D0)        XX 04550
      PARAMETER(W25=24.985837D0,W26=25.982593D0,W27=26.98153D0)         XX 04560
      PARAMETER(W28=27.976927D0,W29=28.976495D0,W30=29.973770D0)        XX 04570
      PARAMETER(W31=30.973762D0,W32=31.972070D0,W33=32.971456D0)        XX 04580
      PARAMETER(W34=33.967866D0,W35=34.968852D0,W36A=35.967080D0)       XX 04590
      PARAMETER(W36B=35.967545D0,W37=36.965903D0,W38=37.962732D0)       XX 04600
      PARAMETER(W40A=39.962384D0)                                       XX 04610
      COMMON /IO/ IIN1,IOUT,IIN2,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25,
     $   I31,I32,I33,I35,I36,I37,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10,
     $   ISCR11,ISCR12,ISCR13,ISCR14
   1  FORMAT(A72)                                                       XX 04660
   2  FORMAT(F12.6)                                                     XX 04670
      IFLAG=0                                                           XX 04680
      MN=(NA-1)/6+1                                                     XX 04690
      MI=-5                                                             XX 04700
      DO 1000 MM=1,MN                                                   XX 04710
      MI=MI+6                                                           XX 04720
      MF=MI+5                                                           XX 04730
      IF(MF.GT.NA) MF=NA                                                XX 04740
      READ(IIN1,1) B                                                    XX 04750
      NI=-11                                                            XX 04760
      DO 500 M=MI,MF                                                    XX 04770
      NI=NI+12                                                          XX 04780
      NF=NI+11                                                          XX 04790
      A(1:12)=B(NI:NF)                                                  XX 04800
      DO 100 I=1,12                                                     XX 04810
         P=LGE(A(I:I),'A')                                              XX 04820
         Q=LLE(A(I:I),'Z')                                              XX 04830
         IF(P.AND.Q) GO TO 110                                          XX 04840
 100  CONTINUE                                                          XX 04850
      READ(A,2) XMASS(M)                                                XX 04860
      GO TO 500                                                         XX 04870
 110  XMASS(M)=-1.0D0                                                   XX 04880
C HYDROGEN                                                              XX 04890
      IF(A(I:I).EQ.'H') THEN                                            XX 04900
         IF(I.EQ.12) THEN                                               XX 04910
            XMASS(M)=W1                                                 XX 04920
         ELSE                                                           XX 04930
            J=I+1                                                       XX 04940
            IF(A(J:J).EQ.'1'.OR.A(J:J).EQ.' ') XMASS(M)=W1              XX 04950
            IF(A(J:J).EQ.'2') XMASS(M)=W2                               XX 04960
            IF(A(J:J).EQ.'3') XMASS(M)=W3A                              XX 04970
         END IF                                                         XX 04980
      END IF                                                            XX 04990
      IF(A(I:I).EQ.'D') THEN                                            XX 05000
         IF(I.EQ.12) THEN                                               XX 05010
            XMASS(M)=W2                                                 XX 05020
         ELSE                                                           XX 05030
            J=I+1                                                       XX 05040
            IF(A(J:J).EQ.' ') XMASS(M)=W2                               XX 05050
         END IF                                                         XX 05060
      END IF                                                            XX 05070
      IF(A(I:I).EQ.'T') THEN                                            XX 05080
         IF(I.EQ.12) THEN                                               XX 05090
            XMASS(M)=W3A                                                XX 05100
         ELSE                                                           XX 05110
            J=I+1                                                       XX 05120
            IF(A(J:J).EQ.' ') XMASS(M)=W3A                              XX 05130
         END IF                                                         XX 05140
      END IF                                                            XX 05150
C HELIUM                                                                XX 05160
      IF(A(I:I+1).EQ.'HE') THEN                                         XX 05170
         IF(I.EQ.11) THEN                                               XX 05180
            XMASS(M)=W4                                                 XX 05190
         ELSE                                                           XX 05200
            J=I+2                                                       XX 05210
            IF(A(J:J).EQ.'4'.OR.A(J:J).EQ.' ') XMASS(M)=W4              XX 05220
            IF(A(J:J).EQ.'3') XMASS(M)=W3B                              XX 05230
         END IF                                                         XX 05240
      END IF                                                            XX 05250
C LITHIUM                                                               XX 05260
      IF(A(I:I+1).EQ.'LI') THEN                                         XX 05270
         IF(I.EQ.11) THEN                                               XX 05280
            XMASS(M)=W7                                                 XX 05290
         ELSE                                                           XX 05300
            J=I+2                                                       XX 05310
            IF(A(J:J).EQ.'7'.OR.A(J:J).EQ.' ') XMASS(M)=W7              XX 05320
            IF(A(J:J).EQ.'6') XMASS(M)=W6                               XX 05330
         END IF                                                         XX 05340
      END IF                                                            XX 05350
C BERYLLIUM                                                             XX 05360
      IF(A(I:I+1).EQ.'BE') THEN                                         XX 05370
         IF(I.EQ.11) THEN                                               XX 05380
            XMASS(M)=W9                                                 XX 05390
         ELSE                                                           XX 05400
            J=I+2                                                       XX 05410
            IF(A(J:J).EQ.'9'.OR.A(J:J).EQ.' ') XMASS(M)=W9              XX 05420
         END IF                                                         XX 05430
      END IF                                                            XX 05440
C BORON                                                                 XX 05450
      IF(A(I:I).EQ.'B') THEN                                            XX 05460
         IF(I.EQ.12) THEN                                               XX 05470
            XMASS(M)=W11                                                XX 05480
         ELSE                                                           XX 05490
            J=I+1                                                       XX 05500
            K=I+2                                                       XX 05510
            IF(A(J:K).EQ.'11'.OR.A(J:J).EQ.' ') XMASS(M)=W11            XX 05520
            IF(A(J:K).EQ.'10') XMASS(M)=W10                             XX 05530
         END IF                                                         XX 05540
      END IF                                                            XX 05550
C CARBON                                                                XX 05560
      IF(A(I:I).EQ.'C') THEN                                            XX 05570
         IF(I.EQ.12) THEN                                               XX 05580
            XMASS(M)=W12                                                XX 05590
         ELSE                                                           XX 05600
            J=I+1                                                       XX 05610
            K=I+2                                                       XX 05620
            IF(A(J:K).EQ.'12'.OR.A(J:J).EQ.' ') XMASS(M)=W12            XX 05630
            IF(A(J:K).EQ.'13') XMASS(M)=W13                             XX 05640
         END IF                                                         XX 05650
      END IF                                                            XX 05660
C NITROGEN                                                              XX 05670
      IF(A(I:I).EQ.'N') THEN                                            XX 05680
         IF(I.EQ.12) THEN                                               XX 05690
            XMASS(M)=W14                                                XX 05700
         ELSE                                                           XX 05710
            J=I+1                                                       XX 05720
            K=I+2                                                       XX 05730
            IF(A(J:K).EQ.'14'.OR.A(J:J).EQ.' ') XMASS(M)=W14            XX 05740
            IF(A(J:K).EQ.'15') XMASS(M)=W15                             XX 05750
         END IF                                                         XX 05760
      END IF                                                            XX 05770
C OXYGEN                                                                XX 05780
      IF(A(I:I).EQ.'O') THEN                                            XX 05790
         IF(I.EQ.12) THEN                                               XX 05800
            XMASS(M)=W16                                                XX 05810
         ELSE                                                           XX 05820
            J=I+1                                                       XX 05830
            K=I+2                                                       XX 05840
            IF(A(J:K).EQ.'16'.OR.A(J:J).EQ.' ') XMASS(M)=W16            XX 05850
            IF(A(J:K).EQ.'17') XMASS(M)=W17                             XX 05860
            IF(A(J:K).EQ.'18') XMASS(M)=W18                             XX 05870
         END IF                                                         XX 05880
      END IF                                                            XX 05890
C FLUORINE                                                              XX 05900
      IF(A(I:I).EQ.'F') THEN                                            XX 05910
         IF(I.EQ.12) THEN                                               XX 05920
            XMASS(M)=W19                                                XX 05930
         ELSE                                                           XX 05940
            J=I+1                                                       XX 05950
            K=I+2                                                       XX 05960
            IF(A(J:K).EQ.'19'.OR.A(J:J).EQ.' ') XMASS(M)=W19            XX 05970
         END IF                                                         XX 05980
      END IF                                                            XX 05990
C NEON                                                                  XX 06000
      IF(A(I:I+1).EQ.'NE') THEN                                         XX 06010
         IF(I.EQ.11) THEN                                               XX 06020
            XMASS(M)=W20                                                XX 06030
         ELSE                                                           XX 06040
            J=I+2                                                       XX 06050
            K=I+3                                                       XX 06060
            IF(A(J:K).EQ.'20'.OR.A(J:J).EQ.' ') XMASS(M)=W20            XX 06070
            IF(A(J:K).EQ.'21') XMASS(M)=W21                             XX 06080
            IF(A(J:K).EQ.'22') XMASS(M)=W22                             XX 06090
         END IF                                                         XX 06100
      END IF                                                            XX 06110
C SODIUM                                                                XX 06120
      IF(A(I:I+1).EQ.'NA') THEN                                         XX 06130
         IF(I.EQ.11) THEN                                               XX 06140
            XMASS(M)=W23                                                XX 06150
         ELSE                                                           XX 06160
            J=I+2                                                       XX 06170
            K=I+3                                                       XX 06180
            IF(A(J:K).EQ.'23'.OR.A(J:J).EQ.' ') XMASS(M)=W23            XX 06190
         END IF                                                         XX 06200
      END IF                                                            XX 06210
C MAGNESIUM                                                             XX 06220
      IF(A(I:I+1).EQ.'MG') THEN                                         XX 06230
         IF(I.EQ.11) THEN                                               XX 06240
            XMASS(M)=W24                                                XX 06250
         ELSE                                                           XX 06260
            J=I+2                                                       XX 06270
            K=I+3                                                       XX 06280
            IF(A(J:K).EQ.'24'.OR.A(J:J).EQ.' ') XMASS(M)=W24            XX 06290
            IF(A(J:K).EQ.'25') XMASS(M)=W25                             XX 06300
            IF(A(J:K).EQ.'26') XMASS(M)=W26                             XX 06310
         END IF                                                         XX 06320
      END IF                                                            XX 06330
C ALUMINUM                                                              XX 06340
      IF(A(I:I+1).EQ.'AL') THEN                                         XX 06350
         IF(I.EQ.11) THEN                                               XX 06360
            XMASS(M)=W27                                                XX 06370
         ELSE                                                           XX 06380
            J=I+2                                                       XX 06390
            K=I+3                                                       XX 06400
            IF(A(J:K).EQ.'27'.OR.A(J:J).EQ.' ') XMASS(M)=W27            XX 06410
         END IF                                                         XX 06420
      END IF                                                            XX 06430
C SILICON                                                               XX 06440
      IF(A(I:I+1).EQ.'SI') THEN                                         XX 06450
         IF(I.EQ.11) THEN                                               XX 06460
            XMASS(M)=W28                                                XX 06470
         ELSE                                                           XX 06480
            J=I+2                                                       XX 06490
            K=I+3                                                       XX 06500
            IF(A(J:K).EQ.'28'.OR.A(J:J).EQ.' ') XMASS(M)=W28            XX 06510
            IF(A(J:K).EQ.'29') XMASS(M)=W29                             XX 06520
            IF(A(J:K).EQ.'30') XMASS(M)=W30                             XX 06530
         END IF                                                         XX 06540
      END IF                                                            XX 06550
C PHOSPHORUS                                                            XX 06560
      IF(A(I:I).EQ.'P') THEN                                            XX 06570
         IF(I.EQ.12) THEN                                               XX 06580
            XMASS(M)=W31                                                XX 06590
         ELSE                                                           XX 06600
            J=I+1                                                       XX 06610
            K=I+2                                                       XX 06620
            IF(A(J:K).EQ.'31'.OR.A(J:J).EQ.' ') XMASS(M)=W31            XX 06630
         END IF                                                         XX 06640
      END IF                                                            XX 06650
C SULFUR                                                                XX 06660
      IF(A(I:I).EQ.'S') THEN                                            XX 06670
         IF(I.EQ.12) THEN                                               XX 06680
            XMASS(M)=W32                                                XX 06690
         ELSE                                                           XX 06700
            J=I+1                                                       XX 06710
            K=I+2                                                       XX 06720
            IF(A(J:K).EQ.'32'.OR.A(J:J).EQ.' ') XMASS(M)=W32            XX 06730
            IF(A(J:K).EQ.'33') XMASS(M)=W33                             XX 06740
            IF(A(J:K).EQ.'34') XMASS(M)=W34                             XX 06750
            IF(A(J:K).EQ.'36') XMASS(M)=W36A                            XX 06760
         END IF                                                         XX 06770
      END IF                                                            XX 06780
C CHLORINE                                                              XX 06790
      IF(A(I:I+1).EQ.'CL') THEN                                         XX 06800
         IF(I.EQ.11) THEN                                               XX 06810
            XMASS(M)=W35                                                XX 06820
         ELSE                                                           XX 06830
            J=I+2                                                       XX 06840
            K=I+3                                                       XX 06850
            IF(A(J:K).EQ.'35'.OR.A(J:J).EQ.' ') XMASS(M)=W35            XX 06860
            IF(A(J:K).EQ.'37') XMASS(M)=W37                             XX 06870
         END IF                                                         XX 06880
      END IF                                                            XX 06890
C ARGON                                                                 XX 06900
      IF(A(I:I+1).EQ.'AR') THEN                                         XX 06910
         IF(I.EQ.11) THEN                                               XX 06920
            XMASS(M)=W40A                                               XX 06930
         ELSE                                                           XX 06940
            J=I+2                                                       XX 06950
            K=I+3                                                       XX 06960
            IF(A(J:K).EQ.'36') XMASS(M)=W36B                            XX 06970
            IF(A(J:K).EQ.'38') XMASS(M)=W38                             XX 06980
            IF(A(J:K).EQ.'40'.OR.A(J:J).EQ.' ') XMASS(M)=W40A           XX 06990
         END IF                                                         XX 07000
      END IF                                                            XX 07010
C STOP                                                                  XX 07020
      IF(XMASS(M).LE.ZERO) IFLAG=1                                      XX 07030
 500  CONTINUE                                                          XX 07040
 1000 CONTINUE                                                          XX 07050
      RETURN                                                            XX 07060
      END                                                               XX 07070
C///////////////////////////////////////////////////////////////////////INT56370
C                                                                       INT56380
       SUBROUTINE TSTART(ITAPE)                                         INT56390
C                                                                       INT56800
      RETURN                                                            INT57100
C                                                                       INT57110
       ENTRY TSTOP(ITAPE)                                               INT57190
C                                                                       INT57200
      RETURN                                                            INT57710
C                                                                       INT57720
       ENTRY TSET(ITAPE)                                                INT57730
C                                                                       INT57740
      RETURN                                                            INT57810
C                                                                       INT57820
       ENTRY TRSET(ITAPE)                                               INT57830
C                                                                       INT58120
      RETURN                                                            INT58130
      END                                                               INT58140
C///////////////////////////////////////////////////////////////////////INT58150
C     ----- SEARCH THROUGH INPUT FILE FOR TOKEN BEGINNING               INT58890
C           WITH # TO LOCATE INPUT FOR PROGRAM.  IERROR IS              INT58900
C           SET TO 0 IF NO ERRORS, 1 IF ANY ERROR OCCURS.               INT58910
      SUBROUTINE LOCATE(INPUT,TOKEN,IERROR)                             INT58920
      CHARACTER*10 TOKEN,LINE                                           INT58930
      REWIND (UNIT=INPUT,ERR=99)                                        INT58940
    1 CONTINUE                                                          INT58950
      READ (UNIT=INPUT,FMT='(A10)',END=99,ERR=99) LINE                  INT58960
      IF (LINE.NE.TOKEN) GO TO 1                                        INT58970
      IERROR=0                                                          INT58980
      RETURN                                                            INT58990
   99 CONTINUE                                                          INT59000
      IERROR=1                                                          INT59010
      RETURN                                                            INT59020
      END                                                               INT59030
C///////////////////////////////////////////////////////////////////////INT59040
      INTEGER FUNCTION INTOWP(N)                                        INT59050
      INTOWP=2*N                                                        INT59060
      RETURN                                                            INT59070
      END                                                               INT59080
CC//////////////////////////////////////////////////////////////////////INT59090

      SUBROUTINE RFILE(ITAPE)
C
C USE THIS TO OPEN ALL DIRECT ACCESS FILES
C
      IMPLICIT INTEGER (A-Z)      
	 
C
      parameter (reclen = 1024, numbuf = 10)
      common /pointr/ wptr(128),tptr(128)
      COMMON /SECT/ SECTOR
      common /iobufs/ nxtbuf, wrtbuf(numbuf), bufunt(numbuf),
     &   untbuf(128), bufrec(numbuf), buffer(reclen,numbuf)
      character*1 digit
      character*71 file
      character*80 prefix
      character*71 scrtch
      data first/1/
C
      ITAP3 = 3
      ITAP6 = 6
 
c Initialize on the first call to rfile.
      if (first .eq. 1) then
        first = 0
        do 101 i = 1, numbuf
          bufunt(i) = 0
          bufrec(i) = 0
          wrtbuf(i) = 0
  101     continue
        do 102 I = 1, 128
          untbuf(i) = 0
  102     continue
        nxtbuf = 1
        endif
c
c  Some programs require SECTOR = 1024
      SECTOR = 1024
 
      call locate(5,'# FILES ##',ierr)
      if (ierr .ne. 0) then
        iprfx = 0
      else
 804    continue
        prefix = ' '
        read(5,'(A80)',err=904) prefix
 
 904    continue
        iprfx = 0
 805    continue
        iprfx = iprfx + 1
        if (iprfx .eq. 81) goto 905
        if (prefix(iprfx:iprfx).ne.' ') goto 805
 905    continue
        iprfx = iprfx - 1
        endif
      file='file'//digit(itape/10)//digit(itape-10*(itape/10))
      if (iprfx .ne. 0) then
        scrtch = file
        file=prefix(1:iprfx)//'.'//scrtch
        endif



      open (unit=itape,access='DIRECT',recl=reclen*4,form='UNFORMATTED',
     &      iostat=ierr,status='UNKNOWN',file=file)
C
      IF (IERR .ne. 0) then
        WRITE(*,*) ' ERROR ENCOUNTERED OPENING FILE ',ITAPE,' IERR =',
     &            ierr
        WRITE(6,*) ' ERROR ENCOUNTERED OPENING FILE ',ITAPE,' IERR =',
     &            ierr
        write(*,*) ' rfile: filename=',file
        write(6,*) ' rfile: filename=',file
        call ioerr(ierr)
        call mabort
        endif
C
C PLACE THE POINTER OF ITAPE AT THE BEGINNING
C
      CALL SREW(ITAPE)
C
      RETURN
      END

      SUBROUTINE SREW(ITAPE)                                            
C                                                                       
      IMPLICIT INTEGER (A-Z)                                            
C                                                                       
      common /pointr/ wptr(128),tptr(128)
C                                                                       
C THIS ROUTINE REWINDS FILES USING THE FORTRAN UTILITIES DIRECT ACCESS  
C                                                                       
C PLACE THE POINTER OF ITAPE AT THE BEGINNING                           
C                                                                       
      wptr(itape)=1
C                                                                       
      RETURN                                                            
      END         

      SUBROUTINE WREADW(file,array,nwords,fword1,nxtwrd)
c
c  NB. The first word of the file is considered to be "1" by R
c
      implicit integer(a-z)
c
      parameter (RECLEN=1024,numbuf=10)
c

      COMMON /IO/ IIN1,IOUT,IIN2,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25,
     $   I31,I32,I33,I35,I36,I37,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10,
     $   ISCR11,ISCR12,ISCR13,ISCR14
      common /pointr/ wptr(128),tptr(128)
      common /iobufs/ nxtbuf, wrtbuf(numbuf), bufunt(numbuf),
     &  untbuf(128), bufrec(numbuf), buffer(reclen,numbuf)
c
      dimension array(nwords)
c
c      write(iout,1000) file,nwords,fword
 1000 format(' wread:    file',i3,', ',i12,' words,    first word',i8)
c
      fword=fword1
c
c     ---- determine which records need to be read ----
c
      frec=(fword-1)/reclen + 1
      lword=fword+nwords-1
      lrec=(lword-1)/reclen + 1
c
c      ---- loop over the records ----
c
      do 10 irec=frec,lrec
c
c       If this file does not have a buffer then grab one.
        if (untbuf(file) .eq. 0) then
c         Deallocate buffer from previous file.
          if (bufunt(nxtbuf) .ne. 0) then
c           if buffer has been written to then flush it
            call bioflu(bufunt(nxtbuf),'WREADW')
            untbuf(bufunt(nxtbuf)) = 0
            endif
          untbuf(file) = nxtbuf
          bufunt(nxtbuf) = file
          bufrec(nxtbuf) = 0
          nxtbuf = nxtbuf + 1
          if (nxtbuf .gt. numbuf) nxtbuf = 1
          endif
 
c        If we don't have irec yet, then read it in.
         if (bufrec(untbuf(file)) .ne. irec) then
           call bioflu(file,'WREADW')
           read(file,rec=irec,iostat=ierr)
     &       (buffer(i,untbuf(file)),i=1,reclen)
           if(ierr.ne.0) then
             write(iout,1000) file,nwords,fword
             write(iout,*) 'frec,lrec,irec=',frec,lrec,irec
             call ioerr(ierr)
             call mabort
             endif
           bufrec(untbuf(file)) = irec
           endif
c
c
c           -- transfer the appropriate portion --
c           --      of buffer to array          --
c
            if(irec.eq.frec) then
c
               offset=fword-(irec-1)*reclen - 1
               count=reclen-offset
               if(count.gt.nwords) count=nwords
               do 20 i=1,count
                  array(i)=buffer(i+offset,untbuf(file))
   20          continue
c
            else if(irec.eq.lrec) then
c
               count=lword-(irec-1)*reclen
               offset=nwords-count
               do 30 i=1,count
                  array(i+offset)=buffer(i,untbuf(file))
   30          continue
c
            else
c
              offset=reclen*(irec-1)+1-fword
              do 40 i=1,reclen
                 array(i+offset)=buffer(i,untbuf(file))
   40         continue
c
            end if
c
c
   10 continue
c
c      ---- update pointer ----
c
c
c     wptr is the next word after the last word transferred
      wptr(file)=lword+1
      nxtwrd=wptr(file)
c
      return
      END

      SUBROUTINE WWRITW(file,array,nwords,fword1,nxtwrd)
c
c   NB.  The file is considered to start at word "1" in R
c
      implicit integer(a-z)
c
      parameter (RECLEN=1024,numbuf=10)
c
      COMMON /IO/ IIN1,IOUT,IIN2,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25, 
     $   I31,I32,I33,I35,I36,I37,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10,
     $   ISCR11,ISCR12,ISCR13,ISCR14
      common /pointr/ wptr(128),tptr(128)
      common /iobufs/ nxtbuf, wrtbuf(numbuf), bufunt(numbuf),
     &  untbuf(128), bufrec(numbuf), buffer(reclen,numbuf)
c
      dimension array(nwords)
c
 1000 format(' wwrit:    file',i3,i12,' words   first word:',i8)
c
      fword=fword1
c
c     ---- determine which records need to be written ----
c
      frec=(fword-1)/reclen + 1
      lword=fword+nwords-1
      lrec=(lword-1)/reclen + 1
c
c      ---- loop over the records ----
c
      do 10 irec=frec,lrec
c
c       If file does not have a buffer, then allocate one.
        if (untbuf(file) .eq. 0) then
c         Deallocate buffer from previous file
          if (bufunt(nxtbuf) .ne. 0) then
c           Flush buffer if written to
            call bioflu(bufunt(nxtbuf),'WWRITW(PURGE)')
            untbuf(bufunt(nxtbuf)) = 0
            endif
          untbuf(file) = nxtbuf
          bufunt(nxtbuf) = file
          bufrec(nxtbuf) = 0
          nxtbuf = nxtbuf + 1
          if (nxtbuf .gt. numbuf) nxtbuf = 1
          endif
c
        if (bufrec(untbuf(file)) .ne. irec) then
          call bioflu(file,'WWRITW(NEW RECORD)')
cwa  Warning!
          read(file,rec=irec,err=300,iostat=ierr) 
     &        (buffer(i,untbuf(file)),i=1,reclen)
          if (ierr.eq.-1) go to 300
csa
          if (ierr .ne. 0) then
            write(iout,8893) ierr
 8893       format(//' error code for read of buffer',i10)
            write(iout,*) 'file =',file,'irec =',irec
            call ioerr(ierr)
            call mabort
            endif
          bufrec(untbuf(file)) = irec
          endif
          go to 310
  300     continue
          bufrec(untbuf(file)) = irec
 1300     format('wwritw: first access of record',i10,' file ',i2,
     &            ' err ',i4)
  310     continue
c
c           -- update buffer --
c
            if(irec.eq.frec) then
c
               offset=fword-(irec-1)*reclen - 1
               count=reclen-offset
               if(count.gt.nwords) count=nwords
               do 20 i=1,count
                  buffer(i+offset,untbuf(file))=array(i)
   20          continue
c
            else if(irec.eq.lrec) then
c
               count=lword-(irec-1)*reclen
               offset=nwords-count
               do 30 i=1,count
                  buffer(i,untbuf(file))=array(i+offset)
   30          continue
c
            else
c
              offset=reclen*(irec-1)+1-fword
              do 40 i=1,reclen
                 buffer(i,untbuf(file))=array(i+offset)
   40         continue
c
            end if
c
c      ---- mark the buffer as written to ----
c
            wrtbuf(untbuf(file)) = 1
c
   10 continue
c
c      ---- update pointers ----
c
c     tptr is the highest record written
      if(lrec.gt.tptr(file)) tptr(file)=lrec
c
c     wptr is the next word after the last word transferred
      wptr(file)=lword+1
      nxtwrd=wptr(file)

c     write(iout,8893) lword,wptr(file),nxtwrd
c
c
      return
      end
      SUBROUTINE bioflu(file,caller)
c
c  NB. The first word of the file is considered to be "1" by R
c
      implicit integer(a-z)
c
      parameter (RECLEN=1024,numbuf=10)
      character*(*) caller
c
      COMMON /IO/ IIN1,IOUT,IIN2,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25, 
     $   I31,I32,I33,I35,I36,I37,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10,
     $   ISCR11,ISCR12,ISCR13,ISCR14
      common /pointr/ wptr(128),tptr(128)
      common /iobufs/ nxtbuf, wrtbuf(numbuf), bufunt(numbuf),
     &  untbuf(128), bufrec(numbuf), buffer(reclen,numbuf)
 
      flubuf = untbuf(file)
c If the file does not have a buffer then return
      if (flubuf .eq. 0) return
c
c
c           if buffer has been written to then flush it
            if (wrtbuf(flubuf) .ne. 0) then
              write(bufunt(flubuf),rec=bufrec(flubuf),iostat=ierr)
     &           (buffer(i,flubuf),i=1,reclen)
              if (ierr .ne. 0) then
              write(iout,*) 'bioflu: error in write '
              write(iout,*) 'bioflu: called by ', caller
              write(iout,*) ' nxtbuf ',nxtbuf,' file ',file,
     &                                         ' flubuf ',flubuf
              write(iout,*) ' untbuf(file) ', untbuf(file)
              write(iout,*) ' record ', bufrec(flubuf)
              write(iout,*) ' byte ', (bufrec(flubuf)-1)*reclen*4 + 1
              write(iout,2) ' buffer ',' wrtbuf ',' bufunt ',' bufrec '
   2            format(4a10)
                do 101 i = 1, numbuf
                  write(iout,1) i,wrtbuf(i),bufunt(i),bufrec(i)
   1              format(4i10)
 101              continue
                call ioerr(ierr)
                call mabort
                endif
              wrtbuf(flubuf) = 0
              endif
c
      return
      END

      character*1 function digit(i)
      character*1 digits(0:9)
      data digits/'0','1','2','3','4','5','6','7','8','9'/
      digit = digits(i)
      return
      end


      subroutine ioerr(ierr)
      write(6,*) 'ioerr: error ', ierr, ' encountered.'
      return
      end


      SUBROUTINE MABORT                                                 
                                                             
      write(6,*) 'mabort: ERROR, aborting'
      write(*,*) 'mabort: ERROR, aborting'

      write(0,*) 'mabort: ERROR, aborting'

      call exit(3)

      stop 876
C                                                                       
      END



C///////////////////////////////////////////////////////////////////////	
      SUBROUTINE RCLOSE(ITAPE,JCODE)	
      IMPLICIT INTEGER (A-Z)	

C JCODE = 4     CLOSE AND DELETE FILE	
C JCODE = 3     CLOSE AND SAVE FILE	
	
      IF(ITAPE.EQ.6) STOP ' YOU CANNOT CLOSE A FILE ON UNIT 6'	

	IF (JCODE.EQ.4) THEN
             CLOSE (UNIT=ITAPE,IOSTAT=IERR,STATUS='DELETE')
                      ELSE
             CLOSE (UNIT=ITAPE,IOSTAT=IERR,STATUS='KEEP')
             END IF
		           
      IF (IERR.NE.0) THEN
             WRITE(*,*) ' ERROR ENCOUNTERED CLOSING FILE',ITAPE	
             WRITE(*,*) ' IERR,JCODE = ',IERR,JCODE	
             END IF
 	
      RETURN	
      END
C///////////////////////////////////////////////////////////////////////
      SUBROUTINE NCLOSE(ITAPE)
C
      IF(ITAPE.EQ.6) STOP ' YOU CANNOT CLOSE A FILE ON UNIT 6'
      JCODE=4
      REWIND(ITAPE)
      READ(ITAPE,*,END=20)
      JCODE=3
 20   CONTINUE
C
      IF (JCODE.EQ.4) THEN
             CLOSE (UNIT=ITAPE,IOSTAT=IERR,STATUS='DELETE')
      ELSE
             CLOSE (UNIT=ITAPE,IOSTAT=IERR,STATUS='KEEP')
      END IF
C
      IF (IERR.NE.0) THEN
             WRITE(*,*) ' ERROR ENCOUNTERED CLOSING FILE',ITAPE
             WRITE(*,*) ' IERR,JCODE = ',IERR,JCODE
      END IF
C
      RETURN
      END
C///////////////////////////////////////////////////////////////////////

      SUBROUTINE OPENFF (ITAPE, IRECL)

 	IMPLICIT INTEGER (A-Z)      
	character*1 digit
      character*71 file
      character*80 prefix
      character*71 scrtch


C Open formatted file.  Set IRECL=0 for sequential access.

      call locate(5,'# FILES ##',ierr)
      if (ierr .ne. 0) then
        iprfx = 0
      else
 804    continue
        prefix = ' '
        read(5,'(A80)',err=904) prefix
 
 904    continue
        iprfx = 0
 805    continue
        iprfx = iprfx + 1
        if (iprfx .eq. 81) goto 905
        if (prefix(iprfx:iprfx).ne.' ') goto 805
 905    continue
        iprfx = iprfx - 1
        endif
      
      if (itape.ne.3) then
         file='file'//digit(itape/10)//digit(itape-10*(itape/10))
         if (iprfx .ne. 0) then
              scrtch = file
              file=prefix(1:iprfx)//'.'//scrtch
              endif
                      else
         file = prefix(1:iprfx)//'.CHECK'
         endif

      if (irecl.eq.0) then
CSA   write (iout,*) 'open seq',itape,irecl 
	open (unit=itape, file=file, form='FORMATTED',
     $            access='SEQUENTIAL', iostat=ierr)

			     else
      write (iout,*) 'open direct',itape,irecl
		open (unit=itape, file=file, form='FORMATTED',
     $            access='DIRECT', recl=irecl, iostat=ierr)
		end if

	if (ierr.ne.0) then

        WRITE(*,*) ' ERROR ENCOUNTERED OPENING FILE ',ITAPE,' IERR =',
     &            ierr
        WRITE(iout,*) ' ERROR ENCOUNTERED OPENING FILE ',ITAPE,
     &           ' IERR =', ierr
        write(*,*) ' rfile: filename=',file
        write(iout,*) ' rfile: filename=',file
        call ioerr(ierr)
        call mabort
        endif

	return

	end

C///////////////////////////////////////////////////////////////////////
      subroutine nounfl
C Set cpu control register so no floating-point underflow exceptions
C Taken out for Mac OS X version.

	return
	end
C     ////////////////////////////////////////////////////////////      INT14360
      SUBROUTINE LINTR(NC,NSX,NSY,IOPT,A,XS,F1,F2,F3,F4,V,              INT14370
     $            NMM,NRUN,IFLAG)                                       INT14380
      IMPLICIT REAL*8 (A-H,O-Z)                                         INT14390
      INTEGER P,Q,QQ                                                    INT14400
      DIMENSION F1(NC),F2(NC,NC),F3(NC,NC,NC),F4(NC,NC,NC,NC)           INT14410
      DIMENSION IOPT(30),A(NC,NC),XS(NC,NC),V(NMM),DIP(3)               INT14420
      PARAMETER(ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0)                         INT14430
      COMMON /IO/ IIN1,IOUT,IIN2,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25,
     $   I31,I32,I33,I35,I36,I37,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10,
     $   ISCR11,ISCR12,ISCR13,ISCR14
      COMMON /PHYCON/ BOHR,DEBYE,HART,WAVE0,CINT
      COMMON /DIPOLE/DIP,QQ                                             INT14480
 1    FORMAT(20X,3F20.10)                                               INT14490
 2    FORMAT(8F10.7)                                                    INT14500
 3    FORMAT(3F20.10)                                                   INT14510
 4    FORMAT(F16.12)                                                    INT14520
 5    FORMAT(6F12.6)                                                    INT14530
 6    FORMAT(2I5)                                                       INT14540
 7    FORMAT(/' LINEAR TRANSFORMATION CONTRIBUTIONS TO INTERNAL'/       INT14550
     $   ' COORDINATE FORCE CONSTANTS:')                                INT14560
 8    FORMAT(/' LINEAR TRANSFORMATION CONTRIBUTIONS TO CARTESIAN'/      INT14570
     $   ' COORDINATE FORCE CONSTANTS:')                                INT14580
 9    FORMAT(/' QUADRATIC FORCE CONSTANTS')                             INT14590
 10   FORMAT(/' CUBIC FORCE CONSTANTS')                                 INT14600
 12   FORMAT(/' QUARTIC FORCE CONSTANTS')                               INT14610
 13   FORMAT(3F20.10)                                                   INT14620
 14   FORMAT(20X,F20.10)                                                INT14630
 15   FORMAT(2I5,3F20.10)                                               INT14640
 16   FORMAT(/1X,'TOTAL CHARGE AND DIPOLE MOMENT')                      INT14650
 17   FORMAT(I5,3F20.10)                                                INT14660
      IFLAG=0                                                           INT14670
      NA=IOPT(1)                                                        INT14680
      NDER=IOPT(4)                                                      INT14690
      NEQ=IOPT(5)                                                       INT14700
      NINV=IOPT(7)                                                      INT14710
CSA
      NNINV=0
      IF (ABS(NINV).NE.3) GO TO 25
         NEQ=0
         NNINV=NINV
         NINV=0
 25   NVEC=IOPT(13)                                                     INT14720
      INP1=I11                                                          INT14730
      INP2=I15                                                          INT14740
      INP3=I20                                                          INT14750
      INP4=I24                                                          INT14760
      INP5=I17                                                          INT14770
      IF(ABS(NNINV).NE.3) THEN
C        CF1=8.238857606D0                                              INT14780
C        CF2=15.56918890D0                                              INT14790
C        CF3=29.42151140D0                                              INT14800
C        CF4=55.59861458D0                                              INT14810
         CF1=HART/BOHR
         CF2=CF1/BOHR
         CF3=CF2/BOHR
         CF4=CF3/BOHR
      ELSE
         CF1=ONE
         CF2=ONE
         CF3=ONE
         CF4=ONE
      END IF
CSA
 30   IF(NINV.GE.1) THEN                                                INT14820
      INP1=I12                                                          INT14830
      INP2=I16                                                          INT14840
      INP3=I21                                                          INT14850
      INP4=I25                                                          INT14860
      INP5=I18                                                          INT14870
      CF1=ONE                                                           INT14880
      CF2=ONE                                                           INT14890
      CF3=ONE                                                           INT14900
      CF4=ONE                                                           INT14910
      END IF                                                            INT14920
      IF(NEQ.EQ.0) GO TO 120                                            INT14930
      IF(NVEC.EQ.1) THEN                                                INT14940
        IF(NRUN.LE.1) THEN                                              INT14950
        READ(INP5,15) M,QQ,(DIP(I),I=1,3)                               INT14960
        WRITE(IOUT,16)                                                  INT14970
        WRITE(IOUT,17) QQ,(DIP(I),I=1,3)                                INT14980
          IF(M.NE.NA) THEN                                              INT14990
             IFLAG=5                                                    INT15000
             RETURN                                                     INT15010
          END IF                                                        INT15020
        END IF                                                          INT15030
        READ(INP5,13) (V(IK),IK=1,NSY)                                  INT15040
        GO TO 122                                                       INT15050
      END IF                                                            INT15060
      IF(NINV.GE.1) THEN                                                INT15070
        READ(INP1,6) M                                                  INT15080
        IF(M.NE.NA) THEN                                                INT15090
             IFLAG=1                                                    INT15100
             RETURN                                                     INT15110
        END IF                                                          INT15120
        READ(INP1,14)(V(IK),IK=1,NSY)                                   INT15130
      ELSE                                                              INT15140
        READ(INP1,1)(V(IK),IK=1,NSY)                                    INT15150
      END IF                                                            INT15160
 122  DO 110  M=1,NSX                                                   INT15170
 110  F1(M)=ZERO                                                        INT15180
      DO 115  I=1,NSY                                                   INT15190
      V(I)=V(I)*CF1                                                     INT15200
      DO 115  M=1,NSX                                                   INT15210
 115  F1(M)=F1(M)+V(I)*A(I,M)                                           INT15220
 120  IF(NDER.LE.1) GO TO 1000                                          INT15230
      NNC=NSY*NSY                                                       INT15240
      NM=NMM                                                            INT15250
      NR=NNC/NMM                                                        INT15260
      NL=NNC-NMM*NR                                                     INT15270
      READ(INP2,6) M,N                                                  INT15280
      IF(M.NE.NA) THEN                                                  INT15290
           IFLAG=2                                                      INT15300
           RETURN                                                       INT15310
      END IF                                                            INT15320
      DO 125  N=1,NSX                                                   INT15330
      DO 125  I=1,NSY                                                   INT15340
 125  XS(I,N)=ZERO                                                      INT15350
      KK=0                                                              INT15360
      DO 128  II=1,NR+1                                                 INT15370
      IF(II.EQ.NR+1) NM=NL                                              INT15380
      READ(INP2,3) (V(IK),IK=1,NM)                                      INT15390
           DO 128  IK=1,NM                                              INT15400
           V(IK)=V(IK)*CF2                                              INT15410
           KK=KK+1                                                      INT15420
           J=(KK-1)/NSY+1                                               INT15430
           I=KK-NSY*(J-1)                                               INT15440
           DO 130  N=1,NSX                                              INT15450
 130       XS(I,N)=XS(I,N)+A(J,N)*V(IK)                                 INT15460
 128  CONTINUE                                                          INT15470
      DO 131  N=1,NSX                                                   INT15480
      DO 131  M=1,NSX                                                   INT15490
      XX=ZERO                                                           INT15500
      DO 132  I=1,NSY                                                   INT15510
 132  XX=XX+A(I,M)*XS(I,N)                                              INT15520
 131  F2(M,N)=XX                                                        INT15530
      DO 134  M=2,NSX                                                   INT15540
      DO 134  N=1,M-1                                                   INT15550
      F2(M,N)=(F2(M,N)+F2(N,M))/TWO                                     INT15560
 134  F2(N,M)=F2(M,N)                                                   INT15570
      IF(NDER.LE.2) GO TO 1000                                          INT15580
      NM=NMM                                                            INT15590
      NNC=(NSY*(NSY+1)*(NSY+2))/6                                       INT15600
      NR=NNC/NMM                                                        INT15610
      NL=NNC-NMM*NR                                                     INT15620
Cwa 12/5/2003
C     READ(INP3,6) M,N                                                  INT15630
C     IF(M.NE.NA.OR.N.NE.NNC) THEN                                      INT15640
      READ(INP3,6) M
      IF(M.NE.NA) THEN
           IFLAG=3                                                      INT15650
           RETURN                                                       INT15660
      END IF                                                            INT15670
      DO 135  P=1,NSX                                                   INT15680
      DO 135  I=1,NSY                                                   INT15690
      DO 135  J=1,I                                                     INT15700
 135  F3(I,J,P)=ZERO                                                    INT15710
      I=0                                                               INT15720
      J=0                                                               INT15730
      K=0                                                               INT15740
      DO 138  II=1,NR+1                                                 INT15750
      IF(II.EQ.NR+1) NM=NL                                              INT15760
      READ(INP3,3) (V(IK),IK=1,NM)                                      INT15770
          DO 138  IK=1,NM                                               INT15780
          V(IK)=V(IK)*CF3                                               INT15790
          IF(K.LT.J) THEN                                               INT15800
             K=K+1                                                      INT15810
             GO TO 139                                                  INT15820
          END IF                                                        INT15830
          IF(J.LT.I) THEN                                               INT15840
             J=J+1                                                      INT15850
             K=1                                                        INT15860
             GO TO 139                                                  INT15870
          END IF                                                        INT15880
          I=I+1                                                         INT15890
          J=1                                                           INT15900
          K=1                                                           INT15910
 139  CONTINUE                                                          INT15920
          IF(I.NE.J) THEN                                               INT15930
             IF(J.NE.K) THEN                                            INT15940
                DO 140  P=1,NSX                                         INT15950
                F3(I,J,P)=V(IK)*A(K,P)+F3(I,J,P)                        INT15960
                F3(I,K,P)=V(IK)*A(J,P)+F3(I,K,P)                        INT15970
 140            F3(J,K,P)=V(IK)*A(I,P)+F3(J,K,P)                        INT15980
             ELSE                                                       INT15990
                DO 141  P=1,NSX                                         INT16000
                F3(I,J,P)=V(IK)*A(J,P)+F3(I,J,P)                        INT16010
 141            F3(J,J,P)=V(IK)*A(I,P)+F3(J,J,P)                        INT16020
             END IF                                                     INT16030
          ELSE                                                          INT16040
             IF(J.NE.K) THEN                                            INT16050
                DO 142  P=1,NSX                                         INT16060
                F3(I,I,P)=V(IK)*A(K,P)+F3(I,I,P)                        INT16070
 142            F3(I,K,P)=V(IK)*A(I,P)+F3(I,K,P)                        INT16080
             ELSE                                                       INT16090
                DO 143  P=1,NSX                                         INT16100
 143            F3(I,I,P)=V(IK)*A(I,P)+F3(I,I,P)                        INT16110
             END IF                                                     INT16120
          END IF                                                        INT16130
 138  CONTINUE                                                          INT16140
      REWIND ISCR5                                                      INT16150
      WRITE(ISCR5,3)(((F3(I,J,P),P=1,NSX),J=1,I),I=1,NSY)               INT16160
      DO 144  P=1,NSX                                                   INT16170
      DO 144  N=1,P                                                     INT16180
      DO 144  I=1,NSY                                                   INT16190
 144  F3(I,N,P)=ZERO                                                    INT16200
      NND=NSY*(NSY+1)*NSX/2                                             INT16210
      NR=NND/NMM                                                        INT16220
      NL=NND-NMM*NR                                                     INT16230
      NM=NMM                                                            INT16240
      I=1                                                               INT16250
      J=1                                                               INT16260
      P=0                                                               INT16270
      REWIND ISCR5                                                      INT16280
      DO 146  II=1,NR+1                                                 INT16290
      IF(II.EQ.NR+1) NM=NL                                              INT16300
      READ(ISCR5,3) (V(IK),IK=1,NM)                                     INT16310
          DO 146  IK=1,NM                                               INT16320
          IF(P.LT.NSX) THEN                                             INT16330
             P=P+1                                                      INT16340
             GO TO 147                                                  INT16350
          END IF                                                        INT16360
          IF(J.LT.I) THEN                                               INT16370
             J=J+1                                                      INT16380
             P=1                                                        INT16390
             GO TO 147                                                  INT16400
          END IF                                                        INT16410
          I=I+1                                                         INT16420
          J=1                                                           INT16430
          P=1                                                           INT16440
 147  CONTINUE                                                          INT16450
          IF(I.NE.J) THEN                                               INT16460
                DO 148  N=1,P                                           INT16470
                F3(I,N,P)=V(IK)*A(J,N)+F3(I,N,P)                        INT16480
 148            F3(J,N,P)=V(IK)*A(I,N)+F3(J,N,P)                        INT16490
          ELSE                                                          INT16500
                DO 149  N=1,P                                           INT16510
 149            F3(I,N,P)=V(IK)*A(I,N)+F3(I,N,P)                        INT16520
          END IF                                                        INT16530
 146  CONTINUE                                                          INT16540
      REWIND ISCR5                                                      INT16550
      WRITE(ISCR5,3)(((F3(I,N,P),I=1,NSY),N=1,P),P=1,NSX)               INT16560
      DO 150  P=1,NSX                                                   INT16570
      DO 150  N=1,P                                                     INT16580
      DO 150  M=1,N                                                     INT16590
 150  F3(M,N,P)=ZERO                                                    INT16600
      NND=NSY*NSX*(NSX+1)/2                                             INT16610
      NR=NND/NMM                                                        INT16620
      NL=NND-NMM*NR                                                     INT16630
      NM=NMM                                                            INT16640
      I=0                                                               INT16650
      N=1                                                               INT16660
      P=1                                                               INT16670
      REWIND ISCR5                                                      INT16680
      DO 152  II=1,NR+1                                                 INT16690
      IF(II.EQ.NR+1) NM=NL                                              INT16700
      READ(ISCR5,3) (V(IK),IK=1,NM)                                     INT16710
          DO 152  IK=1,NM                                               INT16720
          IF(I.LT.NSY) THEN                                             INT16730
             I=I+1                                                      INT16740
             GO TO 153                                                  INT16750
          END IF                                                        INT16760
          IF(N.LT.P) THEN                                               INT16770
             N=N+1                                                      INT16780
             I=1                                                        INT16790
             GO TO 153                                                  INT16800
          END IF                                                        INT16810
          P=P+1                                                         INT16820
          N=1                                                           INT16830
          I=1                                                           INT16840
 153  CONTINUE                                                          INT16850
          DO 154  M=1,N                                                 INT16860
 154      F3(M,N,P)=V(IK)*A(I,M)+F3(M,N,P)                              INT16870
 152  CONTINUE                                                          INT16880
      REWIND ISCR5                                                      INT16890
      CALL FILL3A(NC,NSX,F3)                                            INT16900
      IF(NDER.LE.3) GO TO 1000                                          INT16910
      NM=NMM                                                            INT16920
      NNC=(NSY*(NSY+1)*(NSY+2)*(NSY+3))/24                              INT16930
      NR=NNC/NMM                                                        INT16940
      NL=NNC-NMM*NR                                                     INT16950
Cwa 12/5/2003
C     READ(INP4,6) M,N                                                  INT16960
      READ(INP4,6) M
C     IF(M.NE.NA.OR.N.NE.NNC) THEN                                      INT16970
      IF(M.NE.NA) THEN
           IFLAG=4                                                      INT16980
           RETURN                                                       INT16990
      END IF                                                            INT17000
      DO 170  Q=1,NSX                                                   INT17010
      DO 170  I=1,NSY                                                   INT17020
      DO 170  J=1,I                                                     INT17030
      DO 170  K=1,J                                                     INT17040
 170  F4(I,J,K,Q)=ZERO                                                  INT17050
      I=0                                                               INT17060
      J=0                                                               INT17070
      K=0                                                               INT17080
      L=0                                                               INT17090
      DO 178  II=1,NR+1                                                 INT17100
      IF(II.EQ.NR+1) NM=NL                                              INT17110
      READ(INP4,3) (V(IK),IK=1,NM)                                      INT17120
          DO 178  IK=1,NM                                               INT17130
          V(IK)=V(IK)*CF4                                               INT17140
          IF(L.LT.K) THEN                                               INT17150
             L=L+1                                                      INT17160
             GO TO 180                                                  INT17170
          END IF                                                        INT17180
          IF(K.LT.J) THEN                                               INT17190
             K=K+1                                                      INT17200
             L=1                                                        INT17210
             GO TO 180                                                  INT17220
          END IF                                                        INT17230
          IF(J.LT.I) THEN                                               INT17240
             J=J+1                                                      INT17250
             K=1                                                        INT17260
             L=1                                                        INT17270
             GO TO 180                                                  INT17280
          END IF                                                        INT17290
             I=I+1                                                      INT17300
             J=1                                                        INT17310
             K=1                                                        INT17320
             L=1                                                        INT17330
 180      CONTINUE                                                      INT17340
          IF(I.NE.J) THEN                                               INT17350
             IF(J.NE.K) THEN                                            INT17360
                IF(K.NE.L) THEN                                         INT17370
                DO 182  Q=1,NSX                                         INT17380
                F4(I,J,K,Q)=V(IK)*A(L,Q)+F4(I,J,K,Q)                    INT17390
                F4(I,J,L,Q)=V(IK)*A(K,Q)+F4(I,J,L,Q)                    INT17400
                F4(I,K,L,Q)=V(IK)*A(J,Q)+F4(I,K,L,Q)                    INT17410
 182            F4(J,K,L,Q)=V(IK)*A(I,Q)+F4(J,K,L,Q)                    INT17420
                ELSE                                                    INT17430
                DO 184  Q=1,NSX                                         INT17440
                F4(I,J,K,Q)=V(IK)*A(K,Q)+F4(I,J,K,Q)                    INT17450
                F4(I,K,K,Q)=V(IK)*A(J,Q)+F4(I,K,K,Q)                    INT17460
 184            F4(J,K,K,Q)=V(IK)*A(I,Q)+F4(J,K,K,Q)                    INT17470
                END IF                                                  INT17480
             ELSE                                                       INT17490
                IF(K.NE.L) THEN                                         INT17500
                DO 186  Q=1,NSX                                         INT17510
                F4(I,J,J,Q)=V(IK)*A(L,Q)+F4(I,J,J,Q)                    INT17520
                F4(I,J,L,Q)=V(IK)*A(J,Q)+F4(I,J,L,Q)                    INT17530
 186            F4(J,J,L,Q)=V(IK)*A(I,Q)+F4(J,J,L,Q)                    INT17540
                ELSE                                                    INT17550
                DO 188  Q=1,NSX                                         INT17560
                F4(I,J,J,Q)=V(IK)*A(J,Q)+F4(I,J,J,Q)                    INT17570
 188            F4(J,J,J,Q)=V(IK)*A(I,Q)+F4(J,J,J,Q)                    INT17580
                END IF                                                  INT17590
             END IF                                                     INT17600
          ELSE                                                          INT17610
             IF(J.NE.K) THEN                                            INT17620
                IF(K.NE.L) THEN                                         INT17630
                DO 190  Q=1,NSX                                         INT17640
                F4(I,I,K,Q)=V(IK)*A(L,Q)+F4(I,I,K,Q)                    INT17650
                F4(I,I,L,Q)=V(IK)*A(K,Q)+F4(I,I,L,Q)                    INT17660
 190            F4(I,K,L,Q)=V(IK)*A(I,Q)+F4(I,K,L,Q)                    INT17670
                ELSE                                                    INT17680
                DO 192  Q=1,NSX                                         INT17690
                F4(I,I,K,Q)=V(IK)*A(K,Q)+F4(I,I,K,Q)                    INT17700
 192            F4(I,K,K,Q)=V(IK)*A(I,Q)+F4(I,K,K,Q)                    INT17710
                END IF                                                  INT17720
             ELSE                                                       INT17730
                IF(K.NE.L) THEN                                         INT17740
                DO 194  Q=1,NSX                                         INT17750
                F4(I,I,I,Q)=V(IK)*A(L,Q)+F4(I,I,I,Q)                    INT17760
 194            F4(I,I,L,Q)=V(IK)*A(I,Q)+F4(I,I,L,Q)                    INT17770
                ELSE                                                    INT17780
                DO 196  Q=1,NSX                                         INT17790
 196            F4(I,I,I,Q)=V(IK)*A(I,Q)+F4(I,I,I,Q)                    INT17800
                END IF                                                  INT17810
             END IF                                                     INT17820
          END IF                                                        INT17830
 178  CONTINUE                                                          INT17840
      REWIND ISCR5                                                      INT17850
      WRITE(ISCR5,3) ((((F4(I,J,K,Q),Q=1,NSX),K=1,J),J=1,I),I=1,NSY)    INT17860
      DO 198  Q=1,NSX                                                   INT17870
      DO 198  P=1,Q                                                     INT17880
      DO 198  I=1,NSY                                                   INT17890
      DO 198  J=1,I                                                     INT17900
 198  F4(I,J,P,Q)=ZERO                                                  INT17910
      NND=NSY*(NSY+1)*(NSY+2)*NSX/6                                     INT17920
      NR=NND/NMM                                                        INT17930
      NL=NND-NMM*NR                                                     INT17940
      NM=NMM                                                            INT17950
      I=1                                                               INT17960
      J=1                                                               INT17970
      K=1                                                               INT17980
      Q=0                                                               INT17990
      REWIND ISCR5                                                      INT18000
      DO 200  II=1,NR+1                                                 INT18010
      IF(II.EQ.NR+1) NM=NL                                              INT18020
      READ(ISCR5,3) (V(IK),IK=1,NM)                                     INT18030
          DO 200  IK=1,NM                                               INT18040
          IF(Q.LT.NSX) THEN                                             INT18050
            Q=Q+1                                                       INT18060
            GO TO 202                                                   INT18070
          END IF                                                        INT18080
          IF(K.LT.J) THEN                                               INT18090
            K=K+1                                                       INT18100
            Q=1                                                         INT18110
            GO TO 202                                                   INT18120
          END IF                                                        INT18130
          IF(J.LT.I) THEN                                               INT18140
            J=J+1                                                       INT18150
            K=1                                                         INT18160
            Q=1                                                         INT18170
            GO TO 202                                                   INT18180
          END IF                                                        INT18190
            I=I+1                                                       INT18200
            J=1                                                         INT18210
            K=1                                                         INT18220
            Q=1                                                         INT18230
 202      CONTINUE                                                      INT18240
          IF(I.NE.J) THEN                                               INT18250
             IF(J.NE.K) THEN                                            INT18260
                DO 204  P=1,Q                                           INT18270
                F4(I,J,P,Q)=V(IK)*A(K,P)+F4(I,J,P,Q)                    INT18280
                F4(I,K,P,Q)=V(IK)*A(J,P)+F4(I,K,P,Q)                    INT18290
 204            F4(J,K,P,Q)=V(IK)*A(I,P)+F4(J,K,P,Q)                    INT18300
             ELSE                                                       INT18310
                DO 206  P=1,Q                                           INT18320
                F4(I,J,P,Q)=V(IK)*A(J,P)+F4(I,J,P,Q)                    INT18330
 206            F4(J,J,P,Q)=V(IK)*A(I,P)+F4(J,J,P,Q)                    INT18340
             END IF                                                     INT18350
          ELSE                                                          INT18360
             IF(J.NE.K) THEN                                            INT18370
                DO 208  P=1,Q                                           INT18380
                F4(I,I,P,Q)=V(IK)*A(K,P)+F4(I,I,P,Q)                    INT18390
 208            F4(I,K,P,Q)=V(IK)*A(I,P)+F4(I,K,P,Q)                    INT18400
             ELSE                                                       INT18410
                DO 210  P=1,Q                                           INT18420
 210            F4(I,I,P,Q)=V(IK)*A(I,P)+F4(I,I,P,Q)                    INT18430
             END IF                                                     INT18440
          END IF                                                        INT18450
 200  CONTINUE                                                          INT18460
      REWIND ISCR5                                                      INT18470
      WRITE(ISCR5,3) ((((F4(I,J,P,Q),P=1,Q),Q=1,NSX),J=1,I),I=1,NSY)    INT18480
      DO 212  Q=1,NSX                                                   INT18490
      DO 212  P=1,Q                                                     INT18500
      DO 212  N=1,P                                                     INT18510
      DO 212  I=1,NSY                                                   INT18520
 212  F4(I,N,P,Q)=ZERO                                                  INT18530
      NND=NSX*(NSX+1)*NSY*(NSY+1)/4                                     INT18540
      NR=NND/NMM                                                        INT18550
      NL=NND-NMM*NR                                                     INT18560
      NM=NMM                                                            INT18570
      P=0                                                               INT18580
      Q=0                                                               INT18590
      I=1                                                               INT18600
      J=1                                                               INT18610
      REWIND ISCR5                                                      INT18620
      DO 214  II=1,NR+1                                                 INT18630
      IF(II.EQ.NR+1) NM=NL                                              INT18640
      READ(ISCR5,3) (V(IK),IK=1,NM)                                     INT18650
          DO 214  IK=1,NM                                               INT18660
          IF(P.LT.Q) THEN                                               INT18670
            P=P+1                                                       INT18680
            GO TO 216                                                   INT18690
          END IF                                                        INT18700
          IF(Q.LT.NSX) THEN                                             INT18710
            Q=Q+1                                                       INT18720
            P=1                                                         INT18730
            GO TO 216                                                   INT18740
          END IF                                                        INT18750
          IF(J.LT.I) THEN                                               INT18760
            J=J+1                                                       INT18770
            P=1                                                         INT18780
            Q=1                                                         INT18790
            GO TO 216                                                   INT18800
          END IF                                                        INT18810
            I=I+1                                                       INT18820
            J=1                                                         INT18830
            P=1                                                         INT18840
            Q=1                                                         INT18850
 216      CONTINUE                                                      INT18860
          IF(I.NE.J) THEN                                               INT18870
                DO 218  N=1,P                                           INT18880
                F4(I,N,P,Q)=V(IK)*A(J,N)+F4(I,N,P,Q)                    INT18890
 218            F4(J,N,P,Q)=V(IK)*A(I,N)+F4(J,N,P,Q)                    INT18900
          ELSE                                                          INT18910
                DO 220  N=1,P                                           INT18920
 220            F4(I,N,P,Q)=V(IK)*A(I,N)+F4(I,N,P,Q)                    INT18930
          END IF                                                        INT18940
 214  CONTINUE                                                          INT18950
      REWIND ISCR5                                                      INT18960
      WRITE(ISCR5,3) ((((F4(I,N,P,Q),I=1,NSY),N=1,P),P=1,Q),Q=1,NSX)    INT18970
      DO 222  Q=1,NSX                                                   INT18980
      DO 222  P=1,Q                                                     INT18990
      DO 222  N=1,P                                                     INT19000
      DO 222  M=1,N                                                     INT19010
 222  F4(M,N,P,Q)=ZERO                                                  INT19020
      NND=NSX*(NSX+1)*(NSX+2)*NSY/6                                     INT19030
      NR=NND/NMM                                                        INT19040
      NL=NND-NMM*NR                                                     INT19050
      NM=NMM                                                            INT19060
      N=1                                                               INT19070
      P=1                                                               INT19080
      Q=1                                                               INT19090
      I=0                                                               INT19100
      REWIND ISCR5                                                      INT19110
      DO 224  II=1,NR+1                                                 INT19120
      IF(II.EQ.NR+1) NM=NL                                              INT19130
      READ(ISCR5,3) (V(IK),IK=1,NM)                                     INT19140
          DO 224  IK=1,NM                                               INT19150
          IF(I.LT.NSY) THEN                                             INT19160
            I=I+1                                                       INT19170
            GO TO 226                                                   INT19180
          END IF                                                        INT19190
          IF(N.LT.P) THEN                                               INT19200
            N=N+1                                                       INT19210
            I=1                                                         INT19220
            GO TO 226                                                   INT19230
          END IF                                                        INT19240
          IF(P.LT.Q) THEN                                               INT19250
            P=P+1                                                       INT19260
            N=1                                                         INT19270
            I=1                                                         INT19280
            GO TO 226                                                   INT19290
          END IF                                                        INT19300
            Q=Q+1                                                       INT19310
            P=1                                                         INT19320
            N=1                                                         INT19330
            I=1                                                         INT19340
 226      CONTINUE                                                      INT19350
          DO 228  M=1,N                                                 INT19360
 228      F4(M,N,P,Q)=V(IK)*A(I,M)+F4(M,N,P,Q)                          INT19370
 224  CONTINUE                                                          INT19380
      REWIND ISCR5                                                      INT19390
      CALL FILL4A(NC,NSX,F4)                                            INT19400
C     ****************************************                          INT19410
 1000 IF(LPRT(1,NPRT).LT.5) RETURN                                      INT19420
      IF(NDER.LE.1) RETURN                                              INT19430
      IF(NINV.LE.0) THEN                                                INT19440
         WRITE(IOUT,7)                                                  INT19450
      ELSE                                                              INT19460
         WRITE(IOUT,8)                                                  INT19470
      END IF                                                            INT19480
      WRITE(IOUT,9)                                                     INT19490
      DO 1015  M=1,NSX                                                  INT19500
 1015 WRITE(IOUT,5) (F2(M,N),N=1,NSX)                                   INT19510
      IF(NDER.LE.2) RETURN                                              INT19520
      WRITE(IOUT,10)                                                    INT19530
      DO 1020  M=1,NSX                                                  INT19540
      WRITE(IOUT,*) 'M=',M                                              INT19550
      DO 1020  N=1,NSX                                                  INT19560
 1020 WRITE(IOUT,5) (F3(M,N,P),P=1,NSX)                                 INT19570
      IF(NDER.LE.3) RETURN                                              INT19580
      WRITE(IOUT,12)                                                    INT19590
      DO 1030  M=1,NSX                                                  INT19600
      DO 1030  N=1,NSX                                                  INT19610
      WRITE(IOUT,*) 'M=',M,'N=',N                                       INT19620
      DO 1030  P=1,NSX                                                  INT19630
 1030 WRITE(IOUT,5) (F4(M,N,P,Q),Q=1,NSX)                               INT19640
      RETURN                                                            INT19650
      END                                                               INT19660
C     ////////////////////////////////////////////////////////////      INT20820
      SUBROUTINE MACHF2(NC,NS,NINV,F1,F2,X,V)                           INT20830
      IMPLICIT REAL*8 (A-H,O-Z)                                         INT20840
      INTEGER R                                                         INT20850
      DIMENSION F1(NC),F2(NC,NC),X(NC,NC),V(NC)                         INT20860
      COMMON /IO/ IIN1,IOUT,IIN2,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25,
     $   I31,I32,I33,I35,I36,I37,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10,
     $   ISCR11,ISCR12,ISCR13,ISCR14
   4  FORMAT(I5)                                                        INT20910
   5  FORMAT(20X,F20.10)                                                INT20920
      IF(NINV.GT.0) GO TO 30                                            INT20930
      DO 10  R=1,NS                                                     INT20940
         CALL XIN(NC,NS,X,R,ISCR1)                                      INT20950
         DO 20  N=1,NS                                                  INT20960
         DO 20  M=1,NS                                                  INT20970
 20      F2(M,N)=F2(M,N)-F1(R)*X(M,N)                                   INT20980
 10   CONTINUE                                                          INT20990
      RETURN                                                            INT21000
 30   REWIND I12                                                        INT21010
      READ(I12,4) R                                                     INT21020
      READ(I12,5)(V(R),R=1,NS)                                          INT21030
 40   DO 50  R=1,NS                                                     INT21040
         CALL XIN(NC,NC,X,-R,ISCR1)                                     INT21050
         DO 60  N=1,NC                                                  INT21060
         DO 60  M=1,NC                                                  INT21070
 60      F2(M,N)=F2(M,N)+V(R)*X(M,N)                                    INT21080
 50   CONTINUE                                                          INT21090
      RETURN                                                            INT21100
      END                                                               INT21110
C     ////////////////////////////////////////////////////////////      INT21120
C             (NUMB)                                                    INT21130
      SUBROUTINE MACHF3(NC,NS,NINV,F1,F3,YR,V)                          INT21140
      IMPLICIT REAL*8 (A-H,O-Z)                                         INT21150
      INTEGER M,N,P,R                                                   INT21160
      DIMENSION F1(NC),F3(NC,NC,NC)                                     INT21170
      DIMENSION YR(NC,NC,NC),V(NC)                                      INT21180
      COMMON /IO/ IIN1,IOUT,IIN2,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25,
     $   I31,I32,I33,I35,I36,I37,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10,
     $   ISCR11,ISCR12,ISCR13,ISCR14
   4  FORMAT(I5)                                                        INT21230
   5  FORMAT(20X,F20.10)                                                INT21240
      IF(NINV.GT.0) GO TO 200                                           INT21250
      DO 100  R=1,NS                                                    INT21260
         CALL YIN(NC,NS,YR,R,ISCR3)                                     INT21270
         DO 110  P=1,NS                                                 INT21280
         DO 110  N=1,P                                                  INT21290
         DO 110  M=1,N                                                  INT21300
 110     F3(M,N,P)=F3(M,N,P)-F1(R)*YR(M,N,P)                            INT21310
 100  CONTINUE                                                          INT21320
      CALL FILL3A(NC,NS,F3)                                             INT21330
      RETURN                                                            INT21340
 200  REWIND I12                                                        INT21350
      READ(I12,4) R                                                     INT21360
      READ(I12,5)(V(R),R=1,NS)                                          INT21370
      DO 210  R=1,NS                                                    INT21380
         CALL YIN(NC,NC,YR,-R,ISCR3)                                    INT21390
         DO 220  P=1,NC                                                 INT21400
         DO 220  N=1,P                                                  INT21410
         DO 220  M=1,N                                                  INT21420
 220     F3(M,N,P)=F3(M,N,P)+V(R)*YR(M,N,P)                             INT21430
 210  CONTINUE                                                          INT21440
      CALL FILL3A(NC,NC,F3)                                             INT21450
      RETURN                                                            INT21460
      END                                                               INT21470
C     //////////////////////////////////////////////////////////////    INT21850
      SUBROUTINE MACHF4(NC,NS,NINV,F1,F4,YR,V)                          YY 01930
      IMPLICIT REAL*8 (A-H,O-Z)                                         YY 01940
      INTEGER M,N,P,Q,R                                                 YY 01950
      DIMENSION F1(NC),F4(NC,NC,NC,NC)                                  YY 01960
      DIMENSION YR(NC,NC,NC),V(NC)                                      YY 01970
      COMMON /IO/ IIN1,IOUT,IIN2,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25,
     $   I31,I32,I33,I35,I36,I37,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10,
     $   ISCR11,ISCR12,ISCR13,ISCR14
   4  FORMAT(I5)                                                        YY 02020
   5  FORMAT(20X,F20.10)                                                YY 02030
      IF(NINV.GT.0) GO TO 200                                           YY 02040
      DO 100  R=1,NS                                                    YY 02050
         DO 110  Q=1,NS                                                 YY 02060
         CALL YIN2(NC,NS,YR,R,Q,ISCR9)                                  YY 02070
         DO 110  P=1,Q                                                  YY 02080
         DO 110  N=1,P                                                  YY 02090
         DO 110  M=1,N                                                  YY 02100
 110     F4(M,N,P,Q)=F4(M,N,P,Q)-F1(R)*YR(M,N,P)                        YY 02110
 100  CONTINUE                                                          YY 02120
      CALL FILL4A(NC,NS,F4)                                             YY 02130
      RETURN                                                            YY 02140
 200  REWIND I12                                                        YY 02150
      READ(I12,4) R                                                     YY 02160
      READ(I12,5)(V(R),R=1,NS)                                          YY 02170
      DO 210  R=1,NS                                                    YY 02180
         DO 220  Q=1,NC                                                 YY 02190
         CALL YIN2(NC,NC,YR,-R,Q,ISCR9)                                 YY 02200
         DO 220  P=1,Q                                                  YY 02210
         DO 220  N=1,P                                                  YY 02220
         DO 220  M=1,N                                                  YY 02230
 220     F4(M,N,P,Q)=F4(M,N,P,Q)+V(R)*YR(M,N,P)                         YY 02240
 210  CONTINUE                                                          YY 02250
      CALL FILL4A(NC,NC,F4)                                             YY 02260
      RETURN                                                            YY 02270
      END                                                               YY 02280
C     ///////////////////////////////////////////////////////////////// YY 02290
      SUBROUTINE XF2(NDER,NC,NS,NINV,BS,F2,F3,F4,V,XR,XS,XT)            INT21860
      IMPLICIT REAL*8 (A-H,O-Z)                                         INT21870
      INTEGER P,Q,R,S                                                   INT21880
      DIMENSION F2(NC,NC),F3(NC,NC,NC),F4(NC,NC,NC,NC)                  INT21890
      DIMENSION XR(NC,NC),XS(NC,NC),XT(NC,NC),BS(NC,NC),V(NC,NC)        INT21900
      PARAMETER(ZERO=0.0D0)                                             INT21910
      COMMON /IO/ IIN1,IOUT,IIN2,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25,
     $   I31,I32,I33,I35,I36,I37,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10,
     $   ISCR11,ISCR12,ISCR13,ISCR14
  1   FORMAT(2I5)                                                       INT21960
  2   FORMAT(3F20.10)                                                   INT21970
      IF(NINV.GT.0) GO TO 100                                           INT21980
      DO 10  R=1,NS                                                     INT21990
         CALL XIN(NC,NS,XR,R,ISCR1)                                     INT22000
         DO 20  P=1,NS                                                  INT22010
         DO 20  N=1,P                                                   INT22020
         DO 20  M=1,N                                                   INT22030
         W=F2(P,R)*XR(M,N)+F2(N,R)*XR(M,P)+F2(M,R)*XR(N,P)              INT22040
 20      F3(M,N,P)=F3(M,N,P)-W                                          INT22050
 10   CONTINUE                                                          INT22060
      CALL FILL3A(NC,NS,F3)                                             INT22070
      IF(NDER.LE.3) RETURN                                              INT22080
      DO 60  R=1,NS                                                     INT22090
         DO 30  N=1,NS                                                  INT22100
         DO 30  M=1,NS                                                  INT22110
 30      XR(M,N)=ZERO                                                   INT22120
           DO 35  S=1,NS                                                INT22130
           CALL XIN(NC,NS,XS,S,ISCR1)                                   INT22140
           IF(R.EQ.S) THEN                                              INT22150
             DO 40  N=1,NS                                              INT22160
             DO 40  M=1,NS                                              INT22170
 40          XT(M,N)=XS(M,N)                                            INT22180
           END IF                                                       INT22190
           DO 45  N=1,NS                                                INT22200
           DO 45  M=1,NS                                                INT22210
 45        XR(M,N)=XR(M,N)+F2(R,S)*XS(M,N)                              INT22220
 35        CONTINUE                                                     INT22230
         DO 55  Q=1,NS                                                  INT22240
         DO 55  P=1,Q                                                   INT22250
         DO 55  N=1,P                                                   INT22260
         DO 55  M=1,N                                                   INT22270
         W=XT(M,N)*XR(P,Q)+XT(M,P)*XR(N,Q)+XT(M,Q)*XR(N,P)              INT22280
 55      F4(M,N,P,Q)=F4(M,N,P,Q)-W                                      INT22290
 60      CONTINUE                                                       INT22300
      CALL FILL4A(NC,NS,F4)                                             INT22310
      RETURN                                                            INT22320
C                                                                       INT22330
 100  REWIND I16                                                        INT22340
      READ(I16,1) M,N                                                   INT22350
      READ(I16,2) ((V(M,N),N=1,NS),M=1,NS)                              INT22360
      DO 110  M=1,NS                                                    INT22370
      DO 110  J=1,NC                                                    INT22380
      XS(M,J)=ZERO                                                      INT22390
      DO 110  N=1,NS                                                    INT22400
 110  XS(M,J)=XS(M,J)+V(M,N)*BS(N,J)                                    INT22410
      DO 120  R=1,NS                                                    INT22420
         CALL XIN(NC,NC,XR,-R,ISCR1)                                    INT22430
         DO 125  K=1,NC                                                 INT22440
         DO 125  J=1,K                                                  INT22450
         DO 125  I=1,J                                                  INT22460
         W=XR(I,J)*XS(R,K)+XR(I,K)*XS(R,J)+XR(J,K)*XS(R,I)              INT22470
 125     F3(I,J,K)=F3(I,J,K)+W                                          INT22480
 120  CONTINUE                                                          INT22490
      CALL FILL3A(NC,NC,F3)                                             INT22500
      IF(NDER.LE.3) RETURN                                              INT22510
      DO 130  R=1,NS                                                    INT22520
         DO 140  J=1,NC                                                 INT22530
         DO 140  I=1,NC                                                 INT22540
 140     XR(I,J)=ZERO                                                   INT22550
         DO 145  S=1,NS                                                 INT22560
         CALL XIN(NC,NC,XS,-S,ISCR1)                                    INT22570
         IF(R.EQ.S) THEN                                                INT22580
            DO 150  J=1,NC                                              INT22590
            DO 150  I=1,NC                                              INT22600
 150        XT(I,J)=XS(I,J)                                             INT22610
         END IF                                                         INT22620
         DO 155  J=1,NC                                                 INT22630
         DO 155  I=1,NC                                                 INT22640
 155     XR(I,J)=XR(I,J)+V(R,S)*XS(I,J)                                 INT22650
 145     CONTINUE                                                       INT22660
         DO 160  L=1,NC                                                 INT22670
         DO 160  K=1,L                                                  INT22680
         DO 160  J=1,K                                                  INT22690
         DO 160  I=1,J                                                  INT22700
         W=XT(I,J)*XR(K,L)+XT(I,K)*XR(J,L)+XT(I,L)*XR(J,K)              INT22710
 160     F4(I,J,K,L)=F4(I,J,K,L)+W                                      INT22720
 130  CONTINUE                                                          INT22730
      CALL FILL4A(NC,NC,F4)                                             INT22740
      RETURN                                                            INT22750
      END                                                               INT22760
C     ////////////////////////////////////////////////////////////////  INT22770
      SUBROUTINE XF3(NC,NS,NINV,BS,XR,XS,XT,F3,F4,YR)                   INT22780
      IMPLICIT REAL*8 (A-H,O-Z)                                         INT22790
      INTEGER P,Q,R,S                                                   INT22800
      DIMENSION F3(NC,NC,NC),F4(NC,NC,NC,NC),YR(NC,NC,NC)               INT22810
      DIMENSION XR(NC,NC),XS(NC,NC),XT(NC,NC),BS(NC,NC)                 INT22820
      PARAMETER(ZERO=0.0D0)                                             INT22830
      COMMON /IO/ IIN1,IOUT,IIN2,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25,
     $   I31,I32,I33,I35,I36,I37,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10,
     $   ISCR11,ISCR12,ISCR13,ISCR14
  1   FORMAT(2I5)                                                       INT22880
  2   FORMAT(3F20.10)                                                   INT22890
      IF(NINV.GT.0) GO TO 100                                           INT22900
      DO 10  R=1,NS                                                     INT22910
         CALL XIN(NC,NS,XR,R,ISCR1)                                     INT22920
         DO 20  Q=1,NS                                                  INT22930
         DO 20  P=1,Q                                                   INT22940
         DO 20  N=1,P                                                   INT22950
         DO 20  M=1,N                                                   INT22960
         W=F3(M,N,R)*XR(P,Q)+F3(M,P,R)*XR(N,Q)+F3(M,Q,R)*XR(N,P)        INT22970
         W=W+F3(N,P,R)*XR(M,Q)+F3(N,Q,R)*XR(M,P)+F3(P,Q,R)*XR(M,N)      INT22980
 20      F4(M,N,P,Q)=F4(M,N,P,Q)-W                                      INT22990
 10   CONTINUE                                                          INT23000
      CALL FILL4A(NC,NS,F4)                                             INT23010
      RETURN                                                            INT23020
C                                                                       INT23030
 100  REWIND I21                                                        INT23040
      READ(I21,1) M,N                                                   INT23050
      READ(I21,2) (((YR(M,N,P),P=1,N),N=1,M),M=1,NS)                    INT23060
      CALL FILL3B(NC,NS,YR)                                             INT23070
      DO 110  R=1,NS                                                    INT23080
         CALL XIN(NC,NC,XR,-R,ISCR1)                                    INT23090
         DO 120  L=1,NC                                                 INT23100
         DO 120  N=1,NC                                                 INT23110
           XT(N,L)=ZERO                                                 INT23120
           DO 120  P=1,NC                                               INT23130
 120       XT(N,L)=XT(N,L)+YR(R,N,P)*BS(P,L)                            INT23140
         DO 125  L=1,NC                                                 INT23150
         DO 125  K=1,NC                                                 INT23160
           XS(K,L)=ZERO                                                 INT23170
           DO 125  N=1,NC                                               INT23180
 125       XS(K,L)=XS(K,L)+XT(N,L)*BS(N,K)                              INT23190
         DO 130  L=1,NC                                                 INT23200
         DO 130  K=1,L                                                  INT23210
         DO 130  J=1,K                                                  INT23220
         DO 130  I=1,J                                                  INT23230
         W=XR(I,J)*XS(K,L)+XR(I,K)*XS(J,L)+XR(J,K)*XS(I,L)              INT23240
         W=W+XR(K,L)*XS(I,J)+XR(J,L)*XS(I,K)+XR(I,L)*XS(J,K)            INT23250
 130     F4(I,J,K,L)=F4(I,J,K,L)+W                                      INT23260
 110  CONTINUE                                                          INT23270
      CALL FILL4A(NC,NC,F4)                                             INT23280
      RETURN                                                            INT23290
      END                                                               INT23300
C     ////////////////////////////////////////////////////////////////  INT23310
C     (OPT1) (NUMB)                                                     INT23320
      SUBROUTINE YF2(NC,NS,NINV,BS,XR,XS,F2,F4,YR)                      INT23330
      IMPLICIT REAL*8 (A-H,O-Z)                                         INT23340
      INTEGER P,Q,R,S                                                   INT23350
      DIMENSION F2(NC,NC),F4(NC,NC,NC,NC)                               INT23360
      DIMENSION YR(NC,NC,NC),BS(NC,NC),XR(NC,NC),XS(NC,NC)              INT23370
      PARAMETER(ZERO=0.0D0)                                             INT23380
      COMMON /IO/ IIN1,IOUT,IIN2,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25,
     $   I31,I32,I33,I35,I36,I37,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10,
     $   ISCR11,ISCR12,ISCR13,ISCR14
  1   FORMAT(2I5)                                                       INT23430
  2   FORMAT(3F20.10)                                                   INT23440
      IF(NINV.GT.0) GO TO 200                                           INT23450
      DO 110  R=1,NS                                                    INT23460
         CALL YIN(NC,NS,YR,R,ISCR3)                                     INT23470
         DO 120  Q=1,NS                                                 INT23480
         DO 120  P=1,Q                                                  INT23490
         DO 120  N=1,P                                                  INT23500
         DO 120  M=1,N                                                  INT23510
         W=F2(M,R)*YR(N,P,Q)+F2(N,R)*YR(M,P,Q)                          INT23520
         W=W+F2(P,R)*YR(M,N,Q)+F2(Q,R)*YR(M,N,P)                        INT23530
 120     F4(M,N,P,Q)=F4(M,N,P,Q)-W                                      INT23540
 110  CONTINUE                                                          INT23550
      CALL FILL4A(NC,NS,F4)                                             INT23560
      RETURN                                                            INT23570
C                                                                       INT23580
 200  REWIND I16                                                        INT23590
      READ(I16,1) M,N                                                   INT23600
      READ(I16,2) ((XS(M,N),N=1,NS),M=1,NS)                             INT23610
      DO 210  M=1,NS                                                    INT23620
      DO 210  J=1,NC                                                    INT23630
      XR(M,J)=ZERO                                                      INT23640
      DO 210  N=1,NS                                                    INT23650
 210  XR(M,J)=XR(M,J)+XS(M,N)*BS(N,J)                                   INT23660
      DO 220  R=1,NS                                                    INT23670
         CALL YIN(NC,NC,YR,-R,ISCR3)                                    INT23680
         DO 230  L=1,NC                                                 INT23690
         DO 230  K=1,L                                                  INT23700
         DO 230  J=1,K                                                  INT23710
         DO 230  I=1,J                                                  INT23720
         W=YR(I,J,K)*XR(R,L)+YR(I,J,L)*XR(R,K)                          INT23730
         W=W+YR(I,K,L)*XR(R,J)+YR(J,K,L)*XR(R,I)                        INT23740
 230     F4(I,J,K,L)=F4(I,J,K,L)+W                                      INT23750
 220  CONTINUE                                                          INT23760
      CALL FILL4A(NC,NC,F4)                                             INT23770
      RETURN                                                            INT23780
      END                                                               INT23790
C     ////////////////////////////////////////////////////////////      INT14360
      SUBROUTINE MACHB(NAD,NC,NS,XA,XMASS,TYPE,IA,B,S)                  INT05310
      IMPLICIT REAL*8 (A-H,O-Z)                                         INT05320
      INTEGER IL(5)                                                     INT05330
      CHARACTER TYPE*5                                                  INT05340
      DIMENSION TYPE(NS),IA(NS,6),B(NS,NC),XA(NAD,3),S(NS)              INT05350
      DIMENSION V(3,5),Q(3,5),Z(3,5),XMASS(1)                           INT05360
      PARAMETER(ZERO=0.0D0,ONE=1.0D0)                                   INT05370
      DO 100  I=1,NS                                                    INT05380
         DO 5  J=1,NC                                                   INT05390
 5       B(I,J)=ZERO                                                    INT05400
           NIA=0                                                        INT05410
           K1=IA(I,1)                                                   INT05420
           K2=IA(I,2)                                                   INT05430
           K3=IA(I,3)                                                   INT05440
           K4=IA(I,4)                                                   INT05450
           K5=IA(I,5)                                                   INT05460
           IL(1)=3*(K1-1)                                               INT05470
           IL(2)=3*(K2-1)                                               INT05480
           IL(3)=3*(K3-1)                                               INT05490
           IL(4)=3*(K4-1)                                               INT05500
           IL(5)=3*(K5-1)                                               INT05510
      IF(TYPE(I).NE.' STRE') GO TO 10                                   INT05520
         CALL VECT1(NAD,K1,K2,V(1,1),XA,W)                              INT05530
         DO 6  K=1,3                                                    INT05540
  6      V(K,2)=-V(K,1)                                                 INT05550
         NIA=2                                                          INT05560
         GO TO 80                                                       INT05570
 10   IF(TYPE(I).NE.' BEND') GO TO 15                                   INT05580
         CALL VECT2(NAD,K1,K2,K3,V(1,1),V(1,2),V(1,3),XA,W)             INT05590
         NIA=3                                                          INT05600
         GO TO 80                                                       INT05610
 15   IF(TYPE(I).NE.' LIN1') GO TO 25                                   INT05620
         CALL VECT3(NAD,K1,K2,K3,K4,V(1,1),V(1,2),V(1,3),XA,W)          INT05630
         NIA=3                                                          INT05640
         GO TO 80                                                       INT05650
 25   IF(TYPE(I).NE.'  OUT') GO TO 30                                   INT05660
         CALL VECT5(NAD,K1,K2,K3,K4,V(1,1),V(1,2),V(1,3),V(1,4),XA,W)   INT05670
         NIA=4                                                          INT05680
         GO TO 80                                                       INT05690
 30   IF(TYPE(I).NE.' TORS') GO TO 35                                   INT05700
         CALL VECT6(NAD,K1,K2,K3,K4,V(1,1),V(1,2),V(1,3),V(1,4),XA,W)   INT05710
         NIA=4                                                          INT05720
         GO TO 80                                                       INT05730
 35   IF(TYPE(I).NE.'  SPF') GO TO 40                                   INT05740
         CALL VECT1(NAD,K1,K2,V(1,1),XA,W)                              INT05750
         DO 36  K=1,3                                                   INT05760
 36      V(K,2)=-V(K,1)                                                 INT05770
         NIA=2                                                          INT05780
         FACT=S(I)/(W*W)                                                INT05790
         W=ONE-S(I)/W                                                   INT05800
         DO 38  L=1,2                                                   INT05810
         DO 38  K=1,3                                                   INT05820
 38      V(K,L)=V(K,L)*FACT                                             INT05830
         GO TO 80                                                       INT05840
 40   IF(TYPE(I).NE.' LINX') GO TO 45
         CALL VECT8(NAD,K1,K2,K3,K4,V(1,1),V(1,2),V(1,3),V(1,4),XA,W)   INT05710
         NIA=4                                                          INT05720
         GO TO 80                                                       INT05730
 45   IF(TYPE(I).NE.' LINY') GO TO 50
         CALL VECT9(NAD,K1,K2,K3,K4,V(1,1),V(1,2),V(1,3),V(1,4),XA,W)   INT05710
         NIA=4                                                          INT05720
         GO TO 80                                                       INT05730
 50   IF(TYPE(I).NE.' RCOM') GO TO 55
         CALL VECT10(NAD,K1,K2,K3,K4,V(1,1),XA,XMASS,XMA,XMB,W) 
         DO 51  K=K1,K2
           L=3*(K-1)
           DO 51  J=1,3
 51        B(I,L+J)=XMASS(K)*V(J,1)/XMA
         DO 52  K=K3,K4
           L=3*(K-1)
           DO 52  J=1,3
 52        B(I,L+J)=-XMASS(K)*V(J,1)/XMB
         GO TO 90
 55   CONTINUE                                                          INT05850
 80   DO 85  K=1,NIA                                                    INT05860
      DO 85  J=1,3                                                      INT05870
 85   B(I,IL(K)+J)=V(J,K)                                               INT05880
 90   S(I)=W                                                            INT05890
 100  CONTINUE                                                          INT05900
      RETURN                                                            INT05910
      END                                                               INT05920
C     ///////////////////////////////////////////////////////////////   INT46820
      SUBROUTINE MACHB2(NAD,NC,XA,B2)                                   INT46830
      IMPLICIT REAL*8 (A-H,O-Z)                                         INT46840
      DIMENSION XA(NAD,3),B2(6,NC)                                      INT46850
      PARAMETER(ZERO=0.0D0,ONE=1.0D0)                                   INT46860
      NA=NC/3                                                           INT46870
      T=DSQRT(DBLE(NA))                                                 INT46880
      DO 10  J=1,NC                                                     INT46890
      DO 10  I=1,6                                                      INT46900
 10   B2(I,J)=ZERO                                                      INT46910
      DO 15  K=1,NA                                                     INT46920
      DO 15  I=1,3                                                      INT46930
      L=3*(K-1)+I                                                       INT46940
 15   B2(I,L)=ONE/T                                                     INT46950
      DO 20  K=1,NA                                                     INT46960
      L=3*K                                                             INT46970
      B2(4,L-1)=-XA(K,3)                                                INT46980
      B2(4,L)=XA(K,2)                                                   INT46990
      B2(5,L)=-XA(K,1)                                                  INT47000
      B2(5,L-2)=XA(K,3)                                                 INT47010
      B2(6,L-1)=XA(K,1)                                                 INT47020
 20   B2(6,L-2)=-XA(K,2)                                                INT47030
      DO 25  I=4,6                                                      INT47040
      T=ZERO                                                            INT47050
      DO 30  J=1,NC                                                     INT47060
 30   T=T+B2(I,J)*B2(I,J)                                               INT47070
      T=DSQRT(T)                                                        INT47080
      DO 35  J=1,NC                                                     INT47090
 35   B2(I,J)=B2(I,J)/T                                                 INT47100
 25   CONTINUE                                                          INT47110
      RETURN                                                            INT47120
      END                                                               INT47130
C     ///////////////////////////////////////////////////////////////// YY 01920
      SUBROUTINE BINVRT(NS,NC,XMASS,B,D,A,AV,V,IFLAG,TOLINV,ITST)       INT13740
      IMPLICIT REAL*8 (A-H,O-Z)                                         INT13750
      DIMENSION B(NC,NC),A(NC,NC),D(NC,1),XMASS(1)                      INT13760
      DIMENSION AV(1),V(NC,3)
C  Note that AV is stored in the second part of D.
      PARAMETER(ZERO=0.0D0,ONE=1.0D0)                                   INT13770
      COMMON /IO/ IIN1,IOUT,IIN2,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25,
     $   I31,I32,I33,I35,I36,I37,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10,
     $   ISCR11,ISCR12,ISCR13,ISCR14
 1    FORMAT(//' B*BT MATRIX FOR (SYMMETRY) INTERNAL COORDINATES',/)    INT13820
 2    FORMAT(8F12.6)                                                    INT13830
 3    FORMAT(/1X,'DETERMINANT OF B*BT MATRIX=',G12.5)                   INT13840
 4    FORMAT(//,' A MATRIX FOR (SYMMETRY) INTERNAL COORDINATES',/)      INT13850
 5    FORMAT(8F12.6)                                                    INT13860
 6    FORMAT(//,' NORMALIZED OVERLAP MATRIX FOR INTERNAL COORDINATES',/)
 7    FORMAT(8F12.6)
 8    FORMAT(//,' EIGENVECTORS AND EIGENVALUES OF',
     $          ' NORMALIZED OVERLAP MATRIX')
 9    FORMAT(//,' DETERMINANT OF OVERLAP MATRIX   ', E12.5)
      IFLAG=0                                                           INT13870
      DO 10  J=1,NS                                                     INT13880
      DO 10  I=1,NS                                                     INT13890
      D(I,J)=ZERO                                                       INT13900
      D(I,J+NS)=ZERO                                                    INT13910
      DO 15  K=1,NC                                                     INT13920
      IK=(K-1)/3+1                                                      INT13930
 15   D(I,J)=D(I,J)+B(I,K)*B(J,K)/XMASS(IK)                             INT13940
 10   CONTINUE                                                          INT13950
      IF(LPRT(1,NPRT).GE.4) THEN                                        INT13960
         WRITE(IOUT,1)                                                  INT13970
         DO 20  I=1,NS                                                  INT13980
         WRITE(IOUT,*) 'I=',I                                           INT13990
 20      WRITE(IOUT,2) (D(I,J),J=1,NS)                                  INT14000
      END IF                                                            INT14010
      DO 25  K=1,NS                                                     INT14020
 25   D(K,K+NS)=ONE                                                     INT14030
      CALL FLIN(D,NC,NS,NS,DET)                                         INT14040
      WRITE(IOUT,3) DET                                                 INT14050
      IF(DET.EQ.ZERO) THEN                                              INT14060
         IFLAG=1                                                        INT14070
         RETURN                                                         INT14080
      END IF                                                            INT14090
      DO 35  J=1,NS                                                     INT14100
      DO 35  I=1,NC                                                     INT14110
      A(I,J)=ZERO                                                       INT14120
      DO 40  K=1,NS                                                     INT14130
 40   A(I,J)=A(I,J)+B(K,I)*D(K,J+NS)                                    INT14140
      IK=(I-1)/3+1                                                      INT14150
      A(I,J)=A(I,J)/XMASS(IK)                                           INT14160
 35   CONTINUE                                                          INT14170
      IF(LPRT(1,NPRT).GE.3) THEN                                        INT14180
         WRITE(IOUT,4)                                                  INT14190
         DO 45  I=1,NC                                                  INT14200
         WRITE(IOUT,*) 'I=',I                                           INT14210
 45      WRITE(IOUT,5) (A(I,J),J=1,NS)                                  INT14220
      END IF                                                            INT14230
      IF(ITST.EQ.0) RETURN                                              INT14240
CSA   IF(ITST.GT.1) GO TO 60
      DO 50  J=1,NS                                                     INT14250
      DO 50  I=1,NS                                                     INT14260
      W=ZERO                                                            INT14270
      DO 55  K=1,NC                                                     INT14280
 55   W=W+B(I,K)*A(K,J)                                                 INT14290
      D(I,J)=W                                                          INT14300
      IF(I.EQ.J) W=W-ONE                                                INT14310
      IF(DABS(W).GT.TOLINV) IFLAG=2                                     INT14320
 50   CONTINUE
 60   IF(ITST.LT.2) GO TO 80
      DO 62  J=1,NS                                                     INT13880
      DO 62  I=1,NS                                                     INT13890
      D(I,J)=ZERO                                                       INT13900
      DO 63  K=1,NC                                                     INT13920
      IK=(K-1)/3+1                                                      INT13930
 63   D(I,J)=D(I,J)+B(I,K)*B(J,K)/XMASS(IK)                             INT13940
 62   CONTINUE                                                          INT13950
      DO 64  I=1,NS
      DO 64  J=1,NS
      IF(I.EQ.J) GO TO 64
      D(I,J)=D(I,J)/DSQRT(D(I,I)*D(J,J))
 64   CONTINUE
      DO 65  I=1,NS
 65   D(I,I)=ONE
      WRITE(IOUT,6)
      DO 66  I=1,NS
 66   WRITE(IOUT,7) (D(I,J),J=1,NS)
      II=0
      DO 67  I=1,NS
      DO 67  J=1,I
      II=II+1
 67   AV(II)=D(I,J)
      NV=NS*(NS+1)/2
      CALL RSP(NC,NS,NV,AV,V(1,1),1,D,V(1,2),V(1,3))
      WRITE(IOUT,8)
      CALL TABLE5(NC,NS,NS,V(1,1),D)
      DETS=ONE
      DO  68  I=1,NS
 68   DETS=DETS*V(I,1)
      WRITE(IOUT,9) DETS
 80   CONTINUE
      RETURN                                                            
      END                                                               
C     ///////////////////////////////////////////////////////////////   INT05930
      SUBROUTINE MACHX(NAD,NC,NS,NSX,IOPT,XA,XMASS,
     $                  TYPE,IA,A,S,U,IU,X,SR) 
      IMPLICIT REAL*8 (A-H,O-Z)                                         YY 00020
      CHARACTER TYPE*5                                                  YY 00030
      INTEGER R,RP,NSX,IOPT(30)                                         YY 00040
      DIMENSION TYPE(NS),IA(NS,6),A(NC,NC),X(NC,NC),SR(NC,NC)           YY 00050
      DIMENSION XA(NAD,3),S(NS),XMASS(1),U(NS,1),IU(NS,0:1)             YY 00060
      DIMENSION H11(3,3),H21(3,3),H22(3,3),H31(3,3),H32(3,3)            YY 00070
      DIMENSION H33(3,3),H41(3,3),H42(3,3),H43(3,3),H44(3,3)            YY 00080
      DIMENSION E21(3)                                                  YY 00090
      PARAMETER(ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0,THREE=3.0D0)             YY 00100
      COMMON /IO/ IIN1,IOUT,IIN2,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25,
     $   I31,I32,I33,I35,I36,I37,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10,
     $   ISCR11,ISCR12,ISCR13,ISCR14
    1 FORMAT(/,1X,'NUMERICAL SR(I,J) AND X(M,N) MATRICES USED FOR',/,   YY 00150
     $  ' SIMPLE INTERNAL COORDINATE',I5)                               YY 00160
    2 FORMAT(/,1X,'SR(I,J) AND X(M,N) MATRICES SET TO ZERO FOR',/,      YY 00170
     $  ' SIMPLE INTERNAL COORDINATE',I5)                               YY 00180
      NSYM=IOPT(3)                                                      YY 00190
      IF(NSYM.NE.0) THEN                                                YY 00200
         ISCR=ISCR2                                                     YY 00210
      ELSE                                                              YY 00220
         ISCR=ISCR1                                                     YY 00230
      END IF                                                            YY 00240
C     ***********************************                               YY 00250
      DO 500  R=1,NS                                                    YY 00260
         DO 10  N=1,NC                                                  YY 00270
         DO 10  M=1,NC                                                  YY 00280
 10      X(M,N)=ZERO                                                    YY 00290
C  OPTION                                                               YY 00300
         DO 49  J=1,NC                                                  YY 00310
         DO 49  I=1,NC                                                  YY 00320
 49      SR(I,J)=ZERO                                                   YY 00330
C  **************                                                       YY 00340
         IF(IA(R,6).EQ.-1) GO TO 275                                    YY 00350
         IF(TYPE(R).NE.' STRE') GO TO 25                                YY 00360
         K1=IA(R,1)                                                     YY 00370
         K2=IA(R,2)                                                     YY 00380
         L1=3*(K1-1)                                                    YY 00390
         L2=3*(K2-1)                                                    YY 00400
         CALL HIJS1(NAD,K1,K2,XA,H11)                                   YY 00410
C   OPTION                                                              YY 00420
              DO 28  J=1,3                                              YY 00430
              DO 28  I=1,3                                              YY 00440
              SR(L1+I,L1+J)=H11(I,J)                                    YY 00450
              SR(L2+I,L2+J)=H11(I,J)                                    YY 00460
              SR(L1+I,L2+J)=-H11(I,J)                                   YY 00470
 28           SR(L2+I,L1+J)=-H11(I,J)                                   YY 00480
C  *************                                                        YY 00490
         CALL AHX2(NC,NSX,L1,L2,H11,A,X)                                YY 00500
         GO TO 300                                                      YY 00510
 25      IF(TYPE(R).NE.' BEND') GO TO 75                                YY 00520
         K1=IA(R,1)                                                     YY 00530
         K2=IA(R,2)                                                     YY 00540
         K3=IA(R,3)                                                     YY 00550
         L1=3*(K1-1)                                                    YY 00560
         L2=3*(K2-1)                                                    YY 00570
         L3=3*(K3-1)                                                    YY 00580
         CALL HIJS2(NAD,K1,K2,K3,XA,H11,H21,H31,H22,H32,H33)            YY 00590
C  OPTION                                                               YY 00600
      DO 48  J=1,3                                                      YY 00610
      DO 48  I=1,3                                                      YY 00620
      SR(L1+I,L1+J)=H11(I,J)                                            YY 00630
      SR(L2+I,L1+J)=H21(I,J)                                            YY 00640
      SR(L3+I,L1+J)=H31(I,J)                                            YY 00650
      SR(L1+I,L2+J)=H21(J,I)                                            YY 00660
      SR(L2+I,L2+J)=H22(I,J)                                            YY 00670
      SR(L3+I,L2+J)=H32(I,J)                                            YY 00680
      SR(L1+I,L3+J)=H31(J,I)                                            YY 00690
      SR(L2+I,L3+J)=H32(J,I)                                            YY 00700
 48   SR(L3+I,L3+J)=H33(I,J)                                            YY 00710
C     ***********************************                               YY 00720
         CALL AHX3(NC,NSX,L1,L2,L3,H11,H21,H31,H22,H32,H33,A,X)         YY 00730
         GO TO 300                                                      YY 00740
 75      IF(TYPE(R).NE.' LIN1') GO TO 100                               YY 00750
         K1=IA(R,1)                                                     YY 00760
         K2=IA(R,2)                                                     YY 00770
         K3=IA(R,3)                                                     YY 00780
         K4=IA(R,4)                                                     YY 00790
         L1=3*(K1-1)                                                    YY 00800
         L2=3*(K2-1)                                                    YY 00810
         L3=3*(K3-1)                                                    YY 00820
         CALL HIJS3(NAD,K1,K2,K3,K4,XA,H11,H21,H31,H22,H32,H33)         YY 00830
C  OPTION                                                               YY 00840
      DO 86  J=1,3                                                      YY 00850
      DO 86  I=1,3                                                      YY 00860
      SR(L1+I,L1+J)=H11(I,J)                                            YY 00870
      SR(L2+I,L1+J)=H21(I,J)                                            YY 00880
      SR(L3+I,L1+J)=H31(I,J)                                            YY 00890
      SR(L1+I,L2+J)=H21(J,I)                                            YY 00900
      SR(L2+I,L2+J)=H22(I,J)                                            YY 00910
      SR(L3+I,L2+J)=H32(I,J)                                            YY 00920
      SR(L1+I,L3+J)=H31(J,I)                                            YY 00930
      SR(L2+I,L3+J)=H32(J,I)                                            YY 00940
 86   SR(L3+I,L3+J)=H33(I,J)                                            YY 00950
C     ***********************************                               YY 00960
         CALL AHX3(NC,NSX,L1,L2,L3,H11,H21,H31,H22,H32,H33,A,X)         YY 00970
         GO TO 300                                                      YY 00980
 100     IF(TYPE(R).NE.'  SPF') GO TO 125                               YY 00990
              K1=IA(R,1)                                                YY 01000
              K2=IA(R,2)                                                YY 01010
              L1=3*(K1-1)                                               YY 01020
              L2=3*(K2-1)                                               YY 01030
              CALL VECT1(NAD,K1,K2,E21,XA,T21)                          YY 01040
              CALL HIJS1(NAD,K1,K2,XA,H11)                              YY 01050
              FACT1=(ONE-S(R))/T21                                      YY 01060
              FACT2=TWO*FACT1/T21                                       YY 01070
              DO 102  J=1,3                                             YY 01080
              DO 104  I=1,3                                             YY 01090
 104          H11(I,J)=THREE*H11(I,J)*FACT1                             YY 01100
 102          H11(J,J)=H11(J,J)-FACT2                                   YY 01110
C   OPTION                                                              YY 01120
              DO 108  J=1,3                                             YY 01130
              DO 108  I=1,3                                             YY 01140
              SR(L1+I,L1+J)=H11(I,J)                                    YY 01150
              SR(L2+I,L2+J)=H11(I,J)                                    YY 01160
              SR(L1+I,L2+J)=-H11(I,J)                                   YY 01170
 108          SR(L2+I,L1+J)=-H11(I,J)                                   YY 01180
C  *************                                                        YY 01190
         CALL AHX2(NC,NSX,L1,L2,H11,A,X)                                YY 01200
         GO TO 300                                                      YY 01210
 125     IF(TYPE(R).NE.' TORS') GO TO 150                               YY 01220
         K1=IA(R,1)                                                     YY 01230
         K2=IA(R,2)                                                     YY 01240
         K3=IA(R,3)                                                     YY 01250
         K4=IA(R,4)                                                     YY 01260
         L1=3*(K1-1)                                                    YY 01270
         L2=3*(K2-1)                                                    YY 01280
         L3=3*(K3-1)                                                    YY 01290
         L4=3*(K4-1)                                                    YY 01300
         CALL HIJS6(NAD,K1,K2,K3,K4,XA,H11,H21,H31,H41,H22,H32,H42,     YY 01310
     $               H33,H43,H44)                                       YY 01320
C  OPTION                                                               YY 01330
      DO 148  J=1,3                                                     YY 01340
      DO 148  I=1,3                                                     YY 01350
      SR(L1+I,L1+J)=H11(I,J)                                            YY 01360
      SR(L2+I,L1+J)=H21(I,J)                                            YY 01370
      SR(L3+I,L1+J)=H31(I,J)                                            YY 01380
      SR(L4+I,L1+J)=H41(I,J)                                            YY 01390
      SR(L1+I,L2+J)=H21(J,I)                                            YY 01400
      SR(L2+I,L2+J)=H22(I,J)                                            YY 01410
      SR(L3+I,L2+J)=H32(I,J)                                            YY 01420
      SR(L4+I,L2+J)=H42(I,J)                                            YY 01430
      SR(L1+I,L3+J)=H31(J,I)                                            YY 01440
      SR(L2+I,L3+J)=H32(J,I)                                            YY 01450
      SR(L3+I,L3+J)=H33(I,J)                                            YY 01460
      SR(L4+I,L3+J)=H43(I,J)                                            YY 01470
      SR(L1+I,L4+J)=H41(J,I)                                            YY 01480
      SR(L2+I,L4+J)=H42(J,I)                                            YY 01490
      SR(L3+I,L4+J)=H43(J,I)                                            YY 01500
 148  SR(L4+I,L4+J)=H44(I,J)                                            YY 01510
C     ***********************************                               YY 01520
         CALL AHX4(NC,NSX,L1,L2,L3,L4,H11,H21,H31,H41,H22,H32,H42,      YY 01530
     $             H33,H43,H44,A,X)                                     YY 01540
         GO TO 300                                                      YY 01550
 150     IF(TYPE(R).NE.'  OUT') GO TO 175                               YY 01570
         K1=IA(R,1)                                                     YY 01580
         K2=IA(R,2)                                                     YY 01590
         K3=IA(R,3)                                                     YY 01600
         K4=IA(R,4)                                                     YY 01610
         L1=3*(K1-1)                                                    YY 01620
         L2=3*(K2-1)                                                    YY 01630
         L3=3*(K3-1)                                                    YY 01640
         L4=3*(K4-1)                                                    YY 01650
         CALL HIJS7(NAD,K1,K2,K3,K4,XA,H11,H21,H31,H41,H22,H32,H42,     YY 01660
     $               H33,H43,H44)                                       YY 01670
C  OPTION                                                               YY 01680
      DO 158  J=1,3                                                     YY 01690
      DO 158  I=1,3                                                     YY 01700
      SR(L1+I,L1+J)=H11(I,J)                                            YY 01710
      SR(L2+I,L1+J)=H21(I,J)                                            YY 01720
      SR(L3+I,L1+J)=H31(I,J)                                            YY 01730
      SR(L4+I,L1+J)=H41(I,J)                                            YY 01740
      SR(L1+I,L2+J)=H21(J,I)                                            YY 01750
      SR(L2+I,L2+J)=H22(I,J)                                            YY 01760
      SR(L3+I,L2+J)=H32(I,J)                                            YY 01770
      SR(L4+I,L2+J)=H42(I,J)                                            YY 01780
      SR(L1+I,L3+J)=H31(J,I)                                            YY 01790
      SR(L2+I,L3+J)=H32(J,I)                                            YY 01800
      SR(L3+I,L3+J)=H33(I,J)                                            YY 01810
      SR(L4+I,L3+J)=H43(I,J)                                            YY 01820
      SR(L1+I,L4+J)=H41(J,I)                                            YY 01830
      SR(L2+I,L4+J)=H42(J,I)                                            YY 01840
      SR(L3+I,L4+J)=H43(J,I)                                            YY 01850
 158  SR(L4+I,L4+J)=H44(I,J)                                            YY 01860
C     ***********************************                               YY 01870
         CALL AHX4(NC,NSX,L1,L2,L3,L4,H11,H21,H31,H41,H22,H32,H42,      YY 01880
     $             H33,H43,H44,A,X)                                     YY 01890
         GO TO 300                                                      YY 01900
 175     IF(TYPE(R).NE.' LINX'.AND.TYPE(R).NE.' LINY') GO TO 200        YY 01220
         K1=IA(R,1)                                                     YY 01230
         K2=IA(R,2)                                                     YY 01240
         K3=IA(R,3)                                                     YY 01250
         K4=IA(R,4)                                                     YY 01260
         L1=3*(K1-1)                                                    YY 01270
         L2=3*(K2-1)                                                    YY 01280
         L3=3*(K3-1)                                                    YY 01290
         L4=3*(K4-1)                                                    YY 01300
         IF(TYPE(R).EQ.' LINX') CALL HIJS8(NAD,K1,K2,K3,K4,XA,H11,
     $       H21,H31,H41,H22,H32,H42,H33,H43,H44)
         IF(TYPE(R).EQ.' LINY') CALL HIJS9(NAD,K1,K2,K3,K4,XA,H11,
     $       H21,H31,H41,H22,H32,H42,H33,H43,H44)
C  OPTION                                                               YY 01330
      DO 178  J=1,3                                                     YY 01340
      DO 178  I=1,3                                                     YY 01350
      SR(L1+I,L1+J)=H11(I,J)                                            YY 01360
      SR(L2+I,L1+J)=H21(I,J)                                            YY 01370
      SR(L3+I,L1+J)=H31(I,J)                                            YY 01380
      SR(L4+I,L1+J)=H41(I,J)                                            YY 01390
      SR(L1+I,L2+J)=H21(J,I)                                            YY 01400
      SR(L2+I,L2+J)=H22(I,J)                                            YY 01410
      SR(L3+I,L2+J)=H32(I,J)                                            YY 01420
      SR(L4+I,L2+J)=H42(I,J)                                            YY 01430
      SR(L1+I,L3+J)=H31(J,I)                                            YY 01440
      SR(L2+I,L3+J)=H32(J,I)                                            YY 01450
      SR(L3+I,L3+J)=H33(I,J)                                            YY 01460
      SR(L4+I,L3+J)=H43(I,J)                                            YY 01470
      SR(L1+I,L4+J)=H41(J,I)                                            YY 01480
      SR(L2+I,L4+J)=H42(J,I)                                            YY 01490
      SR(L3+I,L4+J)=H43(J,I)                                            YY 01500
 178  SR(L4+I,L4+J)=H44(I,J)                                            YY 01510
C     ***********************************                               YY 01520
         CALL AHX4(NC,NSX,L1,L2,L3,L4,H11,H21,H31,H41,H22,H32,H42,      YY 01530
     $             H33,H43,H44,A,X)                                     YY 01540
         GO TO 300                                                      YY 01550
 200     IF(TYPE(R).NE.' RCOM') GO TO 225
         K1=IA(R,1)                                                     YY 01230
         K2=IA(R,2)                                                     YY 01240
         K3=IA(R,3)                                                     YY 01250
         K4=IA(R,4)                                                     YY 01260
         CALL HIJS10(NAD,K1,K2,K3,K4,XA,XMASS,XMA,XMB,H11)
         DO 202  K = K1,K2
         DO 202  L = K1,K2
         FF=XMASS(K)*XMASS(L)/(XMA*XMA)
            DO 202  J=1,3
            DO 202  I=1,3
 202     SR(3*(K-1)+I,3*(L-1)+J)=H11(I,J)*FF
         DO 204  K = K3,K4
         DO 204  L = K3,K4
         FF=XMASS(K)*XMASS(L)/(XMB*XMB)
            DO 204  J=1,3
            DO 204  I=1,3
 204     SR(3*(K-1)+I,3*(L-1)+J)=H11(I,J)*FF
         DO 206  K = K1,K2
         DO 206  L = K3,K4
         FF=XMASS(K)*XMASS(L)/(XMA*XMB)
            DO 206  J=1,3
            DO 206  I=1,3
         SR(3*(K-1)+I,3*(L-1)+J)=-H11(I,J)*FF
 206     SR(3*(L-1)+I,3*(K-1)+J)=-H11(I,J)*FF
         CALL AHX(NC,NSX,SR,A,X)
         GO TO 300
C
 225     IF(IABS(IA(R,6)).NE.1) THEN                                    YY 01920
             WRITE(IOUT,2) R                                            YY 01930
             GO TO 300                                                  YY 01940
         END IF                                                         YY 01950
 275     CALL XIN(NC,NC,SR,-R,ISCR6)                                    YY 01960
         CALL XIN(NC,NSX,X,R,ISCR6)                                     YY 01970
         WRITE(IOUT,1) R                                                YY 01980
 300     CONTINUE                                                       YY 01990
      CALL XOUT(NC,NC,SR,-R,ISCR)                                       YY 02000
      CALL XOUT(NC,NSX,X,R,ISCR)                                        YY 02010
 500  CONTINUE                                                          YY 02020
 600  IF(NSYM.EQ.0) RETURN                                              YY 02030
      DO 650  R=1,NSYM                                                  YY 02040
         DO 610  N=1,NSYM                                               YY 02050
         DO 610  M=1,NSYM                                               YY 02060
 610     X(M,N)=ZERO                                                    YY 02070
         L=IU(R,0)                                                      YY 02080
         DO 620  I=1,L                                                  YY 02090
         RP=IU(R,I)                                                     YY 02100
         CALL XIN(NC,NSYM,SR,RP,ISCR2)                                  YY 02110
         DO 630  N=1,NSYM                                               YY 02120
         DO 630  M=1,NSYM                                               YY 02130
 630     X(M,N)=X(M,N)+U(R,I)*SR(M,N)                                   YY 02140
 620     CONTINUE                                                       YY 02150
         CALL XOUT(NC,NSYM,X,R,ISCR1)                                   YY 02160
 650  CONTINUE                                                          YY 02170
C   OPTION                                                              YY 02180
      DO 750  R=1,NSYM                                                  YY 02190
         DO 710  N=1,NC                                                 YY 02200
         DO 710  M=1,NC                                                 YY 02210
 710     SR(M,N)=ZERO                                                   YY 02220
         L=IU(R,0)                                                      YY 02230
         DO 720  I=1,L                                                  YY 02240
         RP=IU(R,I)                                                     YY 02250
         CALL XIN(NC,NC,X,-RP,ISCR2)                                    YY 02260
         DO 730  N=1,NC                                                 YY 02270
         DO 730  M=1,NC                                                 YY 02280
 730     SR(M,N)=SR(M,N)+U(R,I)*X(M,N)                                  YY 02290
 720     CONTINUE                                                       YY 02300
         CALL XOUT(NC,NC,SR,-R,ISCR1)                                   YY 02310
 750  CONTINUE                                                          YY 02320
      END                                                               YY 02330
C     ///////////////////////////////////////////////////////////////// YY 04480
      SUBROUTINE MACHY(NAD,NC,NS,NSX,IOPT,XA,TYPE,IA,A,S,U,IU,Y,SR)     YY 02840
      IMPLICIT REAL*8 (A-H,O-Z)                                         YY 02850
      INTEGER R,P,Q,NSX,RP,IOPT(30)                                     YY 02860
      CHARACTER TYPE*5                                                  YY 02870
      DIMENSION TYPE(NS),IA(NS,6),S(NS),U(NS,1),IU(NS,0:1)              YY 02880
      DIMENSION XA(NAD,3),A(NC,NC),Y(NC,NC,NC),SR(NC,NC,NC)             YY 02890
      DIMENSION H111(3,3,3),H112(3,3,3),H221(3,3,3),H222(3,3,3)         YY 02900
      DIMENSION H113(3,3,3),H123(3,3,3),H223(3,3,3),H331(3,3,3)         YY 02910
      DIMENSION H332(3,3,3),H333(3,3,3),H411(3,3,3),H421(3,3,3)         YY 02920
      DIMENSION H422(3,3,3),H431(3,3,3),H432(3,3,3),H433(3,3,3)         YY 02930
      DIMENSION H441(3,3,3),H442(3,3,3),H443(3,3,3),H444(3,3,3)         YY 02940
      DIMENSION E21(3)                                                  YY 02950
      PARAMETER(ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0,THREE=3.0D0)             YY 02960
      COMMON /IO/ IIN1,IOUT,IIN2,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25,
     $   I31,I32,I33,I35,I36,I37,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10,
     $   ISCR11,ISCR12,ISCR13,ISCR14
    1 FORMAT(/,1X,'NUMERICAL SR(I,J,K) AND Y(M,N,P) MATRICES USED FOR', YY 03010
     $ /, ' SIMPLE INTERNAL COORDINATE',I5)                             YY 03020
    2 FORMAT(/,1X,'SR(I,J,K) AND Y(M,N,P) MATRICES SET TO ZERO FOR',    YY 03030
     $ /, ' SIMPLE INTERNAL COORDINATE',I5)                             YY 03040
      NSYM=IOPT(3)                                                      YY 03050
      IF(NSYM.NE.0) THEN                                                YY 03060
         ISCR=ISCR4                                                     YY 03070
      ELSE                                                              YY 03080
         ISCR=ISCR3                                                     YY 03090
      END IF                                                            YY 03100
      DO 500  R=1,NS                                                    YY 03110
         DO 10  P=1,NC                                                  YY 03120
         DO 10  N=1,NC                                                  YY 03130
         DO 10  M=1,NC                                                  YY 03140
 10      Y(M,N,P)=ZERO                                                  YY 03150
         DO 11  I=1,NC                                                  YY 03160
         DO 11  J=1,NC                                                  YY 03170
         DO 11  K=1,NC                                                  YY 03180
 11      SR(I,J,K)=ZERO                                                 YY 03190
      IF(IA(R,6).EQ.-2) GO TO 275                                       YY 03200
      K1=IA(R,1)                                                        YY 03210
      K2=IA(R,2)                                                        YY 03220
      K3=IA(R,3)                                                        YY 03230
      K4=IA(R,4)                                                        YY 03240
      L1=3*(K1-1)                                                       YY 03250
      L2=3*(K2-1)                                                       YY 03260
      L3=3*(K3-1)                                                       YY 03270
      L4=3*(K4-1)                                                       YY 03280
C                                                                       YY 03290
      IF(TYPE(R).NE.' STRE') GO TO 25                                   YY 03300
      CALL HIJKS1(NAD,K1,K2,XA,H111)                                    YY 03310
C  OPTION                                                               YY 03320
      CALL HSRY2(NC,L1,L2,H111,SR)                                      YY 03330
C  OPTION                                                               YY 03340
      CALL AHY2(NC,NSX,L1,L2,H111,A,Y)                                  YY 03350
      GO TO 300                                                         YY 03360
 25   IF(TYPE(R).NE.' BEND') GO TO 75                                   YY 03370
      CALL HIJKS2(NAD,K1,K2,K3,XA,H111,H112,H113,H123,H221,             YY 03380
     $    H222,H223,H331,H332,H333)                                     YY 03390
C  OPTION                                                               YY 03400
      CALL HSRY3(NC,L1,L2,L3,H111,H112,H113,H123,H221,                  YY 03410
     $    H222,H223,H331,H332,H333,SR)                                  YY 03420
C  OPTION                                                               YY 03430
      CALL AHY3(NC,NSX,L1,L2,L3,H111,H112,H113,H123,H221,H222,          YY 03440
     $      H223,H331,H332,H333,A,Y)                                    YY 03450
         GO TO 300                                                      YY 03460
 75      IF(TYPE(R).NE.' LIN1') GO TO 100                               YY 03470
      CALL HIJKS3(NAD,K1,K2,K3,K4,XA,H111,H112,H113,H123,               YY 03480
     $            H221,H222,H223,H331,H332,H333)                        YY 03490
C  OPTION                                                               YY 03500
      CALL HSRY3(NC,L1,L2,L3,H111,H112,H113,H123,H221,                  YY 03510
     $    H222,H223,H331,H332,H333,SR)                                  YY 03520
C  OPTION                                                               YY 03530
      CALL AHY3(NC,NSX,L1,L2,L3,H111,H112,H113,H123,H221,H222,          YY 03540
     $      H223,H331,H332,H333,A,Y)                                    YY 03550
      GO TO 300                                                         YY 03560
 100  IF(TYPE(R).NE.'  SPF') GO TO 125                                  YY 03570
           CALL VECT1(NAD,K1,K2,E21,XA,T21)                             YY 03580
           CALL HIJKS1(NAD,K1,K2,XA,H111)                               YY 03590
           FACT1=(ONE-S(R))*THREE/T21                                   YY 03600
           FACT2=TWO/(T21*T21)                                          YY 03610
           DO 102  K=1,3                                                YY 03620
           DO 102  J=1,3                                                YY 03630
           DO 102  I=1,3                                                YY 03640
 102       H111(I,J,K)=FACT1*(H111(I,J,K)+FACT2*E21(I)*E21(J)*E21(K))   YY 03650
C  OPTION                                                               YY 03660
      CALL HSRY2(NC,L1,L2,H111,SR)                                      YY 03670
C  OPTION                                                               YY 03680
      CALL AHY2(NC,NSX,L1,L2,H111,A,Y)                                  YY 03690
      GO TO 300                                                         YY 03700
 125  IF(TYPE(R).NE.' TORS') GO TO 150                                  YY 03710
      CALL HIJKS6(NAD,K1,K2,K3,K4,XA,H111,H112,H221,H222,H113,H123,     YY 03720
     $      H223,H331,H332,H333,H411,H421,H422,H431,H432,H433,H441,     YY 03730
     $      H442,H443,H444)                                             YY 03740
C  OPTION                                                               YY 03750
      CALL HSRY4(NC,L1,L2,L3,L4,H111,H112,H221,H222,H113,H123,          YY 03760
     $     H223,H331,H332,H333,H411,H421,H422,H431,H432,H433,H441,      YY 03770
     $     H442,H443,H444,SR)                                           YY 03780
C  OPTION                                                               YY 03790
      CALL AHY4(NC,NSX,L1,L2,L3,L4,H111,H112,H221,H222,H113,H123,       YY 03800
     $     H223,H331,H332,H333,H411,H421,H422,H431,H432,H433,H441,      YY 03810
     $     H442,H443,H444,A,Y)                                          YY 03820
      GO TO 300                                                         YY 03830
 150  IF(TYPE(R).NE.'  OUT') GO TO 175                                  YY 03840
      CALL HIJKS7(NAD,K1,K2,K3,K4,XA,H111,H112,H221,H222,H113,H123,     YY 03850
     $      H223,H331,H332,H333,H411,H421,H422,H431,H432,H433,H441,     YY 03860
     $      H442,H443,H444)                                             YY 03870
C  OPTION                                                               YY 03880
      CALL HSRY4(NC,L1,L2,L3,L4,H111,H112,H221,H222,H113,H123,          YY 03890
     $     H223,H331,H332,H333,H411,H421,H422,H431,H432,H433,H441,      YY 03900
     $     H442,H443,H444,SR)                                           YY 03910
C  OPTION                                                               YY 03920
      CALL AHY4(NC,NSX,L1,L2,L3,L4,H111,H112,H221,H222,H113,H123,       YY 03930
     $     H223,H331,H332,H333,H411,H421,H422,H431,H432,H433,H441,      YY 03940
     $     H442,H443,H444,A,Y)                                          YY 03950
      GO TO 300                                                         YY 03960
 175  IF(TYPE(R).NE.' LINX'.AND.TYPE(R).NE.' LINY') GO TO 200           YY 03710
      IF(TYPE(R).EQ.' LINX') CALL HIJKS8(NAD,K1,K2,K3,K4,XA,H111,
     $      H112,H221,H222,H113,H123,H223,H331,H332,H333,H411,
     $      H421,H422,H431,H432,H433,H441,H442,H443,H444)
      IF(TYPE(R).EQ.' LINY') CALL HIJKS9(NAD,K1,K2,K3,K4,XA,H111,
     $      H112,H221,H222,H113,H123,H223,H331,H332,H333,H411,
     $      H421,H422,H431,H432,H433,H441,H442,H443,H444)
C  OPTION                                                               YY 03750
      CALL HSRY4(NC,L1,L2,L3,L4,H111,H112,H221,H222,H113,H123,          YY 03760
     $     H223,H331,H332,H333,H411,H421,H422,H431,H432,H433,H441,      YY 03770
     $     H442,H443,H444,SR)                                           YY 03780
C  OPTION                                                               YY 03790
      CALL AHY4(NC,NSX,L1,L2,L3,L4,H111,H112,H221,H222,H113,H123,       YY 03800
     $     H223,H331,H332,H333,H411,H421,H422,H431,H432,H433,H441,      YY 03810
     $     H442,H443,H444,A,Y)                                          YY 03820
      GO TO 300                                                         YY 03830
 200  IF(IABS(IA(R,6)).NE.2) THEN                                       YY 03970
         WRITE(IOUT,2) R                                                YY 03980
         GO TO 300                                                      YY 03990
      END IF                                                            YY 04000
 275  CALL YIN(NC,NC,SR,-R,ISCR7)                                       YY 04010
      CALL YIN(NC,NSX,Y,R,ISCR7)                                        YY 04020
      WRITE(IOUT,1) R                                                   YY 04030
 300  CONTINUE                                                          YY 04040
      CALL FILL3A(NC,NSX,Y)                                             YY 04050
      CALL YOUT(NC,NC,SR,-R,ISCR)                                       YY 04060
      CALL YOUT(NC,NSX,Y,R,ISCR)                                        YY 04070
 500  CONTINUE                                                          YY 04080
 600  IF(NSYM.EQ.0) RETURN                                              YY 04090
      DO 650  R=1,NSYM                                                  YY 04100
         DO 610  P=1,NSYM                                               YY 04110
         DO 610  N=1,P                                                  YY 04120
         DO 610  M=1,N                                                  YY 04130
 610     Y(M,N,P)=ZERO                                                  YY 04140
         L=IU(R,0)                                                      YY 04150
         DO 620  I=1,L                                                  YY 04160
         RP=IU(R,I)                                                     YY 04170
         CALL YIN(NC,NSYM,SR,RP,ISCR4)                                  YY 04180
         W1=U(R,I)                                                      YY 04190
         DO 630  P=1,NSYM                                               YY 04200
         DO 630  N=1,P                                                  YY 04210
         DO 630  M=1,N                                                  YY 04220
 630     Y(M,N,P)=Y(M,N,P)+W1*SR(M,N,P)                                 YY 04230
 620     CONTINUE                                                       YY 04240
         CALL FILL3A(NC,NSYM,Y)                                         YY 04250
         CALL YOUT(NC,NSYM,Y,R,ISCR3)                                   YY 04260
 650  CONTINUE                                                          YY 04270
C   OPTION                                                              YY 04280
      DO 750  R=1,NSYM                                                  YY 04290
         DO 710  P=1,NC                                                 YY 04300
         DO 710  N=1,P                                                  YY 04310
         DO 710  M=1,N                                                  YY 04320
 710     SR(M,N,P)=ZERO                                                 YY 04330
         L=IU(R,0)                                                      YY 04340
         DO 720  I=1,L                                                  YY 04350
         RP=IU(R,I)                                                     YY 04360
         CALL YIN(NC,NC,Y,-RP,ISCR4)                                    YY 04370
         W1=U(R,I)                                                      YY 04380
         DO 730  P=1,NC                                                 YY 04390
         DO 730  N=1,P                                                  YY 04400
         DO 730  M=1,N                                                  YY 04410
 730     SR(M,N,P)=SR(M,N,P)+W1*Y(M,N,P)                                YY 04420
 720     CONTINUE                                                       YY 04430
         CALL FILL3A(NC,NC,SR)                                          YY 04440
         CALL YOUT(NC,NC,SR,-R,ISCR3)                                   YY 04450
 750  CONTINUE                                                          YY 04460
      END                                                               YY 04470
C     //////////////////////////////////////////////////////////////    INT39050
      SUBROUTINE BROW(NAD,NC,NS,XA,XMASS,TYPE,IA,RB,S,IR)               INT39060
      IMPLICIT REAL*8 (A-H,O-Z)                                         INT39070
      INTEGER IL(5)                                                     INT39080
      CHARACTER TYPE*5                                                  INT39090
      DIMENSION TYPE(NS),IA(NS,6),RB(NC),XA(NAD,3),XMASS(1),S(NS)       INT39100
      DIMENSION V(3,5),Q(3,5)                                           INT39110
      PARAMETER(ZERO=0.0D0,ONE=1.0D0)                                   INT39120
         DO 5  J=1,NC                                                   INT39130
 5       RB(J)=ZERO                                                     INT39140
           NIA=0                                                        INT39150
           K1=IA(IR,1)                                                  INT39160
           K2=IA(IR,2)                                                  INT39170
           K3=IA(IR,3)                                                  INT39180
           K4=IA(IR,4)                                                  INT39190
           K5=IA(IR,5)                                                  INT39200
           IL(1)=3*(K1-1)                                               INT39210
           IL(2)=3*(K2-1)                                               INT39220
           IL(3)=3*(K3-1)                                               INT39230
           IL(4)=3*(K4-1)                                               INT39240
           IL(5)=3*(K5-1)                                               INT39250
      IF(TYPE(IR).NE.' STRE') GO TO 10                                  INT39260
         CALL VECT1(NAD,K1,K2,V(1,1),XA,W)                              INT39270
         DO 6  K=1,3                                                    INT39280
  6      V(K,2)=-V(K,1)                                                 INT39290
         NIA=2                                                          INT39300
         GO TO 80                                                       INT39310
 10   IF(TYPE(IR).NE.' BEND') GO TO 15                                  INT39320
         CALL VECT2(NAD,K1,K2,K3,V(1,1),V(1,2),V(1,3),XA,W)             INT39330
         NIA=3                                                          INT39340
         GO TO 80                                                       INT39350
 15   IF(TYPE(IR).NE.' LIN1') GO TO 25                                  INT39360
           CALL VECT3(NAD,K1,K2,K3,K4,V(1,1),V(1,2),V(1,3),XA,W)        INT39370
           NIA=3                                                        INT39380
           GO TO 80                                                     INT39390
 25   IF(TYPE(IR).NE.'  OUT') GO TO 30                                  INT39400
           CALL VECT5(NAD,K1,K2,K3,K4,V(1,1),V(1,2),V(1,3),V(1,4),XA,W) INT39410
           NIA=4                                                        INT39420
           GO TO 80                                                     INT39430
 30   IF(TYPE(IR).NE.' TORS') GO TO 35                                  INT39440
           CALL VECT6(NAD,K1,K2,K3,K4,V(1,1),V(1,2),V(1,3),V(1,4),XA,W) INT39450
           NIA=4                                                        INT39460
           GO TO 80                                                     INT39470
 35   IF(TYPE(IR).NE.'  SPF') GO TO 40                                  INT39480
           CALL VECT1(NAD,K1,K2,V(1,1),XA,W)                            INT39490
           DO 36  K=1,3                                                 INT39500
 36        V(K,2)=-V(K,1)                                               INT39510
           NIA=2                                                        INT39520
           FACT=S(IR)/(W*W)                                             INT39530
           DO 38  L=1,2                                                 INT39540
           DO 38  K=1,3                                                 INT39550
 38        V(K,L)=V(K,L)*FACT                                           INT39560
           GO TO 80                                                     INT39570
 40   IF(TYPE(IR).NE.' LINX') GO TO 45
         CALL VECT8(NAD,K1,K2,K3,K4,V(1,1),V(1,2),V(1,3),V(1,4),XA,W)   INT05710
         NIA=4                                                          INT05720
         GO TO 80                                                       INT05730
 45   IF(TYPE(IR).NE.' LINY') GO TO 50
         CALL VECT9(NAD,K1,K2,K3,K4,V(1,1),V(1,2),V(1,3),V(1,4),XA,W)   INT05710
         NIA=4                                                          INT05720
         GO TO 80
 50   IF(TYPE(IR).NE.' RCOM') GO TO 55
         CALL VECT10(NAD,K1,K2,K3,K4,V(1,1),XA,XMASS,XMA,XMB,W)
         DO 51  K=K1,K2
           L=3*(K-1)
           DO 51  J=1,3
 51        RB(L+J)=XMASS(K)*V(J,1)/XMA
         DO 52  K=K3,K4
           L=3*(K-1)
           DO 52  J=1,3
 52        RB(L+J)=-XMASS(K)*V(J,1)/XMB
         GO TO 90
 55   CONTINUE
 80   DO 85  K=1,NIA                                                    INT39590
      DO 85  J=1,3                                                      INT39600
 85   RB(IL(K)+J)=V(J,K)                                                INT39610
 90   CONTINUE
      RETURN                                                            INT39620
      END                                                               INT39630
C     ///////////////////////////////////////////////////////////////   INT39640
      SUBROUTINE XROW(NAD,NC,NS,XA,TYPE,IA,S,SR,R)                      YY 03150
      IMPLICIT REAL*8 (A-H,O-Z)                                         YY 03160
      CHARACTER TYPE*5                                                  YY 03170
      INTEGER R,RP                                                      YY 03180
      DIMENSION TYPE(NS),IA(NS,6),S(NS),XA(NAD,3),SR(NC,NC)             YY 03190
      DIMENSION H11(3,3),H21(3,3),H22(3,3),H31(3,3),H32(3,3)            YY 03200
      DIMENSION H33(3,3),H41(3,3),H42(3,3),H43(3,3),H44(3,3)            YY 03210
      DIMENSION H11A(3,3),H33A(3,3),EA(3),E21(3)                        YY 03220
      PARAMETER(ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0,THREE=3.0D0)             YY 03230
      COMMON /IO/ IIN1,IOUT,IIN2,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25,
     $   I31,I32,I33,I35,I36,I37,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10,
     $   ISCR11,ISCR12,ISCR13,ISCR14
         DO 49  J=1,NC                                                  YY 03280
         DO 49  I=1,NC                                                  YY 03290
 49      SR(I,J)=ZERO                                                   YY 03300
         IF(TYPE(R).NE.' STRE') GO TO 25                                YY 03310
         K1=IA(R,1)                                                     YY 03320
         K2=IA(R,2)                                                     YY 03330
         L1=3*(K1-1)                                                    YY 03340
         L2=3*(K2-1)                                                    YY 03350
         CALL HIJS1(NAD,K1,K2,XA,H11)                                   YY 03360
C   OPTION                                                              YY 03370
              DO 28  J=1,3                                              YY 03380
              DO 28  I=1,3                                              YY 03390
              SR(L1+I,L1+J)=H11(I,J)                                    YY 03400
              SR(L2+I,L2+J)=H11(I,J)                                    YY 03410
              SR(L1+I,L2+J)=-H11(I,J)                                   YY 03420
 28           SR(L2+I,L1+J)=-H11(I,J)                                   YY 03430
C  *************                                                        YY 03440
         GO TO 300                                                      YY 03450
 25      IF(TYPE(R).NE.' BEND') GO TO 75                                YY 03460
         K1=IA(R,1)                                                     YY 03470
         K2=IA(R,2)                                                     YY 03480
         K3=IA(R,3)                                                     YY 03490
         L1=3*(K1-1)                                                    YY 03500
         L2=3*(K2-1)                                                    YY 03510
         L3=3*(K3-1)                                                    YY 03520
         CALL HIJS2(NAD,K1,K2,K3,XA,H11,H21,H31,H22,H32,H33)            YY 03530
C  OPTION                                                               YY 03540
      DO 48  J=1,3                                                      YY 03550
      DO 48  I=1,3                                                      YY 03560
      SR(L1+I,L1+J)=H11(I,J)                                            YY 03570
      SR(L2+I,L1+J)=H21(I,J)                                            YY 03580
      SR(L3+I,L1+J)=H31(I,J)                                            YY 03590
      SR(L1+I,L2+J)=H21(J,I)                                            YY 03600
      SR(L2+I,L2+J)=H22(I,J)                                            YY 03610
      SR(L3+I,L2+J)=H32(I,J)                                            YY 03620
      SR(L1+I,L3+J)=H31(J,I)                                            YY 03630
      SR(L2+I,L3+J)=H32(J,I)                                            YY 03640
 48   SR(L3+I,L3+J)=H33(I,J)                                            YY 03650
         GO TO 300                                                      YY 03660
 75      IF(TYPE(R).NE.' LIN1') GO TO 100                               YY 03670
         K1=IA(R,1)                                                     YY 03680
         K2=IA(R,2)                                                     YY 03690
         K3=IA(R,3)                                                     YY 03700
         K4=IA(R,4)                                                     YY 03710
         L1=3*(K1-1)                                                    YY 03720
         L2=3*(K2-1)                                                    YY 03730
         L3=3*(K3-1)                                                    YY 03740
         CALL HIJS3(NAD,K1,K2,K3,K4,XA,H11,H21,H31,H22,H32,H33)         YY 03750
C  OPTION                                                               YY 03760
      DO 86  J=1,3                                                      YY 03770
      DO 86  I=1,3                                                      YY 03780
      SR(L1+I,L1+J)=H11(I,J)                                            YY 03790
      SR(L2+I,L1+J)=H21(I,J)                                            YY 03800
      SR(L3+I,L1+J)=H31(I,J)                                            YY 03810
      SR(L1+I,L2+J)=H21(J,I)                                            YY 03820
      SR(L2+I,L2+J)=H22(I,J)                                            YY 03830
      SR(L3+I,L2+J)=H32(I,J)                                            YY 03840
      SR(L1+I,L3+J)=H31(J,I)                                            YY 03850
      SR(L2+I,L3+J)=H32(J,I)                                            YY 03860
 86   SR(L3+I,L3+J)=H33(I,J)                                            YY 03870
         GO TO 300                                                      YY 03880
 100  IF(TYPE(R).NE.'  SPF') GO TO 125                                  YY 03890
         K1=IA(R,1)                                                     YY 03900
         K2=IA(R,2)                                                     YY 03910
         L1=3*(K1-1)                                                    YY 03920
         L2=3*(K2-1)                                                    YY 03930
         CALL VECT1(NAD,K1,K2,E21,XA,T21)                               YY 03940
         CALL HIJS1(NAD,K1,K2,XA,H11)                                   YY 03950
         FACT1=S(R)/(T21*T21)                                           YY 03960
         FACT2=TWO*FACT1/T21                                            YY 03970
         DO 102  J=1,3                                                  YY 03980
         DO 104  I=1,3                                                  YY 03990
 104     H11(I,J)=THREE*H11(I,J)*FACT1                                  YY 04000
 102     H11(J,J)=H11(J,J)-FACT2                                        YY 04010
         DO 106  J=1,3                                                  YY 04020
         DO 106  I=1,3                                                  YY 04030
         SR(L1+I,L1+J)=H11(I,J)                                         YY 04040
         SR(L2+I,L2+J)=H11(I,J)                                         YY 04050
         SR(L1+I,L2+J)=-H11(I,J)                                        YY 04060
 106     SR(L2+I,L1+J)=-H11(I,J)                                        YY 04070
         GO TO 300                                                      YY 04080
 125     IF(TYPE(R).NE.' TORS') GO TO 150                               YY 04090
         K1=IA(R,1)                                                     YY 04100
         K2=IA(R,2)                                                     YY 04110
         K3=IA(R,3)                                                     YY 04120
         K4=IA(R,4)                                                     YY 04130
         L1=3*(K1-1)                                                    YY 04140
         L2=3*(K2-1)                                                    YY 04150
         L3=3*(K3-1)                                                    YY 04160
         L4=3*(K4-1)                                                    YY 04170
         CALL HIJS6(NAD,K1,K2,K3,K4,XA,H11,H21,H31,H41,H22,H32,H42,     YY 04180
     $               H33,H43,H44)                                       YY 04190
C  OPTION                                                               YY 04200
      DO 148  J=1,3                                                     YY 04210
      DO 148  I=1,3                                                     YY 04220
      SR(L1+I,L1+J)=H11(I,J)                                            YY 04230
      SR(L2+I,L1+J)=H21(I,J)                                            YY 04240
      SR(L3+I,L1+J)=H31(I,J)                                            YY 04250
      SR(L4+I,L1+J)=H41(I,J)                                            YY 04260
      SR(L1+I,L2+J)=H21(J,I)                                            YY 04270
      SR(L2+I,L2+J)=H22(I,J)                                            YY 04280
      SR(L3+I,L2+J)=H32(I,J)                                            YY 04290
      SR(L4+I,L2+J)=H42(I,J)                                            YY 04300
      SR(L1+I,L3+J)=H31(J,I)                                            YY 04310
      SR(L2+I,L3+J)=H32(J,I)                                            YY 04320
      SR(L3+I,L3+J)=H33(I,J)                                            YY 04330
      SR(L4+I,L3+J)=H43(I,J)                                            YY 04340
      SR(L1+I,L4+J)=H41(J,I)                                            YY 04350
      SR(L2+I,L4+J)=H42(J,I)                                            YY 04360
      SR(L3+I,L4+J)=H43(J,I)                                            YY 04370
 148  SR(L4+I,L4+J)=H44(I,J)                                            YY 04380
 150     IF(TYPE(R).NE.'  OUT') GO TO 175                               YY 04390
         K1=IA(R,1)                                                     YY 04400
         K2=IA(R,2)                                                     YY 04410
         K3=IA(R,3)                                                     YY 04420
         K4=IA(R,4)                                                     YY 04430
         L1=3*(K1-1)                                                    YY 04440
         L2=3*(K2-1)                                                    YY 04450
         L3=3*(K3-1)                                                    YY 04460
         L4=3*(K4-1)                                                    YY 04470
         CALL HIJS7(NAD,K1,K2,K3,K4,XA,H11,H21,H31,H41,H22,H32,H42,     YY 04480
     $               H33,H43,H44)                                       YY 04490
C  OPTION                                                               YY 04500
      DO 158  J=1,3                                                     YY 04510
      DO 158  I=1,3                                                     YY 04520
      SR(L1+I,L1+J)=H11(I,J)                                            YY 04530
      SR(L2+I,L1+J)=H21(I,J)                                            YY 04540
      SR(L3+I,L1+J)=H31(I,J)                                            YY 04550
      SR(L4+I,L1+J)=H41(I,J)                                            YY 04560
      SR(L1+I,L2+J)=H21(J,I)                                            YY 04570
      SR(L2+I,L2+J)=H22(I,J)                                            YY 04580
      SR(L3+I,L2+J)=H32(I,J)                                            YY 04590
      SR(L4+I,L2+J)=H42(I,J)                                            YY 04600
      SR(L1+I,L3+J)=H31(J,I)                                            YY 04610
      SR(L2+I,L3+J)=H32(J,I)                                            YY 04620
      SR(L3+I,L3+J)=H33(I,J)                                            YY 04630
      SR(L4+I,L3+J)=H43(I,J)                                            YY 04640
      SR(L1+I,L4+J)=H41(J,I)                                            YY 04650
      SR(L2+I,L4+J)=H42(J,I)                                            YY 04660
      SR(L3+I,L4+J)=H43(J,I)                                            YY 04670
 158  SR(L4+I,L4+J)=H44(I,J)                                            YY 04680
         GO TO 300                                                      YY 04690
 175     IF(TYPE(R).NE.' LINX'.AND.TYPE(R).NE.' LINY') GO TO 200        YY 01220
         K1=IA(R,1)                                                     YY 01230
         K2=IA(R,2)                                                     YY 01240
         K3=IA(R,3)                                                     YY 01250
         K4=IA(R,4)                                                     YY 01260
         L1=3*(K1-1)                                                    YY 01270
         L2=3*(K2-1)                                                    YY 01280
         L3=3*(K3-1)                                                    YY 01290
         L4=3*(K4-1)                                                    YY 01300
         IF(TYPE(R).EQ.' LINX') CALL HIJS8(NAD,K1,K2,K3,K4,XA,H11,
     $       H21,H31,H41,H22,H32,H42,H33,H43,H44)
         IF(TYPE(R).EQ.' LINY') CALL HIJS9(NAD,K1,K2,K3,K4,XA,H11,
     $       H21,H31,H41,H22,H32,H42,H33,H43,H44)
C  OPTION
      DO 178  J=1,3                                                     YY 01340
      DO 178  I=1,3                                                     YY 01350
      SR(L1+I,L1+J)=H11(I,J)                                            YY 01360
      SR(L2+I,L1+J)=H21(I,J)                                            YY 01370
      SR(L3+I,L1+J)=H31(I,J)                                            YY 01380
      SR(L4+I,L1+J)=H41(I,J)                                            YY 01390
      SR(L1+I,L2+J)=H21(J,I)                                            YY 01400
      SR(L2+I,L2+J)=H22(I,J)                                            YY 01410
      SR(L3+I,L2+J)=H32(I,J)                                            YY 01420
      SR(L4+I,L2+J)=H42(I,J)                                            YY 01430
      SR(L1+I,L3+J)=H31(J,I)                                            YY 01440
      SR(L2+I,L3+J)=H32(J,I)                                            YY 01450
      SR(L3+I,L3+J)=H33(I,J)                                            YY 01460
      SR(L4+I,L3+J)=H43(I,J)                                            YY 01470
      SR(L1+I,L4+J)=H41(J,I)                                            YY 01480
      SR(L2+I,L4+J)=H42(J,I)                                            YY 01490
      SR(L3+I,L4+J)=H43(J,I)                                            YY 01500
 178  SR(L4+I,L4+J)=H44(I,J)                                            YY 01510
C     ***********************************                               YY 01520
         CALL AHX4(NC,NSX,L1,L2,L3,L4,H11,H21,H31,H41,H22,H32,H42,      YY 01530
     $             H33,H43,H44,A,X)                                     YY 01540
         GO TO 300                                                      YY 01550
 200  CONTINUE                                                          YY 04700
 300  CONTINUE                                                          YY 04710
      END                                                               YY 04720
C     /////////////////////////////////////////////////////////////////
      SUBROUTINE YROW(NAD,NC,NS,XA,TYPE,IA,S,SR,R)                      YY 07470
      IMPLICIT REAL*8 (A-H,O-Z)                                         YY 07480
      INTEGER R,P,Q                                                     YY 07490
      CHARACTER TYPE*5                                                  YY 07500
      DIMENSION TYPE(NS),S(NS),IA(NS,6)                                 YY 07510
      DIMENSION XA(NAD,3),SR(NC,NC,NC)                                  YY 07520
      DIMENSION H111(3,3,3),H112(3,3,3),H221(3,3,3),H222(3,3,3)         YY 07530
      DIMENSION H113(3,3,3),H123(3,3,3),H223(3,3,3),H331(3,3,3)         YY 07540
      DIMENSION H332(3,3,3),H333(3,3,3),H411(3,3,3),H421(3,3,3)         YY 07550
      DIMENSION H422(3,3,3),H431(3,3,3),H432(3,3,3),H433(3,3,3)         YY 07560
      DIMENSION H441(3,3,3),H442(3,3,3),H443(3,3,3),H444(3,3,3)         YY 07570
      DIMENSION E21(3)                                                  YY 07580
      PARAMETER(ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0,THREE=3.0D0)             YY 07590
      COMMON /IO/ IIN1,IOUT,IIN2,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25,
     $   I31,I32,I33,I35,I36,I37,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10,
     $   ISCR11,ISCR12,ISCR13,ISCR14
         DO 10  I=1,NC                                                  YY 07640
         DO 10  J=1,NC                                                  YY 07650
         DO 10  K=1,NC                                                  YY 07660
 10      SR(I,J,K)=ZERO                                                 YY 07670
      IF(TYPE(R).NE.' STRE') GO TO 25                                   YY 07680
      K1=IA(R,1)                                                        YY 07690
      K2=IA(R,2)                                                        YY 07700
      L1=3*(K1-1)                                                       YY 07710
      L2=3*(K2-1)                                                       YY 07720
      CALL HIJKS1(NAD,K1,K2,XA,H111)                                    YY 07730
C  OPTION                                                               YY 07740
      DO 17  K=1,3                                                      YY 07750
      DO 17  J=1,3                                                      YY 07760
      DO 17  I=1,3                                                      YY 07770
      Z=H111(I,J,K)                                                     YY 07780
      SR(L1+I,L1+J,L1+K)=Z                                              YY 07790
      SR(L1+I,L1+J,L2+K)=-Z                                             YY 07800
      SR(L1+I,L2+J,L1+K)=-Z                                             YY 07810
      SR(L1+I,L2+J,L2+K)=Z                                              YY 07820
      SR(L2+I,L1+J,L1+K)=-Z                                             YY 07830
      SR(L2+I,L1+J,L2+K)=Z                                              YY 07840
      SR(L2+I,L2+J,L1+K)=Z                                              YY 07850
 17   SR(L2+I,L2+J,L2+K)=-Z                                             YY 07860
C  OPTION                                                               YY 07870
      GO TO 300                                                         YY 07880
 25   IF(TYPE(R).NE.' BEND') GO TO 75                                   YY 07890
      K1=IA(R,1)                                                        YY 07900
      K2=IA(R,2)                                                        YY 07910
      K3=IA(R,3)                                                        YY 07920
      L1=3*(K1-1)                                                       YY 07930
      L2=3*(K2-1)                                                       YY 07940
      L3=3*(K3-1)                                                       YY 07950
      CALL HIJKS2(NAD,K1,K2,K3,XA,H111,H112,H113,H123,H221,             YY 07960
     $    H222,H223,H331,H332,H333)                                     YY 07970
C  OPTION                                                               YY 07980
      DO 37  K=1,3                                                      YY 07990
      DO 37  J=1,3                                                      YY 08000
      DO 37  I=1,3                                                      YY 08010
      SR(L1+I,L1+J,L1+K)=H111(I,J,K)                                    YY 08020
      SR(L1+I,L1+J,L2+K)=H112(I,J,K)                                    YY 08030
      SR(L1+I,L1+J,L3+K)=H113(I,J,K)                                    YY 08040
      SR(L1+I,L2+J,L1+K)=H112(I,K,J)                                    YY 08050
      SR(L1+I,L2+J,L2+K)=H221(J,K,I)                                    YY 08060
      SR(L1+I,L2+J,L3+K)=H123(I,J,K)                                    YY 08070
      SR(L2+I,L2+J,L1+K)=H221(I,J,K)                                    YY 08080
      SR(L2+I,L2+J,L2+K)=H222(I,J,K)                                    YY 08090
      SR(L2+I,L2+J,L3+K)=H223(I,J,K)                                    YY 08100
      SR(L2+I,L1+J,L1+K)=H112(J,K,I)                                    YY 08110
      SR(L2+I,L1+J,L2+K)=H221(I,K,J)                                    YY 08120
      SR(L2+I,L1+J,L3+K)=H123(J,I,K)                                    YY 08130
      SR(L1+I,L3+J,L1+K)=H113(I,K,J)                                    YY 08140
      SR(L1+I,L3+J,L2+K)=H123(I,K,J)                                    YY 08150
      SR(L1+I,L3+J,L3+K)=H331(J,K,I)                                    YY 08160
      SR(L2+I,L3+J,L1+K)=H123(K,I,J)                                    YY 08170
      SR(L2+I,L3+J,L2+K)=H223(I,K,J)                                    YY 08180
      SR(L2+I,L3+J,L3+K)=H332(J,K,I)                                    YY 08190
      SR(L3+I,L1+J,L1+K)=H113(J,K,I)                                    YY 08200
      SR(L3+I,L1+J,L2+K)=H123(J,K,I)                                    YY 08210
      SR(L3+I,L1+J,L3+K)=H331(I,K,J)                                    YY 08220
      SR(L3+I,L2+J,L1+K)=H123(K,J,I)                                    YY 08230
      SR(L3+I,L2+J,L2+K)=H223(J,K,I)                                    YY 08240
      SR(L3+I,L2+J,L3+K)=H332(I,K,J)                                    YY 08250
      SR(L3+I,L3+J,L1+K)=H331(I,J,K)                                    YY 08260
      SR(L3+I,L3+J,L2+K)=H332(I,J,K)                                    YY 08270
 37   SR(L3+I,L3+J,L3+K)=H333(I,J,K)                                    YY 08280
C  OPTION                                                               YY 08290
         GO TO 300                                                      YY 08300
 75      IF(TYPE(R).NE.' LIN1') GO TO 100                               YY 08310
         K1=IA(R,1)                                                     YY 08320
         K2=IA(R,2)                                                     YY 08330
         K3=IA(R,3)                                                     YY 08340
         K4=IA(R,4)                                                     YY 08350
         L1=3*(K1-1)                                                    YY 08360
         L2=3*(K2-1)                                                    YY 08370
         L3=3*(K3-1)                                                    YY 08380
      CALL HIJKS3(NAD,K1,K2,K3,K4,XA,H111,H112,H113,H123,               YY 08390
     $            H221,H222,H223,H331,H332,H333)                        YY 08400
C  OPTION                                                               YY 08410
      DO 88  K=1,3                                                      YY 08420
      DO 88  J=1,3                                                      YY 08430
      DO 88  I=1,3                                                      YY 08440
      SR(L1+I,L1+J,L1+K)=H111(I,J,K)                                    YY 08450
      SR(L1+I,L1+J,L2+K)=H112(I,J,K)                                    YY 08460
      SR(L1+I,L1+J,L3+K)=H113(I,J,K)                                    YY 08470
      SR(L1+I,L2+J,L1+K)=H112(I,K,J)                                    YY 08480
      SR(L1+I,L2+J,L2+K)=H221(J,K,I)                                    YY 08490
      SR(L1+I,L2+J,L3+K)=H123(I,J,K)                                    YY 08500
      SR(L2+I,L2+J,L1+K)=H221(I,J,K)                                    YY 08510
      SR(L2+I,L2+J,L2+K)=H222(I,J,K)                                    YY 08520
      SR(L2+I,L2+J,L3+K)=H223(I,J,K)                                    YY 08530
      SR(L2+I,L1+J,L1+K)=H112(J,K,I)                                    YY 08540
      SR(L2+I,L1+J,L2+K)=H221(I,K,J)                                    YY 08550
      SR(L2+I,L1+J,L3+K)=H123(J,I,K)                                    YY 08560
      SR(L1+I,L3+J,L1+K)=H113(I,K,J)                                    YY 08570
      SR(L1+I,L3+J,L2+K)=H123(I,K,J)                                    YY 08580
      SR(L1+I,L3+J,L3+K)=H331(J,K,I)                                    YY 08590
      SR(L2+I,L3+J,L1+K)=H123(K,I,J)                                    YY 08600
      SR(L2+I,L3+J,L2+K)=H223(I,K,J)                                    YY 08610
      SR(L2+I,L3+J,L3+K)=H332(J,K,I)                                    YY 08620
      SR(L3+I,L1+J,L1+K)=H113(J,K,I)                                    YY 08630
      SR(L3+I,L1+J,L2+K)=H123(J,K,I)                                    YY 08640
      SR(L3+I,L1+J,L3+K)=H331(I,K,J)                                    YY 08650
      SR(L3+I,L2+J,L1+K)=H123(K,J,I)                                    YY 08660
      SR(L3+I,L2+J,L2+K)=H223(J,K,I)                                    YY 08670
      SR(L3+I,L2+J,L3+K)=H332(I,K,J)                                    YY 08680
      SR(L3+I,L3+J,L1+K)=H331(I,J,K)                                    YY 08690
      SR(L3+I,L3+J,L2+K)=H332(I,J,K)                                    YY 08700
 88   SR(L3+I,L3+J,L3+K)=H333(I,J,K)                                    YY 08710
C  OPTION                                                               YY 08720
      GO TO 300                                                         YY 08730
 100  IF(TYPE(R).NE.'  SPF') GO TO 125                                  YY 08740
           K1=IA(R,1)                                                   YY 08750
           K2=IA(R,2)                                                   YY 08760
           L1=3*(K1-1)                                                  YY 08770
           L2=3*(K2-1)                                                  YY 08780
           CALL VECT1(NAD,K1,K2,E21,XA,T21)                             YY 08790
           CALL HIJKS1(NAD,K1,K2,XA,H111)                               YY 08800
           FACT1=THREE*S(R)/(T21*T21)                                   YY 08810
           FACT2=TWO/(T21*T21)                                          YY 08820
           DO 102  K=1,3                                                YY 08830
           DO 102  J=1,3                                                YY 08840
           DO 102  I=1,3                                                YY 08850
 102       H111(I,J,K)=FACT1*(H111(I,J,K)+FACT2*E21(I)*E21(J)*E21(K))   YY 08860
C  OPTION                                                               YY 08870
      DO 107  K=1,3                                                     YY 08880
      DO 107  J=1,3                                                     YY 08890
      DO 107  I=1,3                                                     YY 08900
      Z=H111(I,J,K)                                                     YY 08910
      SR(L1+I,L1+J,L1+K)=Z                                              YY 08920
      SR(L1+I,L1+J,L2+K)=-Z                                             YY 08930
      SR(L1+I,L2+J,L1+K)=-Z                                             YY 08940
      SR(L1+I,L2+J,L2+K)=Z                                              YY 08950
      SR(L2+I,L1+J,L1+K)=-Z                                             YY 08960
      SR(L2+I,L1+J,L2+K)=Z                                              YY 08970
      SR(L2+I,L2+J,L1+K)=Z                                              YY 08980
 107  SR(L2+I,L2+J,L2+K)=-Z                                             YY 08990
      GO TO 300                                                         YY 09000
 125  IF(TYPE(R).NE.' TORS') GO TO 150                                  YY 09010
      K1=IA(R,1)                                                        YY 09020
      K2=IA(R,2)                                                        YY 09030
      K3=IA(R,3)                                                        YY 09040
      K4=IA(R,4)                                                        YY 09050
      L1=3*(K1-1)                                                       YY 09060
      L2=3*(K2-1)                                                       YY 09070
      L3=3*(K3-1)                                                       YY 09080
      L4=3*(K4-1)                                                       YY 09090
      CALL HIJKS6(NAD,K1,K2,K3,K4,XA,H111,H112,H221,H222,H113,H123,     YY 09100
     $      H223,H331,H332,H333,H411,H421,H422,H431,H432,H433,H441,     YY 09110
     $      H442,H443,H444)                                             YY 09120
C  OPTION                                                               YY 09130
      DO 137  K=1,3                                                     YY 09140
      DO 137  J=1,3                                                     YY 09150
      DO 137  I=1,3                                                     YY 09160
      SR(L1+I,L1+J,L1+K)=H111(I,J,K)                                    YY 09170
      SR(L1+I,L1+J,L2+K)=H112(I,J,K)                                    YY 09180
      SR(L1+I,L1+J,L3+K)=H113(I,J,K)                                    YY 09190
      SR(L1+I,L1+J,L4+K)=H411(K,J,I)                                    YY 09200
      SR(L1+I,L2+J,L1+K)=H112(I,K,J)                                    YY 09210
      SR(L1+I,L2+J,L2+K)=H221(J,K,I)                                    YY 09220
      SR(L1+I,L2+J,L3+K)=H123(I,J,K)                                    YY 09230
      SR(L1+I,L2+J,L4+K)=H421(K,J,I)                                    YY 09240
      SR(L1+I,L3+J,L1+K)=H113(I,K,J)                                    YY 09250
      SR(L1+I,L3+J,L2+K)=H123(I,K,J)                                    YY 09260
      SR(L1+I,L3+J,L3+K)=H331(J,K,I)                                    YY 09270
      SR(L1+I,L3+J,L4+K)=H431(K,J,I)                                    YY 09280
      SR(L1+I,L4+J,L1+K)=H411(J,K,I)                                    YY 09290
      SR(L1+I,L4+J,L2+K)=H421(J,K,I)                                    YY 09300
      SR(L1+I,L4+J,L3+K)=H431(J,K,I)                                    YY 09310
      SR(L1+I,L4+J,L4+K)=H441(J,K,I)                                    YY 09320
      SR(L2+I,L1+J,L1+K)=H112(J,K,I)                                    YY 09330
      SR(L2+I,L1+J,L2+K)=H221(I,K,J)                                    YY 09340
      SR(L2+I,L1+J,L3+K)=H123(J,I,K)                                    YY 09350
      SR(L2+I,L1+J,L4+K)=H421(K,I,J)                                    YY 09360
      SR(L2+I,L2+J,L1+K)=H221(I,J,K)                                    YY 09370
      SR(L2+I,L2+J,L2+K)=H222(I,J,K)                                    YY 09380
      SR(L2+I,L2+J,L3+K)=H223(I,J,K)                                    YY 09390
      SR(L2+I,L2+J,L4+K)=H422(K,J,I)                                    YY 09400
      SR(L2+I,L3+J,L1+K)=H123(K,I,J)                                    YY 09410
      SR(L2+I,L3+J,L2+K)=H223(I,K,J)                                    YY 09420
      SR(L2+I,L3+J,L3+K)=H332(J,K,I)                                    YY 09430
      SR(L2+I,L3+J,L4+K)=H432(K,J,I)                                    YY 09440
      SR(L2+I,L4+J,L1+K)=H421(J,I,K)                                    YY 09450
      SR(L2+I,L4+J,L2+K)=H422(J,I,K)                                    YY 09460
      SR(L2+I,L4+J,L3+K)=H432(J,K,I)                                    YY 09470
      SR(L2+I,L4+J,L4+K)=H442(J,K,I)                                    YY 09480
      SR(L3+I,L1+J,L1+K)=H113(J,K,I)                                    YY 09490
      SR(L3+I,L1+J,L2+K)=H123(J,K,I)                                    YY 09500
      SR(L3+I,L1+J,L3+K)=H331(I,K,J)                                    YY 09510
      SR(L3+I,L1+J,L4+K)=H431(K,I,J)                                    YY 09520
      SR(L3+I,L2+J,L1+K)=H123(K,J,I)                                    YY 09530
      SR(L3+I,L2+J,L2+K)=H223(J,K,I)                                    YY 09540
      SR(L3+I,L2+J,L3+K)=H332(I,K,J)                                    YY 09550
      SR(L3+I,L2+J,L4+K)=H432(K,I,J)                                    YY 09560
      SR(L3+I,L3+J,L1+K)=H331(I,J,K)                                    YY 09570
      SR(L3+I,L3+J,L2+K)=H332(I,J,K)                                    YY 09580
      SR(L3+I,L3+J,L3+K)=H333(I,J,K)                                    YY 09590
      SR(L3+I,L3+J,L4+K)=H433(K,I,J)                                    YY 09600
      SR(L3+I,L4+J,L1+K)=H431(J,I,K)                                    YY 09610
      SR(L3+I,L4+J,L2+K)=H432(J,I,K)                                    YY 09620
      SR(L3+I,L4+J,L3+K)=H433(J,I,K)                                    YY 09630
      SR(L3+I,L4+J,L4+K)=H443(J,K,I)                                    YY 09640
      SR(L4+I,L1+J,L1+K)=H411(I,J,K)                                    YY 09650
      SR(L4+I,L1+J,L2+K)=H421(I,K,J)                                    YY 09660
      SR(L4+I,L1+J,L3+K)=H431(I,K,J)                                    YY 09670
      SR(L4+I,L1+J,L4+K)=H441(I,K,J)                                    YY 09680
      SR(L4+I,L2+J,L1+K)=H421(I,J,K)                                    YY 09690
      SR(L4+I,L2+J,L2+K)=H422(I,J,K)                                    YY 09700
      SR(L4+I,L2+J,L3+K)=H432(I,K,J)                                    YY 09710
      SR(L4+I,L2+J,L4+K)=H442(I,K,J)                                    YY 09720
      SR(L4+I,L3+J,L1+K)=H431(I,J,K)                                    YY 09730
      SR(L4+I,L3+J,L2+K)=H432(I,J,K)                                    YY 09740
      SR(L4+I,L3+J,L3+K)=H433(I,J,K)                                    YY 09750
      SR(L4+I,L3+J,L4+K)=H443(I,K,J)                                    YY 09760
      SR(L4+I,L4+J,L1+K)=H441(I,J,K)                                    YY 09770
      SR(L4+I,L4+J,L2+K)=H442(I,J,K)                                    YY 09780
      SR(L4+I,L4+J,L3+K)=H443(I,J,K)                                    YY 09790
 137  SR(L4+I,L4+J,L4+K)=H444(I,J,K)                                    YY 09800
      GO TO 300                                                         YY 09810
 150  IF(TYPE(R).NE.'  OUT') GO TO 175                                  YY 09820
      K1=IA(R,1)                                                        YY 09830
      K2=IA(R,2)                                                        YY 09840
      K3=IA(R,3)                                                        YY 09850
      K4=IA(R,4)                                                        YY 09860
      L1=3*(K1-1)                                                       YY 09870
      L2=3*(K2-1)                                                       YY 09880
      L3=3*(K3-1)                                                       YY 09890
      L4=3*(K4-1)                                                       YY 09900
      CALL HIJKS7(NAD,K1,K2,K3,K4,XA,H111,H112,H221,H222,H113,H123,     YY 09910
     $      H223,H331,H332,H333,H411,H421,H422,H431,H432,H433,H441,     YY 09920
     $      H442,H443,H444)                                             YY 09930
C  OPTION                                                               YY 09940
      DO 160  K=1,3                                                     YY 09950
      DO 160  J=1,3                                                     YY 09960
      DO 160  I=1,3                                                     YY 09970
      SR(L1+I,L1+J,L1+K)=H111(I,J,K)                                    YY 09980
      SR(L1+I,L1+J,L2+K)=H112(I,J,K)                                    YY 09990
      SR(L1+I,L1+J,L3+K)=H113(I,J,K)                                    YY 10000
      SR(L1+I,L1+J,L4+K)=H411(K,J,I)                                    YY 10010
      SR(L1+I,L2+J,L1+K)=H112(I,K,J)                                    YY 10020
      SR(L1+I,L2+J,L2+K)=H221(J,K,I)                                    YY 10030
      SR(L1+I,L2+J,L3+K)=H123(I,J,K)                                    YY 10040
      SR(L1+I,L2+J,L4+K)=H421(K,J,I)                                    YY 10050
      SR(L1+I,L3+J,L1+K)=H113(I,K,J)                                    YY 10060
      SR(L1+I,L3+J,L2+K)=H123(I,K,J)                                    YY 10070
      SR(L1+I,L3+J,L3+K)=H331(J,K,I)                                    YY 10080
      SR(L1+I,L3+J,L4+K)=H431(K,J,I)                                    YY 10090
      SR(L1+I,L4+J,L1+K)=H411(J,K,I)                                    YY 10100
      SR(L1+I,L4+J,L2+K)=H421(J,K,I)                                    YY 10110
      SR(L1+I,L4+J,L3+K)=H431(J,K,I)                                    YY 10120
      SR(L1+I,L4+J,L4+K)=H441(J,K,I)                                    YY 10130
      SR(L2+I,L1+J,L1+K)=H112(J,K,I)                                    YY 10140
      SR(L2+I,L1+J,L2+K)=H221(I,K,J)                                    YY 10150
      SR(L2+I,L1+J,L3+K)=H123(J,I,K)                                    YY 10160
      SR(L2+I,L1+J,L4+K)=H421(K,I,J)                                    YY 10170
      SR(L2+I,L2+J,L1+K)=H221(I,J,K)                                    YY 10180
      SR(L2+I,L2+J,L2+K)=H222(I,J,K)                                    YY 10190
      SR(L2+I,L2+J,L3+K)=H223(I,J,K)                                    YY 10200
      SR(L2+I,L2+J,L4+K)=H422(K,J,I)                                    YY 10210
      SR(L2+I,L3+J,L1+K)=H123(K,I,J)                                    YY 10220
      SR(L2+I,L3+J,L2+K)=H223(I,K,J)                                    YY 10230
      SR(L2+I,L3+J,L3+K)=H332(J,K,I)                                    YY 10240
      SR(L2+I,L3+J,L4+K)=H432(K,J,I)                                    YY 10250
      SR(L2+I,L4+J,L1+K)=H421(J,I,K)                                    YY 10260
      SR(L2+I,L4+J,L2+K)=H422(J,I,K)                                    YY 10270
      SR(L2+I,L4+J,L3+K)=H432(J,K,I)                                    YY 10280
      SR(L2+I,L4+J,L4+K)=H442(J,K,I)                                    YY 10290
      SR(L3+I,L1+J,L1+K)=H113(J,K,I)                                    YY 10300
      SR(L3+I,L1+J,L2+K)=H123(J,K,I)                                    YY 10310
      SR(L3+I,L1+J,L3+K)=H331(I,K,J)                                    YY 10320
      SR(L3+I,L1+J,L4+K)=H431(K,I,J)                                    YY 10330
      SR(L3+I,L2+J,L1+K)=H123(K,J,I)                                    YY 10340
      SR(L3+I,L2+J,L2+K)=H223(J,K,I)                                    YY 10350
      SR(L3+I,L2+J,L3+K)=H332(I,K,J)                                    YY 10360
      SR(L3+I,L2+J,L4+K)=H432(K,I,J)                                    YY 10370
      SR(L3+I,L3+J,L1+K)=H331(I,J,K)                                    YY 10380
      SR(L3+I,L3+J,L2+K)=H332(I,J,K)                                    YY 10390
      SR(L3+I,L3+J,L3+K)=H333(I,J,K)                                    YY 10400
      SR(L3+I,L3+J,L4+K)=H433(K,I,J)                                    YY 10410
      SR(L3+I,L4+J,L1+K)=H431(J,I,K)                                    YY 10420
      SR(L3+I,L4+J,L2+K)=H432(J,I,K)                                    YY 10430
      SR(L3+I,L4+J,L3+K)=H433(J,I,K)                                    YY 10440
      SR(L3+I,L4+J,L4+K)=H443(J,K,I)                                    YY 10450
      SR(L4+I,L1+J,L1+K)=H411(I,J,K)                                    YY 10460
      SR(L4+I,L1+J,L2+K)=H421(I,K,J)                                    YY 10470
      SR(L4+I,L1+J,L3+K)=H431(I,K,J)                                    YY 10480
      SR(L4+I,L1+J,L4+K)=H441(I,K,J)                                    YY 10490
      SR(L4+I,L2+J,L1+K)=H421(I,J,K)                                    YY 10500
      SR(L4+I,L2+J,L2+K)=H422(I,J,K)                                    YY 10510
      SR(L4+I,L2+J,L3+K)=H432(I,K,J)                                    YY 10520
      SR(L4+I,L2+J,L4+K)=H442(I,K,J)                                    YY 10530
      SR(L4+I,L3+J,L1+K)=H431(I,J,K)                                    YY 10540
      SR(L4+I,L3+J,L2+K)=H432(I,J,K)                                    YY 10550
      SR(L4+I,L3+J,L3+K)=H433(I,J,K)                                    YY 10560
      SR(L4+I,L3+J,L4+K)=H443(I,K,J)                                    YY 10570
      SR(L4+I,L4+J,L1+K)=H441(I,J,K)                                    YY 10580
      SR(L4+I,L4+J,L2+K)=H442(I,J,K)                                    YY 10590
      SR(L4+I,L4+J,L3+K)=H443(I,J,K)                                    YY 10600
 160  SR(L4+I,L4+J,L4+K)=H444(I,J,K)                                    YY 10610
      GO TO 300                                                         YY 10620
 175  IF(TYPE(R).NE.' LINX'.AND.TYPE(R).NE.' LINY') GO TO 200           YY 03710
      K1=IA(R,1)                                                        YY 09020
      K2=IA(R,2)                                                        YY 09030
      K3=IA(R,3)                                                        YY 09040
      K4=IA(R,4)                                                        YY 09050
      L1=3*(K1-1)                                                       YY 09060
      L2=3*(K2-1)                                                       YY 09070
      L3=3*(K3-1)                                                       YY 09080
      L4=3*(K4-1)                                                       YY 09090
      IF(TYPE(R).EQ.' LINX') CALL HIJKS8(NAD,K1,K2,K3,K4,XA,
     $      H111,H112,H221,H222,H113,H123,H223,H331,H332,H333,
     $      H411,H421,H422,H431,H432,H433,H441,H442,H443,H444)
      IF(TYPE(R).EQ.' LINY') CALL HIJKS9(NAD,K1,K2,K3,K4,XA,
     $      H111,H112,H221,H222,H113,H123,H223,H331,H332,H333,
     $      H411,H421,H422,H431,H432,H433,H441,H442,H443,H444)
C  OPTION                                                               YY 09130
      DO 177  K=1,3                                                     YY 09140
      DO 177  J=1,3                                                     YY 09150
      DO 177  I=1,3                                                     YY 09160
      SR(L1+I,L1+J,L1+K)=H111(I,J,K)                                    YY 09170
      SR(L1+I,L1+J,L2+K)=H112(I,J,K)                                    YY 09180
      SR(L1+I,L1+J,L3+K)=H113(I,J,K)                                    YY 09190
      SR(L1+I,L1+J,L4+K)=H411(K,J,I)                                    YY 09200
      SR(L1+I,L2+J,L1+K)=H112(I,K,J)                                    YY 09210
      SR(L1+I,L2+J,L2+K)=H221(J,K,I)                                    YY 09220
      SR(L1+I,L2+J,L3+K)=H123(I,J,K)                                    YY 09230
      SR(L1+I,L2+J,L4+K)=H421(K,J,I)                                    YY 09240
      SR(L1+I,L3+J,L1+K)=H113(I,K,J)                                    YY 09250
      SR(L1+I,L3+J,L2+K)=H123(I,K,J)                                    YY 09260
      SR(L1+I,L3+J,L3+K)=H331(J,K,I)                                    YY 09270
      SR(L1+I,L3+J,L4+K)=H431(K,J,I)                                    YY 09280
      SR(L1+I,L4+J,L1+K)=H411(J,K,I)                                    YY 09290
      SR(L1+I,L4+J,L2+K)=H421(J,K,I)                                    YY 09300
      SR(L1+I,L4+J,L3+K)=H431(J,K,I)                                    YY 09310
      SR(L1+I,L4+J,L4+K)=H441(J,K,I)                                    YY 09320
      SR(L2+I,L1+J,L1+K)=H112(J,K,I)                                    YY 09330
      SR(L2+I,L1+J,L2+K)=H221(I,K,J)                                    YY 09340
      SR(L2+I,L1+J,L3+K)=H123(J,I,K)                                    YY 09350
      SR(L2+I,L1+J,L4+K)=H421(K,I,J)                                    YY 09360
      SR(L2+I,L2+J,L1+K)=H221(I,J,K)                                    YY 09370
      SR(L2+I,L2+J,L2+K)=H222(I,J,K)                                    YY 09380
      SR(L2+I,L2+J,L3+K)=H223(I,J,K)                                    YY 09390
      SR(L2+I,L2+J,L4+K)=H422(K,J,I)                                    YY 09400
      SR(L2+I,L3+J,L1+K)=H123(K,I,J)                                    YY 09410
      SR(L2+I,L3+J,L2+K)=H223(I,K,J)                                    YY 09420
      SR(L2+I,L3+J,L3+K)=H332(J,K,I)                                    YY 09430
      SR(L2+I,L3+J,L4+K)=H432(K,J,I)                                    YY 09440
      SR(L2+I,L4+J,L1+K)=H421(J,I,K)                                    YY 09450
      SR(L2+I,L4+J,L2+K)=H422(J,I,K)                                    YY 09460
      SR(L2+I,L4+J,L3+K)=H432(J,K,I)                                    YY 09470
      SR(L2+I,L4+J,L4+K)=H442(J,K,I)                                    YY 09480
      SR(L3+I,L1+J,L1+K)=H113(J,K,I)                                    YY 09490
      SR(L3+I,L1+J,L2+K)=H123(J,K,I)                                    YY 09500
      SR(L3+I,L1+J,L3+K)=H331(I,K,J)                                    YY 09510
      SR(L3+I,L1+J,L4+K)=H431(K,I,J)                                    YY 09520
      SR(L3+I,L2+J,L1+K)=H123(K,J,I)                                    YY 09530
      SR(L3+I,L2+J,L2+K)=H223(J,K,I)                                    YY 09540
      SR(L3+I,L2+J,L3+K)=H332(I,K,J)                                    YY 09550
      SR(L3+I,L2+J,L4+K)=H432(K,I,J)                                    YY 09560
      SR(L3+I,L3+J,L1+K)=H331(I,J,K)                                    YY 09570
      SR(L3+I,L3+J,L2+K)=H332(I,J,K)                                    YY 09580
      SR(L3+I,L3+J,L3+K)=H333(I,J,K)                                    YY 09590
      SR(L3+I,L3+J,L4+K)=H433(K,I,J)                                    YY 09600
      SR(L3+I,L4+J,L1+K)=H431(J,I,K)                                    YY 09610
      SR(L3+I,L4+J,L2+K)=H432(J,I,K)                                    YY 09620
      SR(L3+I,L4+J,L3+K)=H433(J,I,K)                                    YY 09630
      SR(L3+I,L4+J,L4+K)=H443(J,K,I)                                    YY 09640
      SR(L4+I,L1+J,L1+K)=H411(I,J,K)                                    YY 09650
      SR(L4+I,L1+J,L2+K)=H421(I,K,J)                                    YY 09660
      SR(L4+I,L1+J,L3+K)=H431(I,K,J)                                    YY 09670
      SR(L4+I,L1+J,L4+K)=H441(I,K,J)                                    YY 09680
      SR(L4+I,L2+J,L1+K)=H421(I,J,K)                                    YY 09690
      SR(L4+I,L2+J,L2+K)=H422(I,J,K)                                    YY 09700
      SR(L4+I,L2+J,L3+K)=H432(I,K,J)                                    YY 09710
      SR(L4+I,L2+J,L4+K)=H442(I,K,J)                                    YY 09720
      SR(L4+I,L3+J,L1+K)=H431(I,J,K)                                    YY 09730
      SR(L4+I,L3+J,L2+K)=H432(I,J,K)                                    YY 09740
      SR(L4+I,L3+J,L3+K)=H433(I,J,K)                                    YY 09750
      SR(L4+I,L3+J,L4+K)=H443(I,K,J)                                    YY 09760
      SR(L4+I,L4+J,L1+K)=H441(I,J,K)                                    YY 09770
      SR(L4+I,L4+J,L2+K)=H442(I,J,K)                                    YY 09780
      SR(L4+I,L4+J,L3+K)=H443(I,J,K)                                    YY 09790
 177  SR(L4+I,L4+J,L4+K)=H444(I,J,K)                                    YY 09800
      GO TO 300                                                         YY 09810
 200  CONTINUE                                                          YY 10630
      GO TO 300                                                         YY 10640
 300  CONTINUE                                                          YY 10650
      RETURN                                                            YY 10660
      END                                                               YY 10670
C     /////////////////////////////////////////////////////////////////
      SUBROUTINE MACHZ(NAD,NC,NS,NSX,IOPT,XA,TYPE,IA,A,S,U,IU,YR1,YR2)  YY 01070
      IMPLICIT REAL*8 (A-H,O-Z)                                         YY 01080
      INTEGER R,P,Q,NSX,RP,IOPT(30)                                     YY 01090
      CHARACTER TYPE*5                                                  YY 01100
      DIMENSION TYPE(NS),IA(NS,6),S(NS),U(NS,1),IU(NS,0:1)              YY 01110
      DIMENSION XA(NAD,3),A(NC,NC),YR1(NC,NC,NC),YR2(NC,NC,NC)          YY 01120
      PARAMETER(ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0,THREE=3.0D0)             YY 01130
      COMMON /IO/ IIN1,IOUT,IIN2,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25,
     $   I31,I32,I33,I35,I36,I37,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10,
     $   ISCR11,ISCR12,ISCR13,ISCR14
    1 FORMAT(/,1X,'NUMERICAL SR(I,J,K,L) AND Z(M,N,P,Q) MATRICES USED', YY 01180
     $ /, ' FOR SIMPLE INTERNAL COORDINATE',I5)                         YY 01190
      NSYM=IOPT(3)                                                      YY 01200
      IF(NSYM.NE.0) THEN                                                YY 01210
         ISCR=ISCR10                                                    YY 01220
      ELSE                                                              YY 01230
         ISCR=ISCR9                                                     YY 01240
      END IF                                                            YY 01250
C                                                                       YY 01260
CCC  BLOCK FOR INSERTION OF ANALYTICAL Z MATRIX SECTION                 YY 01270
C                                                                       YY 01280
C     DO 500  R=1,NS                                                    YY 01290
C READ NUMERICAL Z MATRIX AND TRANSFER IF DESIRED.                      YY 01300
C     CALL ZIN(NC,NC,SR,-R,ISCR8)                                       YY 01310
C     CALL ZIN(NC,NSX,Z,R,ISCR8)                                        YY 01320
C     WRITE(IOUT,1) R                                                   YY 01330
C300  CONTINUE                                                          YY 01340
C     CALL FILL4A(NC,NC,SR)                                             YY 01350
C     CALL FILL4A(NC,NSX,Z)                                             YY 01360
C     CALL ZOUT(NC,NC,SR,-R,ISCR)                                       YY 01370
C     CALL ZOUT(NC,NSX,Z,R,ISCR)                                        YY 01380
C500  CONTINUE                                                          YY 01390
C                                                                       YY 01400
 600  IF(NSYM.EQ.0) RETURN                                              YY 01410
      DO 650  R=1,NSYM                                                  YY 01450
      DO 660  Q=1,NSYM                                                  YY 01490
         DO 610  P=1,Q                                                  YY 01500
         DO 610  N=1,P                                                  YY 01510
         DO 610  M=1,N                                                  YY 01520
 610     YR1(M,N,P)=ZERO                                                YY 01530
         L=IU(R,0)                                                      YY 01540
         DO 620  I=1,L                                                  YY 01550
         RP=IU(R,I)                                                     YY 01560
         CALL YIN2(NC,NSYM,YR2,RP,Q,ISCR10)                             YY 01570
         W1=U(R,I)                                                      YY 01580
         DO 630  P=1,Q                                                  YY 01590
         DO 630  N=1,P                                                  YY 01600
         DO 630  M=1,N                                                  YY 01610
 630     YR1(M,N,P)=YR1(M,N,P)+W1*YR2(M,N,P)                            YY 01620
 620     CONTINUE                                                       YY 01630
         CALL FILL3A(NC,NSYM,YR1)                                       YY 01640
 660     CALL YOUT2(NC,NSYM,YR1,R,Q,ISCR9)                              YY 01650
 650  CONTINUE                                                          YY 01660
C   OPTION                                                              YY 01670
      DO 750  R=1,NSYM                                                  YY 01680
      DO 760  Q=1,NC                                                    YY 01720
         DO 710  P=1,Q                                                  YY 01730
         DO 710  N=1,P                                                  YY 01740
         DO 710  M=1,N                                                  YY 01750
 710     YR1(M,N,P)=ZERO                                                YY 01760
         L=IU(R,0)                                                      YY 01770
         DO 720  I=1,L                                                  YY 01780
         RP=IU(R,I)                                                     YY 01790
         CALL YIN2(NC,NC,YR2,-RP,Q,ISCR10)                              YY 01800
         W1=U(R,I)                                                      YY 01810
         DO 730  P=1,Q                                                  YY 01820
         DO 730  N=1,P                                                  YY 01830
         DO 730  M=1,N                                                  YY 01840
 730     YR1(M,N,P)=YR1(M,N,P)+W1*YR2(M,N,P)                            YY 01850
 720     CONTINUE                                                       YY 01860
         CALL FILL3A(NC,NC,YR1)                                         YY 01870
 760     CALL YOUT2(NC,NC,YR1,-R,Q,ISCR9)                               YY 01880
 750  CONTINUE                                                          YY 01890
      RETURN                                                            YY 01900
      END                                                               YY 01910
C     ///////////////////////////////////////////////////////////////// YY 05960
      SUBROUTINE HSRY2(NC,L1,L2,H111,SR)                                YY 05970
      IMPLICIT REAL*8 (A-H,O-Z)                                         YY 05980
      DIMENSION SR(NC,NC,NC)                                            YY 05990
      DIMENSION H111(3,3,3)                                             YY 06000
      DO 10  K=1,3                                                      YY 06010
      DO 10  J=1,3                                                      YY 06020
      DO 10  I=1,3                                                      YY 06030
      Z=H111(I,J,K)                                                     YY 06040
      SR(L1+I,L1+J,L1+K)=Z                                              YY 06050
      SR(L1+I,L1+J,L2+K)=-Z                                             YY 06060
      SR(L1+I,L2+J,L1+K)=-Z                                             YY 06070
      SR(L1+I,L2+J,L2+K)=Z                                              YY 06080
      SR(L2+I,L1+J,L1+K)=-Z                                             YY 06090
      SR(L2+I,L1+J,L2+K)=Z                                              YY 06100
      SR(L2+I,L2+J,L1+K)=Z                                              YY 06110
 10   SR(L2+I,L2+J,L2+K)=-Z                                             YY 06120
      RETURN                                                            YY 06130
      END                                                               YY 06140
C     ///////////////////////////////////////////////////////////////// YY 06150
      SUBROUTINE HSRY3(NC,L1,L2,L3,H111,H112,H113,H123,H221,            YY 06160
     $    H222,H223,H331,H332,H333,SR)                                  YY 06170
      IMPLICIT REAL*8 (A-H,O-Z)                                         YY 06180
      DIMENSION SR(NC,NC,NC)                                            YY 06190
      DIMENSION H111(3,3,3),H112(3,3,3),H221(3,3,3),H222(3,3,3)         YY 06200
      DIMENSION H113(3,3,3),H123(3,3,3),H223(3,3,3),H331(3,3,3)         YY 06210
      DIMENSION H332(3,3,3),H333(3,3,3)                                 YY 06220
      DO 10  K=1,3                                                      YY 06230
      DO 10  J=1,3                                                      YY 06240
      DO 10  I=1,3                                                      YY 06250
      SR(L1+I,L1+J,L1+K)=H111(I,J,K)                                    YY 06260
      SR(L1+I,L1+J,L2+K)=H112(I,J,K)                                    YY 06270
      SR(L1+I,L1+J,L3+K)=H113(I,J,K)                                    YY 06280
      SR(L1+I,L2+J,L1+K)=H112(I,K,J)                                    YY 06290
      SR(L1+I,L2+J,L2+K)=H221(J,K,I)                                    YY 06300
      SR(L1+I,L2+J,L3+K)=H123(I,J,K)                                    YY 06310
      SR(L2+I,L2+J,L1+K)=H221(I,J,K)                                    YY 06320
      SR(L2+I,L2+J,L2+K)=H222(I,J,K)                                    YY 06330
      SR(L2+I,L2+J,L3+K)=H223(I,J,K)                                    YY 06340
      SR(L2+I,L1+J,L1+K)=H112(J,K,I)                                    YY 06350
      SR(L2+I,L1+J,L2+K)=H221(I,K,J)                                    YY 06360
      SR(L2+I,L1+J,L3+K)=H123(J,I,K)                                    YY 06370
      SR(L1+I,L3+J,L1+K)=H113(I,K,J)                                    YY 06380
      SR(L1+I,L3+J,L2+K)=H123(I,K,J)                                    YY 06390
      SR(L1+I,L3+J,L3+K)=H331(J,K,I)                                    YY 06400
      SR(L2+I,L3+J,L1+K)=H123(K,I,J)                                    YY 06410
      SR(L2+I,L3+J,L2+K)=H223(I,K,J)                                    YY 06420
      SR(L2+I,L3+J,L3+K)=H332(J,K,I)                                    YY 06430
      SR(L3+I,L1+J,L1+K)=H113(J,K,I)                                    YY 06440
      SR(L3+I,L1+J,L2+K)=H123(J,K,I)                                    YY 06450
      SR(L3+I,L1+J,L3+K)=H331(I,K,J)                                    YY 06460
      SR(L3+I,L2+J,L1+K)=H123(K,J,I)                                    YY 06470
      SR(L3+I,L2+J,L2+K)=H223(J,K,I)                                    YY 06480
      SR(L3+I,L2+J,L3+K)=H332(I,K,J)                                    YY 06490
      SR(L3+I,L3+J,L1+K)=H331(I,J,K)                                    YY 06500
      SR(L3+I,L3+J,L2+K)=H332(I,J,K)                                    YY 06510
 10   SR(L3+I,L3+J,L3+K)=H333(I,J,K)                                    YY 06520
      RETURN                                                            YY 06530
      END                                                               YY 06540
C     ////////////////////////////////////////////////////////////////  YY 06550
      SUBROUTINE HSRY4(NC,L1,L2,L3,L4,H111,H112,H221,H222,H113,H123,    YY 06560
     $     H223,H331,H332,H333,H411,H421,H422,H431,H432,H433,H441,      YY 06570
     $     H442,H443,H444,SR)                                           YY 06580
      IMPLICIT REAL*8 (A-H,O-Z)                                         YY 06590
      DIMENSION SR(NC,NC,NC)                                            YY 06600
      DIMENSION H111(3,3,3),H112(3,3,3),H221(3,3,3),H222(3,3,3)         YY 06610
      DIMENSION H113(3,3,3),H123(3,3,3),H223(3,3,3),H331(3,3,3)         YY 06620
      DIMENSION H332(3,3,3),H333(3,3,3),H411(3,3,3),H421(3,3,3)         YY 06630
      DIMENSION H422(3,3,3),H431(3,3,3),H432(3,3,3),H433(3,3,3)         YY 06640
      DIMENSION H441(3,3,3),H442(3,3,3),H443(3,3,3),H444(3,3,3)         YY 06650
      DO 10  K=1,3                                                      YY 06660
      DO 10  J=1,3                                                      YY 06670
      DO 10  I=1,3                                                      YY 06680
      SR(L1+I,L1+J,L1+K)=H111(I,J,K)                                    YY 06690
      SR(L1+I,L1+J,L2+K)=H112(I,J,K)                                    YY 06700
      SR(L1+I,L1+J,L3+K)=H113(I,J,K)                                    YY 06710
      SR(L1+I,L1+J,L4+K)=H411(K,J,I)                                    YY 06720
      SR(L1+I,L2+J,L1+K)=H112(I,K,J)                                    YY 06730
      SR(L1+I,L2+J,L2+K)=H221(J,K,I)                                    YY 06740
      SR(L1+I,L2+J,L3+K)=H123(I,J,K)                                    YY 06750
      SR(L1+I,L2+J,L4+K)=H421(K,J,I)                                    YY 06760
      SR(L1+I,L3+J,L1+K)=H113(I,K,J)                                    YY 06770
      SR(L1+I,L3+J,L2+K)=H123(I,K,J)                                    YY 06780
      SR(L1+I,L3+J,L3+K)=H331(J,K,I)                                    YY 06790
      SR(L1+I,L3+J,L4+K)=H431(K,J,I)                                    YY 06800
      SR(L1+I,L4+J,L1+K)=H411(J,K,I)                                    YY 06810
      SR(L1+I,L4+J,L2+K)=H421(J,K,I)                                    YY 06820
      SR(L1+I,L4+J,L3+K)=H431(J,K,I)                                    YY 06830
      SR(L1+I,L4+J,L4+K)=H441(J,K,I)                                    YY 06840
      SR(L2+I,L1+J,L1+K)=H112(J,K,I)                                    YY 06850
      SR(L2+I,L1+J,L2+K)=H221(I,K,J)                                    YY 06860
      SR(L2+I,L1+J,L3+K)=H123(J,I,K)                                    YY 06870
      SR(L2+I,L1+J,L4+K)=H421(K,I,J)                                    YY 06880
      SR(L2+I,L2+J,L1+K)=H221(I,J,K)                                    YY 06890
      SR(L2+I,L2+J,L2+K)=H222(I,J,K)                                    YY 06900
      SR(L2+I,L2+J,L3+K)=H223(I,J,K)                                    YY 06910
      SR(L2+I,L2+J,L4+K)=H422(K,J,I)                                    YY 06920
      SR(L2+I,L3+J,L1+K)=H123(K,I,J)                                    YY 06930
      SR(L2+I,L3+J,L2+K)=H223(I,K,J)                                    YY 06940
      SR(L2+I,L3+J,L3+K)=H332(J,K,I)                                    YY 06950
      SR(L2+I,L3+J,L4+K)=H432(K,J,I)                                    YY 06960
      SR(L2+I,L4+J,L1+K)=H421(J,I,K)                                    YY 06970
      SR(L2+I,L4+J,L2+K)=H422(J,I,K)                                    YY 06980
      SR(L2+I,L4+J,L3+K)=H432(J,K,I)                                    YY 06990
      SR(L2+I,L4+J,L4+K)=H442(J,K,I)                                    YY 07000
      SR(L3+I,L1+J,L1+K)=H113(J,K,I)                                    YY 07010
      SR(L3+I,L1+J,L2+K)=H123(J,K,I)                                    YY 07020
      SR(L3+I,L1+J,L3+K)=H331(I,K,J)                                    YY 07030
      SR(L3+I,L1+J,L4+K)=H431(K,I,J)                                    YY 07040
      SR(L3+I,L2+J,L1+K)=H123(K,J,I)                                    YY 07050
      SR(L3+I,L2+J,L2+K)=H223(J,K,I)                                    YY 07060
      SR(L3+I,L2+J,L3+K)=H332(I,K,J)                                    YY 07070
      SR(L3+I,L2+J,L4+K)=H432(K,I,J)                                    YY 07080
      SR(L3+I,L3+J,L1+K)=H331(I,J,K)                                    YY 07090
      SR(L3+I,L3+J,L2+K)=H332(I,J,K)                                    YY 07100
      SR(L3+I,L3+J,L3+K)=H333(I,J,K)                                    YY 07110
      SR(L3+I,L3+J,L4+K)=H433(K,I,J)                                    YY 07120
      SR(L3+I,L4+J,L1+K)=H431(J,I,K)                                    YY 07130
      SR(L3+I,L4+J,L2+K)=H432(J,I,K)                                    YY 07140
      SR(L3+I,L4+J,L3+K)=H433(J,I,K)                                    YY 07150
      SR(L3+I,L4+J,L4+K)=H443(J,K,I)                                    YY 07160
      SR(L4+I,L1+J,L1+K)=H411(I,J,K)                                    YY 07170
      SR(L4+I,L1+J,L2+K)=H421(I,K,J)                                    YY 07180
      SR(L4+I,L1+J,L3+K)=H431(I,K,J)                                    YY 07190
      SR(L4+I,L1+J,L4+K)=H441(I,K,J)                                    YY 07200
      SR(L4+I,L2+J,L1+K)=H421(I,J,K)                                    YY 07210
      SR(L4+I,L2+J,L2+K)=H422(I,J,K)                                    YY 07220
      SR(L4+I,L2+J,L3+K)=H432(I,K,J)                                    YY 07230
      SR(L4+I,L2+J,L4+K)=H442(I,K,J)                                    YY 07240
      SR(L4+I,L3+J,L1+K)=H431(I,J,K)                                    YY 07250
      SR(L4+I,L3+J,L2+K)=H432(I,J,K)                                    YY 07260
      SR(L4+I,L3+J,L3+K)=H433(I,J,K)                                    YY 07270
      SR(L4+I,L3+J,L4+K)=H443(I,K,J)                                    YY 07280
      SR(L4+I,L4+J,L1+K)=H441(I,J,K)                                    YY 07290
      SR(L4+I,L4+J,L2+K)=H442(I,J,K)                                    YY 07300
      SR(L4+I,L4+J,L3+K)=H443(I,J,K)                                    YY 07310
 10   SR(L4+I,L4+J,L4+K)=H444(I,J,K)                                    YY 07320
      RETURN                                                            YY 07330
      END                                                               YY 07340
C     ////////////////////////////////////////////////////////////
      SUBROUTINE AHX(NC,NSX,SR,A,X)
      IMPLICIT REAL*8 (A-H,O-Z)                                         INT07930
      DIMENSION A(NC,NC),X(NC,NC),SR(NC,NC)
      PARAMETER(ZERO=0.0D0)
      DO 5  L=1,NSX
      DO 5  M=1,L
         X(L,M)=ZERO
         DO 5  I=1,NC
         DO 5  J=1,NC
 5       X(L,M)=X(L,M)+SR(I,J)*A(I,L)*A(J,M) 
      DO 10  L=2,NSX
      DO 10  M=1,L-1
 10      X(M,L)=X(L,M)
      RETURN
      END
C     ////////////////////////////////////////////////////////////      INT14360
      SUBROUTINE AHX2(NC,NSX,L1,L2,H11,A,X)                             INT07920
      IMPLICIT REAL*8 (A-H,O-Z)                                         INT07930
      DIMENSION A(NC,NC),X(NC,NC),H11(3,3)                              INT07940
      DO 10  N=1,NSX                                                    INT07950
      DO 10  M=1,N                                                      INT07960
          DO 15  I=1,3                                                  INT07970
          DO 15  J=1,3                                                  INT07980
          W1=(A(L1+I,M)-A(L2+I,M))*(A(L1+J,N)-A(L2+J,N))                INT07990
  15      X(M,N)=X(M,N)+W1*H11(I,J)                                     INT08000
  10  CONTINUE                                                          INT08010
      DO 20  N=1,NSX                                                    INT08020
      DO 20  M=1,N                                                      INT08030
  20  X(N,M)=X(M,N)                                                     INT08040
      RETURN                                                            INT08050
      END                                                               INT08060
C     /////////////////////////////////////////////////////////////     INT08070
      SUBROUTINE AHX3(NC,NSX,L1,L2,L3,H11,H21,H31,H22,H32,H33,A,X)      INT08080
      IMPLICIT REAL*8 (A-H,O-Z)                                         INT08090
      DIMENSION A(NC,NC),X(NC,NC)                                       INT08100
      DIMENSION H11(3,3),H21(3,3),H31(3,3),H22(3,3),H32(3,3),H33(3,3)   INT08110
      DO 10  N=1,NSX                                                    INT08120
      DO 10  M=1,N                                                      INT08130
          DO 15  I=1,3                                                  INT08140
          DO 15  J=1,3                                                  INT08150
          W1=A(L1+I,M)*A(L1+J,N)                                        INT08160
          W2=A(L2+I,M)*A(L2+J,N)                                        INT08170
          W3=A(L3+I,M)*A(L3+J,N)                                        INT08180
          X(M,N)=X(M,N)+W1*H11(I,J)+W2*H22(I,J)+W3*H33(I,J)             INT08190
          W1=A(L2+I,M)*A(L1+J,N)+A(L1+J,M)*A(L2+I,N)                    INT08200
          W2=A(L3+I,M)*A(L1+J,N)+A(L1+J,M)*A(L3+I,N)                    INT08210
          W3=A(L3+I,M)*A(L2+J,N)+A(L2+J,M)*A(L3+I,N)                    INT08220
  15      X(M,N)=X(M,N)+W1*H21(I,J)+W2*H31(I,J)+W3*H32(I,J)             INT08230
  10  CONTINUE                                                          INT08240
      DO 20  N=1,NSX                                                    INT08250
      DO 20  M=1,N                                                      INT08260
  20  X(N,M)=X(M,N)                                                     INT08270
      RETURN                                                            INT08280
      END                                                               INT08290
C     /////////////////////////////////////////////////////////////     INT08300
      SUBROUTINE AHX4(NC,NSX,L1,L2,L3,L4,H11,H21,H31,H41,               INT08310
     $            H22,H32,H42,H33,H43,H44,A,X)                          INT08320
      IMPLICIT REAL*8 (A-H,O-Z)                                         INT08330
      DIMENSION A(NC,NC),X(NC,NC)                                       INT08340
      DIMENSION H11(3,3),H21(3,3),H31(3,3),H41(3,3),H22(3,3),H32(3,3)   INT08350
      DIMENSION H42(3,3),H33(3,3),H43(3,3),H44(3,3)                     INT08360
      DO 10  N=1,NSX                                                    INT08370
      DO 10  M=1,N                                                      INT08380
          DO 15  I=1,3                                                  INT08390
          DO 15  J=1,3                                                  INT08400
          W1=A(L1+I,M)*A(L1+J,N)                                        INT08410
          W2=A(L2+I,M)*A(L2+J,N)                                        INT08420
          W3=A(L3+I,M)*A(L3+J,N)                                        INT08430
          W4=A(L4+I,M)*A(L4+J,N)                                        INT08440
          X(M,N)=X(M,N)+W1*H11(I,J)+W2*H22(I,J)+W3*H33(I,J)+W4*H44(I,J) INT08450
          W1=A(L2+I,M)*A(L1+J,N)+A(L1+J,M)*A(L2+I,N)                    INT08460
          W2=A(L3+I,M)*A(L1+J,N)+A(L1+J,M)*A(L3+I,N)                    INT08470
          W3=A(L4+I,M)*A(L1+J,N)+A(L1+J,M)*A(L4+I,N)                    INT08480
          X(M,N)=X(M,N)+W1*H21(I,J)+W2*H31(I,J)+W3*H41(I,J)             INT08490
          W1=A(L3+I,M)*A(L2+J,N)+A(L2+J,M)*A(L3+I,N)                    INT08500
          W2=A(L4+I,M)*A(L2+J,N)+A(L2+J,M)*A(L4+I,N)                    INT08510
          W3=A(L4+I,M)*A(L3+J,N)+A(L3+J,M)*A(L4+I,N)                    INT08520
  15      X(M,N)=X(M,N)+W1*H32(I,J)+W2*H42(I,J)+W3*H43(I,J)             INT08530
  10  CONTINUE                                                          INT08540
      DO 20  N=1,NSX                                                    INT08550
      DO 20  M=1,N                                                      INT08560
  20  X(N,M)=X(M,N)                                                     INT08570
      RETURN                                                            INT08580
      END                                                               INT08590
C     /////////////////////////////////////////////////////////////     INT08600
      SUBROUTINE AHY2(NC,NSX,L1,L2,H111,A,Y)                            YY 04490
      IMPLICIT REAL*8 (A-H,O-Z)                                         YY 04500
      INTEGER P,NSX                                                     YY 04510
      DIMENSION A(NC,NC),Y(NC,NC,NC)                                    YY 04520
      DIMENSION H111(3,3,3)                                             YY 04530
      DO 10  P=1,NSX                                                    YY 04540
      DO 10  N=1,P                                                      YY 04550
      DO 10  M=1,N                                                      YY 04560
         DO 12  I=1,3                                                   YY 04570
         DO 12  J=1,3                                                   YY 04580
         DO 12  K=1,3                                                   YY 04590
      W1=A(L1+J,N)*(A(L1+K,P)-A(L2+K,P))-A(L2+J,N)*(A(L1+K,P)-A(L2+K,P))YY 04600
      W1=(A(L1+I,M)-A(L2+I,M))*W1                                       YY 04610
 12      Y(M,N,P)=Y(M,N,P)+W1*H111(I,J,K)                               YY 04620
 10   CONTINUE                                                          YY 04630
      RETURN                                                            YY 04640
      END                                                               YY 04650
C     ///////////////////////////////////////////////////////////////// YY 04660
      SUBROUTINE AHY3(NC,NSX,L1,L2,L3,H111,H112,H113,H123,H221,H222,    YY 04670
     $      H223,H331,H332,H333,A,Y)                                    YY 04680
      IMPLICIT REAL*8 (A-H,O-Z)                                         YY 04690
      INTEGER P,NSX                                                     YY 04700
      DIMENSION A(NC,NC),Y(NC,NC,NC)                                    YY 04710
      DIMENSION H111(3,3,3),H112(3,3,3),H221(3,3,3),H222(3,3,3)         YY 04720
      DIMENSION H113(3,3,3),H123(3,3,3),H223(3,3,3),H331(3,3,3)         YY 04730
      DIMENSION H332(3,3,3),H333(3,3,3)                                 YY 04740
      DO 10  P=1,NSX                                                    YY 04750
      DO 10  N=1,P                                                      YY 04760
      DO 10  M=1,N                                                      YY 04770
         DO 12  I=1,3                                                   YY 04780
         DO 12  J=1,3                                                   YY 04790
         V1=A(L1+I,M)*A(L1+J,N)                                         YY 04800
         V2=A(L2+I,M)*A(L2+J,N)                                         YY 04810
         V3=A(L3+I,M)*A(L3+J,N)                                         YY 04820
         V4=A(L1+I,M)*A(L1+J,P)                                         YY 04830
         V5=A(L1+I,N)*A(L1+J,P)                                         YY 04840
         V6=A(L2+I,M)*A(L2+J,P)                                         YY 04850
         V7=A(L2+I,N)*A(L2+J,P)                                         YY 04860
         V8=A(L3+I,M)*A(L3+J,P)                                         YY 04870
         V9=A(L3+I,N)*A(L3+J,P)                                         YY 04880
         V10=A(L1+I,M)*A(L2+J,N)                                        YY 04890
         V11=A(L1+I,M)*A(L2+J,P)                                        YY 04900
         V12=A(L1+I,N)*A(L2+J,M)                                        YY 04910
         V13=A(L1+I,N)*A(L2+J,P)                                        YY 04920
         V14=A(L1+I,P)*A(L2+J,M)                                        YY 04930
         V15=A(L1+I,P)*A(L2+J,N)                                        YY 04940
C                                                                       YY 04950
         DO 14  K=1,3                                                   YY 04960
         W1=V1*A(L1+K,P)                                                YY 04970
         W2=V2*A(L2+K,P)                                                YY 04980
         W3=V3*A(L3+K,P)                                                YY 04990
         Y(M,N,P)=Y(M,N,P)+W1*H111(I,J,K)+W2*H222(I,J,K)+W3*H333(I,J,K) YY 05000
         W1=V1*A(L2+K,P)+V4*A(L2+K,N)+V5*A(L2+K,M)                      YY 05010
         W2=V1*A(L3+K,P)+V4*A(L3+K,N)+V5*A(L3+K,M)                      YY 05020
         W3=V3*A(L2+K,P)+V8*A(L2+K,N)+V9*A(L2+K,M)                      YY 05030
         W4=V3*A(L1+K,P)+V8*A(L1+K,N)+V9*A(L1+K,M)                      YY 05040
         W5=V2*A(L1+K,P)+V6*A(L1+K,N)+V7*A(L1+K,M)                      YY 05050
         W6=V2*A(L3+K,P)+V6*A(L3+K,N)+V7*A(L3+K,M)                      YY 05060
         Y(M,N,P)=Y(M,N,P)+W1*H112(I,J,K)+W2*H113(I,J,K)+W3*H332(I,J,K) YY 05070
         Y(M,N,P)=Y(M,N,P)+W4*H331(I,J,K)+W5*H221(I,J,K)+W6*H223(I,J,K) YY 05080
         W1=A(L3+K,P)*(V10+V12)+A(L3+K,N)*(V11+V14)+A(L3+K,M)*(V13+V15) YY 05090
 14      Y(M,N,P)=Y(M,N,P)+W1*H123(I,J,K)                               YY 05100
 12   CONTINUE                                                          YY 05110
 10   CONTINUE                                                          YY 05120
      RETURN                                                            YY 05130
      END                                                               YY 05140
C     ///////////////////////////////////////////////////////////////// YY 05150
      SUBROUTINE AHY4(NC,NSX,L1,L2,L3,L4,H111,H112,H221,H222,H113,H123, YY 05160
     $      H223,H331,H332,H333,H411,H421,H422,H431,H432,H433,H441,     YY 05170
     $      H442,H443,H444,A,Y)                                         YY 05180
      IMPLICIT REAL*8 (A-H,O-Z)                                         YY 05190
      INTEGER P,NSX                                                     YY 05200
      DIMENSION A(NC,NC),Y(NC,NC,NC)                                    YY 05210
      DIMENSION H111(3,3,3),H112(3,3,3),H221(3,3,3),H222(3,3,3)         YY 05220
      DIMENSION H113(3,3,3),H123(3,3,3),H223(3,3,3),H331(3,3,3)         YY 05230
      DIMENSION H332(3,3,3),H333(3,3,3),H411(3,3,3),H421(3,3,3)         YY 05240
      DIMENSION H422(3,3,3),H431(3,3,3),H432(3,3,3),H433(3,3,3)         YY 05250
      DIMENSION H441(3,3,3),H442(3,3,3),H443(3,3,3),H444(3,3,3)         YY 05260
      DO 10  P=1,NSX                                                    YY 05270
      DO 10  N=1,P                                                      YY 05280
      DO 10  M=1,N                                                      YY 05290
         DO 12  I=1,3                                                   YY 05300
         DO 12  J=1,3                                                   YY 05310
         V1=A(L1+I,M)*A(L1+J,N)                                         YY 05320
         V2=A(L2+I,M)*A(L2+J,N)                                         YY 05330
         V3=A(L3+I,M)*A(L3+J,N)                                         YY 05340
         V4=A(L4+I,M)*A(L4+J,N)                                         YY 05350
         V5=A(L1+I,M)*A(L1+J,P)                                         YY 05360
         V6=A(L1+I,N)*A(L1+J,P)                                         YY 05370
         V7=A(L2+I,M)*A(L2+J,P)                                         YY 05380
         V8=A(L2+I,N)*A(L2+J,P)                                         YY 05390
         V9=A(L3+I,M)*A(L3+J,P)                                         YY 05400
         V10=A(L3+I,N)*A(L3+J,P)                                        YY 05410
         V11=A(L4+I,M)*A(L4+J,P)                                        YY 05420
         V12=A(L4+I,N)*A(L4+J,P)                                        YY 05430
         V13=A(L1+I,M)*A(L2+J,N)                                        YY 05440
         V14=A(L1+I,M)*A(L2+J,P)                                        YY 05450
         V15=A(L1+I,N)*A(L2+J,M)                                        YY 05460
         V16=A(L1+I,N)*A(L2+J,P)                                        YY 05470
         V17=A(L1+I,P)*A(L2+J,M)                                        YY 05480
         V18=A(L1+I,P)*A(L2+J,N)                                        YY 05490
         V19=A(L4+I,M)*A(L2+J,N)                                        YY 05500
         V20=A(L4+I,M)*A(L2+J,P)                                        YY 05510
         V21=A(L4+I,N)*A(L2+J,M)                                        YY 05520
         V22=A(L4+I,N)*A(L2+J,P)                                        YY 05530
         V23=A(L4+I,P)*A(L2+J,M)                                        YY 05540
         V24=A(L4+I,P)*A(L2+J,N)                                        YY 05550
         V25=A(L4+I,M)*A(L3+J,N)                                        YY 05560
         V26=A(L4+I,M)*A(L3+J,P)                                        YY 05570
         V27=A(L4+I,N)*A(L3+J,M)                                        YY 05580
         V28=A(L4+I,N)*A(L3+J,P)                                        YY 05590
         V29=A(L4+I,P)*A(L3+J,M)                                        YY 05600
         V30=A(L4+I,P)*A(L3+J,N)                                        YY 05610
C                                                                       YY 05620
         DO 14  K=1,3                                                   YY 05630
         W1=V1*A(L1+K,P)                                                YY 05640
         W2=V2*A(L2+K,P)                                                YY 05650
         W3=V3*A(L3+K,P)                                                YY 05660
         W4=V4*A(L4+K,P)                                                YY 05670
         Y(M,N,P)=Y(M,N,P)+W1*H111(I,J,K)+W2*H222(I,J,K)                YY 05680
         Y(M,N,P)=Y(M,N,P)+W3*H333(I,J,K)+W4*H444(I,J,K)                YY 05690
         W1=V1*A(L2+K,P)+V5*A(L2+K,N)+V6*A(L2+K,M)                      YY 05700
         W2=V1*A(L3+K,P)+V5*A(L3+K,N)+V6*A(L3+K,M)                      YY 05710
         W3=V3*A(L2+K,P)+V9*A(L2+K,N)+V10*A(L2+K,M)                     YY 05720
         W4=V3*A(L1+K,P)+V9*A(L1+K,N)+V10*A(L1+K,M)                     YY 05730
         W5=V2*A(L1+K,P)+V7*A(L1+K,N)+V8*A(L1+K,M)                      YY 05740
         W6=V2*A(L3+K,P)+V7*A(L3+K,N)+V8*A(L3+K,M)                      YY 05750
         W7=V1*A(L4+K,P)+V5*A(L4+K,N)+V6*A(L4+K,M)                      YY 05760
         W8=V2*A(L4+K,P)+V7*A(L4+K,N)+V8*A(L4+K,M)                      YY 05770
         W9=V3*A(L4+K,P)+V9*A(L4+K,N)+V10*A(L4+K,M)                     YY 05780
         W10=V4*A(L1+K,P)+V11*A(L1+K,N)+V12*A(L1+K,M)                   YY 05790
         W11=V4*A(L2+K,P)+V11*A(L2+K,N)+V12*A(L2+K,M)                   YY 05800
         W12=V4*A(L3+K,P)+V11*A(L3+K,N)+V12*A(L3+K,M)                   YY 05810
         Y(M,N,P)=Y(M,N,P)+W1*H112(I,J,K)+W2*H113(I,J,K)+W3*H332(I,J,K) YY 05820
         Y(M,N,P)=Y(M,N,P)+W4*H331(I,J,K)+W5*H221(I,J,K)+W6*H223(I,J,K) YY 05830
         Y(M,N,P)=Y(M,N,P)+W7*H411(K,I,J)+W8*H422(K,I,J)+W9*H433(K,I,J) YY 05840
      Y(M,N,P)=Y(M,N,P)+W10*H441(I,J,K)+W11*H442(I,J,K)+W12*H443(I,J,K) YY 05850
         W1=A(L3+K,P)*(V13+V15)+A(L3+K,N)*(V14+V17)+A(L3+K,M)*(V16+V18) YY 05860
         W2=A(L1+K,P)*(V19+V21)+A(L1+K,N)*(V20+V23)+A(L1+K,M)*(V22+V24) YY 05870
         W3=A(L1+K,P)*(V25+V27)+A(L1+K,N)*(V26+V29)+A(L1+K,M)*(V28+V30) YY 05880
         W4=A(L2+K,P)*(V25+V27)+A(L2+K,N)*(V26+V29)+A(L2+K,M)*(V28+V30) YY 05890
         Y(M,N,P)=Y(M,N,P)+W1*H123(I,J,K)+W2*H421(I,J,K)                YY 05900
 14      Y(M,N,P)=Y(M,N,P)+W3*H431(I,J,K)+W4*H432(I,J,K)                YY 05910
 12   CONTINUE                                                          YY 05920
 10   CONTINUE                                                          YY 05930
      RETURN                                                            YY 05940
      END                                                                      
C     ///////////////////////////////////////////////////////////////   INT32920
      SUBROUTINE FLIN(A,IDIM,IN,IM,DET)                                 INT32930
      IMPLICIT REAL*8 (A-H,O-Z)                                         INT32940
C                                                                       INT32950
C     LINEAR SIMULTANEOUS EQUATION                                      INT32960
C                                                                       INT32970
C     A(IN*IN) * X(IN*IM) = B(IN*IM)                                    INT32980
C                                                                       INT32990
C     A & B SHOULD BE STORED ON A(IN*(IN+IM))                           INT33000
C     SOLUTION X WILL BE STORED ON B PART IN DIMENSION A.               INT33010
C                                                                       INT33020
      DIMENSION A(IDIM,1)                                               INT33030
      PARAMETER (ZERO=0.0D0,ONE=1.0D0)                                  INT33040
C                                                                       INT33050
      N=IN                                                              INT33060
      NR=IM                                                             INT33070
      JMAX=N+NR                                                         INT33080
      SIGN=ONE                                                          INT33090
C M IS THE STAGE OF ELIMINATION                                         INT33100
      DO 49 M=1,N                                                       INT33110
      TEMP=ZERO                                                         INT33120
      DO 41 I=M,N                                                       INT33130
      IF(M.GT.1)A(I,M)=A(I,M)-DOTX(A(I,1),IDIM,A(1,M),1,M-1)            INT33140
      AVAL=A(I,M)                                                       INT33150
      IF(DABS(AVAL).LE.TEMP)GOTO 41                                     INT33160
      TEMP=DABS(AVAL)                                                   INT33170
      IMAX=I                                                            INT33180
 41   CONTINUE                                                          INT33190
      IF(TEMP.LE.ZERO)GOTO 999                                          INT33200
      IF(IMAX.EQ.M)GOTO 45                                              INT33210
      SIGN=-SIGN                                                        INT33220
      DO 44 J=1,JMAX                                                    INT33230
      STOR=A(M,J)                                                       INT33240
      A(M,J)=A(IMAX,J)                                                  INT33250
      A(IMAX,J)=STOR                                                    INT33260
 44   CONTINUE                                                          INT33270
 45   CONTINUE                                                          INT33280
      JJ=M+1                                                            INT33290
      IF(JJ.GT.JMAX)GOTO 49                                             INT33300
      IF(M.GT.1)GOTO 47                                                 INT33310
      DO 46 J=JJ,JMAX                                                   INT33320
      A(1,J)=A(1,J)/A(1,1)                                              INT33330
 46   CONTINUE                                                          INT33340
      D=A(1,1)                                                          INT33350
      GOTO 49                                                           INT33360
 47   CONTINUE                                                          INT33370
      DO 48 J=JJ,JMAX                                                   INT33380
      A(M,J)=(A(M,J)-DOTX(A(M,1),IDIM,A(1,J),1,M-1))/A(M,M)             INT33390
 48   CONTINUE                                                          INT33400
      D=D*A(M,M)                                                        INT33410
 49   CONTINUE                                                          INT33420
      IF(NR.EQ.0) RETURN                                                INT33430
      DO 59 I=1,NR                                                      INT33440
      NPI=N+I                                                           INT33450
      DO 58 K=2,N                                                       INT33460
      J=N+1-K                                                           INT33470
      A(J,NPI)=A(J,NPI)-DOTX(A(J,J+1),IDIM,A(J+1,NPI),1,K-1)            INT33480
 58   CONTINUE                                                          INT33490
 59   CONTINUE                                                          INT33500
      DET=D*SIGN                                                        INT33510
      RETURN                                                            INT33520
C ON ZERO PIVOT, SET DET=0.AND RETURN TO CALLING PROGRAM NOV 1972       INT33530
 999  DET=ZERO                                                          INT33540
      RETURN                                                            INT33550
      END                                                               INT33560
C     //////////////////////////////////////////////////////////////    INT33570
      FUNCTION DOTX(A,NA,B,NB,N)                                        INT33580
      IMPLICIT REAL*8 (A-H,O-Z)                                         INT33590
      DIMENSION A(1),B(1)                                               INT33600
      PARAMETER(ZERO=0.0D0)                                             INT33610
      IAPT=1                                                            INT33620
      IBPT=1                                                            INT33630
      D=ZERO                                                            INT33640
      DO 10 I=1,N                                                       INT33650
      D=D+A(IAPT)*B(IBPT)                                               INT33660
      IAPT=IAPT+NA                                                      INT33670
      IBPT=IBPT+NB                                                      INT33680
 10   CONTINUE                                                          INT33690
      DOTX=D                                                            INT33700
      RETURN                                                            INT33710
      END                                                               INT33720
C     ///////////////////////////////////////////////////////////////   INT33730
      SUBROUTINE VECPRO(U,V,W)                                          INT33740
      REAL*8 U(3),V(3),W(3)                                             INT33750
      W(1)=U(2)*V(3)-V(2)*U(3)                                          INT33760
      W(2)=U(3)*V(1)-U(1)*V(3)                                          INT33770
      W(3)=U(1)*V(2)-U(2)*V(1)                                          INT33780
      RETURN                                                            INT33790
      END                                                               INT33800
C     //////////////////////////////////////////////////////////////    INT33810
      SUBROUTINE SCAPRO(U,V,D)                                          INT33820
      REAL*8 D,U(3),V(3),ZERO                                           INT33830
      PARAMETER(ZERO=0.0D0)                                             INT33840
      D=ZERO                                                            INT33850
      DO 5  I=1,3                                                       INT33860
 5    D=D+U(I)*V(I)                                                     INT33870
      RETURN                                                            INT33880
      END                                                               INT33890
C     //////////////////////////////////////////////////////////////    INT33900
      SUBROUTINE MAT1(EM,V)                                             INT33910
C       EM(I,J)=(EI X EJ)*V                                             INT33920
      IMPLICIT REAL*8 (A-H,O-Z)                                         INT33930
      DIMENSION EM(3,3),V(3)                                            INT33940
      PARAMETER(ZERO=0.0D0)                                             INT33950
      EM(2,1)=-V(3)                                                     INT33960
      EM(3,1)=V(2)                                                      INT33970
      EM(3,2)=-V(1)                                                     INT33980
      EM(1,2)=-EM(2,1)                                                  INT33990
      EM(1,3)=-EM(3,1)                                                  INT34000
      EM(2,3)=-EM(3,2)                                                  INT34010
      DO 5  I=1,3                                                       INT34020
 5    EM(I,I)=ZERO                                                      INT34030
      RETURN                                                            INT34040
      END                                                               INT34050
C     //////////////////////////////////////////////////////////////    INT34060
      SUBROUTINE MAT2(AA,VK)                                            INT34070
      IMPLICIT REAL*8 (A-H,O-Z)                                         INT34080
      DIMENSION AA(3,3),VK(3)                                           INT34090
      PARAMETER(ZERO=0.0D0)                                             INT34100
      DO 10  I=1,3                                                      INT34110
  10  AA(I,I)=ZERO                                                      INT34120
      AA(1,2)=VK(1)                                                     INT34130
      AA(2,1)=-VK(1)                                                    INT34140
      AA(1,3)=VK(2)                                                     INT34150
      AA(3,1)=-VK(2)                                                    INT34160
      AA(2,3)=VK(3)                                                     INT34170
      AA(3,2)=-VK(3)                                                    INT34180
      END                                                               INT34190
C     //////////////////////////////////////////////////////////////    INT34200
      SUBROUTINE TRIPRO(PROD)                                           YY 01550
C        PROD(I,J,K)=(EI X EJ) * EK                                     YY 01560
      IMPLICIT REAL*8 (A-H,O-Z)                                         YY 01570
      DIMENSION PROD(3,3,3),RMAT(3,3),VECT(3)                           YY 01580
      PARAMETER(ZERO=0.0D0,ONE=1.0D0)                                   YY 01590
      DO 10 K=1,3                                                       YY 01600
      DO 5 N=1,3                                                        YY 01610
5     VECT(N)=ZERO                                                      YY 01620
      VECT(K)=ONE                                                       YY 01630
      CALL MAT1(RMAT,VECT)                                              YY 01640
      DO 10 J=1,3                                                       YY 01650
      DO 10 I=1,3                                                       YY 01660
10    PROD(I,J,K)=RMAT(I,J)                                             YY 01670
      RETURN                                                            YY 01680
      END                                                               YY 01690
C     //////////////////////////////////////////////////////////////////DIS03230
      SUBROUTINE EXPMAT(RK,W)                                           DIS03240
      IMPLICIT REAL*8 (A-H,O-Z)                                         DIS03250
      REAL*8 K(3,3)                                                     DIS03260
      DIMENSION RK(3),W(3,3),C(3,3),D(3,3)                              DIS03270
      PARAMETER(ZERO=0.0D0,ONE=1.0D0)                                   DIS03280
      PARAMETER(XTOL=1.0D-20,MAXIT=50)                                  DIS03290
      COMMON /IO/ IIN1,IOUT,IIN2,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25,
     $   I31,I32,I33,I35,I36,I37,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10,
     $   ISCR11,ISCR12,ISCR13,ISCR14
  1   FORMAT(/'CONVERGENCE NOT REACHED IN EXPMAT AFTER',I5,             DIS03340
     $   ' ITERATIONS.'/' CNORM= ',E10.3,' XTOL= ',E10.3)               DIS03350
      DO 10  I=1,3                                                      DIS03360
 10   K(I,I)=ZERO                                                       DIS03370
      K(1,2)=RK(1)                                                      DIS03380
      K(2,1)=-RK(1)                                                     DIS03390
      K(1,3)=RK(2)                                                      DIS03400
      K(3,1)=-RK(2)                                                     DIS03410
      K(2,3)=RK(3)                                                      DIS03420
      K(3,2)=-RK(3)                                                     DIS03430
      DO 20  I=1,3                                                      DIS03440
      DO 30  J=1,3                                                      DIS03450
      C(I,J)=ZERO                                                       DIS03460
 30   W(I,J)=ZERO                                                       DIS03470
      C(I,I)=ONE                                                        DIS03480
 20   W(I,I)=ONE                                                        DIS03490
      ITER=0                                                            DIS03500
C                                                                       DIS03510
 100  CNORM=ZERO                                                        DIS03520
      ITER=ITER+1
      DO 50  I=1,3                                                      DIS03530
      DO 50  J=1,3                                                      DIS03540
      D(I,J)=ZERO                                                       DIS03550
      DO 60  L=1,3                                                      DIS03560
 60   D(I,J)=D(I,J)+C(I,L)*K(L,J)                                       DIS03570
      D(I,J)=D(I,J)/DFLOAT(ITER)                                        DIS03580
      IF(DABS(D(I,J)).GT.CNORM) CNORM=DABS(D(I,J))                      DIS03590
 50   CONTINUE                                                          DIS03600
      DO 70  I=1,3                                                      DIS03610
      DO 70  J=1,3                                                      DIS03620
      W(I,J)=W(I,J)+D(I,J)                                              DIS03630
 70   C(I,J)=D(I,J)                                                     DIS03640
      IF(ITER.LT.MAXIT.AND.CNORM.GT.XTOL) GO TO 100                     DIS03650
      IF(CNORM.GT.XTOL) THEN                                            DIS03660
      WRITE(IOUT,1) ITER,CNORM,XTOL                                     DIS03670
      END IF                                                            DIS03680
      RETURN                                                            DIS03690
      END                                                               DIS03700
C////////////////////////////////////////////////////////////////////// INT51150
      SUBROUTINE RSP(NM,N,NV,A,W,MATZ,Z,FV1,FV2)                        INT51740
      IMPLICIT REAL*8 (A-H,O-Z)                                         INT51750
      INTEGER I,J,N,NM,NV,IERR,MATZ                                     INT51760
      DIMENSION A(NV),W(N),Z(NM,N),FV1(N),FV2(N)                        INT51770
      DATA ZERO,ONE / 0.0D+00 , 1.0D+00 /                               INT51780
    1 FORMAT(//,2X,' IERR = ',I5//)                                     INT51790
C     THIS SUBROUTINE CALLS THE RECOMMENDED SEQUENCE OF                 INT51800
C     SUBROUTINES FROM THE EIGENSYSTEM SUBROUTINE PACKAGE (EISPACK)     INT51810
C     TO FIND THE EIGENVALUES AND EIGENVECTORS (IF DESIRED)             INT51820
C     OF A REAL SYMMETRIC PACKED MATRIX.                                INT51830
C     ON INPUT-                                                         INT51840
C        NM  MUST BE SET TO THE ROW DIMENSION OF THE TWO-DIMENSIONAL    INT51850
C        ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM            INT51860
C        DIMENSION STATEMENT,                                           INT51870
C        N  IS THE ORDER OF THE MATRIX  A,                              INT51880
C        NV  IS AN INTEGER VARIABLE SET EQUAL TO THE                    INT51890
C        DIMENSION OF THE ARRAY  A  AS SPECIFIED FOR                    INT51900
C        A  IN THE CALLING PROGRAM.  NV  MUST NOT BE                    INT51910
C        LESS THAN  N*(N+1)/2,                                          INT51920
C        A  CONTAINS THE LOWER TRIANGLE OF THE REAL SYMMETRIC           INT51930
C        PACKED MATRIX STORED ROW-WISE,                                 INT51940
C        MATZ  IS AN INTEGER VARIABLE SET EQUAL TO ZERO IF              INT51950
C        ONLY EIGENVALUES ARE DESIRED,  OTHERWISE IT IS SET TO          INT51960
C        ANY NON-ZERO INTEGER FOR BOTH EIGENVALUES AND EIGENVECTORS.    INT51970
C     ON OUTPUT-                                                        INT51980
C        W  CONTAINS THE EIGENVALUES IN ASCENDING ORDER,                INT51990
C        Z  CONTAINS THE EIGENVECTORS IF MATZ IS NOT ZERO,              INT52000
C        IERR  IS AN INTEGER OUTPUT VARIABLE SET EQUAL TO AN            INT52010
C        ERROR COMPLETION CODE DESCRIBED IN SECTION 2B OF THE           INT52020
C        DOCUMENTATION.  THE NORMAL COMPLETION CODE IS ZERO,            INT52030
C        FV1  AND  FV2  ARE TEMPORARY STORAGE ARRAYS.                   INT52040
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,        INT52050
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY         INT52060
C     ------------------------------------------------------------------INT52070
      IF (N .LE. NM) GO TO 5                                            INT52080
      IERR = 10 * N                                                     INT52090
      GO TO 50                                                          INT52100
    5 IF (NV .GE. (N * (N + 1)) / 2) GO TO 10                           INT52110
      IERR = 20 * N                                                     INT52120
      GO TO 50                                                          INT52130
   10 CALL  TRED3(N,NV,A,W,FV1,FV2)                                     INT52140
      IF (MATZ .NE. 0) GO TO 20                                         INT52150
C     ********** FIND EIGENVALUES ONLY **********                       INT52160
      CALL  TQLRAT(N,W,FV2,IERR)                                        INT52170
      GO TO 50                                                          INT52180
C     ********** FIND BOTH EIGENVALUES AND EIGENVECTORS **********      INT52190
   20 DO 40 I = 1, N                                                    INT52200
         DO 30 J = 1, N                                                 INT52210
            Z(J,I) = ZERO                                               INT52220
   30    CONTINUE                                                       INT52230
         Z(I,I) = ONE                                                   INT52240
   40 CONTINUE                                                          INT52250
      CALL  TQL2(NM,N,W,FV1,Z,IERR)                                     INT52260
      IF (IERR .NE. 0) GO TO 50                                         INT52270
      CALL  TRBAK3(NM,N,NV,A,N,Z)                                       INT52280
   50 IF(IERR.NE.0) GO TO 60                                            INT52290
      RETURN                                                            INT52300
   60 WRITE(6,1) IERR                                                   INT52310
      RETURN                                                            INT52320
C     ********** LAST CARD OF RSP **********                            INT52330
      END                                                               INT52340
C///////////////////////////////////////////////////////////////////////INT52350
      SUBROUTINE TQLRAT(N,D,E2,IERR)                                    INT52360
      IMPLICIT REAL*8 (A-H,O-Z)                                         INT52370
      INTEGER I,J,L,M,N,II,L1,MML,IERR                                  INT52380
      REAL*8 MACHEP                                                     INT52390
      DIMENSION D(N),E2(N)                                              INT52400
      DATA ZERO,ONE,TWO / 0.0D+00 , 1.0D+00 , 2.0D+00 /                 INT52410
C     REAL SQRT,ABS,SIGN                                                INT52420
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TQLRAT,   INT52430
C     ALGORITHM 464, COMM. ACM 16, 689(1973) BY REINSCH.                INT52440
C     THIS SUBROUTINE FINDS THE EIGENVALUES OF A SYMMETRIC              INT52450
C     TRIDIAGONAL MATRIX BY THE RATIONAL QL METHOD.                     INT52460
C     ON INPUT-                                                         INT52470
C        N IS THE ORDER OF THE MATRIX,                                  INT52480
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX,          INT52490
C        E2 CONTAINS THE SQUARES OF THE SUBDIAGONAL ELEMENTS OF THE     INT52500
C      ON OUTPUT-                                                       INT52510
C        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN          INT52520
C          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT AND          INT52530
C          ORDERED FOR INDICES 1,2,...IERR-1, BUT MAY NOT BE            INT52540
C          THE SMALLEST EIGENVALUES,                                    INT52550
C        E2 HAS BEEN DESTROYED,                                         INT52560
C        IERR IS SET TO                                                 INT52570
C          ZERO       FOR NORMAL RETURN,                                INT52580
C          J          IF THE J-TH EIGENVALUE HAS NOT BEEN               INT52590
C                     DETERMINED AFTER 30 ITERATIONS.                   INT52600
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,        INT52610
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY         INT52620
C     ------------------------------------------------------------------INT52630
C     ********** MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING     INT52640
C                THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.   INT52650
C                **********                                             INT52660
      MACHEP = 2.**(-47)                                                INT52670
      IERR = 0                                                          INT52680
      IF (N .EQ. 1) GO TO 1001                                          INT52690
      DO 100 I = 2, N                                                   INT52700
  100 E2(I-1) = E2(I)                                                   INT52710
      F = ZERO                                                          INT52720
      B = ZERO                                                          INT52730
      E2(N) = ZERO                                                      INT52740
      DO 290 L = 1, N                                                   INT52750
         J = 0                                                          INT52760
         H = MACHEP * (DABS(D(L)) + DSQRT(E2(L)))                       INT52770
         IF (B .GT. H) GO TO 105                                        INT52780
         B = H                                                          INT52790
         C = B * B                                                      INT52800
C     ********** LOOK FOR SMALL SQUARED SUB-DIAGONAL ELEMENT ********** INT52810
  105    DO 110 M = L, N                                                INT52820
            IF (E2(M) .LE. C) GO TO 120                                 INT52830
C     ********** E2(N) IS ALWAYS ZERO, SO THERE IS NO EXIT              INT52840
C                THROUGH THE BOTTOM OF THE LOOP **********              INT52850
  110    CONTINUE                                                       INT52860
  120    IF (M .EQ. L) GO TO 210                                        INT52870
  130    IF (J .EQ. 30) GO TO 1000                                      INT52880
         J = J + 1                                                      INT52890
C     ********** FORM SHIFT **********                                  INT52900
         L1 = L + 1                                                     INT52910
         S = DSQRT(E2(L))                                               INT52920
         G = D(L)                                                       INT52930
         P = (D(L1) - G) / (TWO * S)                                    INT52940
         R = DSQRT(P*P+ONE)                                             INT52950
         D(L) = S / (P + DSIGN(R,P))                                    INT52960
         H = G - D(L)                                                   INT52970
         DO 140 I = L1, N                                               INT52980
  140    D(I) = D(I) - H                                                INT52990
         F = F + H                                                      INT53000
C     ********** RATIONAL QL TRANSFORMATION **********                  INT53010
         G = D(M)                                                       INT53020
         IF (G .EQ. ZERO) G = B                                         INT53030
         H = G                                                          INT53040
         S = ZERO                                                       INT53050
         MML = M - L                                                    INT53060
C     ********** FOR I=M-1 STEP -1 UNTIL L DO -- **********             INT53070
         DO 200 II = 1, MML                                             INT53080
            I = M - II                                                  INT53090
            P = G * H                                                   INT53100
            R = P + E2(I)                                               INT53110
            E2(I+1) = S * R                                             INT53120
            S = E2(I) / R                                               INT53130
            D(I+1) = H + S * (H + D(I))                                 INT53140
            G = D(I) - E2(I) / G                                        INT53150
            IF (G .EQ. ZERO) G = B                                      INT53160
            H = G * P / R                                               INT53170
  200    CONTINUE                                                       INT53180
         E2(L) = S * G                                                  INT53190
         D(L) = H                                                       INT53200
C     ********** GUARD AGAINST UNDERFLOW IN CONVERGENCE TEST ********** INT53210
         IF (H .EQ. ZERO) GO TO 210                                     INT53220
         IF (DABS(E2(L)) .LE. DABS(C/H)) GO TO 210                      INT53230
         E2(L) = H * E2(L)                                              INT53240
         IF (E2(L) .NE. ZERO) GO TO 130                                 INT53250
  210    P = D(L) + F                                                   INT53260
C     ********** ORDER EIGENVALUES **********                           INT53270
         IF (L .EQ. 1) GO TO 250                                        INT53280
C     ********** FOR I=L STEP -1 UNTIL 2 DO -- **********               INT53290
         DO 230 II = 2, L                                               INT53300
            I = L + 2 - II                                              INT53310
            IF (P .GE. D(I-1)) GO TO 270                                INT53320
            D(I) = D(I-1)                                               INT53330
  230    CONTINUE                                                       INT53340
  250    I = 1                                                          INT53350
  270    D(I) = P                                                       INT53360
  290 CONTINUE                                                          INT53370
      GO TO 1001                                                        INT53380
C     ********** SET ERROR -- NO CONVERGENCE TO AN                      INT53390
C                EIGENVALUE AFTER 30 ITERATIONS **********              INT53400
 1000 IERR = L                                                          INT53410
 1001 RETURN                                                            INT53420
C     ********** LAST CARD OF TQLRAT **********                         INT53430
      END                                                               INT53440
C///////////////////////////////////////////////////////////////////////INT53450
      SUBROUTINE TQL2(NM,N,D,E,Z,IERR)                                  INT53460
      IMPLICIT REAL*8 (A-H,O-Z)                                         INT53470
      INTEGER I,J,K,L,M,N,II,L1,NM,MML,IERR                             INT53480
      REAL*8 MACHEP                                                     INT53490
      DIMENSION D(N),E(N),Z(NM,N)                                       INT53500
      DATA ZERO,ONE,TWO / 0.0D+00 , 1.0D+00 , 2.0D+00 /                 INT53510
C     REAL SQRT,ABS,SIGN                                                INT53520
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TQL2,     INT53530
C     NUM. MATH. 11, 293-306(1968) BY BOWDLER, MARTIN, REINSCH, AND     INT53540
C     WILKINSON.                                                        INT53550
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 227-240(1971).   INT53560
C     THIS SUBROUTINE FINDS THE EIGENVALUES AND EIGENVECTORS            INT53570
C     OF A SYMMETRIC TRIDIAGONAL MATRIX BY THE QL METHOD.               INT53580
C     THE EIGENVECTORS OF A FULL SYMMETRIC MATRIX CAN ALSO              INT53590
C     BE FOUND IF  TRED2  HAS BEEN USED TO REDUCE THIS                  INT53600
C     FULL MATRIX TO TRIDIAGONAL FORM.                                  INT53610
C     ON INPUT-                                                         INT53620
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL         INT53630
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM          INT53640
C          DIMENSION STATEMENT,                                         INT53650
C        N IS THE ORDER OF THE MATRIX,                                  INT53660
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX,          INT53670
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX        INT53680
C          IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY,               INT53690
C        Z CONTAINS THE TRANSFORMATION MATRIX PRODUCED IN THE           INT53700
C          REDUCTION BY  TRED2, IF PERFORMED.  IF THE EIGENVECTORS      INT53710
C          OF THE TRIDIAGONAL MATRIX ARE DESIRED, Z MUST CONTAIN        INT53720
C          THE IDENTITY MATRIX.                                         INT53730
C      ON OUTPUT-                                                       INT53740
C        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN          INT53750
C          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT BUT          INT53760
C          UNORDERED FOR INDICES 1,2,...,IERR-1,                        INT53770
C        E HAS BEEN DESTROYED,                                          INT53780
C        Z CONTAINS ORTHONORMAL EIGENVECTORS OF THE SYMMETRIC           INT53790
C          TRIDIAGONAL (OR FULL) MATRIX.  IF AN ERROR EXIT IS MADE,     INT53800
C          Z CONTAINS THE EIGENVECTORS ASSOCIATED WITH THE STORED       INT53810
C          EIGENVALUES,                                                 INT53820
C        IERR IS SET TO                                                 INT53830
C          ZERO       FOR NORMAL RETURN,                                INT53840
C          J          IF THE J-TH EIGENVALUE HAS NOT BEEN               INT53850
C                     DETERMINED AFTER 30 ITERATIONS.                   INT53860
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,        INT53870
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY         INT53880
C     ------------------------------------------------------------------INT53890
C     ********** MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING     INT53900
C                THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.   INT53910
C                **********                                             INT53920
      MACHEP = 2.**(-47)                                                INT53930
      IERR = 0                                                          INT53940
      IF (N .EQ. 1) GO TO 1001                                          INT53950
      DO 100 I = 2, N                                                   INT53960
  100 E(I-1) = E(I)                                                     INT53970
      F = ZERO                                                          INT53980
      B = ZERO                                                          INT53990
      E(N) = ZERO                                                       INT54000
      DO 240 L = 1, N                                                   INT54010
         J = 0                                                          INT54020
         H = MACHEP * (DABS(D(L)) + DABS(E(L)))                         INT54030
         IF (B .LT. H) B = H                                            INT54040
C     ********** LOOK FOR SMALL SUB-DIAGONAL ELEMENT **********         INT54050
         DO 110 M = L, N                                                INT54060
            IF (DABS(E(M)) .LE. B) GO TO 120                            INT54070
C     ********** E(N) IS ALWAYS ZERO, SO THERE IS NO EXIT               INT54080
C                THROUGH THE BOTTOM OF THE LOOP **********              INT54090
  110    CONTINUE                                                       INT54100
  120    IF (M .EQ. L) GO TO 220                                        INT54110
  130    IF (J .EQ. 30) GO TO 1000                                      INT54120
         J = J + 1                                                      INT54130
C     ********** FORM SHIFT **********                                  INT54140
         L1 = L + 1                                                     INT54150
         G = D(L)                                                       INT54160
         P = (D(L1) - G) / (TWO * E(L))                                 INT54170
         R = DSQRT(P*P+ONE)                                             INT54180
         D(L) = E(L) / (P + DSIGN(R,P))                                 INT54190
         H = G - D(L)                                                   INT54200
         DO 140 I = L1, N                                               INT54210
  140    D(I) = D(I) - H                                                INT54220
         F = F + H                                                      INT54230
C     ********** QL TRANSFORMATION **********                           INT54240
         P = D(M)                                                       INT54250
         C = ONE                                                        INT54260
         S = ZERO                                                       INT54270
         MML = M - L                                                    INT54280
C     ********** FOR I=M-1 STEP -1 UNTIL L DO -- **********             INT54290
         DO 200 II = 1, MML                                             INT54300
            I = M - II                                                  INT54310
            G = C * E(I)                                                INT54320
            H = C * P                                                   INT54330
            IF (DABS(P) .LT. DABS(E(I))) GO TO 150                      INT54340
            C = E(I) / P                                                INT54350
            R = DSQRT(C*C+ONE)                                          INT54360
            E(I+1) = S * P * R                                          INT54370
            S = C / R                                                   INT54380
            C = ONE / R                                                 INT54390
            GO TO 160                                                   INT54400
  150       C = P / E(I)                                                INT54410
            R = DSQRT(C*C+ONE)                                          INT54420
            E(I+1) = S * E(I) * R                                       INT54430
            S = ONE / R                                                 INT54440
            C = C * S                                                   INT54450
  160       P = C * D(I) - S * G                                        INT54460
            D(I+1) = H + S * (C * G + S * D(I))                         INT54470
C     ********** FORM VECTOR **********                                 INT54480
            DO 180 K = 1, N                                             INT54490
               H = Z(K,I+1)                                             INT54500
               Z(K,I+1) = S * Z(K,I) + C * H                            INT54510
               Z(K,I) = C * Z(K,I) - S * H                              INT54520
  180       CONTINUE                                                    INT54530
  200    CONTINUE                                                       INT54540
         E(L) = S * P                                                   INT54550
         D(L) = C * P                                                   INT54560
         IF (DABS(E(L)) .GT. B) GO TO 130                               INT54570
  220    D(L) = D(L) + F                                                INT54580
  240 CONTINUE                                                          INT54590
C     ********** ORDER EIGENVALUES AND EIGENVECTORS **********          INT54600
      DO 300 II = 2, N                                                  INT54610
         I = II - 1                                                     INT54620
         K = I                                                          INT54630
         P = D(I)                                                       INT54640
         DO 260 J = II, N                                               INT54650
            IF (D(J) .GE. P) GO TO 260                                  INT54660
            K = J                                                       INT54670
            P = D(J)                                                    INT54680
  260    CONTINUE                                                       INT54690
         IF (K .EQ. I) GO TO 300                                        INT54700
         D(K) = D(I)                                                    INT54710
         D(I) = P                                                       INT54720
         DO 280 J = 1, N                                                INT54730
            P = Z(J,I)                                                  INT54740
            Z(J,I) = Z(J,K)                                             INT54750
            Z(J,K) = P                                                  INT54760
  280    CONTINUE                                                       INT54770
  300 CONTINUE                                                          INT54780
      GO TO 1001                                                        INT54790
C     ********** SET ERROR -- NO CONVERGENCE TO AN                      INT54800
C                EIGENVALUE AFTER 30 ITERATIONS **********              INT54810
 1000 IERR = L                                                          INT54820
 1001 RETURN                                                            INT54830
C     ********** LAST CARD OF TQL2 **********                           INT54840
      END                                                               INT54850
C///////////////////////////////////////////////////////////////////////INT54860
      SUBROUTINE TRBAK3(NM,N,NV,A,M,Z)                                  INT54870
      IMPLICIT REAL*8 (A-H,O-Z)                                         INT54880
      INTEGER I,J,K,L,M,N,IK,IZ,NM,NV                                   INT54890
      DIMENSION A(NV),Z(NM,N)                                           INT54900
      DATA ZERO / 0.0D+00 /                                             INT54910
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TRBAK3,   INT54920
C     NUM. MATH. 11, 181-195(1968) BY MARTIN, REINSCH, AND WILKINSON.   INT54930
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).   INT54940
C     THIS SUBROUTINE FORMS THE EIGENVECTORS OF A REAL SYMMETRIC        INT54950
C     MATRIX BY BACK TRANSFORMING THOSE OF THE CORRESPONDING            INT54960
C     SYMMETRIC TRIDIAGONAL MATRIX DETERMINED BY  TRED3.                INT54970
C     ON INPUT-                                                         INT54980
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL         INT54990
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM          INT55000
C          DIMENSION STATEMENT,                                         INT55010
C        N IS THE ORDER OF THE MATRIX,                                  INT55020
C        NV MUST BE SET TO THE DIMENSION OF THE ARRAY PARAMETER A       INT55030
C          AS DECLARED IN THE CALLING PROGRAM DIMENSION STATEMENT,      INT55040
C        A CONTAINS INFORMATION ABOUT THE ORTHOGONAL TRANSFORMATIONS    INT55050
C          USED IN THE REDUCTION BY  TRED3  IN ITS FIRST                INT55060
C          N*(N+1)/2 POSITIONS,                                         INT55070
C        M IS THE NUMBER OF EIGENVECTORS TO BE BACK TRANSFORMED,        INT55080
C        Z CONTAINS THE EIGENVECTORS TO BE BACK TRANSFORMED             INT55090
C          IN ITS FIRST M COLUMNS.                                      INT55100
C     ON OUTPUT-                                                        INT55110
C        Z CONTAINS THE TRANSFORMED EIGENVECTORS                        INT55120
C          IN ITS FIRST M COLUMNS.                                      INT55130
C     NOTE THAT TRBAK3 PRESERVES VECTOR EUCLIDEAN NORMS.                INT55140
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,        INT55150
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY         INT55160
C     ------------------------------------------------------------------INT55170
      IF (M .EQ. 0) GO TO 200                                           INT55180
      IF (N .EQ. 1) GO TO 200                                           INT55190
      DO 140 I = 2, N                                                   INT55200
         L = I - 1                                                      INT55210
         IZ = (I * L) / 2                                               INT55220
         IK = IZ + I                                                    INT55230
         H = A(IK)                                                      INT55240
         IF (H .EQ. ZERO) GO TO 140                                     INT55250
         DO 130 J = 1, M                                                INT55260
            S = ZERO                                                    INT55270
            IK = IZ                                                     INT55280
            DO 110 K = 1, L                                             INT55290
               IK = IK + 1                                              INT55300
               S = S + A(IK) * Z(K,J)                                   INT55310
  110       CONTINUE                                                    INT55320
C     ********** DOUBLE DIVISION AVOIDS POSSIBLE UNDERFLOW **********   INT55330
            S = (S / H) / H                                             INT55340
            IK = IZ                                                     INT55350
            DO 120 K = 1, L                                             INT55360
               IK = IK + 1                                              INT55370
               Z(K,J) = Z(K,J) - S * A(IK)                              INT55380
  120       CONTINUE                                                    INT55390
  130    CONTINUE                                                       INT55400
  140 CONTINUE                                                          INT55410
  200 RETURN                                                            INT55420
C     ********** LAST CARD OF TRBAK3 **********                         INT55430
      END                                                               INT55440
C///////////////////////////////////////////////////////////////////////INT55450
      SUBROUTINE TRED3(N,NV,A,D,E,E2)                                   INT55460
      IMPLICIT REAL*8 (A-H,O-Z)                                         INT55470
      INTEGER I,J,K,L,N,II,IZ,JK,NV                                     INT55480
      DIMENSION A(NV),D(N),E(N),E2(N)                                   INT55490
      DATA ZERO / 0.0D+00 /                                             INT55500
C     REAL SQRT,ABS,SIGN                                                INT55510
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TRED3,    INT55520
C     NUM. MATH. 11, 181-195(1968) BY MARTIN, REINSCH, AND WILKINSON.   INT55530
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).   INT55540
C     THIS SUBROUTINE REDUCES A REAL SYMMETRIC MATRIX, STORED AS        INT55550
C     A ONE-DIMENSIONAL ARRAY, TO A SYMMETRIC TRIDIAGONAL MATRIX        INT55560
C     USING ORTHOGONAL SIMILARITY TRANSFORMATIONS.                      INT55570
C     ON INPUT-                                                         INT55580
C        N IS THE ORDER OF THE MATRIX,                                  INT55590
C        NV MUST BE SET TO THE DIMENSION OF THE ARRAY PARAMETER A       INT55600
C          AS DECLARED IN THE CALLING PROGRAM DIMENSION STATEMENT,      INT55610
C        A CONTAINS THE LOWER TRIANGLE OF THE REAL SYMMETRIC            INT55620
C          INPUT MATRIX, STORED ROW-WISE AS A ONE-DIMENSIONAL           INT55630
C          ARRAY, IN ITS FIRST N*(N+1)/2 POSITIONS.                     INT55640
C     ON OUTPUT-                                                        INT55650
C        A CONTAINS INFORMATION ABOUT THE ORTHOGONAL                    INT55660
C          TRANSFORMATIONS USED IN THE REDUCTION,                       INT55670
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE TRIDIAGONAL MATRIX,    INT55680
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE TRIDIAGONAL         INT55690
C          MATRIX IN ITS LAST N-1 POSITIONS.  E(1) IS SET TO ZERO,      INT55700
C        E2 CONTAINS THE SQUARES OF THE CORRESPONDING ELEMENTS OF E.    INT55710
C          E2 MAY COINCIDE WITH E IF THE SQUARES ARE NOT NEEDED.        INT55720
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,        INT55730
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY         INT55740
C     ------------------------------------------------------------------INT55750
C     ********** FOR I=N STEP -1 UNTIL 1 DO -- **********               INT55760
      DO  300 II = 1, N                                                 INT55770
         I = N + 1 - II                                                 INT55780
         L = I - 1                                                      INT55790
         IZ = (I * L) / 2                                               INT55800
         H = ZERO                                                       INT55810
         SCALE = ZERO                                                   INT55820
         IF (L .LT. 1) GO TO 130                                        INT55830
C     ********** SCALE ROW (ALGOL TOL THEN NOT NEEDED) **********       INT55840
         DO 120 K = 1, L                                                INT55850
            IZ = IZ + 1                                                 INT55860
            D(K) = A(IZ)                                                INT55870
            SCALE = SCALE + DABS(D(K))                                  INT55880
  120    CONTINUE                                                       INT55890
         IF (SCALE .NE. ZERO) GO TO 140                                 INT55900
  130    E(I) = ZERO                                                    INT55910
         E2(I) = ZERO                                                   INT55920
         GO TO 290                                                      INT55930
  140    DO 150 K = 1, L                                                INT55940
            D(K) = D(K) / SCALE                                         INT55950
            H = H + D(K) * D(K)                                         INT55960
  150    CONTINUE                                                       INT55970
         E2(I) = SCALE * SCALE * H                                      INT55980
         F = D(L)                                                       INT55990
         G = -DSIGN(DSQRT(H),F)                                         INT56000
         E(I) = SCALE * G                                               INT56010
         H = H - F * G                                                  INT56020
         D(L) = F - G                                                   INT56030
         A(IZ) = SCALE * D(L)                                           INT56040
         IF (L .EQ. 1) GO TO 290                                        INT56050
         F = ZERO                                                       INT56060
         DO 240 J = 1, L                                                INT56070
            G = ZERO                                                    INT56080
            JK = (J * (J-1)) / 2                                        INT56090
C     ********** FORM ELEMENT OF A*U **********                         INT56100
            DO 180 K = 1, L                                             INT56110
               JK = JK + 1                                              INT56120
               IF (K .GT. J) JK = JK + K - 2                            INT56130
               G = G + A(JK) * D(K)                                     INT56140
  180       CONTINUE                                                    INT56150
C     ********** FORM ELEMENT OF P **********                           INT56160
            E(J) = G / H                                                INT56170
            F = F + E(J) * D(J)                                         INT56180
  240    CONTINUE                                                       INT56190
         HH = F / (H + H)                                               INT56200
         JK = 0                                                         INT56210
C     ********** FORM REDUCED A **********                              INT56220
         DO 260 J = 1, L                                                INT56230
            F = D(J)                                                    INT56240
            G = E(J) - HH * F                                           INT56250
            E(J) = G                                                    INT56260
            DO 260 K = 1, J                                             INT56270
               JK = JK + 1                                              INT56280
               A(JK) = A(JK) - F * E(K) - G * D(K)                      INT56290
  260    CONTINUE                                                       INT56300
  290    D(I) = A(IZ+1)                                                 INT56310
         A(IZ+1) = SCALE * DSQRT(H)                                     INT56320
  300 CONTINUE                                                          INT56330
      RETURN                                                            INT56340
C     ********** LAST CARD OF TRED3 **********                          INT56350
      END                                                               INT56360
C//////////////////////////////////////////////////////////////////////
      SUBROUTINE NUMX(NAD,NC,NS,NSX,XA,XMASS,TYPE,IA,RB,A,S,X,SR)       INT38510
      IMPLICIT REAL*8 (A-H,O-Z)                                         INT38520
      CHARACTER TYPE*5                                                  INT38530
      INTEGER R                                                         INT38540
      DIMENSION XA(NAD,3),XMASS(1),TYPE(NS),IA(NS,6)                    INT38550
      DIMENSION RB(NC),A(NC,NC),X(NC,NC),S(NS),SR(NC,NC)                INT38560
      COMMON /IO/ IIN1,IOUT,IIN2,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25,
     $   I31,I32,I33,I35,I36,I37,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10,
     $   ISCR11,ISCR12,ISCR13,ISCR14
      PARAMETER(ZERO=0.0D0,ONE=1.0D0,DELTA=1.0D-4)                      INT38610
      TDELTA=DELTA+DELTA                                                INT38620
      DO 10  R=1,NS                                                     INT38630
              DO 700  I=1,NC                                            INT38640
              DO 700  J=1,NC                                            INT38650
              X(I,J)=ZERO                                               INT38660
 700          SR(I,J)=ZERO                                              INT38670
         IF(TYPE(R).EQ.'  SPF') THEN
              CALL VECT1(NAD,IA(R,1),IA(R,2),RB,XA,T21)
              S(R)=(ONE-S(R))*T21
         END IF
         IF(IA(R,6).EQ.0) GO TO 30                                      INT38720
         DO 25  I=1,NC                                                  INT38730
         IP=(I-1)/3+1                                                   INT38740
         JP=I-3*(IP-1)                                                  INT38750
         XA(IP,JP)=XA(IP,JP)+DELTA                                      INT38760
         CALL BROW(NAD,NC,NS,XA,XMASS,TYPE,IA,RB,S,R)                   INT38770
         DO 27  J=1,NC                                                  INT38780
 27           SR(I,J)=RB(J)                                             INT38790
         XA(IP,JP)=XA(IP,JP)-TDELTA                                     INT38800
         CALL BROW(NAD,NC,NS,XA,XMASS,TYPE,IA,RB,S,R)                   INT38810
         DO 28  J=1,NC                                                  INT38820
 28           SR(I,J)=(SR(I,J)-RB(J))/TDELTA                            INT38830
 25      XA(IP,JP)=XA(IP,JP)+DELTA                                      INT38840
 30   CALL XOUT(NC,NC,SR,-R,ISCR6)                                      INT38850
      IF(IABS(IA(R,6)).NE.1) GO TO 100                                  INT38860
      DO 40  N=1,NSX                                                    INT38870
      DO 40  I=1,NC                                                     INT38880
      DO 40  J=1,NC                                                     INT38890
 40   X(I,N)=X(I,N)+SR(I,J)*A(J,N)                                      INT38900
      DO 45  N=1,NSX                                                    INT38910
      DO 45  I=1,NC                                                     INT38920
      SR(I,N)=X(I,N)                                                    INT38930
 45   X(I,N)=ZERO                                                       INT38940
      DO 50  M=1,NSX                                                    INT38950
      DO 50  N=1,NSX                                                    INT38960
      DO 50  I=1,NC                                                     INT38970
 50   X(M,N)=X(M,N)+SR(I,N)*A(I,M)                                      INT38980
100   CALL XOUT(NC,NSX,X,R,ISCR6)                                       INT38990
         IF(TYPE(R).EQ.'  SPF') THEN
              S(R)=ONE-S(R)/T21
         END IF
 10   CONTINUE                                                          INT39030
      END                                                               INT39040
C     //////////////////////////////////////////////////////////////    INT39050
      SUBROUTINE NUMY(NAD,NC,NS,NSX,XA,TYPE,IA,A,S,Y,SR,XT)             INT39650
      IMPLICIT REAL*8 (A-H,O-Z)                                         INT39660
      CHARACTER TYPE*5                                                  INT39670
      INTEGER R,P,Q                                                     INT39680
      DIMENSION XA(NAD,3),TYPE(NS),IA(NS,6),XT(NC,NC)                   INT39690
      DIMENSION A(NC,NC),S(NS),SR(NC,NC,NC),Y(NC,NC,NC)                 INT39700
      COMMON /IO/ IIN1,IOUT,IIN2,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25,
     $   I31,I32,I33,I35,I36,I37,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10,
     $   ISCR11,ISCR12,ISCR13,ISCR14
      PARAMETER(ZERO=0.0D0,ONE=1.0D0,DELTA=1.0D-4)                      INT39750
      TDELTA=DELTA+DELTA                                                INT39760
      DO 10  R=1,NS                                                     INT39770
              DO 700  K=1,NC                                            INT39780
              DO 700  J=1,NC                                            INT39790
              DO 700  I=1,NC                                            INT39800
              Y(I,J,K)=ZERO                                             INT39810
 700          SR(I,J,K)=ZERO                                            INT39820
         IF(TYPE(R).EQ.'  SPF') THEN
              CALL VECT1(NAD,IA(R,1),IA(R,2),XT(1,1),XA,T21)
              S(R)=(ONE-S(R))*T21
         END IF
         IF(IA(R,6).EQ.0) GO TO 30                                      INT39870
         DO 25  I=1,NC                                                  INT39880
         IP=(I-1)/3+1                                                   INT39890
         JP=I-3*(IP-1)                                                  INT39900
         XA(IP,JP)=XA(IP,JP)+DELTA                                      INT39910
         CALL XROW(NAD,NC,NS,XA,TYPE,IA,S,XT,R)                         INT39920
         DO 27  K=1,NC                                                  INT39930
         DO 27  J=1,NC                                                  INT39940
 27           SR(I,J,K)=XT(J,K)                                         INT39950
         XA(IP,JP)=XA(IP,JP)-TDELTA                                     INT39960
         CALL XROW(NAD,NC,NS,XA,TYPE,IA,S,XT,R)                         INT39970
         DO 28  K=1,NC                                                  INT39980
         DO 28  J=1,NC                                                  INT39990
 28           SR(I,J,K)=(SR(I,J,K)-XT(J,K))/TDELTA                      INT40000
 25      XA(IP,JP)=XA(IP,JP)+DELTA                                      INT40010
 30   CALL YOUT(NC,NC,SR,-R,ISCR7)                                      INT40020
      IF(IABS(IA(R,6)).NE.2) GO TO 100                                  INT40030
      DO 35  P=1,NSX                                                    INT40040
      DO 35  J=1,NC                                                     INT40050
      DO 35  I=1,NC                                                     INT40060
      DO 35  K=1,NC                                                     INT40070
 35   Y(I,J,P)=Y(I,J,P)+SR(I,J,K)*A(K,P)                                INT40080
      DO 40  P=1,NSX                                                    INT40090
      DO 40  N=1,NSX                                                    INT40100
      DO 40  I=1,NC                                                     INT40110
      SR(I,N,P)=ZERO                                                    INT40120
      DO 40  J=1,NC                                                     INT40130
 40   SR(I,N,P)=SR(I,N,P)+Y(I,J,P)*A(J,N)                               INT40140
      DO 45  P=1,NSX                                                    INT40150
      DO 45  N=1,NSX                                                    INT40160
      DO 45  M=1,NSX                                                    INT40170
      Y(M,N,P)=ZERO                                                     INT40180
      DO 45  I=1,NC                                                     INT40190
 45   Y(M,N,P)=Y(M,N,P)+SR(I,N,P)*A(I,M)                                INT40200
 100  CALL YOUT(NC,NSX,Y,R,ISCR7)                                       INT40210
         IF(TYPE(R).EQ.'  SPF') THEN
              S(R)=ONE-S(R)/T21
         END IF
 10   CONTINUE                                                          INT40250
      END                                                               INT40260
C     /////////////////////////////////////////////////////////////     INT40270
      SUBROUTINE NUMZ(NAD,NC,NS,NSX,IOPT,XA,TYPE,IA,A,S,YR1,YR2)        YY 00010
      IMPLICIT REAL*8 (A-H,O-Z)                                         YY 00020
      CHARACTER TYPE*5                                                  YY 00030
      INTEGER R,P,Q,IOPT(30)                                            YY 00040
      DIMENSION XA(NAD,3),TYPE(NS),S(NS),A(NC,NC),IA(NS,6),V1(3)        YY 00050
      DIMENSION YR1(NC,NC,NC),YR2(NC,NC,NC)                             YY 00060
      PARAMETER(ZERO=0.0D0,ONE=1.0D0,DELTA=1.0D-4)                      YY 00070
      COMMON /IO/ IIN1,IOUT,IIN2,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25,
     $   I31,I32,I33,I35,I36,I37,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10,
     $   ISCR11,ISCR12,ISCR13,ISCR14
      TDELTA=DELTA+DELTA                                                YY 00120
      NSYM=IOPT(3)                                                      YY 00130
      IF(NSYM.NE.0) THEN                                                YY 00140
         ISCR=ISCR10                                                    YY 00150
      ELSE                                                              YY 00160
         ISCR=ISCR9                                                     YY 00170
      END IF                                                            YY 00180
      DO 10  R=1,NS                                                     YY 00200
         IF(TYPE(R).EQ.'  SPF') THEN
              CALL VECT1(NAD,IA(R,1),IA(R,2),V1,XA,T21)
              S(R)=(ONE-S(R))*T21
         END IF
         DO 25  L=1,NC                                                  YY 00280
         IP=(L-1)/3+1                                                   YY 00290
         JP=L-3*(IP-1)                                                  YY 00300
         XA(IP,JP)=XA(IP,JP)+DELTA                                      YY 00310
         CALL YROW(NAD,NC,NS,XA,TYPE,IA,S,YR1,R)                        YY 00320
         XA(IP,JP)=XA(IP,JP)-TDELTA                                     YY 00330
         CALL YROW(NAD,NC,NS,XA,TYPE,IA,S,YR2,R)                        YY 00340
            DO 28  K=1,NC                                               YY 00350
            DO 28  J=1,NC                                               YY 00360
            DO 28  I=1,NC                                               YY 00370
  28        YR1(I,J,K)=(YR1(I,J,K)-YR2(I,J,K))/TDELTA                   YY 00380
         XA(IP,JP)=XA(IP,JP)+DELTA                                      YY 00390
         CALL YOUT2(NC,NC,YR1,-R,L,ISCR)                                YY 00400
  25     CONTINUE                                                       YY 00410
C                                                                       YY 00420
         DO 40  L=1,NC                                                  YY 00450
         CALL YIN2(NC,NC,YR2,-R,L,ISCR)                                 YY 00460
            DO 32  K=1,NC                                               YY 00470
            DO 32  J=1,NC                                               YY 00480
            DO 32  M=1,NSX                                              YY 00490
  32        YR1(M,J,K)=ZERO                                             YY 00500
         DO 42  K=1,NC                                                  YY 00510
         DO 42  J=1,NC                                                  YY 00520
         DO 42  M=1,NSX                                                 YY 00530
         DO 42  I=1,NC                                                  YY 00540
 42      YR1(M,J,K)=YR1(M,J,K)+YR2(I,J,K)*A(I,M)                        YY 00550
 40      CALL YOUT2(NC,NC,YR1,R,L,ISCR)                                 YY 00560
C                                                                       YY 00570
         DO 45  L=1,NC                                                  YY 00580
         CALL YIN2(NC,NC,YR2,R,L,ISCR)                                  YY 00590
            DO 43  K=1,NC                                               YY 00600
            DO 43  N=1,NSX                                              YY 00610
            DO 43  M=1,NSX                                              YY 00620
 43         YR1(M,N,K)=ZERO                                             YY 00630
         DO 47  K=1,NC                                                  YY 00640
         DO 47  N=1,NSX                                                 YY 00650
         DO 47  M=1,NSX                                                 YY 00660
         DO 47  J=1,NC                                                  YY 00670
 47      YR1(M,N,K)=YR1(M,N,K)+YR2(M,J,K)*A(J,N)                        YY 00680
 45      CALL YOUT2(NC,NC,YR1,R,L,ISCR)                                 YY 00690
C                                                                       YY 00700
         DO 55  L=1,NC                                                  YY 00710
         CALL YIN2(NC,NC,YR2,R,L,ISCR)                                  YY 00720
              DO 56  P=1,NSX                                            YY 00730
              DO 56  N=1,NSX                                            YY 00740
              DO 56  M=1,NSX                                            YY 00750
 56           YR1(M,N,P)=ZERO                                           YY 00760
         DO 58  P=1,NSX                                                 YY 00770
         DO 58  N=1,NSX                                                 YY 00780
         DO 58  M=1,NSX                                                 YY 00790
         DO 58  K=1,NC                                                  YY 00800
 58      YR1(M,N,P)=YR1(M,N,P)+YR2(M,N,K)*A(K,P)                        YY 00810
         LR=NS+1                                                        YY 00820
 55      CALL YOUT2(NC,NC,YR1,-LR,L,ISCR)                               YY 00830
C                                                                       YY 00840
         DO 70  Q=1,NSX                                                 YY 00850
            DO 62  P=1,NSX                                              YY 00860
            DO 62  N=1,NSX                                              YY 00870
            DO 62  M=1,NSX                                              YY 00880
 62         YR1(M,N,P)=ZERO                                             YY 00890
         DO 65  L=1,NC                                                  YY 00900
         LR=NS+1                                                        YY 00910
         CALL YIN2(NC,NC,YR2,-LR,L,ISCR)                                YY 00920
            DO 67  P=1,NSX                                              YY 00930
            DO 67  N=1,NSX                                              YY 00940
            DO 67  M=1,NSX                                              YY 00950
 67         YR1(M,N,P)=YR1(M,N,P)+YR2(M,N,P)*A(L,Q)                     YY 00960
 65      CONTINUE                                                       YY 00970
 70      CALL YOUT2(NC,NC,YR1,R,Q,ISCR)                                 YY 00980
C                                                                       YY 00990
         IF(TYPE(R).EQ.'  SPF') THEN
              S(R)=ONE-S(R)/T21
         END IF
 10   CONTINUE                                                          YY 01030
      RETURN                                                            YY 01040
      END                                                               YY 01050
C     ///////////////////////////////////////////////////////////////   YY 07460
      SUBROUTINE SRTST1(NC,NS,NSX,NSYM,IA,SRA,SRN)                      INT44720
      IMPLICIT REAL*8 (A-H,O-Z)                                         INT44730
      INTEGER R                                                         INT44740
      DIMENSION IA(NS,6),SRA(NC,NC),SRN(NC,NC)                          INT44750
      PARAMETER(ZERO=0.0D0,AVCUT=1.0D-6)                                INT44760
      COMMON /IO/ IIN1,IOUT,IIN2,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25,
     $   I31,I32,I33,I35,I36,I37,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10,
     $   ISCR11,ISCR12,ISCR13,ISCR14
    1 FORMAT(//1X,'TESTING OF BR(I,J) AND CR(M,N) MATRICES')            INT44810
    2 FORMAT(/1X,'R=',I3,' BR:  AVRG=',F12.8,'  RMSDIF=',F12.8)         INT44820
    3 FORMAT(//1X,'BR(I,J) MATRIX FOR COORDINATE ',I5/)                 INT44830
    4 FORMAT(6F12.6)                                                    INT44840
    5 FORMAT(//1X,'BR(I,J) ANALYTIC - BR(I,J) NUMERICAL'/)              INT44850
    6 FORMAT(6F12.6)                                                    INT44860
    7 FORMAT(/1X,'R=',I3,' CR:  AVRG=',F12.8,'  RMSDIF=',F12.8)         INT44870
    8 FORMAT(//1X,'CR(M,N) MATRIX FOR COORDINATE ',I5/)                 INT44880
    9 FORMAT(//1X,'CR(M,N) ANALYTIC - CR(M,N) NUMERICAL'/)              INT44890
      IF(NSYM.EQ.0) THEN                                                INT44900
         ISCR=ISCR1                                                     INT44910
      ELSE                                                              INT44920
         ISCR=ISCR2                                                     INT44930
      END IF                                                            INT44940
      WRITE(IOUT,1)                                                     INT44950
      DO 100  R=1,NS                                                    INT44960
         IF(IA(R,6).NE.1) GO TO 100                                     INT44970
         CALL XIN(NC,NC,SRA,-R,ISCR)                                    INT44980
         CALL XIN(NC,NC,SRN,-R,ISCR6)                                   INT44990
         RMSDIF=ZERO                                                    INT45000
         AVRG=ZERO                                                      INT45010
         NNZ=0                                                          INT45020
         DO 20  I=1,NC                                                  INT45030
         DO 20  J=1,NC                                                  INT45040
         SRN(I,J)=SRA(I,J)-SRN(I,J)                                     INT45050
         IF(DABS(SRA(I,J)).GT.AVCUT) THEN                               INT45060
            NNZ=NNZ+1                                                   INT45070
            AVRG=AVRG+DABS(SRA(I,J))                                    INT45080
            RMSDIF=RMSDIF+SRN(I,J)**2                                   INT45090
         END IF                                                         INT45100
   20    CONTINUE                                                       INT45110
         IF(NNZ.EQ.0) NNZ=1                                             INT45120
         AVRG=AVRG/NNZ                                                  INT45130
         RMSDIF=DSQRT(RMSDIF)/NNZ                                       INT45140
         WRITE(IOUT,2) R,AVRG,RMSDIF                                    INT45150
         IF(LPRT(2,NPRT).GE.1) THEN                                     INT45160
           WRITE(IOUT,3) R                                              INT45170
           DO 30  I=1,NC                                                INT45180
   30      WRITE(IOUT,4) (SRA(I,J),J=1,NC)                              INT45190
         END IF                                                         INT45200
         IF(LPRT(2,NPRT).GE.2) THEN                                     INT45210
           WRITE(IOUT,5)                                                INT45220
           DO 35  I=1,NC                                                INT45230
   35      WRITE(IOUT,6) (SRN(I,J),J=1,NC)                              INT45240
         END IF                                                         INT45250
         CALL XIN(NC,NC,SRA,R,ISCR)                                     INT45260
         CALL XIN(NC,NC,SRN,R,ISCR6)                                    INT45270
         RMSDIF=ZERO                                                    INT45280
         AVRG=ZERO                                                      INT45290
         NNZ=0                                                          INT45300
         DO 40  I=1,NSX                                                 INT45310
         DO 40  J=1,NSX                                                 INT45320
         SRN(I,J)=SRA(I,J)-SRN(I,J)                                     INT45330
         IF(DABS(SRA(I,J)).GT.AVCUT) THEN                               INT45340
            NNZ=NNZ+1                                                   INT45350
            AVRG=AVRG+DABS(SRA(I,J))                                    INT45360
            RMSDIF=RMSDIF+SRN(I,J)**2                                   INT45370
         END IF                                                         INT45380
   40    CONTINUE                                                       INT45390
         IF(NNZ.EQ.0) NNZ=1                                             INT45400
         AVRG=AVRG/NNZ                                                  INT45410
         RMSDIF=DSQRT(RMSDIF)/NNZ                                       INT45420
         WRITE(IOUT,7) R,AVRG,RMSDIF                                    INT45430
         IF(LPRT(2,NPRT).GE.1) THEN                                     INT45440
           WRITE(IOUT,8) R                                              INT45450
           DO 50  I=1,NSX                                               INT45460
   50      WRITE(IOUT,4) (SRA(I,J),J=1,NSX)                             INT45470
         END IF                                                         INT45480
         IF(LPRT(2,NPRT).GE.2) THEN                                     INT45490
           WRITE(IOUT,9)                                                INT45500
           DO 55  I=1,NSX                                               INT45510
   55      WRITE(IOUT,6) (SRN(I,J),J=1,NSX)                             INT45520
         END IF                                                         INT45530
  100 CONTINUE                                                          INT45540
      END                                                               INT45550
C     ///////////////////////////////////////////////////////////////   INT45560
      SUBROUTINE SRTST2(NC,NS,NSX,NSYM,IA,SRA,SRN)                      INT45570
      IMPLICIT REAL*8 (A-H,O-Z)                                         INT45580
      INTEGER R                                                         INT45590
      DIMENSION IA(NS,6),SRA(NC,NC,NC),SRN(NC,NC,NC)                    INT45600
      PARAMETER(ZERO=0.0D0,ONE=1.0D0,AVCUT=1.0D-5,ERCUT=1.0D-5)         INT45610
      COMMON /IO/ IIN1,IOUT,IIN2,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25,
     $   I31,I32,I33,I35,I36,I37,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10,
     $   ISCR11,ISCR12,ISCR13,ISCR14
    1 FORMAT(//1X,'TESTING OF BR(I,J,K) AND CR(M,N,P) MATRICES',/)      INT45660
    2 FORMAT(/1X,'R=',I3,' BR:  AVRG=',F12.8,'  RMSDIF=',F12.8)         INT45670
    3 FORMAT(//1X,'BR(I,J,K) MATRIX FOR COORDINATE ',I5/)               INT45680
    4 FORMAT(6F12.6)                                                    INT45690
    5 FORMAT(//1X,'BR(I,J,K) ANALYTIC - BR(I,J,K) NUMERICAL'/)          INT45700
    6 FORMAT(6F12.6)                                                    INT45710
    7 FORMAT(1X,' I= ',I5)                                              INT45720
    8 FORMAT(1X,' ERROR IN H',3I1)                                      INT45730
    9 FORMAT(/1X,'R=',I3,' CR:  AVRG=',F12.8,'  RMSDIF=',F12.8)         INT45740
   10 FORMAT(//1X,'CR(M,N,P) MATRIX FOR COORDINATE ',I5/)               INT45750
   11 FORMAT(//1X,'CR(M,N,P) ANALYTIC - CR(M,N,P) NUMERICAL'/)          INT45760
      IF(NSYM.EQ.0) THEN                                                INT45770
         ISCR=ISCR3                                                     INT45780
      ELSE                                                              INT45790
         ISCR=ISCR4                                                     INT45800
      END IF                                                            INT45810
      WRITE(IOUT,1)                                                     INT45820
      NA=NC/3                                                           INT45830
      DO 100  R=1,NS                                                    INT45840
         IF(IA(R,6).NE.2) GO TO 100                                     INT45850
         CALL YIN(NC,NC,SRA,-R,ISCR)                                    INT45860
         CALL YIN(NC,NC,SRN,-R,ISCR7)                                   INT45870
         RMSDIF=ZERO                                                    INT45880
         AVRG=ZERO                                                      INT45890
         NNZ=0                                                          INT45900
         DO 20  I=1,NC                                                  INT45910
         DO 20  J=1,NC                                                  INT45920
         DO 20  K=1,NC                                                  INT45930
         SRN(I,J,K)=SRA(I,J,K)-SRN(I,J,K)                               INT45940
         IF(DABS(SRA(I,J,K)).GT.AVCUT) THEN                             INT45950
            NNZ=NNZ+1                                                   INT45960
            AVRG=AVRG+DABS(SRA(I,J,K))                                  INT45970
            RMSDIF=RMSDIF+SRN(I,J,K)**2                                 INT45980
         END IF                                                         INT45990
   20    CONTINUE                                                       INT46000
         IF(NNZ.EQ.0) NNZ=1                                             INT46010
         AVRG=AVRG/NNZ                                                  INT46020
         RMSDIF=DSQRT(RMSDIF)/NNZ                                       INT46030
         WRITE(IOUT,2) R,AVRG,RMSDIF                                    INT46040
         IF(LPRT(2,NPRT).GE.1) THEN                                     INT46050
           WRITE(IOUT,3) R                                              INT46060
           DO 30  I=1,NC                                                INT46070
           WRITE(IOUT,7) I                                              INT46080
           DO 30  J=1,NC                                                INT46090
   30      WRITE(IOUT,4) (SRA(I,J,K),K=1,NC)                            INT46100
         END IF                                                         INT46110
         IF(LPRT(2,NPRT).GE.2) THEN                                     INT46120
           WRITE(IOUT,5)                                                INT46130
           DO 35  I=1,NC                                                INT46140
           WRITE(IOUT,7) I                                              INT46150
           DO 35  J=1,NC                                                INT46160
   35      WRITE(IOUT,6) (SRN(I,J,K),K=1,NC)                            INT46170
         END IF                                                         INT46180
         DO 40  I=1,NC                                                  INT46190
         DO 40  J=1,NC                                                  INT46200
         DO 40  K=1,NC                                                  INT46210
   40    SRA(I,J,K)=ZERO                                                INT46220
         DO 50  I=1,NC                                                  INT46230
         DO 50  J=1,NC                                                  INT46240
         DO 50  K=1,NC                                                  INT46250
         IF(DABS(SRN(I,J,K)).LT.ERCUT) GO TO 50                         INT46260
         INA=(I-1)/3+1                                                  INT46270
         JNA=(J-1)/3+1                                                  INT46280
         KNA=(K-1)/3+1                                                  INT46290
           DO 45  MM=1,5                                                INT46300
         IF(INA.EQ.IA(R,MM)) INA=-MM                                    INT46310
         IF(JNA.EQ.IA(R,MM)) JNA=-MM                                    INT46320
         IF(KNA.EQ.IA(R,MM)) KNA=-MM                                    INT46330
   45      CONTINUE                                                     INT46340
         INA=-INA                                                       INT46350
         JNA=-JNA                                                       INT46360
         KNA=-KNA                                                       INT46370
         SRA(INA,JNA,KNA)=ONE                                           INT46380
   50    CONTINUE                                                       INT46390
         DO 55  I=1,NA                                                  INT46400
         DO 55  J=1,I                                                   INT46410
         DO 55  K=1,J                                                   INT46420
         IF(SRA(I,J,K).GT.ZERO) THEN                                    INT46430
         WRITE(IOUT,8) I,J,K                                            INT46440
         END IF                                                         INT46450
   55    CONTINUE                                                       INT46460
         CALL YIN(NC,NSX,SRA,R,ISCR)                                    INT46470
         CALL YIN(NC,NSX,SRN,R,ISCR7)                                   INT46480
         RMSDIF=ZERO                                                    INT46490
         AVRG=ZERO                                                      INT46500
         NNZ=0                                                          INT46510
         DO 60  I=1,NSX                                                 INT46520
         DO 60  J=1,NSX                                                 INT46530
         DO 60  K=1,NSX                                                 INT46540
         SRN(I,J,K)=SRA(I,J,K)-SRN(I,J,K)                               INT46550
         IF(DABS(SRA(I,J,K)).GT.AVCUT) THEN                             INT46560
            NNZ=NNZ+1                                                   INT46570
            AVRG=AVRG+DABS(SRA(I,J,K))                                  INT46580
            RMSDIF=RMSDIF+SRN(I,J,K)**2                                 INT46590
         END IF                                                         INT46600
   60    CONTINUE                                                       INT46610
         IF(NNZ.EQ.0) NNZ=1                                             INT46620
         AVRG=AVRG/NNZ                                                  INT46630
         RMSDIF=DSQRT(RMSDIF)/NNZ                                       INT46640
         WRITE(IOUT,9) R,AVRG,RMSDIF                                    INT46650
         IF(LPRT(2,NPRT).GE.1) THEN                                     INT46660
           WRITE(IOUT,10) R                                             INT46670
           DO 65  I=1,NSX                                               INT46680
           WRITE(IOUT,7) I                                              INT46690
           DO 65  J=1,NSX                                               INT46700
   65      WRITE(IOUT,4) (SRA(I,J,K),K=1,NSX)                           INT46710
         END IF                                                         INT46720
         IF(LPRT(2,NPRT).GE.2) THEN                                     INT46730
           WRITE(IOUT,11)                                               INT46740
           DO 70  I=1,NSX                                               INT46750
           WRITE(IOUT,7) I                                              INT46760
           DO 70  J=1,NSX                                               INT46770
   70      WRITE(IOUT,6) (SRN(I,J,K),K=1,NSX)                           INT46780
         END IF                                                         INT46790
  100 CONTINUE                                                          INT46800
      END                                                               INT46810
C  \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        SUBROUTINE FORMP (IOPT,XMASS,NA,NC,XAR,DK1,DK2,DK3,
     $			  G3,D3,X2,X4,P2,P3,P4)
C  \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C  FORM CART PROJECTION MATRICES,
C  \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C  ON ENTRY:
C       NA                      NUMBER OF ATOMS
C       NC                      NA*3
C       XAR     (NA,3)	        REFERENCE GEOMETERY
C       DK1     (6,NC)          1ST DERIV'S
C       DK2     (3,NC,NC)       2ND DERIV'S   (USED ONLY IF NDER>1)
C       DK3     (3,NC,NC,NC)    3RD DERIV'S   (USED ONLY IF NDER>2)
C	G3	(NC,NC,NC)	SCRATCH SPACE (USED ONLY IF NDER>1)
C	D3	(NC,NC,NC)	SCRATCH SPACE (USED ONLY IF NDER>1)
C	X2	(NC,NC)		SCRATCH SPACE (USED ONLY IF NDER>1)
C	X4	(NC,NC,NC,NC)	SCRATCH SPACE (USED ONLY IF NDER>2)
C  ON RETURN:
C       P2      (NC,NC)         1ST ORDER PROJECTION MATRIX
C       P3      (NC,NC,NC)      2ND ORDER PROJECTION MATRIX
C       P4      (NC,NC,NC,NC)   3RD ORDER PROJECTION MATRIX
C  \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        IMPLICIT REAL*8 (A-H,O-Z)
        PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)
      COMMON /IO/ IIN1,IOUT,IIN2,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25,
     $   I31,I32,I33,I35,I36,I37,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10,
     $   ISCR11,ISCR12,ISCR13,ISCR14
        DIMENSION XAR(NA,3), P2(NC,NC), P3(NC,NC,NC), P4(NC,NC,NC,NC)
        DIMENSION DK1(6,NC), DK2(3,NC,NC), DK3(3,NC,NC,NC), IOPT(30)
	DIMENSION XMASS(NA)
        DIMENSION G3(NC,NC,NC), D3(NC,NC,NC), X2(NC,NC), X4(NC,NC,NC,NC)
  1	FORMAT (3F20.10)
  2     FORMAT (/' P(I,J) MATRIX'/) 
  3     FORMAT (/' P(I,J,K) MATRIX')
  4     FORMAT (/' P(I,J,K,L) MATRIX')
  5     FORMAT (/' I=',I5/)
  6     FORMAT (6F12.6)
  7     FORMAT (/' I=',I5,' J=',I5/)

        NDER=IOPT(4)
        NSTOP=IOPT(14)

C FORM P2
        DO 10 I = 1, NA
	IX = (I-1)*3+1
        IY = IX + 1
        IZ = IY + 1
        DO 15 J = 1, NC
        P2(IX,J) = -DK1(1,J)-XAR(I,3)*DK1(5,J)+XAR(I,2)*DK1(6,J)
        P2(IY,J) = -DK1(2,J)+XAR(I,3)*DK1(4,J)-XAR(I,1)*DK1(6,J)
  15    P2(IZ,J) = -DK1(3,J)-XAR(I,2)*DK1(4,J)+XAR(I,1)*DK1(5,J)
        P2(IX,IX) = P2(IX,IX) + ONE
        P2(IY,IY) = P2(IY,IY) + ONE
  10    P2(IZ,IZ) = P2(IZ,IZ) + ONE
C SAVE TO DISK IF NSTOP SET
        IF (NSTOP.EQ.3) THEN
                REWIND I35
                WRITE (I35,1) ((P2(I,J), J=1,NC), I=1,NC)
                END IF
C PRINT IF SPECIFIED BY PRINT OPTION
        IF (LPRT(1,IOPT(6)).GE.2) THEN
                WRITE (IOUT,2)
		DO 11 I = 1, NC
  11		WRITE (IOUT,6) (P2(I,J), J=1,NC)
		END IF
        IF (NDER.LE.2) RETURN
C FORM G
        DO 110 J1 = 1, NA
	JX = (J1-1)*3+1
        JY = JX + 1
        JZ = JY + 1
        DO 110 J2 = 1, NC
        DO 110 J3 = 1, J2
        G3(JX,J2,J3) = XAR(J1,3)*DK2(2,J2,J3)-XAR(J1,2)*DK2(3,J2,J3)
        G3(JY,J2,J3) = XAR(J1,1)*DK2(3,J2,J3)-XAR(J1,3)*DK2(1,J2,J3)
        G3(JZ,J2,J3) = XAR(J1,2)*DK2(1,J2,J3)-XAR(J1,1)*DK2(2,J2,J3)
        G3(JX,J3,J2) = G3(JX,J2,J3)
        G3(JY,J3,J2) = G3(JY,J2,J3)
 110    G3(JZ,J3,J2) = G3(JZ,J2,J3)

C FORM D3
	DO 150 JX = 1, NC-1, 3
	JY = JX + 1
	JZ = JY + 1
	DO 150 I1 = 1, NC
	DO 150 I2 = 1, NC
	D3(I1,I2,JX)= P2(JY,I2)*DK1(6,I1)-P2(JZ,I2)*DK1(5,I1)
	D3(I1,I2,JY)=-P2(JX,I2)*DK1(6,I1)+P2(JZ,I2)*DK1(4,I1)
 150	D3(I1,I2,JZ)= P2(JX,I2)*DK1(5,I1)-P2(JY,I2)*DK1(4,I1)
	DO 155 IX = 1, NC-1, 3
	IY = IX + 1
	IZ = IY + 1
	DO 155 J = 1, NC
	D3(IX,J,IY) = D3(IX,J,IY) - DK1(6,J)
	D3(IX,J,IZ) = D3(IX,J,IZ) + DK1(5,J)
	D3(IY,J,IX) = D3(IY,J,IX) + DK1(6,J)
	D3(IY,J,IZ) = D3(IY,J,IZ) - DK1(4,J)
	D3(IZ,J,IX) = D3(IZ,J,IX) - DK1(5,J)
 155	D3(IZ,J,IY) = D3(IZ,J,IY) + DK1(4,J)

C FORM P3
        DO 170 L = 1, NC
        DO 170 I2 = 1, NC
        DO 160 J1 = 1, NC
        DO 160 J2 = 1, NC
 160    X2(J1,J2) = P2(J2,I2)*G3(L,J1,J2)
        DO 170 I1 = 1, NC
        XXX = ZERO
        DO 165 J1 = 1, NC
        XXX = XXX + P2(L,J1)*D3(I1,I2,J1)
        DO 165 J2 = 1, NC
 165    XXX = XXX - P2(J1,I1)*X2(J1,J2)
 170    P3(L,I1,I2) = XXX


C SAVE TO DISK IF NSTOP SET
        IF (IOPT(14).EQ.3) THEN
                REWIND I36
                WRITE (I36,1) (((P3(I,J,K), K=1,NC),J=1,NC),I=1,NC)
                END IF
C PRINT IF SPECIFIED BY PRINT OPTION
        IF (LPRT(1,IOPT(6)).GE.2) THEN
                WRITE (IOUT,3)
		DO 180 I = 1, NC
		WRITE (IOUT,5) I
                DO 180 J = 1, NC
 180          	WRITE (IOUT,6) (P3(I,J,K), K=1,NC)
                END IF
	IF (NDER.LE.3) RETURN

C FORM -G PART OF SUM IN TERM 1 OF P4
        DO 220 J1 = 1, NA
	JX = (J1-1)*3+1
        JY = JX + 1
        JZ = JY + 1
        DO 220 J2 = 1, NC
        DO 220 J3 = 1, J2
	DO 220 J4 = 1, J3
        X4(JX,J2,J3,J4)=-XAR(J1,3)*DK3(2,J2,J3,J4)
     $			+XAR(J1,2)*DK3(3,J2,J3,J4)
        X4(JY,J2,J3,J4)=-XAR(J1,1)*DK3(3,J2,J3,J4)
     $			+XAR(J1,3)*DK3(1,J2,J3,J4)
        X4(JZ,J2,J3,J4)=-XAR(J1,2)*DK3(1,J2,J3,J4)
     $ 			+XAR(J1,1)*DK3(2,J2,J3,J4)
	X4(JX,J2,J4,J3) = X4(JX,J2,J3,J4)
	X4(JX,J3,J2,J4) = X4(JX,J2,J3,J4)
	X4(JX,J3,J4,J2) = X4(JX,J2,J3,J4)
	X4(JX,J4,J2,J3) = X4(JX,J2,J3,J4)
	X4(JX,J4,J3,J2) = X4(JX,J2,J3,J4)
	X4(JY,J2,J4,J3) = X4(JY,J2,J3,J4)
	X4(JY,J3,J2,J4) = X4(JY,J2,J3,J4)
	X4(JY,J3,J4,J2) = X4(JY,J2,J3,J4)
	X4(JY,J4,J2,J3) = X4(JY,J2,J3,J4)
	X4(JY,J4,J3,J2) = X4(JY,J2,J3,J4)
	X4(JZ,J2,J4,J3) = X4(JZ,J2,J3,J4)
	X4(JZ,J3,J2,J4) = X4(JZ,J2,J3,J4)
	X4(JZ,J3,J4,J2) = X4(JZ,J2,J3,J4)
	X4(JZ,J4,J2,J3) = X4(JZ,J2,J3,J4)
 220	X4(JZ,J4,J3,J2) = X4(JZ,J2,J3,J4)

C FORM R PARTS OF SUM IN TERM 1 OF P4
	DO 240 I = 1, NC
	DO 240 JX = 1, NC-1, 3
	JY = JX + 1
	JZ = JY + 1
	DO 240 K = 1, NC
	RX = -P2(I,JY)*DK2(3,K,K)+P2(I,JZ)*DK2(2,K,K)
	RY =  P2(I,JX)*DK2(3,K,K)-P2(I,JZ)*DK2(1,K,K)
	RZ = -P2(I,JX)*DK2(2,K,K)+P2(I,JY)*DK2(1,K,K)
	X4(I,JX,K,K) = X4(I,JX,K,K) + RX
	X4(I,JY,K,K) = X4(I,JY,K,K) + RY
	X4(I,JZ,K,K) = X4(I,JZ,K,K) + RZ
	X4(I,K,JX,K) = X4(I,K,JX,K) + RX
	X4(I,K,JY,K) = X4(I,K,JY,K) + RY
	X4(I,K,JZ,K) = X4(I,K,JZ,K) + RZ
	X4(I,K,K,JX) = X4(I,K,K,JX) + RX
	X4(I,K,K,JY) = X4(I,K,K,JY) + RY
	X4(I,K,K,JZ) = X4(I,K,K,JZ) + RZ
	DO 240 L = 1, K-1
	RX = -P2(I,JY)*DK2(3,K,L)+P2(I,JZ)*DK2(2,K,L)
	RY =  P2(I,JX)*DK2(3,K,L)-P2(I,JZ)*DK2(1,K,L)
	RZ = -P2(I,JX)*DK2(2,K,L)+P2(I,JY)*DK2(1,K,L)
	X4(I,JX,K,L) = X4(I,JX,K,L) + RX
	X4(I,JX,L,K) = X4(I,JX,L,K) + RX
	X4(I,JY,K,L) = X4(I,JY,K,L) + RY
	X4(I,JY,L,K) = X4(I,JY,L,K) + RY
	X4(I,JZ,K,L) = X4(I,JZ,K,L) + RZ
 	X4(I,JZ,L,K) = X4(I,JZ,L,K) + RZ
	X4(I,K,L,JX) = X4(I,K,L,JX) + RX
	X4(I,K,L,JY) = X4(I,K,L,JY) + RY
	X4(I,K,L,JZ) = X4(I,K,L,JZ) + RZ
	X4(I,L,K,JX) = X4(I,L,K,JX) + RX
	X4(I,L,K,JY) = X4(I,L,K,JY) + RY
	X4(I,L,K,JZ) = X4(I,L,K,JZ) + RZ
	X4(I,K,JX,L) = X4(I,K,JX,L) + RX
	X4(I,L,JX,K) = X4(I,L,JX,K) + RX
	X4(I,K,JY,L) = X4(I,K,JY,L) + RY
	X4(I,L,JY,K) = X4(I,L,JY,K) + RY
	X4(I,K,JZ,L) = X4(I,K,JZ,L) + RZ
 240	X4(I,L,JZ,K) = X4(I,L,JZ,K) + RZ

C DO THE SUMMATIONS OF THE FIRST TERMS IN P4

        DO 250  L = 1, NC
        DO 250  I3 = 1, NC
        DO 250  J1 = 1, NC
        DO 250  J2 = 1, NC
        P4(L,J1,J2,I3) = ZERO
        DO 250  J3 = 1, NC
 250    P4(L,J1,J2,I3)=P4(L,J1,J2,I3)+P2(J3,I3)*X4(L,J1,J2,I3)

        DO 252  L = 1, NC
        DO 252  I3 = 1, NC
        DO 252  I2 = 1, NC
        DO 252  J1 = 1, NC
        X4(L,J1,I2,I3) = ZERO
        DO 252  J2 = 1, NC
 252    X4(L,J1,I2,I3)=P4(L,J1,I2,I3)+P2(J2,I2)*P4(L,J1,J2,I3)

        DO 254  L = 1, NC
        DO 254  I3 = 1, NC
        DO 254  I2 = 1, NC
        DO 254  I1 = 1, NC
        P4(L,I1,I2,I3) = ZERO
        DO 254  J1 = 1, NC
 254    P4(L,I1,I2,I3)=P4(L,I1,I2,I3)+P2(J1,I1)*X4(L,J1,I2,I3)

C PUT TERM 3 INTO P4
        DO 230 I1 = 1, NC
        DO 230 I2 = 1, NC
        X44 = DK1(4,I2)*DK1(4,I1)
        X55 = DK1(5,I2)*DK1(5,I1)
        X45 = DK1(4,I2)*DK1(5,I1)
        X46 = DK1(4,I2)*DK1(6,I1)
        X54 = DK1(5,I2)*DK1(4,I1)
        X66 = DK1(6,I2)*DK1(6,I1)
        X56 = DK1(5,I2)*DK1(6,I1)
        X64 = DK1(6,I2)*DK1(4,I1)
        X65 = DK1(6,I2)*DK1(5,I1)
        X6655 = X66 + X55
        X6644 = X66 + X44
        X5544 = X55 + X44
        DO 230 L = 1, NC
        DO 230 I3X = 1, NC-1, 3
        I3Y = I3X + 1
        I3Z = I3Y + 1
        SX = -P2(L,I3X)*X6655 + P2(L,I3Y)*X45 + P2(L,I3Z)*X46
        SY = P2(L,I3X)*X54 - P2(L,I3Y)*X6644 + P2(L,I3Z)*X56
        SZ = P2(L,I3X)*X64 + P2(L,I3Y)*X65 - P2(L,I3Z)*X5544
        DO 235 JX = 1, NC-1, 3
        JY = JX + 1
        JZ = JY + 1
        SX = SX +
     $  P2(L,JX)*(DK1(6,I3X)*D3(I2,I1,JY)-DK1(5,I3X)*D3(I2,I1,JZ)) +
     $  P2(L,JY)*(DK1(4,I3X)*D3(I2,I1,JZ)-DK1(6,I3X)*D3(I2,I1,JX)) +
     $  P2(L,JZ)*(DK1(5,I3X)*D3(I2,I1,JX)-DK1(4,I3X)*D3(I2,I1,JY))
        SY = SY +
     $  P2(L,JX)*(DK1(6,I3Y)*D3(I2,I1,JY)-DK1(5,I3Y)*D3(I2,I1,JZ)) +
     $  P2(L,JY)*(DK1(4,I3Y)*D3(I2,I1,JZ)-DK1(6,I3Y)*D3(I2,I1,JX)) +
     $  P2(L,JZ)*(DK1(5,I3Y)*D3(I2,I1,JX)-DK1(4,I3Y)*D3(I2,I1,JY))
 235    SZ = SZ +
     $  P2(L,JX)*(DK1(6,I3Z)*D3(I2,I1,JY)-DK1(5,I3Z)*D3(I2,I1,JZ)) +
     $  P2(L,JY)*(DK1(4,I3Z)*D3(I2,I1,JZ)-DK1(6,I3Z)*D3(I2,I1,JX)) +
     $  P2(L,JZ)*(DK1(5,I3Z)*D3(I2,I1,JX)-DK1(4,I3Z)*D3(I2,I1,JY))
        P4(L,I1,I2,I3X) = SX
        P4(L,I1,I2,I3Y) = SY
 230    P4(L,I1,I2,I3Z) = SZ
C
C FORM DELTA AND PUT INTO M (HELD IN X4)
	DO 210 I = 1, NC
	DO 211 J = 1, NC
	DO 211 K = 1, NC
 211	X4(J,K,I,I) = ZERO
	DO 210 J = 1, I-1
	VX = DK1(5,I)*DK1(6,J) - DK1(6,I)*DK1(5,J)
	VY = DK1(6,I)*DK1(4,J) - DK1(4,I)*DK1(6,J)
	VZ = DK1(4,I)*DK1(5,J) - DK1(5,I)*DK1(4,J)
	DO 210 IX = 1, NC-1, 3
	IY = IX + 1
	IZ = IY + 1
	DO 210 L = 1, NC
	X4(L,IX,I,J) = P2(L,IY)*VZ - P2(L,IZ)*VY
	X4(L,IY,I,J) = P2(L,IZ)*VX - P2(L,IX)*VZ
	X4(L,IZ,I,J) = P2(L,IX)*VY - P2(L,IY)*VX
	X4(L,IX,J,I) = -X4(L,IX,I,J)
	X4(L,IY,J,I) = -X4(L,IY,I,J)
 210   	X4(L,IZ,J,I) = -X4(L,IZ,I,J)

C SUM 2ND AND 3RD TERMS INTO M
	DO 270 L = 1, NC
	DO 270 I1X = 1, NC-1, 3
	I1Y = I1X + 1
	I1Z = I1Y + 1
	DO 270 I2 = 1, NC
	DO 270 I3 = 1, NC
	XXX = ZERO
	YYY = ZERO
	ZZZ = ZERO
	DO 271 J3 = 1, NC
 	X4(L,I1X,I2,I3) = X4(L,I1X,I2,I3)-
     $		P2(L,I1Y)*DK1(6,J3)*D3(I2,I3,J3)+
     $		P2(L,I1Z)*DK1(5,J3)*D3(I2,I3,J3)
 	X4(L,I1Y,I2,I3) = X4(L,I1Y,I2,I3)+
     $		P2(L,I1X)*DK1(6,J3)*D3(I2,I3,J3)-
     $		P2(L,I1Z)*DK1(4,J3)*D3(I2,I3,J3)
 	X4(L,I1Z,I2,I3) = X4(L,I1Z,I2,I3)-
     $		P2(L,I1X)*DK1(5,J3)*D3(I2,I3,J3)+
     $		P2(L,I1Y)*DK1(4,J3)*D3(I2,I3,J3)
 	XXX = XXX + G3(L,I1X,J3)*P3(J3,I2,I3)
 	YYY = YYY + G3(L,I1Y,J3)*P3(J3,I2,I3)
 271	ZZZ = ZZZ + G3(L,I1Z,J3)*P3(J3,I2,I3)
 	X4(L,I1X,I2,I3) = X4(L,I1X,I2,I3) - XXX
 	X4(L,I1Y,I2,I3) = X4(L,I1Y,I2,I3) - YYY
 270	X4(L,I1Z,I2,I3) = X4(L,I1Z,I2,I3) - ZZZ

C SUM IN TERM 2 INTO P4
	DO 260  L = 1, NC
	DO 260 I1 = 1, NC
	DO 260 I2 = 1, NC
	DO 260 I3 = 1, NC
	XXX = ZERO
	DO 265 J = 1, NC
 265	XXX = XXX + P2(J,I1)*X4(L,J,I2,I3)
 	P4(L,I1,I2,I3) = P4(L,I1,I2,I3) + XXX
  	P4(L,I2,I1,I3) = P4(L,I2,I1,I3) + XXX
 260   	P4(L,I2,I3,I1) = P4(L,I2,I3,I1) + XXX

C SAVE TO DISK IF NSTOP SET
        IF (IOPT(14).EQ.3) THEN
                REWIND I37
                WRITE (I37,1)
     $			((((P4(I,J,K,L),L=1,NC),K=1,NC),J=1,NC),I=1,NC)
                END IF

C PRINT IF SPECIFIED BY PRINT OPTION
        IF (LPRT(1,IOPT(6)).GE.2) THEN
                WRITE (IOUT,4)
                DO 280 I = 1, NC
		DO 280 J = 1, NC
                WRITE (IOUT,7) I, J
                DO 280 K = 1, NC
 280            WRITE (IOUT,6) (P4(I,J,K,L), L=1,NC)
                END IF

	RETURN
        END
C \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
	SUBROUTINE PROJV (NC,IOPT,P2,P3,P4,F1,F2,F3,X3,X4,
     $			  F2P,F3P,F4P)
C \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C PROJECT CARTESIAN FORCE CONSTANTS AND TEST ROTATIONAL INVARIANCE
C RELATIONS IF TEST OPTION SPECIFIED.
C THIS ROUTINE SUMS IN THE NONLINEAR TERMS.  LINEAR TERMS ARE COMPUTED
C IN LINTR WHICH IS CALLED BEFORE PROJV.
C \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C ON ENTRY:
C	NC			NA*3
C	P2	(NC,NC)		1ST ORDER PROJECTION MATRIX
C	P3	(NC,NC,NC)	2ND	"	     "
C	P4	(NC,NC,NC,NC)	3RD	"	     "
C	F1	(NC)		SCRATCH SPACE
C	F2	(NC,NC)		SCRATCH SPACE
C	F3	(NC,NC,NC)	SCRATCH SPACE
C	X3	(NC,NC,NC)	SCRATCH SPACE
C	X4	(NC,NC,NC,NC)	SCRATCH SPACE
C	F2P	(NC,NC)		PROJECTED 2ND DERIVATIVES (COMPUTED PREV.
C				BY LINTR)
C ON RETURN:
C	F3P	(NC,NC,NC)	PROJECTED THIRD DERIVS
C	F4P	(NC,NC,NC,NC)	PROJECTED SECOND DERIVS
C NOTES:
C	1. HIGHEST ORDER PROJECTED DETERMINED BY IOPT(4)
C	2. PROJECTED FORCE CONSTANTS ARE SAVED AND PRINTED IN FCOUT
C	   WHICH IS CALLED ELSEWHERE
C \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

	IMPLICIT REAL*8 (A-H,O-Z)

	PARAMETER (TOL=1.0D-10)
	PARAMETER (ZERO = 0.0D0)

	LOGICAL NOERR, NOF1
	CHARACTER STRING*10
	DIMENSION IOPT(30),P2(NC,NC),P3(NC,NC,NC),P4(NC,NC,NC,NC)
	DIMENSION F1(NC), F2(NC,NC), F3(NC,NC,NC)
	DIMENSION X3(NC,NC,NC), X4(NC,NC,NC,NC)
	DIMENSION F2P(NC,NC), F3P(NC,NC,NC), F4P(NC,NC,NC,NC)
      COMMON /IO/ IIN1,IOUT,IIN2,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25,
     $   I31,I32,I33,I35,I36,I37,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10,
     $   ISCR11,ISCR12,ISCR13,ISCR14
      COMMON /PHYCON/ BOHR,DEBYE,HART,WAVE0,CINT
 2	FORMAT (/' TESTING ROTATIONAL INVARIANCE OF FORCE CONSTANTS:')
 3	FORMAT (/'     NO TESTS PERFORMED: NDER<2')
 4	FORMAT (/'     2ND ORDER RELATION HOLDS FOR ALL CASES.')
 5	FORMAT (/'     3RD ORDER RELATION HOLDS FOR ALL CASES.')
 6	FORMAT (/'     2ND ORDER RELATION FAILS FOR:')
 7	FORMAT ('         I1=',I3,'  I2=',I3)
 8	FORMAT (/'     3RD ORDER RELATION FAILS FOR:')
 9	FORMAT ('         I1=',I3,'  I2=',I3,'  I3=',I3)
 10	FORMAT ('     TOLERANCE IS SET TO:',E10.2)
 11	FORMAT (A10)
 12	FORMAT (20X,3F20.10)
 13	FORMAT (2I15)
 14	FORMAT (3F20.10)
 15	FORMAT (/'     NO RELATIONS WILL BE TESTED: ',
     $		 'NO GRADIENTS FOUND IN FILE 11.')

	BOHR2 = BOHR*BOHR
	NOF1 = .FALSE.
        NTEST=IOPT(9)
        NDER=IOPT(4)

	IF (NTEST.EQ.4) THEN
		WRITE (IOUT,2)
		WRITE (IOUT,10) TOL
		END IF

	IF (NDER.LT.2) THEN
	        IF (NTEST.EQ.4) WRITE (IOUT,3)
		RETURN
		END IF
        IF (NTEST.EQ.4.OR.NDER.GT.2) THEN
	  REWIND (I15)
	  READ(I15,13) I,J
	  READ(I15,14) ((F2(I,J),J=1,NC),I=1,NC)
        END IF

	IF (NDER.GT.3) THEN
        REWIND I20
 	READ (I20,13) M,N
 	READ (I20,14) (((F3(I,J,K), K = 1, J), J = 1, I), I = 1, NC)
 	CALL FILL3B (NC,NC,F3)
	END IF


C TEST 1ST ROTATIONAL INVARIANCE RELATION
	IF (NTEST.EQ.4) THEN
		NA = NC/3
		REWIND I11
 128		READ (I11,11,END=129) STRING
		GOTO 128
 129		DO 127 J = 1, NA+1
 127		BACKSPACE I11
      		DO 130 J = 1, NA
 130		READ (I11,12,END=132,ERR=132) (F1((J-1)*3+K),K=1,3)
		GOTO 133

 132		WRITE (IOUT, 15)
		NOF1 = .TRUE.
		GOTO 101	

 133		NOERR = .TRUE.
		DO 100 I1 = 1, NC
		DO 100 I2 = 1, NC
		XXX = ZERO
		DO 110 J1 = 1, NC
 110			XXX = XXX + F1(J1)*P3(J1,I1,I2)*BOHR
	 	IF (ABS(F2(I1,I2)-F2P(I1,I2)-XXX).GT.TOL) THEN
			IF (NOERR) WRITE (IOUT,6)
			NOERR = .FALSE.
			WRITE (IOUT,7) I1,I2
			END IF
 100		CONTINUE
        	IF (NOERR) WRITE (IOUT,4)
	END IF
 101	CONTINUE
      IF (NDER.LT.3) RETURN
C SUM IN ALL TERMS CONTAINING A P2,P3 PRODUCT
	DO 200 J1 = 1, NC
	DO 200 J2 = 1, NC

	DO 210 I = 1, NC
	DO 210 J = 1, NC
	DO 210 K = 1, NC
 210		X3(I,J,K) = P2(J1,I)*P3(J2,J,K)*BOHR

	DO 220 I1 = 1, NC
	DO 220 I2 = 1, I1
	DO 220 I3 = 1, I2
 220	F3P(I1,I2,I3) = F3P(I1,I2,I3) +
     $		F2(J1,J2)* ( X3(I1,I2,I3)+X3(I2,I1,I3)+X3(I3,I1,I2) )
 
 200    CONTINUE

	IF (NDER.GE.4) THEN

	DO 250 I1 = 1, NC

        DO 252 J2 = 1, NC
        DO 252 J3 = 1, NC
        X4(I1,J2,J3,1) = ZERO
        DO 252 J1 = 1, NC
 252    X4(I1,J2,J3,1) = X4(I1,J2,J3,1) + F3(J1,J2,J3)*P2(J1,I1)

        DO 254 I2 = 1, NC
        DO 254 J3 = 1, NC
        X3(I1,I2,J3) = ZERO
        DO 255 J2 = 1, NC
 255    X3(I1,I2,J3) = X3(I1,I2,J3) + X4(I1,J2,J3,1)*P2(J2,I2)
 254    X3(I1,I2,J3)=X3(I1,I2,J3)*BOHR
 
 250    CONTINUE

        DO 256  I1 = 1, NC
        DO 256  I2 = 1, I1
        DO 256  I3 = 1, I2
        DO 256  I4 = 1, I3
        DO 256  J3 = 1, NC
 256    F4P(I1,I2,I3,I4) = F4P(I1,I2,I3,I4) + 
     $         X3(I1,I2,J3)*P3(J3,I3,I4) + X3(I1,I3,J3)*P3(J3,I2,I4) +
     $         X3(I1,I4,J3)*P3(J3,I2,I3) + X3(I2,I3,J3)*P3(J3,I1,I4) +
     $         X3(I2,I4,J3)*P3(J3,I1,I3) + X3(I3,I4,J3)*P3(J3,I1,I2) 

	END IF

        CALL FILL3B(NC,NC,F3P)

	IF (NDER.LT.4) RETURN

C TEST SECOND INVARIACE RELATION (EQ 66)
	IF (NTEST.EQ.4.AND..NOT.NOF1) THEN
		NOERR = .TRUE.
		DO 290 I1 = 1, NC
		DO 290 I2 = 1, NC
		DO 290 I3 = 1, NC
		XXX = ZERO
		DO 295 J = 1, NC
 295		XXX = XXX + F1(J)*P4(J,I1,I2,I3)*BOHR2
		IF (ABS(F3(I1,I2,I3)-F3P(I1,I2,I3)-XXX).GT.TOL) THEN
			IF (NOERR) WRITE (IOUT,8)
			NOERR = .FALSE.
			WRITE (IOUT,9) I1,I2,I3
			END IF
 290		CONTINUE
		IF (NOERR) WRITE (IOUT,5)
	END IF

C COMPLETE FOURTH ORDER FORCE CONSTANTS

        DO 300 I1 = 1, NC
        DO 300 I2 = 1, NC
        DO 300 J2 = 1, NC
        X3(J2,I1,I2) = ZERO
        DO 301 J1 = 1, NC
 301    X3(J2,I1,I2) = X3(J2,I1,I2)+F2(J1,J2)*P3(J1,I1,I2)
 300    X3(J2,I1,I2)=X3(J2,I1,I2)*BOHR2

        DO 310 I1 = 1, NC
        DO 310 I2 = 1, I1
        DO 310 I3 = 1, I2
        DO 310 I4 = 1, I3
        DO 310 J2 = 1, NC
 310    F4P(I1,I2,I3,I4) = F4P(I1,I2,I3,I4) +
     $     X3(J2,I1,I2)*P3(J2,I3,I4) + X3(J2,I1,I3)*P3(J2,I2,I4) +
     $     X3(J2,I1,I4)*P3(J2,I2,I3) 

	DO 330 I1 = 1, NC
	DO 330 J2 = 1, NC
        X3(J2,I1,1) = ZERO
	DO 331 J1 = 1, NC
 331	X3(J2,I1,1) = X3(J2,I1,1)+F2(J1,J2)*P2(J1,I1)
 330    X3(J2,I1,1)=X3(J2,I1,1)*BOHR2

	DO 340 I1 = 1, NC
	DO 340 I2 = 1, I1
	DO 340 I3 = 1, I2
	DO 340 I4 = 1, I3
        DO 340 J2 = 1, NC
 340	F4P(I1,I2,I3,I4) = F4P(I1,I2,I3,I4) +
     $     X3(J2,I1,1)*P4(J2,I2,I3,I4) + X3(J2,I2,1)*P4(J2,I1,I3,I4) + 
     $     X3(J2,I3,1)*P4(J2,I1,I2,I4) + X3(J2,I4,1)*P4(J2,I1,I2,I3) 

	CALL FILL4B(NC,F4P)
	RETURN
	END
C /////////////////////////////////////////////////////////////////////
      SUBROUTINE PROJK1(NAD,NA,NC,IOPT,XA,B2,A ,P)
C /////////////////////////////////////////////////////////////////////
C COMPUTES THE PROJECTION MATRIX (P(NC,NC)) WHICH REMOVES TRANSLATIONAL
C AND ROTATIONAL COMPONENTS IF THE FORCE CONSTANT MATRICES WERE
C OBTAINED AT A NON-STATIONARY POINT
C /////////////////////////////////////////////////////////////////////
C ON ENTRY:
C NAD         INTEGER        NA+NDUM
C NA          INTEGER        NUMBER OF ATOMS
C NC          INTEGER        3*NA
C XA          (NAD,3)        CARTESIAN COORDINATES (IN ANGSTROMS)
C ON RETURN:
C B2          ( 6,NC)        Bij DERIVATIVES
C A           (NC,6 )        B2 "INVERSE"
C P           (NC,NC)        PROJECTION MATRIX
C /////////////////////////////////////////////////////////////////////
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XA(NAD,3),B2( 6,NC), A(NC,6 ),P(NC,NC)
      DIMENSION W(6,12),IOPT(30),T0(3,6),V(3)
      PARAMETER(ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0,THREE=3.0D0,FOUR=4.0D0)
      COMMON /IO/ IIN1,IOUT,IIN2,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25,
     $   I31,I32,I33,I35,I36,I37,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10,
     $   ISCR11,ISCR12,ISCR13,ISCR14
C
    1 FORMAT(//1X,'B2 MATRIX',/1X,'(ROWS 1-3: TRANSLATION, ROWS 4-6:',
     $            ' ROTATIONS)',/)
    2 FORMAT(12F10.6)
    3 FORMAT(/1X,'PROJECTION MATRIX')
 
C COMPUTE CENTER OF MOLECULE
      CX=ZERO
      CY=ZERO
      CZ=ZERO
      TM=DBLE(NA)
      DO 40 I=1,NA
         CX=CX+XA(I,1)
         CY=CY+XA(I,2)
  40     CZ=CZ+XA(I,3)
      CX=CX/TM
      CY=CY/TM
      CZ=CZ/TM
 
C TRANSLATIONAL PART OF THE B2 MATRIX
      DO 60  K=1,NA
      DO 60  I=1,3
         L=3*(K-1)+I
  60     B2(I,L)=ONE/TM
C FORM AND INVERT THE INERTIA TENSOR
      DO 62  I=1,3
      DO 62  J=1,6
  62     T0(I,J)=ZERO
      DO 63  I=1,3
  63     T0(I,I+3)=ONE
C
      DO 70  K=1,NA
         T0(1,1)=T0(1,1)+(XA(K,2)-CY)**2+(XA(K,3)-CZ)**2
         T0(2,2)=T0(2,2)+(XA(K,1)-CX)**2+(XA(K,3)-CZ)**2
         T0(3,3)=T0(3,3)+(XA(K,1)-CX)**2+(XA(K,2)-CY)**2
         T0(2,1)=T0(2,1)-(XA(K,1)-CX)*(XA(K,2)-CY)
         T0(3,1)=T0(3,1)-(XA(K,1)-CX)*(XA(K,3)-CZ)
  70     T0(3,2)=T0(3,2)-(XA(K,2)-CY)*(XA(K,3)-CZ)
      T0(1,2)=T0(2,1)
      T0(1,3)=T0(3,1)
      T0(2,3)=T0(3,2)
C
      CALL FLIN(T0,3,3,3,DET)
C
C  FORM THE B2 MATRIX FOR ROTATIONS
      DO 72  I=1,NA
        DO 74  K=1,3
        DO 74  L=1,3
  74    T0(K,L)=ZERO
C
        V(1)=ZERO
        V(2)=(XA(I,3)-CZ)
        V(3)=-(XA(I,2)-CY)
        DO 76  K=1,3
        DO 76  L=1,3
  76    T0(K,1)=T0(K,1)+T0(K,L+3)*V(L)
C
        V(1)=-(XA(I,3)-CZ)
        V(2)=ZERO
        V(3)=(XA(I,1)-CX)
        DO 78  K=1,3
        DO 78  L=1,3
  78    T0(K,2)=T0(K,2)+T0(K,L+3)*V(L)
C
        V(1)=(XA(I,2)-CY)
        V(2)=-(XA(I,1)-CX)
        V(3)=ZERO
        DO 80  K=1,3
        DO 80  L=1,3
  80    T0(K,3)=T0(K,3)+T0(K,L+3)*V(L)
C
        L=3*(I-1)
        DO 82  K=1,3
        B2(K+3,L+1)=T0(K,1)
        B2(K+3,L+2)=T0(K,2)
  82    B2(K+3,L+3)=T0(K,3)
  72  CONTINUE
C
C FORM (BT*B)-1
      DO 90  I=1,6
      DO 90  J=1,12
  90    W(I,J)=ZERO
      DO 92  I=1,6
  92    W(I,I+6)=ONE
      DO 94 I=1,6
      DO 94 J=1,6
      DO 94 K=1,NC
  94    W(I,J)=W(I,J)+B2(I,K)*B2(J,K)
      CALL FLIN(W,6,6,6,DET)
C
      DO 96  I=1,6
      DO 96  J=1,NC
          A(J,I)=ZERO
      DO 96  K=1,6
 96       A(J,I)= A(J,I)+W(I,K+6)*B2(K,J)
      DO 98  I=1,NC
      DO 98  J=1,NC
         P(I,J)=ZERO
      DO 98  K=1,6
 98      P(I,J)=P(I,J)+B2(K,I)* A(J,K)
      DO 100 I=2,NC
      DO 100 J=1,I-1
         P(I,J)=-P(I,J)
 100     P(J,I)=-P(J,I)
      DO 102 I=1,NC
 102     P(I,I)=ONE-P(I,I)
C
      WRITE(IOUT,*)'P2:'
      DO 200 I = 1, NC
      DO 200 J = 1, NC
 200    WRITE(IOUT,201) I,J,P(I,J)
 201    FORMAT (I5,I5,'   ',F12.10)
 
      RETURN
      END
C /////////////////////////////////////////////////////////////////////
      SUBROUTINE PROJK2(NAD,NA,NC,IOPT,XA,B2,BB2,DK,A,P,P3,P3A)
C /////////////////////////////////////////////////////////////////////
C COMPUTE SECOND ORDER PROJECTION MATRIX TO REMOVE ROTATIONAL AND
C TRANSLATIONAL COMPONENTS IF THE FORCE CONSTANT MATRICES WERE
C COMPUTED AT A NON-STATIONARY POINT
C /////////////////////////////////////////////////////////////////////
C ON ENTRY:
C NAD      INTEGER      NA+NDUM
C NA       INTEGER      NUMBER OF ATOMS
C NC       INTEGER      NA*3
C XA       (NAD,3)      CARTESIAN COORDINATES (IN ANGSTROMS)
C B2       (6,NC)       B2 DERIVATIVES
C A        (NC,6)       B2 "INVERSE"
C DK       (NC,NC)      SCRATCH SPACE
C P        (NC,NC)      SECOND ORDER PROJECTION MATRIX
C ON RETURN:
C BB2      (3,NC,NC)    B2 2ND DERIVATIVES ROTATIONS ONLY
C P3       (NC,NC,NC)   3RD ORDER PROJECTION MATRIX
C P3A      (NC,NC,NC)   Gijk DERIVATIVES
C /////////////////////////////////////////////////////////////////////
 
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XA(NAD,3),B2(6, NC),A(NC,6 ),P(NC,NC),DK(NC,NC)
      DIMENSION P3(NC,NC,NC), P3A(NC,NC,NC), BB2(3,NC,NC)
      DIMENSION W(6,12),IOPT(30),T0(3,6),T1(3,3),V(3),CM(3)
      DIMENSION V1(3),V2(3),V3(3),V4(3),V5(3),V6(3)
      PARAMETER(ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0,THREE=3.0D0,FOUR=4.0D0)
      COMMON /IO/ IIN1,IOUT,IIN2,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25,
     $   I31,I32,I33,I35,I36,I37,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10,
     $   ISCR11,ISCR12,ISCR13,ISCR14
      COMMON /PHYCON/ BOHR,DEBYE,HART,WAVE0,CINT
C COMPUTE CENTER OF MOLECULE
      TM=DBLE(NA)
      DO 40 K=1,3
         CM(K)=ZERO
         DO 42  I=1,NA
  42     CM(K)=CM(K)+XA(I,K)
  40     CM(K)=CM(K)/TM
 
C FORM AND INVERT THE INERTIA TENSOR                                    PRO00720
      DO 62  I=1,3
      DO 62  J=1,6
  62     T0(I,J)=ZERO
      DO 63  I=1,3
  63     T0(I,I+3)=ONE
C
      DO 65  K=1,NA                                                     PRO00730
         DO 67  L=1,3
  67     V1(L)=XA(K,L)-CM(L)
         DO 69  L=2,3
         DO 69  M=1,L-1
  69     T0(L,M)=T0(L,M)-V1(L)*V1(M)
         T0(1,1)=T0(1,1)+V1(2)*V1(2)+V1(3)*V1(3)
         T0(2,2)=T0(2,2)+V1(1)*V1(1)+V1(3)*V1(3)
  65     T0(3,3)=T0(3,3)+V1(1)*V1(1)+V1(2)*V1(2)
      T0(1,2)=T0(2,1)
      T0(1,3)=T0(3,1)
      T0(2,3)=T0(3,2)
C
      CALL FLIN(T0,3,3,3,DET)
      DO 68  I=1,3
      DO 68  J=1,3
  68  T1(I,J)=T0(I,J+3)
C
C     FORM Gijk DERIVATIVES
C
      DO 104  I=1,NA
      DO 104  J=1,NA
      LI=3*(I-1)
      LJ=3*(J-1)
      DO 104  IB=1,3
      DO 104  JG=1,3
        LL=LI+IB
        LM=LJ+JG
        DO 108  K=1,3
        V1(K)=B2(K+3,LL)
        V2(K)=B2(K+3,LM)
        V3(K)=XA(I,K)-CM(K)
 108    V4(K)=XA(J,K)-CM(K)
        CALL SCAPRO(V3,V2,D1)
        CALL SCAPRO(V4,V1,D2)
        DO 109  K=1,3
  109   V(K)=-V2(K)*V3(IB)-V1(K)*V4(JG)
        V(IB)=V(IB)+D1
        V(JG)=V(JG)+D2
        DO 111  IA=1,NA
           DO 115  K=1,3
  115      V3(K)=XA(IA,K)-CM(K)
           CALL SCAPRO(V1,V3,D1)
           CALL SCAPRO(V2,V3,D2)
           CALL VECPRO(V3,V2,V4)
           CALL VECPRO(V3,V1,V5)
           DO 113  K=1,3
  113      V(K)=V(K)+(V4(K)*D1+V5(K)*D2)/TWO
  111   CONTINUE
        DO 110  K=1,3
        V1(K)=ZERO
        DO 110  L=1,3
        V1(K)=V1(K)+T1(K,L)*V(L)
  110   BB2(L,I,J)=V1(L)
        DO 112  K=1,NC
        XX=ZERO
        DO 114  N=1,3
  114   XX=XX+A(K,N+3)*V1(N)
  112   P3A(K,LL,LM)=XX
  104   CONTINUE
C
C   FORM Pijk
C
      DO 120  I=1,NC
        DO 122  K=1,NC
        DO 122  L=1,NC
        DK(K,L)=ZERO
        DO 122  M=1,NC
  122   DK(K,L)=DK(K,L)+P3A(I,L,M)*P(M,K)
        DO 124  J=1,NC
        DO 124  K=1,NC
        P3(I,J,K)=ZERO
        DO 124  L=1,NC
  124   P3(I,J,K)=P3(I,J,K)-DK(K,L)*P(L,J)
  120 CONTINUE
      DO  130  K=1,NC
        DO 132  J=1,NC
        DO 132  L=1,NC
        DK(J,L)=ZERO
        DO 132  M=1,NC
  132   DK(J,L)=DK(J,L)+P3A(K,L,M)*P(M,J)
        DO 134  I=1,NC
        DO 134  J=1,NC
        DO 134  L=1,NC
  134   P3(I,J,K)=P3(I,J,K)-DK(J,L)*P(L,I)
  130 CONTINUE
      DO 140  I=1,NC
      DO 140  J=1,NC
      DO 140  K=1,NC
      DO 142  L=1,NC
  142 P3(I,J,K)=P3(I,J,K)-P3A(J,L,K)*P(L,I)
  140 P3(I,J,K)=P3(I,J,K)*BOHR
      WRITE (IOUT,*) 'P3:'
      DO 3000 I = 1, NC
 3001 FORMAT (I5,I5,I5,'    ',F13.10)
      WRITE (IOUT,3001) I,7,7,P3(I,7,7)
      WRITE (IOUT,3001) I,8,8,P3(I,8,8)
      WRITE (IOUT,3001) I,9,9,P3(I,9,9)
      WRITE (IOUT,3001) I,7,8,P3(I,7,8)
      WRITE (IOUT,3001) I,7,9,P3(I,7,9)
 3000 WRITE (IOUT,3001) I,8,9,P3(I,9,8)
      RETURN
      END
C /////////////////////////////////////////////////////////////////////
      SUBROUTINE PDER2(NC,IOPT,F2,XS,P,P3,F3)
C /////////////////////////////////////////////////////////////////////
C ADD REMAINING TERMS TO V*ijk WHICH CONTAIN Pijk
C /////////////////////////////////////////////////////////////////////
C ON ENTRY:
C NC        INTEGER      NUMBER OF ATOMS * 3
C F2        (NC,NC)      SCRATCH SPACE
C XS        (NC,NC)      SCRATCH SPACE
C P         (NC,NC)      2ND ORDER PROJECTION MATRIX
C P3        (NC,NC,NC)   3RD ORDER PROJECTION MARTIX
C F3        (NC,NC,NC)   V*ijk WITH TERMS CONTAINING Pij ONLY
C ON RETURN:
C F3        (NC,NC,NC)   COMPLETE 3RD ORDER FORCE CONSTANTS
C /////////////////////////////////////////////////////////////////////
 
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XS(NC,NC),P(NC,NC),F2(NC,NC),P3(NC,NC,NC),F3(NC,NC,NC)
      PARAMETER(ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0)
      COMMON /IO/ IIN1,IOUT,IIN2,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25,
     $   I31,I32,I33,I35,I36,I37,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10,
     $   ISCR11,ISCR12,ISCR13,ISCR14
  1   FORMAT(2I5)
  2   FORMAT(3F20.10)
      REWIND I15
      READ(I15,1) M,N
      READ(I15,2) ((F2(M,N),N=1,NC),M=1,NC)
      DO 10  M=2,NC
      DO 10  N=1,M-1
      XX=(F2(M,N)+F2(N,M))/TWO
      F2(M,N)=XX
 10   F2(N,M)=XX
 
      DO 25  I=1,NC
      DO 25  M=1,NC
      XS(I,M)=ZERO
      DO 25  L=1,NC
  25  XS(I,M)=XS(I,M)+F2(M,L)*P(L,I)
      DO 30  I=1,NC
      DO 30  J=1,NC
      DO 30  K=1,NC
      DO 30  M=1,NC
  30  F3(I,J,K)=F3(I,J,K)+XS(I,M)*P3(M,J,K)+XS(J,M)*P3(M,I,K)
     $          +XS(K,M)*P3(M,I,J)
CWA
      DO 150  I=1,NC
      WRITE(IOUT,*) 'I=',I
        DO 152  J=1,NC
  152   WRITE(IOUT,888) (P3(I,J,K),K=1,NC)
 150  CONTINUE
 888  FORMAT(9F10.6)
      RETURN
      END
C /////////////////////////////////////////////////////////////////////
      SUBROUTINE PROJK3(NAD,NA,NC,IOPT,XA,BETA,M2,B2,A,BB2,P2,P3,
     $                  G3,D,P4)
C /////////////////////////////////////////////////////////////////////
C COMPUTE THIRD  ORDER PROJECTION MATRIX TO REMOVE ROTATIONAL AND
C TRANSLATIONAL COMPONENTS IF THE FORCE CONSTANT MATRICES WERE
C COMPUTED AT A NON-STATIONARY POINT
C /////////////////////////////////////////////////////////////////////
C ON ENTRY:
C    NAD    INTEGER      NA+NDUM
C    NA     INTEGER      NUMBER OF ATOMS
C    NC     INTEGER      NA*3
C    XA     INTEGER      CARTESIAN COORS OF ATOMS
C    BETA   (3,NC,NC)    SCRATCH ARRAY
C    M2     SCRATCH ARRAY PASS IN SAME ADDRESS AS PASSED FOR BETA
C    B2     (6,NC)       B 1ST DERIVATIVES
C    A      (NC,6)       B2 "INVERSE"
C    BB2    (3,NC,NC)    B 2ND DERIVATIVES ROTATIONS ONLY
C    P2     (NC,NC)      2ND ORDER PROJECTION MATRIX
C    P3     (NC,NC,NC)   3RD ORDER PROJECTION MATRIX
C    G3     (NC,NC,NC)   Gijk
C    D      (NC,NC,NC)   SCRATCH SPACE
C ON RETURN:
C    P4     (NC,NC,NC,NC) 4TH ORDER PROJECTION MATRIX
C /////////////////////////////////////////////////////////////////////
 
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      DIMENSION XA(NAD,3),B2( 6,NC),A(NC,6),P2(NC,NC),P4(NC,NC,NC,NC)
      DIMENSION P3(NC,NC,NC), G3(NC,NC,NC),M2(NC,NC), BETA(3,NC,NC)
      DIMENSION BB2(3,NC,NC), W(6,12),IOPT(30), CM(3), D(NC,NC,NC)
      DIMENSION V1(3),V2(3),V3(3),T1(3),T2(3),T3(3), AI(3), AJ(3), AK(3)
      PARAMETER(ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0,THREE=3.0D0,FOUR=4.0D0)
      COMMON /IO/ IIN1,IOUT,IIN2,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25,
     $   I31,I32,I33,I35,I36,I37,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10,
     $   ISCR11,ISCR12,ISCR13,ISCR14  
      COMMON /PHYCON/ BOHR,DEBYE,HART,WAVE0,CINT
       DO 100 I = 1, NC
       DO 100 J = 1, NC
       DO 100 K = 1, NC
       DO 100 L = 1, NC
 100     P4(I,J,K,L) = ZERO
 
C
C 7TH TERM OF P4
C
 
C FORM BETAi,j'S
      DO 300 I = 1, NC
      DO 300 J = 1, NC
      DO 300 L = 1, 3
          BETA(L,I,J) = ZERO
          DO 300 M = 1, NC
          DO 300 N = 1, NC
 300          BETA(L,I,J) = BETA(L,I,J)+P2(M,I)*P2(N,J)*BB2(L,M,N)
 
C FORM (B*BT)-1
      DO 305  I=1,6
      DO 305  J=1,12
 305    W(I,J)=ZERO
      DO 307  I=1,6
 307    W(I,I+6)=ONE
      DO 309 I=1,6
      DO 309 J=1,6
      DO 309 K=1,NC
 309    W(I,J)=W(I,J)+B2(I,K)*B2(J,K)
      CALL FLIN(W,6,6,6,DET)
 
C MOVE IN 7TH TERM OF P4
       DO 310 L = 1, NC
       DO 310 I = 1, NC
       DO 310 J = 1, NC
       DO 310 K = 1, NC
 
       DO 310 N = 1, 3
       DO 310 M = 1, 3
 
           XX = -W(N,M+3)*BETA(N,I,L)*BETA(M,J,K)
           P4(L,I,J,K) = P4(L,I,J,K) + XX
           P4(L,J,I,K) = P4(L,J,I,K) + XX                               
 310       P4(L,K,I,J) = P4(L,K,I,J) + XX
 
C
C 5TH AND 6TH TERMS OF P4
C
 
C CALCULTE Dijk'S
       DO 312 I = 1, NC
       DO 312 J = 1, NC
       DO 312 K = 1, NC
           D(I,J,K) = G3(I,J,K)
           DO 312 L = 1, NC
  312          D(I,J,K) = D(I,J,K) + P2(L,I)*G3(K,L,J)
 
C SUM IN TERMS 5 AND 6 DOING NC**2 M'S AT A TIME
       DO 325 L = 1, NC
 
       DO 325 LL = 1, NC
       DO 325 MM = 1, NC
 
       DO 320 II = 1, NC
       DO 320 JJ = 1, NC
 320      M2(II,JJ) = D(L,LL,MM)*D(II,MM,JJ) -
     $                (P3(MM,II,JJ)+D(II,MM,JJ))*G3(L,LL,MM)
 
       DO 325 I = 1, NC
       DO 325 J = 1, NC
       DO 325 K = 1, NC
C TERM 5
           XX = P2(LL,I)*M2(J,K)
           P4(L,I,J,K) = P4(L,I,J,K) + XX                                
           P4(L,J,I,K) = P4(L,J,I,K) + XX                                
           P4(L,K,I,J) = P4(L,K,I,J) + XX                                
C TERM 6
           XX = -P2(LL,L)*G3(I,MM,K)*G3(J,LL,M)
           P4(L,I,J,K) = P4(L,I,J,K) + XX                                
           P4(L,I,K,J) = P4(L,I,K,J) + XX                                
 
 325   CONTINUE
 
 
C
C TERMS 1-4 OF P4
C
 
C FIND MOLECULE'S CENTER OF MASS
      TM=DBLE(NA)
      DO 332 K=1,3
         CM(K)=ZERO
         DO 333 I=1,NA
 333     CM(K)=CM(K)+XA(I,K)
 332     CM(K)=CM(K)/TM
 
C FORM AND INVERT THE INERTIA TENSOR
 
      DO 326  I=1,3
      DO 326  J=1,6
 326     W(I,J)=ZERO
      DO 327  I=1,3
 327     W(I,I+3)=ONE
 
      DO 328  K=1,NA
         DO 329  L=1,3
 329     V1(L)=XA(K,L)-CM(L)
         DO 330  L=2,3
         DO 330  M=1,L-1
 330     W(L,M)=W(L,M)-V1(L)*V1(M)
         W(1,1)=W(1,1)+V1(2)*V1(2)+V1(3)*V1(3)
         W(2,2)=W(2,2)+V1(1)*V1(1)+V1(3)*V1(3)
 328     W(3,3)=W(3,3)+V1(1)*V1(1)+V1(2)*V1(2)
      W(1,2)=W(2,1)
      W(1,3)=W(3,1)
      W(2,3)=W(3,2)
      CALL FLIN(W,3,3,3,DET)
 
C  ONE Gl,i,j,k WILL BE FORMED AT A TIME AND THEN APPLIED WHERE REQUIRED.
       DO 400 LL = 1, NC
       DO 400 II = 1, NC
       DO 400 JJ = 1, NC
       DO 400 KK = 1, NC
 
C CALC BBB2 DERIVATIVES IN RESPECT TO II, JJ, KK AND PUT IN V2
 
C      CARTESIAN COORD
       IIX = 1 +MOD(II-1,3)
       JJX = 1 +MOD(JJ-1,3)
       KKX = 1 +MOD(KK-1,3)
C      WHICH ATOM
       IIN = (II-IIX)/3 + 1
       JJN = (JJ-JJX)/3 + 1
       KKN = (KK-KKX)/3 + 1
 
       DO 335 K = 1, 3
           AI(K) = XA(IIN,K) - CM(K)
           AJ(K) = XA(JJN,K) - CM(K)
           AK(K) = XA(KKN,K) - CM(K)
           V1(K) = ZERO
           V2(K) = ZERO
 335       V3(K) = ZERO
 
           V1(IIX) = BB2(IIX,JJ,KK)
           V2(JJX) = BB2(JJX,II,KK)
           V3(KKX) = BB2(KKX,II,JJ)
 
       CALL SCAPRO(AI, V1, T1)
       CALL SCAPRO(AJ, V2, T2)
       CALL SCAPRO(AK, V3, T3)
 
       DO 337 I = 1, 3
 337       V1(I) = -T1(I) - T2(I) - T3(I)
 
       DO 339 I = 1, 3
           V2(I) = ZERO
       DO 339 J = 1, 3
 339       V2(I) = V2(I) + W(I,J+3)*V1(J)
 
 
C CALC G(LL,II,JJ,KK)
       G4 = ZERO
       DO 340 I = 1, 3
 340       G4 = G4 + A(LL,I+3)*V2(I)
 
 
C SUM IN GIVEN G4(LL,II,JJ,KK) IN TERMS AS REQUIRED:
 
C TERM 1
       DO 350 I = 1, NC
       DO 350 J = 1, NC
       DO 350 K = 1, NC
 350       P4(LL,I,J,K) = P4(LL,I,J,K) - P2(II,I)*P2(JJ,J)*P2(KK,K)*G4
 
C TERM 2
       DO 360 L = 1, NC
       DO 360 I = 1, NC
       DO 360 J = 1, NC
 360       P4(L,I,J,LL) = P4(L,I,J,LL) - P2(II,L)*P2(JJ,I)*P2(KK,J)*G4
 
C TERM 3
       DO 370 L = 1, NC
       DO 370 I = 1, NC
 370       P4(L,I,LL,II) = P4(L,I,LL,II) - P2(JJ,L)*P2(KK,I)*G4
 
C TERM 4
       DO 380 L = 1, NC
 380       P4(L,LL,II,JJ) = P4(L,LL,II,JJ) - P2(KK,L)*G4
 
 
C CYCLE TO NEXT G
 400   CONTINUE
 
 
 501   FORMAT(3F20.10)
       WRITE (IOUT,*) '4TH ORDER PROJECTION MATRIX'
       DO 500 I = 1, NC
       DO 500 J = 1, NC
            WRITE (IOUT,*)  'I=',I,'J=',J
       DO 499 K = 1, NC
       DO 499 L = 1, NC
 499        P4(I,J,K,L) =  P4(I,J,K,L)*BOHR
 500        WRITE (IOUT, 501)((P4(I,J,K,L),K=1,NC),L=1,NC)
       RETURN
       END
 
 
C /////////////////////////////////////////////////////////////////////
      SUBROUTINE PDER3(NC,IOPT,P2,P3,P4,F2,F3,F4P)
C /////////////////////////////////////////////////////////////////////
C ADD REMAINING TERMS TO V*ijkl WHICH CONTAIN Pijk AND Pijkl
C /////////////////////////////////////////////////////////////////////
C ON ENTRY:
C NC        INTEGER         NUM ATOMS * 3
C P2        (NC,NC)         2ND ORDER PROJECTION MATRIX
C P3        (NC,NC,NC)      3RD   "               "
C P4        (NC,NC,NC,NC)   4TH   "               "
C F2        (NC,NC)         SCRATCH SPACE
C F3        (NC,NC)         SCRATCH SPACE
C F4P       (NC,NC,NC,NC)   V*ijkl WITH Vijkl CONTRIBUTIONS ONLY
C ON EXIT:
C F4P       (NC,NC,NC,NC)   V*ijkl COMPLETE
C /////////////////////////////////////////////////////////////////////
 
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      INTEGER*4 P
      DIMENSION F2(NC,NC), P2(NC,NC),F3(NC,NC,NC),P3(NC,NC,NC)
      DIMENSION F4P(NC,NC,NC,NC), P4(NC,NC,NC,NC)
      PARAMETER(ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0)
      COMMON /IO/ IIN1,IOUT,IIN2,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25,
     $   I31,I32,I33,I35,I36,I37,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10,
     $   ISCR11,ISCR12,ISCR13,ISCR14 
  1   FORMAT (2I5)
  2   FORMAT (3F20.10)
 
C READ IN UNPROJECTED FORCE CONSTANTS
      REWIND I15
      READ (I15,1) M,N
      READ (I15,2) ((F2(I,J), J=1,NC), I=1,NC)
 
      REWIND I20
      READ (I20,1) L,M
      READ (I20,2) (((F3(I,J,K), K=1,J),  J=1,I),  I=1,NC)
 
C FILL OUT F3 MATRIX
      DO 5 I = 1, NC
      DO 5 J = 1, I-1
        F3(J,I,J)=F3(J,J,I)
        F3(I,J,J)=F3(J,J,I)
        F3(I,J,I)=F3(I,I,J)
        F3(J,I,I)=F3(I,I,J)
        DO 5 K = 1,J-1
          F3(I,K,J)=F3(I,J,K)
          F3(J,I,K)=F3(I,J,K)
          F3(J,K,I)=F3(I,J,K)
          F3(K,I,J)=F3(I,J,K)
          F3(K,J,I)=F3(I,J,K)
  5   CONTINUE
 
C ADD THE LAST TWO TERMS INTO V*
      DO 10 I= 1, NC
      DO 10 J= 1, NC
      DO 10 K= 1, NC
      DO 10 L= 1, NC
 
      DO 10 M= 1, NC
      DO 10 N= 1, NC
 
C         TERMS FROM 3RD SUMMATION
          XX = F2(M,N)*P2(M,I)*P4(N,J,K,L)
          F4P(I,J,K,L)=F4P(I,J,K,L)+XX
          F4P(J,I,K,L)=F4P(J,I,K,L)+XX
          F4P(K,I,J,L)=F4P(K,I,J,L)+XX
          F4P(L,I,J,K)=F4P(L,I,J,K)+XX
 
          XX = F2(M,N)*P3(M,I,J)*P3(N,K,L)
          F4P(I,J,K,L)=F4P(I,J,K,L)+XX
          F4P(I,K,J,L)=F4P(I,K,J,L)+XX
          F4P(I,L,J,K)=F4P(I,L,J,K)+XX
 
      DO 10  P = 1, NC
 
C         TERMS FROM 2ND SUMMATION
          XX =  F3(M,N,P)*P2(M,I)*P2(N,J)*P3(P,K,L)
          F4P(I,J,K,L)=F4P(I,J,K,L)+XX
          F4P(I,K,J,L)=F4P(I,K,J,L)+XX
          F4P(I,L,J,K)=F4P(I,L,J,K)+XX
          F4P(J,K,I,L)=F4P(J,K,I,L)+XX
          F4P(J,L,I,K)=F4P(J,L,I,K)+XX
          F4P(K,L,I,J)=F4P(K,L,I,J)+XX
 
   10 CONTINUE
      RETURN
      END
C ////////////////////////////////////////////////////////////////////////
       SUBROUTINE PRJTST (NC,IOPT,P2,P3,P4,F1,F2,F3,VP2,VP3,X)
C ////////////////////////////////////////////////////////////////////////
C VERIFY PROJECTION MATRICIES USING TRANSLATIONAL AND ROTATIONAL
C INVARIANCE RELTIONSHIPS.
C ////////////////////////////////////////////////////////////////////////
C ON ENTRY:
C NC       INTEGER          3 * NUMBER OF ATOMS
C P2       (NC,NC)          1ST ORDER PROJECTION MATRIX
C P3       (NC,NC,NC)       2ND    "            "
C P4       (NC,NC,NC,NC)    3RD    "            "
C F1       (NC)             SCRATCH SPACE
C F2       (NC,NC)             "      "
C F3       (NC,NC,NC)          "      "
C VP2      (NC,NC)             "      "
C VP3      (NC,NC,NC)          "      "
C X        (NC,NC)             "      "
C ////////////////////////////////////////////////////////////////////////
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      DIMENSION P2(NC,NC), P3(NC,NC,NC), P4(NC,NC,NC,NC)
      DIMENSION F2(NC,NC), F3(NC,NC,NC),X(NC,NC),VP2(NC,NC),F1(NC)
      DIMENSION VP3(NC,NC,NC),IOPT(30)
      LOGICAL FLAG1, FLAG2
      CHARACTER STRING*10
      COMMON /IO/ IIN1,IOUT,IIN2,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25,
     $   I31,I32,I33,I35,I36,I37,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10,
     $   ISCR11,ISCR12,ISCR13,ISCR14 
      PARAMETER (ERTOL=1.0D-5)
      PARAMETER (ZERO=0.0D0,TWO=2.0D0)
  1   FORMAT (2I5)
  2   FORMAT (3F20.10)
 
      WRITE (IOUT,4)
  4   FORMAT(/'TESTING OF PROJECTION MATRICIES USING ROTATIONAL ',
     $        'INVARIANCE RELATIONS')
 
      IF (IOPT(4).LT.2 .OR.IOPT(7).NE.3.OR.IOPT(5).EQ.0) THEN
            WRITE (IOUT,5)
  5         FORMAT ('   NO TESTING WAS PERFORMED: INCONSISTENT OPTIONS')
            RETURN
            END IF
 
      NA = NC/3
C
C TEST 2ND ORDER PROJECTION MATRIX
C THE FOLLOWING ASSUMES THAT THE 1ST ORDER PROJECTION MATRIX IS CORRECT
C
      WRITE(IOUT,9)
  9   FORMAT('   TESTING 2ND ORDER PROJECTION MATRIX, P(I,J,K):')
 
C READ IN F1 AND F2
      REWIND I11
 1000 READ (I11,1001,END=1002) STRING
 1001 FORMAT (A10)
      GOTO 1000
 1002 DO 1010 K = 1, NA+1
 1010    BACKSPACE I11
      DO 1030 J = 1, NA
 1030   READ (I11,1031) (F1((J-1)*3+K),K=1,3)
 1031 FORMAT (20X,3F20.10)
 
      REWIND I15
      READ(I15,1) M,N
      READ(I15,2) ((F2(M,N), N=1,NC), M=1,NC)
      DO 20 M=2,NC
      DO 20 N=1,M-1
        XX = (F2(M,N)+F2(N,M))/TWO
        F2(M,N)=XX
  20    F2(N,M)=XX
 
C COMPUTE Vij*'s
 
      DO 30 I = 1, NC
      DO 30 J = 1, I
        DO 27 K = 1, NC
        DO 27 L = 1, NC
   27     VP2(I,J) = VP2(I,J) + P2(I,K)*F2(K,L)*P2(L,J)
   30     VP2(J,I) = VP2(I,J)
 
      FLAG1 = .TRUE.
      DO 50 I = 1, NC
      DO 50 J = 1, NC
 
        XX = ZERO
        DO 40 K = 1, NC
   40        XX = XX + F1(K)*P3(K,I,J)
 
        ERR = F2(I,J) - VP2(I,J) - XX
        IF (ABS(ERR).GT.ERRTOL) THEN
             FLAG1 = .FALSE.
             WRITE (IOUT,42) I,J,ERR
   42        FORMAT ('      I= ?  J= ',I3,' K= ',I3,' ERROR= ',F12.8)
             END IF
 
   50   CONTINUE
 
        IF (FLAG1) WRITE (IOUT,54)
   54              FORMAT ('      NO ERRORS DETECTED.')
 
 
 
      IF (IOPT(4).LE.2) RETURN
C
C TEST 3RD ORDER PROJECTION MATRIX
C RESULTS WILL BE MEANINGFUL DEPENDING UPON ERROR IN THE
C 2ND ORDER PROJECTION MATRIX.
C
      WRITE (IOUT,109)
 109  FORMAT ('   TESTING 3ND ORDER PROJECTION MATRIX, P(I,J,K,L):')
      IF (FLAG1) WRITE (IOUT,110)
 110  FORMAT ('      RESULTS MAY BE UNRELIABLE DUE TO ERRORS IN '/      .')
     $        'P(I,J,K).')
 
C READ IN F3
 
      REWIND I20
      READ(I20, 1) M,N
      READ(I20, 2) (((F3(I,J,K), K=1,J), J=1,I), I=1,NC)
      DO 120 I = 1, NC
      DO 120 J = 1, I-1
        F3(J,I,J) = F3(J,J,I)
        F3(I,I,J) = F3(J,J,I)
        F3(I,J,I) = F3(I,I,J)
        F3(J,I,I) = F3(I,I,J)
        DO 120 K = 1, J-1
          F3(I,K,J)=F3(I,J,K)
          F3(J,I,J)=F3(I,J,K)
          F3(J,K,I)=F3(I,J,K)
          F3(K,I,J)=F3(I,J,K)
          F3(K,J,I)=F3(I,J,K)
  120     CONTINUE
 
C COMPUTE Vijk*'s
 
      DO 150 I = 1, NC
      DO 150 M = 1, NC
      X(I,M)=ZERO
      DO 150 L = 1, NC
 150    X(I,M) = X(I,M) + F2(M,L)*P2(L,I)
 
      DO 160 I = 1, NC
      DO 160 J = 1, NC
      DO 160 K = 1, NC
      DO 160 M = 1, NC
 160    VP3(I,J,K) = VP3(I,J,K) + X(I,M)*P3(M,J,K) + X(J,M)*P3(M,I,K) +
     $              X(K,M)*P3(M,I,J)
 
 
      FLAG2 = .TRUE.
      DO 180 I = 1, NC
      DO 180 J = 1, NC
      DO 180 K = 1, NC
 
        XX = ZERO
        DO 170 L = 1, NC
  170        XX = XX + F1(K)*P4(L,I,J,K)
 
        ERR = F3(I,J,K) - VP3(I,J,K) - XX
        IF (ABS(ERR).GT.ERRTOL) THEN
             FLAG2= .FALSE.
             WRITE (IOUT,172) I,J,K,ERR
  172        FORMAT ('      I= ?  J= ',I3,' K= ',I3,' L= ',I3,          ,F12.8)
     $               'ERROR=',F12.8)
             END IF
 
  180   CONTINUE
 
C       MESSAGE IF NO ERRORS
        IF (FLAG2) WRITE (IOUT,54)
 
        RETURN
 
        END
C/////////////////////////////////////////////////////////////////////
      SUBROUTINE EQ42 (NC,NA,XA,P2,P3,P4,XAS)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      DIMENSION P2(NC,NC),P3(NC,NC,NC),P4(NC,NC,NC,NC)
      DIMENSION XA(NA,3), XAS(NC), X(3)
      COMMON /IO/ IIN1,IOUT,IIN2,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25,
     $   I31,I32,I33,I35,I36,I37,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10,
     $   ISCR11,ISCR12,ISCR13,ISCR14 
      PARAMETER (ZERO=0.0D0)
      WRITE (IOUT,1)
   1  FORMAT ('APPROX PROJECTED CART COORD DISP VIA EQ42')
 
      DO 2000 I = 1,NC
 2000    XAS(I) = ZERO
 
      XAS(7)=.1091033790
      XAS(8)=.1091033790
      XAS(9)=.1091033790
 
C     DO 10 I = 1, NA
C     DO 10 J = 1, 3
C          K = (I-1)*3+J
C 10       XAS(K) = XA(I,J)
 
      DO 100 I = 1, NA
      DO  90 J = 1, 3
        K = (I-1)*3+J
        T1 = ZERO
        T2 = ZERO
        T3 = ZERO
 
        DO 50 L = 1, NC
          T1 = T1 + P2(K,L)*XAS(L)
 
          DO 50 M = 1, NC
            T2 = T2 + P3(K,L,M)*XAS(L)*XAS(M)
 
             DO 50 N = 1, NC
               T3 = T3 + P4(K,L,M,N)*XAS(L)*XAS(M)*XAS(N)
 
  50    CONTINUE
 
  90    X(J) = T1 + T2/2 + T3/6
 
 100  WRITE (IOUT,101) I,X(1),X(2),X(3)
 101  FORMAT (I3,3F20.10)
 
      RETURN
      END
C//////////////////////////////////////////////////////////////////////
	SUBROUTINE ORTHOG (NA,NC,NS,IOPT,XAR,BS1,BS2,BS3,FLAG)
C  \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C  TEST ORTHOGINALITY CONDITIONS CONCERNING EXTERNAL COORDS AND CART. 
C  DERIVATIVES OF INTERNAL COORDS. 
C  \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C  ON ENTRY:
C	NA			NUMBER OF ATOMS
C	NC			NA*3
C	NS			NUMBER OF SYMETRIZED INTERNAL COORDS
C	XAR	(NA,3)		REFERENCE GEOMETRY
C	BS1	(NC,NC)		1ST CART DERIVATIVES OF INTERNAL COORDS
C	BS2	(NC,NC)		SCRATCH SPACE (USED IFF NDER>1)
C	BS3	(NC,NC,NC)	SCRATCH SPACE (USED IFF NDER>2)
C  ON RETURN
C	FLAG	LOGICAL	SET TRUE IF ANY ORTHOG. CONDITIONS FAILS
C  NOTES:
C	1) HIGHEST CONDITION TESTED DETERMINED BY IOPT(4)
C  \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (TOL = 1.0D-8 )
      PARAMETER (ZERO = 0.0D0,ONE=1.0D0)
	LOGICAL FLAG, NOERR, RFLAG
	DIMENSION XAR(NA,3)
	DIMENSION BS1(NC,NC), BS2(NC,NC), BS3(NC,NC,NC), IOPT(30)
	DIMENSION SUMV(3), V(3), V2(3), V3(3), V5(3), E(3)
      COMMON /IO/ IIN1,IOUT,IIN2,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25,
     $   I31,I32,I33,I35,I36,I37,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10,
     $   ISCR11,ISCR12,ISCR13,ISCR14
  1	FORMAT (//' TESTING ORTHOGONALITY CONDITIONS:')
  2	FORMAT ('     TOLERANCE IS SET TO', E10.2)
  3	FORMAT (/'     1ST CONDITION HOLDS FOR ALL CASES.')
  4	FORMAT (/'     2ND CONDITION HOLDS FOR ALL CASES.')
  5	FORMAT (/'     3RD CONDITION HOLDS FOR ALL CASES.')
  6	FORMAT (/'     1ST CONDITION FAILS FOR:')
  7	FORMAT (/'     2ND CONDITION FAILS FOR:')
  8	FORMAT (/'     3RD CONDITION FAILS FOR:')
  9	FORMAT (12X,'P = ',I5,2X,'ETA = ',I5)
 10	FORMAT (12X,'P = ',I5,2X,'ETA = ',I5,2X,'J2 = ',I5)
 11	FORMAT (12X,'P = ',I5,2X,'ETA = ',I5,2X,'J2 = ',I5,2X,'J3 = ',I5)

	FLAG = .FALSE.
	IF (IOPT(4).LT.1) RETURN
	WRITE (IOUT, 1)
	WRITE (IOUT, 2) TOL

        MDER = IOPT(4)+IOPT(5)

C FIRST CONDITION TRANSLATIONS
	NOERR = .TRUE.
	DO 100 IP = 1, NS
	DO 100 I = 1, 3
	SUM = ZERO
	DO 110 J = I, NC, 3
 110	SUM = SUM - BS1(IP,J)
	IF (ABS(SUM).GT.TOL) THEN
		IF (NOERR) WRITE (IOUT,6)
		NOERR = .FALSE.
		WRITE (IOUT,9) IP, I
		END IF
 100	CONTINUE

C FIRST CONDITION ROTATIONS
	DO 160 IP = 1, NS
	DO 163 I = 1, 3  
 163	SUMV(I) = ZERO
	
	DO 150 JA = 1, NA
	V(1) = XAR(JA,1)
	V(2) = XAR(JA,2)
	V(3) = XAR(JA,3)
	DO 150 JB = 1, 3
	J = (JA-1)*3 + JB
	E(1) = ZERO
	E(2) = ZERO
	E(3) = ZERO
	E(JB) = ONE
	CALL VECPRO (E, V, V2)
	DO 150 I = 1, 3 
 150	SUMV(I) = SUMV(I) + V2(I)*BS1(IP,J)

	DO 160 I = 1, 3 
	IF (ABS(SUMV(I)).GT.TOL) THEN
		IF (NOERR) WRITE (IOUT,6)
		NOERR = .FALSE.
		WRITE (IOUT,9) IP, I+3
		END IF
 160	CONTINUE
	IF (NOERR) WRITE (IOUT,3)
	IF (MDER.LE.2) RETURN
C SECOND CONDITION TRANSLATIONS
	NOERR = .TRUE.
	DO 200 IP = 1, NS
	CALL XIN (NC,0,BS2,-IP,ISCR1)
	DO 200 K = 1, NC
	DO 200 JB = 1, 3
	SUM = ZERO
	DO 210 J = JB, NC, 3
 210	SUM = SUM - BS2(J,K)
	IF (ABS(SUM).GT.TOL) THEN
		IF (NOERR) WRITE (IOUT,7)
		NOERR = .FALSE.
		WRITE (IOUT,10) IP, JB, K
		END IF
 200	CONTINUE

C SECOND CONDITION ROTATIONS
	DO 260 IP = 1, NS
	CALL XIN (NC,0,BS2,-IP,ISCR1)
	DO 260 K = 1, NC
	KN = (K-1)/3+1
	KB = K-(KN-1)*3
	V3(1) = ZERO
	V3(2) = ZERO
	V3(3) = ZERO
	V3(KB) = ONE

	SUMV(1) = ZERO
	SUMV(2) = ZERO
	SUMV(3) = ZERO
	DO 250 JA = 1, NA
	V(1) = XAR(JA,1)
	V(2) = XAR(JA,2)
	V(3) = XAR(JA,3)
	DO 250 JB = 1, 3
	J = (JA-1)*3+JB
	E(1) = ZERO
	E(2) = ZERO
	E(3) = ZERO
	E(JB) = ONE
	CALL VECPRO (E, V, V2)
	DO 252 I = 1, 3 
 252	SUMV(I) = SUMV(I) + V2(I)*BS2(J,K)
	IF (JA.EQ.KN) THEN
		CALL VECPRO (E, V3, V2)
		DO 254 I = 1, 3 
 254		SUMV(I) = SUMV(I) + V2(I)*BS1(IP,J)
		END IF
 250	CONTINUE

	DO 260 I = 1, 3 
	IF (ABS(SUMV(I)).GT.TOL) THEN
		IF (NOERR) WRITE (IOUT,7)
		NOERR = .FALSE.
		WRITE (IOUT,10) IP, I+3, K
		END IF
 260	CONTINUE
	IF (NOERR) WRITE (IOUT,4)
	IF (MDER.LE.3) RETURN
C THIRD CONDITION TRANSLATIONS
	NOERR = .TRUE.
	DO 300 IP = 1, NS
	CALL YIN (NC,0,BS3,-IP,ISCR3)

	DO 300 K = 1, NC
	DO 300 L = 1, NC
	DO 300 JB = 1, 3
	SUM = ZERO
	DO 310 J = JB, NC, 3
 310	SUM = SUM - BS3(J,K,L)
	IF (ABS(SUM).GT.TOL) THEN
		IF (NOERR) WRITE (IOUT,8)
		NOERR = .FALSE.
		WRITE (IOUT,11) IP, JB, L, K
		END IF
 300	CONTINUE

C THIRD CONDITION ROTATIONS

	DO 360 IP = 1, NS
	CALL XIN (NC,0,BS2,-IP,ISCR1)
	CALL YIN (NC,0,BS3,-IP,ISCR3)

	DO 360 K = 1, NC
	KN = (K-1)/3+1
	KB = K-(KN-1)*3
	V3(1) = ZERO
	V3(2) = ZERO
	V3(3) = ZERO
	V3(KB) = ONE
	DO 360 L = 1, NC
	LN = (L-1)/3+1
	LB = L-(LN-1)*3
	V5(1) = ZERO
	V5(2) = ZERO
	V5(3) = ZERO
	V5(LB) = ONE

	SUMV(1) = ZERO
	SUMV(2) = ZERO
	SUMV(3) = ZERO
	DO 350 JA = 1, NA
	V(1) = XAR(JA,1)
	V(2) = XAR(JA,2)
	V(3) = XAR(JA,3)
	DO 350 JB = 1, 3
	J = (JA-1)*3+JB
	E(1) = ZERO
	E(2) = ZERO
	E(3) = ZERO
	E(JB) = ONE
	CALL VECPRO (E, V, V2)
	DO 352 I = 1, 3  
 352	SUMV(I) = SUMV(I) + V2(I)*BS3(J,K,L)
	IF (JA.EQ.KN) THEN
		CALL VECPRO (E, V3, V2)
		DO 354 I = 1, 3 
 354		SUMV(I) = SUMV(I) + V2(I)*BS2(J,L)
		END IF
	IF (JA.EQ.LN) THEN
		CALL VECPRO (E, V5, V2)
		DO 356 I = 1, 3  
 356		SUMV(I) = SUMV(I) + V2(I)*BS2(J,K)
		END IF
 350	CONTINUE
	DO 360 I = 1, 3 
	IF (ABS(SUMV(I)).GT.TOL) THEN
		IF (NOERR) WRITE (IOUT,8)
		NOERR = .FALSE.
		WRITE (IOUT,11) IP, I+3, L, K
		END IF	
 360	CONTINUE
	IF (NOERR) WRITE (IOUT,5)

	RETURN
	END
C//////////////////////////////////////////////////////////////////////////
        SUBROUTINE FORMDK0 (NA,NC,XA,XMASS,RFLAG,IOPT, DOT,CP1,DOT2,CP2,
     $                     EFLAG,DK1,DK2,DK3)
C  \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C  COMPUTE THE DERIVATIVES OF THE EXTERNAL COORDS OF THE MOL IN RESPECT
C  TO THE CARTESIAN COORDINATES OF THE ATOMS AT K=0 (IE ALL EULERIAN ANGLES
C  ZERO AND THE CENTER OF MASS AT THE ORIGIN.
C  \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C  ON ENTRY:
C       NA                      NUMBER OF ATOMS
C       NC                      NA*3
C       XA      (NC)            REFERENCE GEOMETRY (ANG.)
C       XMASS   (NA)            MASS OF ATOMS
C	RFLAG	LOGICAL		SET TRUE IF LINEAR MOLECULE
C       DOT     (NA,NC)         SCRATCH SPACE (USED ONLY IF NDER > 1)
C       CP1     (3,NA,NC)       SCRATCH SPACE (	"		"   )
C	DOT2	(NA,NC,NC)	SCRATCH SPACE (USED ONLY IF NDER > 2)
C	CP2	(3,NA,NC,NC)	SCRATCH SPACE (	"		"   )
C  ON RETURN:
C       EFLAG   LOGICAL         SET TRUE ON FATAL ERROR
C       DK1     (6,NC)          1ST DERVIVATIVES
C       DK2     (3,NC,NC)       2ND DERVIVATIVES ROTATIONAL VARIABLES ONLY
C       DK3     (3,NC,NC,NC)    3RD DERVIVATIVES ROTATIONAL VARIABLES ONLY
C  ERROR CONDITIONS:
C	INERTIA TENSOR IS SINGULAR THUS REF GEO OR MASSES NOT VALID
C  NOTES:
C       1) HIGHEST DERIV FORMED DETERMINED BY IOPT(4)
C	2) DERIVAITVE PRINTED GIVEN IOPT(6)
C	3) DERIVATIVES SAVED TO DISK IF IOPT(14)=2
C  \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        IMPLICIT REAL*8 (A-H,O-Z)
        LOGICAL EFLAG, RFLAG
        DIMENSION XA(NA,3),XMASS(NA),DOT(NA,NC),CP1(3,NA,NC),IOPT(30)
        DIMENSION DK1(6,NC), DK2(3,NC,NC), DK3(3,NC,NC,NC)
        DIMENSION V(3), V1(3), T1(3,6), IPRMN(3), IPRMX(3)
	DIMENSION DOT2(NA,NC,NC), CP2(3,NA,NC,NC), EU(3), V2(3)
	DIMENSION IDX1(6), IDX2(6), IDX3(6)

        EQUIVALENCE (IL, IPRMN(1)), (IM, IPRMN(2)), (IN, IPRMN(3))
        EQUIVALENCE (IB, IPRMX(1)), (IG, IPRMX(2)), (ID, IPRMX(3))
	DATA (IDX1(J), J=1,6) /1, 1, 2, 2, 3, 3 /
        DATA (IDX2(J), J=1,6) /2, 3, 1, 3, 1, 2 /
        DATA (IDX3(J), J=1,6) /3, 2, 3, 1, 2, 1 /

	PARAMETER (ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0, THREE=3.0D0)
	PARAMETER (FOUR=4.0D0)
      COMMON /IO/ IIN1,IOUT,IIN2,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25,
     $   I31,I32,I33,I35,I36,I37,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10,
     $   ISCR11,ISCR12,ISCR13,ISCR14
  1     FORMAT(///'ERROR IN FORMDK0: INERTIA TENSOR SINGULAR:'
     $           /'           REFERENCE GEOMETRY OR MASSES NOT VALID.')
  2	FORMAT (3F20.10)
  3	FORMAT (6F12.6)
  4	FORMAT (/' DK(I) MATRIX FOR TRANSLATIONAL COORDINATE',I5/)
  5	FORMAT (/' DK(I) MATRIX FOR ROTATIONAL COORDINATE',I5/)
  6	FORMAT (/' DK(I,J) MATRIX FOR ROTATIONAL COORDINATE',I5/)
  7	FORMAT (/' DK(I,J,K) MATRIX FOR ROTATIONAL COORDINATE',I5)
  8	FORMAT (/' I=',I5/)
 
        EFLAG = .FALSE.

C FORM AND INVERT THE INERTIA TENSOR
        DO 62 I = 1, 3
        DO 62 J = 1, 6
  62    T1(I,J) = ZERO

	IF (RFLAG) THEN
	TM = ZERO
	SUM = ZERO
	DO 40 I = 1, NA
	TM = TM + XMASS(I)
  40	SUM = SUM + XMASS(I)*XA(I,3)**2
	SUM = ONE / SUM
	T1(1,4) = SUM
	T1(2,5) = SUM

	ELSE
        DO 63  I = 1, 3
  63    T1(I,I+3) = ONE

	TM = ZERO
        DO 70  K=1,NA
	TM = TM + XMASS(K)
        T1(1,1)=T1(1,1)+(XA(K,2)**2+XA(K,3)**2)*XMASS(K)
        T1(2,2)=T1(2,2)+(XA(K,1)**2+XA(K,3)**2)*XMASS(K)
        T1(3,3)=T1(3,3)+(XA(K,1)**2+XA(K,2)**2)*XMASS(K)
        T1(2,1)=T1(2,1)-XA(K,1)*XA(K,2)*XMASS(K)
        T1(3,1)=T1(3,1)-XA(K,1)*XA(K,3)*XMASS(K)
  70    T1(3,2)=T1(3,2)-XA(K,2)*XA(K,3)*XMASS(K)
        T1(1,2)=T1(2,1)
        T1(1,3)=T1(3,1)
        T1(2,3)=T1(3,2)
        CALL FLIN(T1,3,3,3,DET)
        IF (DET.EQ.ZERO) THEN
                WRITE (IOUT,1)
                EFLAG = .TRUE.
                RETURN
                END IF
	END IF

C  FIRST DERIVATIVES (TRANSLATIONAL VARIABLES)
	DO 49 I = 1, 3
	DO 49 K = 1, NC
  49	DK1(I,K) = ZERO
        DO 50 K = 1, NA
        XXX = XMASS(K) / TM
        DO 50 I = 1, 3
        L = 3*(K-1)+I
  50    DK1(I,L) = XXX

C  FIRST DERIVATIVES (ROTATIONAL VARIABLES)
        DO 10 I = 1, NA
        IX = (I-1)*3 + 1
        IY = IX + 1
        IZ = IY + 1
        DK1(4,IX) = ( T1(1,5)*XA(I,3) - T1(1,6)*XA(I,2)) * XMASS(I)
        DK1(5,IX) = ( T1(2,5)*XA(I,3) - T1(2,6)*XA(I,2)) * XMASS(I)
        DK1(6,IX) = ( T1(3,5)*XA(I,3) - T1(3,6)*XA(I,2)) * XMASS(I)
        DK1(4,IY) = (-T1(1,4)*XA(I,3) + T1(1,6)*XA(I,1)) * XMASS(I)
        DK1(5,IY) = (-T1(2,4)*XA(I,3) + T1(2,6)*XA(I,1)) * XMASS(I)
        DK1(6,IY) = (-T1(3,4)*XA(I,3) + T1(3,6)*XA(I,1)) * XMASS(I)
        DK1(4,IZ) = ( T1(1,4)*XA(I,2) - T1(1,5)*XA(I,1)) * XMASS(I)
        DK1(5,IZ) = ( T1(2,4)*XA(I,2) - T1(2,5)*XA(I,1)) * XMASS(I)
        DK1(6,IZ) = ( T1(3,4)*XA(I,2) - T1(3,5)*XA(I,1)) * XMASS(I)
10      CONTINUE

C	SAVE TO DISK IF NSTOP SET
	IF (IOPT(14).EQ.2) THEN
		REWIND I31
		WRITE (I31,2) ((DK1(I,J), J=1,NC), I=1,6)
		END IF

C	PRINT IF SPECIFIED BY PRINT OPTION
	IF (LPRT(1,IOPT(6)).GE.3) THEN
		DO 20 I = 1, 3
		WRITE (IOUT,4) I
  20		WRITE (IOUT,3) (DK1(I,J), J=1,NC)
		DO 21 I = 4, 6
		WRITE (IOUT,5) I-3
  21		WRITE (IOUT,3) (DK1(I,J), J=1,NC)
		END IF
        
	IF (IOPT(4).LE.2) RETURN

C  FORM (A . KIB)*MI AND (A X KIB) QUANTITIES
        DO 100 I = 1, NA
	DO 100 J = 1, NC
        DOT(I,J) = (XA(I,1)*DK1(4,J)+XA(I,2)*DK1(5,J)+
     $                     XA(I,3)*DK1(6,J)) * XMASS(I)
        CP1(1,I,J) = XA(I,2)*DK1(6,J) - XA(I,3)*DK1(5,J)
        CP1(2,I,J) = XA(I,3)*DK1(4,J) - XA(I,1)*DK1(6,J)
        CP1(3,I,J) = XA(I,1)*DK1(5,J) - XA(I,2)*DK1(4,J)
  100   CONTINUE

C  SECOND DERIVATIVES
        DO 125 IN = 1, NA
        DO 125 IB = 1, 3
        I = (IN-1)*3+IB
        DO 125 JN = 1, IN
        DO 125 JG = 1, 3
        J = (JN-1)*3+JG

	V(1) = ZERO
	V(2) = ZERO
	V(3) = ZERO
        DO 120 K = 1, NA
        V(1) = V(1) + DOT(K,I)*CP1(1,K,J)+DOT(K,J)*CP1(1,K,I)
        V(2) = V(2) + DOT(K,I)*CP1(2,K,J)+DOT(K,J)*CP1(2,K,I)
        V(3) = V(3) + DOT(K,I)*CP1(3,K,J)+DOT(K,J)*CP1(3,K,I)
  120   CONTINUE

        XXX = -XMASS(IN)*XA(IN,IB)
        YYY = -XMASS(JN)*XA(JN,JG)
        V(1) = V(1)/TWO + XXX*DK1(4,J)+YYY*DK1(4,I)
        V(2) = V(2)/TWO + XXX*DK1(5,J)+YYY*DK1(5,I)
        V(3) = V(3)/TWO + XXX*DK1(6,J)+YYY*DK1(6,I)

        V(IB) = V(IB) + DOT(IN,J)
        V(JG) = V(JG) + DOT(JN,I)

        DK2(1,I,J) = V(1)*T1(1,4)+V(2)*T1(1,5)+V(3)*T1(1,6)
        DK2(2,I,J) = V(1)*T1(2,4)+V(2)*T1(2,5)+V(3)*T1(2,6)
        DK2(3,I,J) = V(1)*T1(3,4)+V(2)*T1(3,5)+V(3)*T1(3,6)
        DK2(1,J,I) = DK2(1,I,J)
        DK2(2,J,I) = DK2(2,I,J)
        DK2(3,J,I) = DK2(3,I,J)

  125   CONTINUE

C       SAVE TO DISK IF NSTOP NONZERO
        IF (IOPT(14).EQ.2) THEN
                REWIND I32
                WRITE (I32,2) (((DK2(I,J,K), K=1,J), J=1,NC), I=1,3)
                END IF

C       PRINT IF SPECIFIED BY THE PRINT OPTION
        IF (LPRT(1,IOPT(6)).GE.3) THEN
                DO 130 I = 1, 3
                WRITE (IOUT,6) I
		DO 130 J = 1, NC
 130            WRITE (IOUT,3) (DK2(I,J,K), K=1,NC)
                END IF

        IF (IOPT(4).LE.3) RETURN

C FORM (A . KQV,RS)*MP AND (A X KQV,RS) TERMS
        DO 400 I = 1, NA
        DO 400 J = 1, NC
        DO 400 K = 1, J

        DOT2(I,J,K) = ( XA(I,1)*DK2(1,J,K) + XA(I,2)*DK2(2,J,K) +
     $                 XA(I,3)*DK2(3,J,K) ) * XMASS(I)
        DOT2(I,K,J) = DOT2(I,J,K)

        CP2(1,I,J,K) = XA(I,2)*DK2(3,J,K) - XA(I,3)*DK2(2,J,K)
        CP2(2,I,J,K) = XA(I,3)*DK2(1,J,K) - XA(I,1)*DK2(3,J,K)
        CP2(3,I,J,K) = XA(I,1)*DK2(2,J,K) - XA(I,2)*DK2(1,J,K)
        CP2(1,I,K,J) = CP2(1,I,J,K)
        CP2(2,I,K,J) = CP2(2,I,J,K)
        CP2(3,I,K,J) = CP2(3,I,J,K)
 400    CONTINUE

C THIRD DERIVATIVES
        DO 500 IL = 1, NA
        DO 500 IB = 1, 3
	ILB = (IL-1)*3+IB
        DO 500 IM = 1, IL
        DO 500 IG = 1, 3
	IMG = (IM-1)*3+IG
        DO 500 IN = 1, IM
        DO 500 ID = 1, 3
        IND = (IN-1)*3+ID

        V1(1) = ZERO
        V1(2) = ZERO
        V1(3) = ZERO

C       LOOP OVER PERMUTATIONS
        DO 600 K = 1, 6
C       NOTE IPRMN(1) AND IL, IPRMX(1) AND IB, ETC HAVE SAME ADDRESS
C	IDX VECTORS ARE PRESET IN A DATA STATEMENT WITH THE PREMUTATIONS
        IP = IPRMN(IDX1(K))
        IU = IPRMX(IDX1(K))
        IPU = (IP-1)*3+IU
        IR = IPRMN(IDX2(K))
        IS = IPRMX(IDX2(K))
        IRS = (IR-1)*3+IS
        IQ = IPRMN(IDX3(K))
        IV = IPRMX(IDX3(K))
        IQV = (IQ-1)*3+IV
       
	XXX = XMASS(IP)/TWO
	V(1)=(-XA(IP,IU)*DK2(1,IQV,IRS)+DK1(IU+3,IQV)*CP1(1,IP,IRS))*XXX
        V(2)=(-XA(IP,IU)*DK2(2,IQV,IRS)+DK1(IU+3,IQV)*CP1(2,IP,IRS))*XXX
        V(3)=(-XA(IP,IU)*DK2(3,IQV,IRS)+DK1(IU+3,IQV)*CP1(3,IP,IRS))*XXX

        V(IU) = V(IU) + DOT2(IP,IQV,IRS)/TWO

        EU(1) = ZERO
        EU(2) = ZERO
        EU(3) = ZERO
        EU(IU) = ONE
        XXX = ( DK1(4,IQV)*DK1(4,IRS) + DK1(5,IQV)*DK1(5,IRS) +
     $          DK1(6,IQV)*DK1(6,IRS) ) * XMASS(IP) / THREE
        V(1) = V(1) - XXX * (XA(IP,2)*EU(3) - XA(IP,3)*EU(2))
        V(2) = V(2) - XXX * (XA(IP,3)*EU(1) - XA(IP,1)*EU(3))
        V(3) = V(3) - XXX * (XA(IP,1)*EU(2) - XA(IP,2)*EU(1))

        V2(1) = ZERO
        V2(2) = ZERO
        V2(3) = ZERO
        DO 550 J = 1, NA
        V2(1) = V2(1) + DOT2(J,IQV,IRS)*CP1(1,J,IPU) +
     $                  DOT(J,IPU)*CP2(1,J,IQV,IRS)
        V2(2) = V2(2) + DOT2(J,IQV,IRS)*CP1(2,J,IPU) +
     $                  DOT(J,IPU)*CP2(2,J,IQV,IRS)
        V2(3) = V2(3) + DOT2(J,IQV,IRS)*CP1(3,J,IPU) +
     $                  DOT(J,IPU)*CP2(3,J,IQV,IRS)
 550    CONTINUE
        V1(1) = V1(1) + V(1) + V2(1)/FOUR
        V1(2) = V1(2) + V(2) + V2(2)/FOUR
        V1(3) = V1(3) + V(3) + V2(3)/FOUR

 600    CONTINUE

        V2(1) = V1(1)*T1(1,4)+V1(2)*T1(1,5)+V1(3)*T1(1,6)
        V2(2) = V1(1)*T1(2,4)+V1(2)*T1(2,5)+V1(3)*T1(2,6)
        V2(3) = V1(1)*T1(3,4)+V1(2)*T1(3,5)+V1(3)*T1(3,6)
        DK3(1,ILB,IMG,IND) = V2(1)
        DK3(2,ILB,IMG,IND) = V2(2)
        DK3(3,ILB,IMG,IND) = V2(3)
        DK3(1,ILB,IND,IMG) = V2(1)
        DK3(2,ILB,IND,IMG) = V2(2)
        DK3(3,ILB,IND,IMG) = V2(3)
        DK3(1,IMG,ILB,IND) = V2(1)
        DK3(2,IMG,ILB,IND) = V2(2)
        DK3(3,IMG,ILB,IND) = V2(3)
        DK3(1,IMG,IND,ILB) = V2(1)
        DK3(2,IMG,IND,ILB) = V2(2)
        DK3(3,IMG,IND,ILB) = V2(3)
        DK3(1,IND,ILB,IMG) = V2(1)
        DK3(2,IND,ILB,IMG) = V2(2)
        DK3(3,IND,ILB,IMG) = V2(3)
        DK3(1,IND,IMG,ILB) = V2(1)
        DK3(2,IND,IMG,ILB) = V2(2)
        DK3(3,IND,IMG,ILB) = V2(3)

 500    CONTINUE

C       SAVE TO DISK IF NSTOP NONZERO
        IF (IOPT(14).EQ.2) THEN
                REWIND I33
                WRITE (I33,2)
     $			 ((((DK3(I,J,K,L),L=1,K),K=1,J),J=1,NC),I=1,3)
                END IF

C       PRINT IF SPECIFIED BY PRINT OPTION
        IF (LPRT(1,IOPT(6)).GE.3) THEN
                DO 510 I = 1, 3
                WRITE (IOUT,7) I
		DO 510 J = 1, NC
		WRITE (IOUT,8) J
		DO 510 K = 1, NC
  510           WRITE (IOUT,3) (DK3(I,J,K,L), L=1,NC)
                END IF
	RETURN
	END
C////////////////////////////////////////////////////////////////////// INT47140
      SUBROUTINE GFMAT(NA,NAD,NC,NSX,NFREQ,IRINT,NMODE,IA,              INT47150
     $                       XMASS,XA,W,BS,F,XR,XS,XT,XU,IFLAG)         INT47160
      IMPLICIT REAL*8 (A-H,O-Z)                                         INT47170
      LOGICAL ITEST                                                     INT47180
      DIMENSION XA(NAD,3),BS(NC,NC),F(NC,NC),XS(NC,NC),XT(NC,NC)        INT47190
      DIMENSION XMASS(NA),W(NSX),XR(NC,1),XU(1),IA(5,NSX)               INT47200
      COMMON /IO/ IIN1,IOUT,IIN2,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25,
     $   I31,I32,I33,I35,I36,I37,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10,
     $   ISCR11,ISCR12,ISCR13,ISCR14
      COMMON /PHYCON/ BOHR,DEBYE,HART,WAVE0,CINT
      PARAMETER(ZERO=0.0D0,ONE=1.0D0)
C     WAVE=SQRT(10*N0)/(2*PI*C)      C IN (M/S), N0=AVOGADRO'S NUMBER   INT47260
C     CINT=                                                             INT47270
      PARAMETER(EDCUT=4.0D0,TWO=2.0D0)                                  INT47280
   1  FORMAT(//,1X,'NUCLEAR CARTESIAN COORDINATES (ANG.) AND MASSES'/)  INT47290
   2  FORMAT(2X,'ATOM',4X,'MASS',13X,'X',13X,'Y',13X,'Z'/)              INT47300
   3  FORMAT(1X,I4,F12.6,3X,3F14.10)                                    INT47310
   4  FORMAT(2I5)                                                       INT47320
   5  FORMAT(/,1X,'FILE16 HEADER INCONSISTENT WITH NA')                 INT47330
   6  FORMAT(3F20.10)                                                   INT47340
   7  FORMAT(//,1X,'QUADRATIC FORCE CONSTANTS IN MDYN/A, MDYN/RAD, ',   INT47350
     $   'OR MDYN*A/RAD**2'/)                                           INT47360
   8  FORMAT(//,1X,'G MATRIX')                                          INT47370
   9  FORMAT(//,1X,'EIGENVALUES OF THE G MATRIX'/)                      INT47380
  10  FORMAT(7F12.6)                                                    INT47390
  11  FORMAT(7F10.6)                                                    INT47400
  12  FORMAT(///,1X,'LOWEST EIGENVALUES OF THE G MATRIX'/)              INT47410
  13  FORMAT(///,1X,'VIBRATIONAL FREQUENCIES (CM-1) AND EIGENVECTORS'//)INT47420
  14  FORMAT(//,1X,'POTENTIAL ENERGY DISTRIBUTIONS (%) ',               INT47430
     $             'AMONG DIAGONAL ELEMENTS'//)                         INT47440
  15  FORMAT(//,1X,'VIBRATIONAL ASSIGNMENTS')                           INT47450
  16  FORMAT(//,3X,'MODE',4X,'CM-1',17X,'DOMINANT COMPONENTS OF PED'/)  INT47460
  17  FORMAT(1X,I5,F10.1,2X,5(2X,I3,1X,'(',F5.1,')'))                   INT47470
  18  FORMAT(3F20.10)                                                   INT47480
  19  FORMAT(//,24X,'IR INTENSITY'/                                     INT47490
     $          3X,'MODE',4X,'CM-1',4X,'(KM/MOLE AND RELATIVE)',        INT47500
     $          15X,'DOMINANT COMPONENTS OF PED'/)                      INT47510
  20  FORMAT(1X,I5,F10.1,1X,F12.5,F10.5,3X,5(2X,I3,1X,'(',F5.1,')'))    INT47520
  21  FORMAT(//,1X,'DIPOLE MOMENT DERIVATIVES IN (SYMMETRY) ',          INT47530
     $  'INTERNAL COORDINATES'//4X,'I',7X,'M(X)',6X,'M(Y)',6X,'M(Z)'/)  INT47540
  22  FORMAT(1X,I4,3X,3F10.6)                                           INT47550
  23  FORMAT(//,1X,'DIPOLE MOMENT DERIVATIVES IN NORMAL COORDINATES'    INT47560
     $      //4X,'I',7X,'M(X)',6X,'M(Y)',6X,'M(Z)'/)                    INT47570
  24  FORMAT(1X,I4,3X,3F10.6)                                           INT47580
  25  FORMAT(2I5)                                                       INT47590
  26  FORMAT(/,1X,'FILE18 HEADER INCONSISTENT WITH NA')                 INT47600
  27  FORMAT(3F20.10)                                                   INT47610
  28  FORMAT(//,1X,'TOTAL ENERGY DISTRIBUTIONS (%)',//)                 INT47430
  29  FORMAT(//,3X,'MODE',4X,'CM-1',17X,'DOMINANT COMPONENTS OF TED'/)  INT47460
  30  FORMAT(//,24X,'IR INTENSITY'/                                     INT47490
     $          3X,'MODE',4X,'CM-1',4X,'(KM/MOLE AND RELATIVE)',        INT47500
     $          15X,'DOMINANT COMPONENTS OF TED'/)                      INT47510
  31  FORMAT(//,3X,'REACTION COORDINATE = ',I5,/)
  32  FORMAT(//,1X,'** PROJECTED VIBRATIONAL ANALYSIS ALONG RXN PATH',
     $             ' ******************************'/)
  33  FORMAT(I5)
  34  FORMAT(20X,F20.10)
  35  FORMAT(//, 3X,'INTERNAL GRADIENT PROJECTION VECTOR (FILE12)'/)
  36  FORMAT(6F12.8)
CSA
      WAVE=WAVE0
CSA
      IRXN=0
      if(NFREQ.gt.10) IRXN=NFREQ-10
      WRITE(IOUT,1)                                                     INT47620
      WRITE(IOUT,2)                                                     INT47630
      DO 100  I=1,NA                                                    INT47640
 100  WRITE(IOUT,3) I,XMASS(I),(XA(I,J),J=1,3)                          INT47650
      IF(NFREQ.EQ.4) THEN                                               INT47660
          READ(IIN1,11) ((F(M,N),N=M,NSX),M=1,NSX)                      INT47670
          DO 101  M=1,NSX                                               INT47680
          DO 101  N=M,NSX                                               INT47690
 101      F(N,M)=F(M,N)                                                 INT47700
          GO TO 102                                                     INT47710
      END IF                                                            INT47720
      REWIND I16                                                        INT47730
      READ(I16,4) M,N                                                   INT47740
      IF(M.NE.NA) THEN                                                  INT47750
          WRITE(IOUT,5)                                                 INT47760
          IFLAG=1                                                       INT47770
          RETURN                                                        INT47780
      END IF                                                            INT47790
      READ(I16,6) ((F(M,N),N=1,NSX),M=1,NSX)                            INT47800
 102  WRITE(IOUT,7)                                                     INT47810
Cwa
      if(IRXN.gt.0) then
        DO 103  M=1,NSX
        if(M.ne.IRXN) then
          F(M,IRXN)=ZERO
          F(IRXN,M)=ZERO
        end if
 103    continue
      end if
C
      CALL TABLE1(NC,NSX,F)                                             INT47820
      NV=NSX*(NSX+1)/2                                                  INT47830
      DO 105  N=1,NSX                                                   INT47840
      DO 105  M=1,NSX                                                   INT47850
 105  XR(M,N)=ZERO                                                      INT47860
      DO 110  I=1,NC                                                    INT47870
      II=(I-1)/3+1                                                      INT47880
      XX=ONE/XMASS(II)                                                  INT47890
         DO 115  N=1,NSX                                                INT47900
         DO 115  M=1,NSX                                                INT47910
 115     XR(M,N)=XR(M,N)+BS(M,I)*BS(N,I)*XX                             INT47920
 110  CONTINUE
Cwa                                                                     INT47930
      if(IRXN.ne.0) then
         DO 116  N=1,NSX
         DO 116  M=1,NSX
         if(M.ne.IRXN.and.N.ne.IRXN) then
           XR(M,N)=XR(M,N)-XR(M,IRXN)*XR(N,IRXN)/XR(IRXN,IRXN)
         end if
 116     continue
         DO 117  M=1,NSX
           if(M.ne.IRXN) then
           XR(M,IRXN)=ZERO
           XR(IRXN,M)=ZERO
           end if
 117     continue
         write(IOUT,31) IRXN
      end if
C
      IF(LPRT(3,NPRT).EQ.1) THEN                                        INT47940
         WRITE(IOUT,8)                                                  INT47950
         CALL TABLE1(NC,NSX,XR)                                         INT47960
      END IF                                                            INT47970
C     NOTE THAT XU AND XS ARE IN THE SAME MEMORY LOCATION.              INT47980
      II=0                                                              INT47990
      DO 120  M=1,NSX                                                   INT48000
      DO 120  N=1,M                                                     INT48010
      II=II+1                                                           INT48020
 120  XU(II)=XR(M,N)                                                    INT48030
      CALL RSP(NC,NSX,NV,XU,W,1,XT,BS(1,1),BS(1,2))                     INT48040
      IF(LPRT(3,NPRT).EQ.1) THEN                                        INT48050
         WRITE(IOUT,9)                                                  INT48060
         WRITE(IOUT,10) (W(I),I=1,NSX)                                  INT48070
         GO TO 126                                                      INT48080
      END IF                                                            INT48090
      WRITE(IOUT,12)                                                    INT48100
      IGE=MIN0(5,NSX)                                                   INT48110
      WRITE(IOUT,10) (W(I),I=1,IGE)                                     INT48120
 126  DO 128  I=1,NSX                                                   INT48130
 128  W(I)=DSQRT(W(I))                                                  INT48140
      DO 130  N=1,NSX                                                   INT48150
      DO 130  M=1,NSX                                                   INT48160
         XR(M,N)=ZERO                                                   INT48170
         DO 140  K=1,NSX                                                INT48180
 140     XR(M,N)=XR(M,N)+XT(M,K)*XT(N,K)*W(K)                           INT48190
 130  CONTINUE                                                          INT48200
      DO 150  N=1,NSX                                                   INT48210
      DO 150  M=1,NSX                                                   INT48220
           XX=ZERO                                                      INT48230
           DO 152  J=1,NSX                                              INT48240
 152       XX=XX+(F(M,J)+F(J,M))*XR(J,N)/TWO                            INT48250
 150  XT(M,N)=XX                                                        INT48260
      II=0                                                              INT48270
      DO 154  M=1,NSX                                                   INT48280
      DO 154  N=1,M                                                     INT48290
      II=II+1                                                           INT48300
           XX=ZERO                                                      INT48310
           DO 156  J=1,NSX                                              INT48320
 156       XX=XX+XR(M,J)*XT(J,N)                                        INT48330
 154  XU(II)=XX                                                         INT48340
Cwa   
Cwa  Projection Section 1999
      IF(NFREQ.NE.6.AND.NFREQ.NE.8) GO TO 157
      WRITE(IOUT,32) 
      REWIND I12
      READ(I12,33) NAX
      READ(I12,34) (W(I),I=1,NSX)
      WRITE(IOUT,35)
      WRITE(IOUT,36) (W(I),I=1,NSX)
Cwa  Form v = G(1/2)g into 1st column of XT.    XR contains G(1/2)
      DO 160 I=1,NSX
      XT(I,1)=ZERO
      DO 160 J=1,NSX
 160  XT(I,1)=XT(I,1)+XR(I,J)*W(J)
Cwa   Form w = (F') v into W vector.      XU contains lt of F' 
      DO 161 I=1,NSX
      W(I)=ZERO
      DO 161 J=1,NSX
      IJ=MAX0(I,J)
      IK=MIN0(I,J)
      II=IJ*(IJ-1)/2+IK
 161  W(I)=W(I)+XU(II)*XT(J,1)
Cwa   Compute scalars R1=v.v and R2=v.w    
      R1=ZERO
      R2=ZERO
      DO 162 I=1,NSX
      R1=R1+XT(I,1)*XT(I,1)
 162  R2=R2+W(I)*XT(I,1)
Cwa   Form (F")ij = (F')ij - (wi vj + wj vi)/R1 + (R2/R1^2) vi vj 
      II=0
      DO 163  M=1,NSX
      DO 163  N=1,M
      II=II+1
 163  XU(II)=XU(II)-(W(M)*XT(N,1)+W(N)*XT(M,1))/R1+
     $               (R2/(R1*R1))*XT(M,1)*XT(N,1) 
Cwa
 157  CALL RSP(NC,NSX,NV,XU,W,1,XT,BS(1,1),BS(1,2))                     INT48350
      DO 164  I=1,NSX                                                   INT48360
           XX=ONE                                                       INT48370
           IF(W(I).LT.ZERO) XX=-XX                                      INT48380
 164       W(I)=DSQRT(DABS(W(I)))*XX*WAVE                               INT48390
      DO 165  N=1,NSX                                                   INT48400
      DO 165  M=1,NSX                                                   INT48410
           XS(M,N)=ZERO                                                 INT48420
           DO 170  I=1,NSX                                              INT48430
 170       XS(M,N)=XS(M,N)+XR(M,I)*XT(I,N)                              INT48440
 165  CONTINUE                                                          INT48450
      WRITE(IOUT,13)                                                    INT48550
      CALL TABLE2(NC,NSX,NSX,W,XS)                                      INT48560
      IF(NMODE.EQ.1) THEN
        DO 175  N=1,NSX                                                 INT48460
             XX=ZERO                                                    INT48470
             DO 180  M=1,NSX                                            INT48480
             XT(M,N)=XS(M,N)*F(M,M)*XS(M,N)                             INT48490
 180         XX=XX+DABS(XT(M,N))                                        INT48500
             IF(XX.EQ.ZERO) GO TO 175                                   INT48510
             DO 185  M=1,NSX                                            INT48520
 185         XT(M,N)=1.0D2*XT(M,N)/XX                                   INT48530
 175    CONTINUE                                                        INT48540
      ELSE
        DO 191  N=1,NSX
        DO 191  M=1,NSX
 191         XR(M,N)=XS(M,N)
C   Note that the following step overwrites the F array in memory.
        DO 192  N=1,NSX
        DO 192  M=1,NSX
 192         XR(M,N+NSX)=ZERO
        DO 193  M=1,NSX
 193         XR(M,M+NSX)=ONE
        CALL FLIN(XR,NC,NSX,NSX,DET)
        DO 194  N=1,NSX
        DO 194  M=1,NSX
 194         XT(M,N)=XS(M,N)*XR(N,M+NSX)*1.0D2
      END IF
      IF(NMODE.EQ.1) THEN
        WRITE(IOUT,14)                                                  INT48570
      ELSE
        WRITE(IOUT,28)                                                  INT48570
      END IF
      CALL TABLE3(NC,NSX,W,XT)                                          INT48580
      DO 200  J=1,NSX                                                   INT48590
      IA(5,J)=4                                                         INT48600
      DO 200  I=1,4                                                     INT48610
      IA(I,J)=0                                                         INT48620
 200  XR(I,J)=ZERO                                                      INT48630
      DO 205  J=1,NSX                                                   INT48640
          DO 210  I=1,NSX                                               INT48650
          IF(DABS(XT(I,J)).GT.DABS(XR(1,J))) THEN                       INT48660
             XR(1,J)=XT(I,J)                                            INT48670
             IA(1,J)=I                                                  INT48680
          END IF                                                        INT48690
 210      CONTINUE                                                      INT48700
          DO 212  I=1,NSX                                               INT48710
          ITEST=I.NE.IA(1,J)                                            INT48720
          IF(DABS(XT(I,J)).GT.DABS(XR(2,J)).AND.ITEST) THEN             INT48730
             XR(2,J)=XT(I,J)                                            INT48740
             IA(2,J)=I                                                  INT48750
          END IF                                                        INT48760
 212      CONTINUE                                                      INT48770
          DO 214  I=1,NSX                                               INT48780
          ITEST=I.NE.IA(1,J).AND.I.NE.IA(2,J)                           INT48790
          IF(DABS(XT(I,J)).GT.DABS(XR(3,J)).AND.ITEST) THEN             INT48800
             XR(3,J)=XT(I,J)                                            INT48810
             IA(3,J)=I                                                  INT48820
          END IF                                                        INT48830
 214      CONTINUE                                                      INT48840
          DO 216  I=1,NSX                                               INT48850
          ITEST=I.NE.IA(1,J).AND.I.NE.IA(2,J).AND.I.NE.IA(3,J)          INT48860
          IF(DABS(XT(I,J)).GT.DABS(XR(4,J)).AND.ITEST) THEN             INT48870
             XR(4,J)=XT(I,J)                                            INT48880
             IA(4,J)=I                                                  INT48890
          END IF                                                        INT48900
 216      CONTINUE                                                      INT48910
          DO 218  K=1,4                                                 INT48920
             IF(DABS(XR(K,J)).LT.EDCUT) THEN                            INT48930
             IA(5,J)=K-1                                                INT48940
             GO TO 205                                                  INT48950
             END IF                                                     INT48960
 218      CONTINUE                                                      INT48970
 205  CONTINUE                                                          INT48980
      DO 230  J=1,NSX                                                   INT48990
          DO 235  K=1,4                                                 INT49000
          XX=XS(IA(K,J),J)                                              INT49010
          IF(XX.LT.ZERO) IA(K,J)=-IA(K,J)                               INT49020
 235      CONTINUE                                                      INT49030
          IF(IA(1,J).LT.0) THEN                                         INT49040
          DO 240  K=1,4                                                 INT49050
 240      IA(K,J)=-IA(K,J)                                              INT49060
          END IF
 230      CONTINUE                                                      INT49080
      IF(IRINT.EQ.0) GO TO 290                                          INT49090
      IF(IRINT.EQ.1) THEN                                               INT49100
         READ(I18,25) I,J                                               INT49110
         IF(I.NE.NA) THEN                                               INT49120
         WRITE(IOUT,26)                                                 INT49130
         IFLAG=2                                                        INT49140
         RETURN                                                         INT49150
         END IF                                                         INT49160
         READ(I18,27) ((BS(I,J),I=1,NSX),J=1,3)                         INT49170
      ELSE                                                              INT49180
      DO 250  I=1,NSX                                                   INT49190
 250  READ(IIN1,18) (BS(I,J),J=1,3)                                     INT49200
      END IF                                                            INT49210
      WRITE(IOUT,21)                                                    INT49220
      DO 255  I=1,NSX                                                   INT49230
 255  WRITE(IOUT,22) I,(BS(I,J),J=1,3)                                  INT49240
      XZ=ZERO                                                           INT49250
      IZ=0                                                              INT49260
      DO 265  I=1,NSX                                                   INT49270
          XY=ZERO                                                       INT49280
          DO 270  J=1,3                                                 INT49290
          XX=ZERO                                                       INT49300
          DO 275  K=1,NSX                                               INT49310
 275      XX=XX+XS(K,I)*BS(K,J)                                         INT49320
          BS(I,J+3)=XX                                                  INT49330
 270      XY=XY+XX*XX                                                   INT49340
          IF(XY.GT.XZ) THEN                                             INT49350
          IZ=I                                                          INT49360
          XZ=XY                                                         INT49370
          END IF                                                        INT49380
 265  XT(I,1)=XY*CINT                                                   INT49390
      XZ=XZ*CINT                                                        INT49400
      IF(XZ.EQ.ZERO) XZ=ONE                                             INT49410
      DO 280  I=1,NSX                                                   INT49420
 280  XT(I,2)=XT(I,1)/XZ                                                INT49430
      IF(LPRT(3,NPRT).EQ.2) THEN                                        INT49440
          WRITE(IOUT,23)                                                INT49450
          DO 282  I=1,NSX                                               INT49460
 282      WRITE(IOUT,24) I,(BS(I,J),J=4,6)                              INT49470
      END IF                                                            INT49480
 290  WRITE(IOUT,15)                                                    INT49490
      IF(IRINT.EQ.0) THEN                                               INT49500
          IF(NMODE.EQ.1) THEN
            WRITE(IOUT,16)                                              INT49510
          ELSE
            WRITE(IOUT,29)                                              INT49510
          END IF
          DO 300  J=1,NSX                                               INT49520
 300      WRITE(IOUT,17) J,W(J),(IA(K,J),XR(K,J),K=1,IA(5,J))           INT49530
      ELSE                                                              INT49540
          IF(NMODE.EQ.1) THEN
            WRITE(IOUT,19)                                              INT49510
          ELSE
            WRITE(IOUT,30)                                              INT49510
          END IF
          DO 310  J=1,NSX                                               INT49560
          WRITE(IOUT,20) J,W(J),XT(J,1),XT(J,2),                        INT49570
     $                  (IA(K,J),XR(K,J),K=1,IA(5,J))                   INT49580
 310      CONTINUE                                                      INT49590
      END IF                                                            INT49600
      RETURN                                                            INT49610
      END                                                               INT49620
C////////////////////////////////////////////////////////////////////// INT49630
      SUBROUTINE NORMCO(NA,NAD,NC,NFREQ,IRINT,                          INT49640
     $                       XMASS,XA,BS,F,XR,XS,XT,XU,IFLAG)           INT49650
      IMPLICIT REAL*8 (A-H,O-Z)                                         INT49660
      DIMENSION XA(NAD,3),BS(NC,NC),F(NC,NC),XS(NC,NC),XT(NC,NC)        INT49670
      DIMENSION XMASS(NA),W(6),XR(NC,NC),XU(1),IZ(6)                    INT49680
      DIMENSION T1(3,3),T2(3,3),V1(3),V2(3),V3(3),V4(3),V5(6)
      DIMENSION IE(2,3)
      CHARACTER LABEL*10
      COMMON /IO/ IIN1,IOUT,IIN2,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25,
     $   I31,I32,I33,I35,I36,I37,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10,
     $   ISCR11,ISCR12,ISCR13,ISCR14
      COMMON /PHYCON/ BOHR,DEBYE,HART,WAVE0,CINT
      PARAMETER(ZERO=0.0D0,ONE=1.0D0)
      PARAMETER(FRQCUT=5.0D0,TWO=2.0D0)                                 INT49740
   1  FORMAT(//,1X,'NUCLEAR CARTESIAN COORDINATES (ANG.) AND MASSES'/)  INT49750
   2  FORMAT(2X,'ATOM',4X,'MASS',13X,'X',13X,'Y',13X,'Z'/)              INT49760
   3  FORMAT(1X,I4,F12.6,3X,3F14.10)                                    INT49770
   4  FORMAT(2I5)                                                       INT49780
   5  FORMAT(/,1X,'FILE15 HEADER INCONSISTENT WITH NA')                 INT49790
   6  FORMAT(3F20.10)                                                   INT49800
   7  FORMAT(//,1X,'QUADRATIC FORCE CONSTANTS IN HARTREE/BOHR**2'/)     INT49810
   8  FORMAT(///,1X,'ZERO FREQUENCIES (CM-1) AND EIGENVECTORS'//)       INT49820
   9  FORMAT(///,1X,'ZERO FREQUENCIES (CM-1)'/)                         INT49830
  10  FORMAT(6F12.4)                                                    INT49840
  11  FORMAT(///,1X,'VIBRATIONAL FREQUENCIES (CM-1) AND EIGENVECTORS'//)INT49850
  12  FORMAT(2I5)                                                       INT49860
  13  FORMAT(/,1X,'FILE17 HEADER INCONSISTENT WITH NA')                 INT49870
  14  FORMAT(3F20.10)                                                   INT49880
  15  FORMAT(//,1X,'DIPOLE MOMENT DERIVATIVES IN NORMAL COORDINATES'    INT49890
     $      //4X,'I',7X,'M(X)',6X,'M(Y)',6X,'M(Z)'/)                    INT49900
  16  FORMAT(1X,I4,3X,3F10.6)                                           INT49910
  17  FORMAT(//,24X,'IR INTENSITY'/                                     INT49920
     $          3X,'MODE',4X,'CM-1',4X,'(KM/MOLE AND RELATIVE)')        INT49930
  18  FORMAT(1X,I5,F10.1,1X,F12.5,F10.5)                                INT49940
  19  FORMAT(//,3X,'MODE',4X,'CM-1'/)                                   INT49950
  20  FORMAT(1X,I5,F10.1)                                               INT49960
  21  FORMAT(A10)
  22  FORMAT(20X,3F20.10)
  23  FORMAT(//,3X,'CARTESIAN GRADIENT VECTOR (FILE11)')
  24  FORMAT(6F12.8)
  25  FORMAT(//,1X,'** PROJECTED VIBRATIONAL ANALYSIS ALONG RXN PATH',
     $             ' ******************************')
CSA
      WAVE=WAVE0*DSQRT(HART)/BOHR
CSA
      WRITE(IOUT,1)                                                     INT49970
      WRITE(IOUT,2)                                                     INT49980
      DO 100  I=1,NA                                                    INT49990
 100  WRITE(IOUT,3) I,XMASS(I),(XA(I,J),J=1,3)                          INT50000
      REWIND I15                                                        INT50010
      READ(I15,4) I,J                                                   INT50020
      IF(I.NE.NA) THEN                                                  INT50030
          WRITE(IOUT,5)                                                 INT50040
          IFLAG=1                                                       INT50050
          RETURN                                                        INT50060
      END IF                                                            INT50070
      READ(I15,6) ((F(I,J),J=1,NC),I=1,NC)                              INT50080
      WRITE(IOUT,7)                                                     INT50090
      CALL TABLE1(NC,NC,F)                                              INT50100
      NV=NC*(NC+1)/2                                                    INT50110
      DO 102  J=1,NC                                                    INT50120
      JK=(J-1)/3+1                                                      INT50130
      DO 102  I=1,NC                                                    INT50140
      IK=(I-1)/3+1                                                      INT50150
 102  XR(I,J)=(F(J,I)+F(I,J))/(TWO*DSQRT(XMASS(IK)*XMASS(JK)))          INT50160
      II=0                                                              INT50170
      DO 104  I=1,NC                                                    INT50180
      DO 104  J=1,I                                                     INT50190
      II=II+1                                                           INT50200
 104  XU(II)=XR(I,J)                                                    INT50210
C     NOTE THAT XU AND XS ARE IN THE SAME MEMORY LOCATION.              INT50220
Cwa
Cwa   Projection section 1999
      IF(NFREQ.NE.7.AND.NFREQ.NE.8) GO TO 132
      WRITE(IOUT,25)
Cwa   Form inertia tensor in c.m. frame and compute I^(-1/2)
      XMT=ZERO
      DO 106  I=1,NA
 106  XMT=XMT+XMASS(I)
      DO 108  J=1,3 
      V1(J)=ZERO
      DO 107  I=1,NA
 107  V1(J)=V1(J)+XA(I,J)*XMASS(I) 
 108  V1(J)=V1(J)/XMT
      DO 110  I=1,6
 110  V5(I)=ZERO    
      DO 112  I=1,NA
      DO 113  J=1,3
 113  V2(J)=XA(I,J)-V1(J)
      V5(1)=V5(1)+XMASS(I)*(V2(2)*V2(2)+V2(3)*V2(3)) 
      V5(3)=V5(3)+XMASS(I)*(V2(1)*V2(1)+V2(3)*V2(3))
      V5(6)=V5(6)+XMASS(I)*(V2(1)*V2(1)+V2(2)*V2(2))
      V5(2)=V5(2)-XMASS(I)*V2(1)*V2(2)
      V5(4)=V5(4)-XMASS(I)*V2(1)*V2(3)
 112  V5(5)=V5(5)-XMASS(I)*V2(2)*V2(3)
      CALL RSP(3,3,6,V5,V2,1,T1,V3,V4)
      DO 114  I=1,3
      DO 114  J=1,3
      T2(I,J)=ZERO
      DO 114  K=1,3
 114  T2(I,J)=T2(I,J)+T1(I,K)*T1(J,K)/DSQRT(V2(K))
Cwa   Read in Cartesian gradient vector, mass weight, normalize.
 116  READ(I11,21,END=117) LABEL
      GO TO 116
 117  IK=NA+1
      DO 118 K=1,IK
 118  BACKSPACE I11
      READ(I11,22) (XT(I,1),I=1,NC)
      WRITE(IOUT,23)
      WRITE(IOUT,24) (XT(I,1),I=1,NC) 
      DO 119  I=1,NA
 119  XT(I,2)=DSQRT(XMASS(I))
      GN=ZERO
      DO 120  I=1,NC
      II=(I-1)/3+1
      XT(I,1)=XT(I,1)/XT(II,2)
 120  GN=GN+XT(I,1)*XT(I,1)
      GN=ONE/DSQRT(GN)
      DO 121  I=1,NC
 121  XT(I,1)=XT(I,1)*GN
Cwa   Compute projection matrix 
      DO 122 I=1,NC
      DO 122 J=1,6
 122  BS(I,J)=ZERO
      XMTR=ONE/DSQRT(XMT)
      DO 123  I=1,NA
      DO 123  J=1,3
      II=3*(I-1) 
 123  BS(II+J,J)=XT(I,2)*XMTR
Cwa   IE(I,J)=MOD(I+J-1,3)+1
      IE(1,1)=2
      IE(2,1)=3
      IE(1,2)=3
      IE(2,2)=1
      IE(1,3)=1
      IE(2,3)=2
      DO 124  I=1,NA
      II=3*(I-1)
      DO 125  K=1,3
 125  V2(K)=(XA(I,K)-V1(K))*XT(I,2)
      DO 126  J=1,3
      DO 126  K=1,3
 126  BS(II+J,K+3)=T2(K,IE(1,J))*V2(IE(2,J))-T2(K,IE(2,J))*V2(IE(1,J))
 124  CONTINUE	
      DO 127  I=1,NC
      DO 127  J=1,NC
      XR(I,J)=-XT(I,1)*XT(J,1)
      DO 127  K=1,6
 127  XR(I,J)=XR(I,J)-BS(I,K)*BS(J,K)
      DO 128  I=1,NC
 128  XR(I,I)=XR(I,I)+ONE
Cwa
      WRITE(ICHECK,2000)
 2000 FORMAT(//,3X,'IDEMPOTENCY CHECK OF CARTESIAN PROJECTION MATRIX')
      DO 2001  I=1,6
      DO 2001  J=1,6
      XX=ZERO
      DO 2002  K=1,NC
 2002 XX=XX+BS(K,I)*BS(K,J)
      WRITE(ICHECK,2003) I,J,XX
 2001 CONTINUE
 2003 FORMAT(2I5,F20.10)
      I=7
      DO 2005  J=1,6
      XX=ZERO
      DO 2006  K=1,NC
 2006 XX=XX+XT(K,1)*BS(K,J)
 2005 WRITE(ICHECK,2003) I,J,XX
      XX=ZERO
      DO 2004  K=1,NC
 2004 XX=XX+XT(K,1)*XT(K,1)
      WRITE(ICHECK,2003) I,I,XX
Cwa
      DO 129  I=1,NC
      DO 129  J=1,NC
      XT(I,J)=ZERO
      DO 129  K=1,NC
      IJ=MAX0(I,K)
      IK=MIN0(I,K)
      II=IJ*(IJ-1)/2+IK
 129  XT(I,J)=XT(I,J)+XU(II)*XR(K,J) 
      II=0
      DO 130  I=1,NC
      DO 130  J=1,I
      II=II+1
      XU(II)=ZERO
      DO 130  K=1,NC
 130  XU(II)=XU(II)+XR(I,K)*XT(K,J)
Cwa End 1999
Cwa
 132  CALL RSP(NC,NC,NV,XU,BS(1,1),1,XT,BS(1,2),BS(1,3))                INT50230
      DO 133  I=1,NC                                                    INT50240
      IK=(I-1)/3+1                                                      INT50250
      XX=ONE/DSQRT(XMASS(IK))                                           INT50260
      DO 133  J=1,NC                                                    INT50270
 133  XR(I,J)=XT(I,J)*XX                                                INT50280
      DO 134  I=1,NC                                                    INT50290
           XX=ONE                                                       INT50300
           IF(BS(I,1).LT.ZERO) XX=-XX                                   INT50310
 134       XU(I)=DSQRT(DABS(BS(I,1)))*XX*WAVE                           INT50320
      II=0                                                              INT50330
      IJ=0                                                              INT50340
      DO 135  J=1,NC                                                    INT50350
           XX=DABS(XU(J))                                               INT50360
           IF(XX.LT.FRQCUT) THEN                                        INT50370
           IJ=IJ+1                                                      INT50380
           IZ(IJ)=J                                                     INT50390
           W(IJ)=XU(J)                                                  INT50400
           DO 138  I=1,NC                                               INT50410
 138       BS(I,IJ)=XR(I,J)                                             INT50420
           GO TO 135                                                    INT50430
           END IF                                                       INT50440
      II=II+1                                                           INT50450
      XU(II)=XU(J)                                                      INT50460
      DO 140  I=1,NC                                                    INT50470
 140  XR(I,II)=XR(I,J)                                                  INT50480
 135  CONTINUE                                                          INT50490
      NZ=IJ                                                             INT50500
      NF=II                                                             INT50510
      IF(LPRT(3,NPRT).EQ.3) THEN                                        INT50520
           WRITE(IOUT,8)                                                INT50530
           CALL TABLE2(NC,NC,NZ,W,BS)                                   INT50540
      ELSE                                                              INT50550
           WRITE(IOUT,9)                                                INT50560
           WRITE(IOUT,10) (W(I),I=1,NZ)                                 INT50570
      END IF                                                            INT50580
      WRITE(IOUT,11)                                                    INT50590
      CALL TABLE2(NC,NC,NF,XU,XR)                                       INT50600
      IF(IRINT.EQ.0) GO TO 184                                          INT50610
      READ(I17,12) I,J                                                  INT50620
           IF(I.NE.NA) THEN                                             INT50630
           WRITE(IOUT,13)                                               INT50640
           IFLAG=2                                                      INT50650
           RETURN                                                       INT50660
           END IF                                                       INT50670
      READ(I17,14) ((BS(I,J),I=1,NC),J=1,3)                             INT50680
      XZ=ZERO                                                           INT50690
      IY=0                                                              INT50700
      DO 165  I=1,NF                                                    INT50710
          XY=ZERO                                                       INT50720
          DO 170  J=1,3                                                 INT50730
          XX=ZERO                                                       INT50740
          DO 175  K=1,NC                                                INT50750
 175      XX=XX+XR(K,I)*BS(K,J)                                         INT50760
          BS(I,J+3)=XX                                                  INT50770
 170      XY=XY+XX*XX                                                   INT50780
          IF(XY.GT.XZ) THEN                                             INT50790
          IY=I                                                          INT50800
          XZ=XY                                                         INT50810
          END IF                                                        INT50820
 165  XT(I,1)=XY*CINT                                                   INT50830
      XZ=XZ*CINT                                                        INT50840
      IF(XZ.EQ.ZERO) XZ=ONE                                             INT50850
      DO 180  I=1,NF                                                    INT50860
 180  XT(I,2)=XT(I,1)/XZ                                                INT50870
      IF(LPRT(3,NPRT).EQ.2) THEN                                        INT50880
          WRITE(IOUT,15)                                                INT50890
          DO 182  I=1,NF                                                INT50900
 182      WRITE(IOUT,16) I,(BS(I,J),J=4,6)                              INT50910
      END IF                                                            INT50920
 184  IF(IRINT.EQ.0) THEN                                               INT50930
      WRITE(IOUT,19)                                                    INT50940
      DO 186  I=1,NF                                                    INT50950
 186  WRITE(IOUT,20) I,XU(I)                                            INT50960
      ELSE                                                              INT50970
      WRITE(IOUT,17)                                                    INT50980
      DO 190  I=1,NF                                                    INT50990
 190  WRITE(IOUT,18) I,XU(I),XT(I,1),XT(I,2)                            INT51000
      END IF                                                            INT51010
      RETURN                                                            INT51020
      END                                                               INT51030
C////////////////////////////////////////////////////////////////////// XX 04410
      SUBROUTINE SQMFC(NA,NAD,NC,NSX,NSF,NISO,NMODE,XA,BS,XM,V,WT,      XX 00010
     $    GH,D,G,G0,F0,DF,HF,FA,W,S1,S2,S3,F,XR,XS,XT,XU,               XX 00020
     $    C,TA,JA,IFLAG)
      IMPLICIT REAL*8 (A-H,O-Z)                                         XX 00030
      REAL*8 JA                                                         XX 00030
      LOGICAL ITEST                                                     XX 00040
      DIMENSION HF(NSF,NSF),G0(NSF),F0(NSF),DF(NSF)                     XX 00050
      DIMENSION NSC(100,0:100),IA(100,100)                              XX 00060
      DIMENSION XA(NAD,3),BS(NC,NC),XM(NA,NISO),V(NSX,NISO),G(NSF)      XX 00070
      DIMENSION GH(NSX,NSX,NISO),FA(NSF),F(NSX,NSX),D(NSX,NSX)          XX 00080
      DIMENSION S1(NSX),S2(NSX),S3(NSX),W(NSX),XR(NC,2),XU(1),XT(NC,NC) XX 00090
      DIMENSION XS(NC,NC),WT(NSX,NISO),TA(NSX,NSX,NSF)                  XX 00100
      DIMENSION C(NSF,NSF),JA(NSX,NSF,NISO)
C     PARAMETER(CINT=42.25472D0)
C     WAVE=SQRT(10*N0)/(2*PI*C)      C IN (M/S), N0=AVOGADRO'S NUMBER   XX 00160
      PARAMETER(ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0,THREE=3.0D0,FOUR=4.0D0)  XX 00170
      PARAMETER(VCUT=1.0D-2,EDCUT=4.0D0,MAXIT=100,XL00=1.0D0)           XX 00180
      PARAMETER(XLMAX=10.0D0,ONEH=1.0D2,FDISP=1.0D-3)
      PARAMETER(SCONV=1.0D-10,GCONV=1.0D-5,XLCUT=1.0D-6,F2CUT=1.0D-6)   XX 00190
      PARAMETER(SCONV2=1.0D-7,GCONV2=1.0D-5)                            XX 00190
      COMMON /IO/ IIN1,IOUT,IIN2,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25,
     $   I31,I32,I33,I35,I36,I37,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10,
     $   ISCR11,ISCR12,ISCR13,ISCR14
      COMMON /PHYCON/ BOHR,DEBYE,HART,WAVE0,CINT
   1  FORMAT(5I5)                                                       XX 00200
   2  FORMAT(F10.5,I5)                                                  XX 00210
   3  FORMAT(16I5)                                                      XX 00220
   4  FORMAT(7F10.7)                                                    XX 00230
   5  FORMAT(I5)                                                        XX 00240
   6  FORMAT(I5,2F10.7)                                                 XX 00250
   7  FORMAT(//,1X,'NUCLEAR CARTESIAN COORDINATES (ANG.)'/)             XX 00260
   8  FORMAT(2X,'ATOM',9X,'X',13X,'Y',13X,'Z'/)                         XX 00270
   9  FORMAT(1X,I4,3X,3F14.10)                                          XX 00280
  10  FORMAT(//,1X,'MASSES FOR ISOTOPOMERS'/)                           XX 00290
  11  FORMAT(1X,'ISO=',5(I10,4X)/)                                      XX 00300
  12  FORMAT(5X,5F14.5)                                                 XX 00310
  13  FORMAT(//,1X,'EXPERIMENTAL FREQUENCIES FOR ISOTOPOMERS',          XX 00320
     $       ' (RELATIVE WEIGHTS GIVEN IN PARENTHESES)')
  14  FORMAT(1X,'ISO=',8I16/)                                           XX 00330
  15  FORMAT(5X,8(F9.2,' (',F7.5,')'))                                  XX 00340
  16  FORMAT(2I5)                                                       XX 00350
  17  FORMAT(/,1X,'FILE16 HEADER INCONSISTENT WITH NA')                 XX 00360
  18  FORMAT(3F20.10)                                                   XX 00370
  19  FORMAT(//,1X,'UNSCALED QUADRATIC FORCE CONSTANTS IN MDYN/A, ',    XX 00380
     $   'MDYN/RAD, OR MDYN*A/RAD**2'/)                                 XX 00390
  20  FORMAT(//' ISOTOPOMER ',I4,/,1X,                                  XX 00400
     $    'LOWEST EIGENVALUES OF THE G MATRIX'/)                        XX 00410
  21  FORMAT(7F12.6)                                                    XX 00420
  22  FORMAT(//,1X,'SCALED QUADRATIC FORCE CONSTANTS IN MDYN/A, ',      XX 00430
     $   'MDYN/RAD, OR MDYN*A/RAD**2'/)                                 XX 00440
  23  FORMAT(//,1X,'OPTIMIZED SCALE FACTORS'/)                          XX 00450
  24  FORMAT(4(I5,1X,F10.7))                                            XX 00460
  25  FORMAT(///,1X,'ISOTOPOMER ',I4,/,1X,                              XX 00470
     $           'VIBRATIONAL FREQUENCIES (CM-1) AND EIGENVECTORS'//)   XX 00480
  26  FORMAT(//,1X,'ISOTOPOMER ',I4,/,1X,'POTENTIAL ENERGY ',           XX 00490
     $         'DISTRIBUTIONS (%) AMONG DIAGONAL ELEMENTS'//)           XX 00500
  27  FORMAT(//,1X,'ISOTOPOMER ',I4,/,1X,'VIBRATIONAL ASSIGNMENTS')     XX 00510
  28  FORMAT(//,3X,'MODE',4X,'EXPT',5X,'THEORY',3X,'RESIDUAL',3X,       XX 00520
     $       '% ERROR',15X,'DOMINANT COMPONENTS OF PED'/)               XX 00530
  29  FORMAT(1X,I5,2F10.1,2F10.2,2X,5(2X,I3,1X,'(',F5.1,')'))           XX 00540
  30  FORMAT(/' ITER','  LEAST-SQUARES SUM ','  DIRECTIONAL DERIVATIVE',XX 00550
     $     '   DISPLACEMENT NORM')
  31  FORMAT(I5,3F20.10)                                                XX 00560
  32  FORMAT(1X,'GRADIENTS (RELATIVE)')                                 XX 00570
  33  FORMAT(5F12.5)                                                    XX 00580
  34  FORMAT(1X,'LAMBDA UPDATE AVERTED   ',E12.5,5X,E12.5)              XX 00590
  35  FORMAT(1X,'LAMBDA=',E12.5)                                        XX 00600
  36  FORMAT(//,1X,'OPTIMIZATION OF SCALE FACTORS'/)                    XX 00610
  38  FORMAT(/,1X,' DIAGONAL ELEMENTS OF INVERSE HESSIAN'/)             XX 00620
  39  FORMAT(5E12.5)                                                    XX 00630
  40  FORMAT(1X,I5,10X,F10.1,22X,5(2X,I3,1X,'(',F5.1,')'))              XX 00640
  41  FORMAT(//,'  NSF',' NISO',' NOPT','   NH','  NWT'/,5I5)           XX 00650
  42  FORMAT(//,1X,'INITIAL DIAGONAL ELEMENTS OF INVERSE HESSIAN'/)     XX 00660
  43  FORMAT(5E10.3)                                                    XX 00670
  44  FORMAT(2I5,10X,F20.10)                                            XX 00680
  45  FORMAT(//,1X,'ISOTOPOMER ',I4,/,1X,'TOTAL ENERGY ',               XX 00490
     $         'DISTRIBUTIONS (%)'//)                                   XX 00500
  46  FORMAT(//,3X,'MODE',4X,'EXPT',5X,'THEORY',3X,'RESIDUAL',3X,       XX 00520
     $       '% ERROR',15X,'DOMINANT COMPONENTS OF TED'/)               XX 00530
  47  FORMAT(//,1X,'SCALE FACTOR',7X,'CONNECTING COORDINATES'/)
  48  FORMAT(3X,I5,8X,16I5)
  49  FORMAT(//,1X,'OPTIMIZED SCALE FACTORS (STD. ERRORS)'/)
  50  FORMAT(2(I5,1X,F10.7,' (',F8.5,')'))
  51  FORMAT(/' STD. DEV. OF OBSERVATION OF UNIT WEIGHT',F12.7,
     $       /' EFFECTIVE DEGREES OF FREEDOM',F12.7,
     $       /' INTEGER DEGREES OF FREEDOM',I5)
  52  FORMAT(/'CORRELATION COEFFICIENTS'/)
  53  FORMAT(//,1X,'EIGENVALUES AND EIGENVECTORS OF SCALE FACTOR ',
     $      'HESSIAN'/)
  56  FORMAT(//,3X,'MODE',4X,'EXPT',5X,'THEORY',3X,'RESIDUAL',3X,       XX 00520
     $  '% ERROR',2X,'STD. DEV.',14X,'DOMINANT COMPONENTS OF PED'/)     XX 00530
  57  FORMAT(//,3X,'MODE',4X,'EXPT',5X,'THEORY',3X,'RESIDUAL',3X,       XX 00520
     $  '% ERROR',2X,'STD. DEV.',14X,'DOMINANT COMPONENTS OF TED'/)     XX 00530
  58  FORMAT(1X,I5,2F10.1,2F10.2,F10.3,2X,5(2X,I3,1X,'(',F5.1,')'))     XX 00540
  59  FORMAT(1X,I5,10X,F10.1,20X,F10.3,2X,5(2X,I3,1X,'(',F5.1,')'))     XX 00640
  60  FORMAT(/,1X,'HESSIAN FOR ITERATION',I5)
  61  FORMAT(3F20.10)
C*****
C   Note that no memory is allocated for JA and TA unless ABS(NH) > 1.
C   Thus, TA and JA cannot be used in first-order algorithms.
C*****
CSA
      WAVE=WAVE0
CSA
      IFLAG=0                                                           XX 00690
      IF(NSF.GT.100) THEN                                               XX 00700
         IFLAG=1                                                        XX 00710
         RETURN                                                         XX 00720
      END IF                                                            XX 00730
      DO 100  N=1,NSF                                                   XX 00740
      DO 100  M=1,NSF                                                   XX 00750
 100  HF(M,N)=ZERO                                                      XX 00760
      CALL LOCATE(IIN1,'# SQMFC ##',IERR)                               XX 00770
      READ(IIN1,1) M,N,NOPT,NH,NWT                                      XX 00780
      DO 102 I=1,NSF                                                    XX 00790
      READ(IIN1,2) FA(I),NN                                             XX 00800
      NSC(I,0)=NN                                                       XX 00810
 102  READ(IIN1,3) (NSC(I,J),J=1,NN)                                    XX 00820
      DO 104 N=1,NISO                                                   XX 00830
      DO 104 I=1,NSX                                                    XX 00840
      V(I,N)=ZERO                                                       XX 00850
 104  WT(I,N)=ZERO                                                      XX 00850
      NFT=0
      DO 106 N=1,NISO                                                   XX 00860
      CALL MASSIN(XM(1,N),NA,IFLAG)                                     XX 00870
      IF(NOPT.EQ.0) GO TO 106                                           XX 00880
      READ(IIN1,5) NFUND                                                XX 00890
      NFT=NFT+NFUND
      DO 108 J=1,NFUND                                                  XX 00900
      READ(IIN1,6) II,V(II,N),WT(II,N)                                  XX 00910
      IF(WT(II,N).EQ.ZERO) THEN
        IF(NWT.EQ.0) THEN
           WT(II,N)=ONE/V(II,N)
        ELSE IF(NWT.EQ.1) THEN
           WT(II,N)=ONE
        ELSE IF(NWT.EQ.2) THEN
           WT(II,N)=ONE/(V(II,N)*V(II,N))
        END IF
      END IF
 108  CONTINUE                                                          XX 00920
 106  CONTINUE                                                          XX 00920
C  Normalize the weights.
      IF(NOPT.EQ.1) THEN
        WTSUM=ZERO
        DO 110  N=1,NISO
        DO 110  I=1,NSX
 110    WTSUM=WTSUM+WT(I,N)
        DO 112  N=1,NISO
        DO 112  I=1,NSX
 112    WT(I,N)=WT(I,N)/WTSUM
      END IF
C
      IF(NH.EQ.1.AND.NOPT.EQ.1) THEN                                    XX 00930
      READ(IIN1,43) (HF(I,I),I=1,NSF)                                   XX 00940
      END IF                                                            XX 00950
C                                                                       XX 00960
      WRITE(IOUT,41) NSF,NISO,NOPT,NH,NWT                               XX 00970
      WRITE(IOUT,7)                                                     XX 00980
      WRITE(IOUT,8)                                                     XX 00990
      DO 114  I=1,NA                                                    XX 01000
 114  WRITE(IOUT,9) I,(XA(I,J),J=1,3)                                   XX 01010
      WRITE(IOUT,47)
      DO 116  I=1,NSF
      NJ=NSC(I,0)
 116  WRITE(IOUT,48) I,(NSC(I,J),J=1,NJ)
      WRITE(IOUT,10)                                                    XX 01020
      WRITE(IOUT,11) (I,I=1,NISO)                                       XX 01030
      DO 120  I=1,NA                                                    XX 01040
 120  WRITE(IOUT,12) (XM(I,J),J=1,NISO)                                 XX 01050
      IF(NOPT.EQ.0) GO TO 124                                           XX 01060
      WRITE(IOUT,13)                                                    XX 01070
      WRITE(IOUT,14) (I,I=1,NISO)                                       XX 01080
      DO 122  I=1,NSX                                                   XX 01090
 122  WRITE(IOUT,15) (V(I,J),WT(I,J),J=1,NISO)                          XX 01100
      IF(NH.EQ.1.AND.NOPT.EQ.1) THEN                                    XX 01110
      WRITE(IOUT,42)                                                    XX 01120
      WRITE(IOUT,43) (HF(I,I),I=1,NSF)                                  XX 01130
      END IF                                                            XX 01140
C                                                                       XX 01150
 124  REWIND I16                                                        XX 01170
      READ(I16,16) M,N                                                  XX 01180
      IF(M.NE.NA) THEN                                                  XX 01190
          WRITE(IOUT,17)                                                XX 01200
          IFLAG=2                                                       XX 01210
          RETURN                                                        XX 01220
      END IF                                                            XX 01230
      READ(I16,18) ((F(M,N),N=1,NSX),M=1,NSX)                           XX 01240
      DO 126  M=1,NSX                                                   XX 01250
      DO 126  N=1,M-1                                                   XX 01260
      XX=(F(M,N)+F(N,M))/TWO                                            XX 01270
      F(M,N)=XX                                                         XX 01280
 126  F(N,M)=XX                                                         XX 01290
      WRITE(IOUT,19)                                                    XX 01300
      CALL TABLE1(NSX,NSX,F)                                            XX 01310
C                                                                       XX 01320
      DO 128  J=1,NSF                                                   XX 01330
      DO 130  I=1,NSX                                                   XX 01340
 130  IA(I,J)=0                                                         XX 01350
      DO 132  I=1,NSC(J,0)                                              XX 01360
 132  IA(NSC(J,I),J)=1                                                  XX 01370
 128  CONTINUE                                                          XX 01380
C  FORM G(1/2) MATRICES FOR ISOTOPOMERS.                                XX 01390
      NV=NSX*(NSX+1)/2                                                  XX 01400
      DO 150  NN=1,NISO                                                 XX 01410
        II=0
        DO 152  M=1,NSX
        DO 152  N=1,M
        II=II+1
        XU(II)=ZERO
          DO 154  I=1,NC
          IJ=(I-1)/3+1
 154      XU(II)=XU(II)+BS(M,I)*BS(N,I)/XM(IJ,NN)
 152    CONTINUE
C    (NOTE THAT XU AND XS ARE IN THE SAME MEMORY LOCATION.)             XX 01520
      CALL RSP(NC,NSX,NV,XU,W,1,XT,S1,S2)                               XX 01580
      WRITE(IOUT,20) NN                                                 XX 01590
      IGE=MIN0(5,NSX)                                                   XX 01600
      WRITE(IOUT,21) (W(I),I=1,IGE)                                     XX 01610
      DO 156  I=1,NSX                                                   XX 01620
 156  W(I)=DSQRT(W(I))                                                  XX 01630
      DO 158  N=1,NSX                                                   XX 01640
      DO 158  M=1,NSX                                                   XX 01650
         XX=ZERO                                                        XX 01660
         DO 160  K=1,NSX                                                XX 01670
 160     XX=XX+XT(M,K)*XT(N,K)*W(K)                                     XX 01680
 158  GH(M,N,NN)=XX                                                     XX 01690
 150  CONTINUE                                                          XX 01730
      IF(NOPT.EQ.0) GO TO 400                                           XX 01740
C  BEGIN ITERATIONS FOR SCALE FACTORS                                   XX 01750
      NTER=1                                                            XX 01760
      IFD=1
      IFS=1
      XL=XL00                                                           XX 01770
      IOP=1                                                             XX 01780
      SUM0=ZERO
      WRITE(IOUT,36)                                                    XX 01790
 200  CONTINUE                                                          XX 01800
      IF(NH.EQ.-1.AND.IFS.NE.0) THEN
         IF(IFS.EQ.1) THEN
            FA(IFD)=FA(IFD)+FDISP
         ELSE IF (IFS.EQ.-1) THEN
            FA(IFD)=FA(IFD)-TWO*FDISP
         END IF
      END IF
      DO 202  I=1,NSF                                                   XX 01810
      DO 202  K=1,NSC(I,0)                                              XX 01820
      J=NSC(I,K)                                                        XX 01830
 202  S1(J)=DSQRT(FA(I))                                                XX 01840
      SUM=ZERO                                                          XX 01890
      SUM1=ZERO                                                         XX 01900
      DO 204  N=1,NSF                                                   XX 01850
      IF(ABS(NH).GT.1) THEN
      DO 206  M=1,NSF                                                   XX 01860
      C(M,N)=ZERO
 206  HF(M,N)=ZERO                                                      XX 01870
      END IF
 204  G(N)=ZERO
C
C Form gradients and Hessian (if NH = 2).
C
      DO 210  NN=1,NISO                                                 XX 01910
      DO 212  N=1,NSX                                                   XX 01920
      DO 212  M=1,NSX                                                   XX 01930
           XX=ZERO                                                      XX 01940
           DO 214  J=1,NSX                                              XX 01950
 214       XX=XX+F(M,J)*GH(J,N,NN)*S1(J)                                XX 01960
 212  XT(M,N)=XX                                                        XX 01970
      II=0                                                              XX 01980
      DO 216  M=1,NSX                                                   XX 01990
      DO 216  N=1,M                                                     XX 02000
      II=II+1                                                           XX 02010
           XX=ZERO                                                      XX 02020
           DO 220  J=1,NSX                                              XX 02030
 220       XX=XX+GH(M,J,NN)*S1(J)*XT(J,N)                               XX 02040
 216  XU(II)=XX                                                         XX 02050
      CALL RSP(NC,NSX,NV,XU,S2,1,XT,W,S3)                               XX 02060
      DO 230  I=1,NSX                                                   XX 02070
           W(I)=DSIGN(DSQRT(DABS(S2(I))),S2(I))*WAVE                    XX 02100
           S3(I)=WT(I,NN)*(ONE-V(I,NN)/W(I))
           SUM=SUM+WT(I,NN)*(W(I)-V(I,NN))**2
 230       SUM1=SUM1+WT(I,NN)*W(I)*(W(I)-V(I,NN)/TWO)                   XX 02200
      DO 236  M=1,NSX                                                   XX 03600
      DO 236  N=1,NSX                                                   XX 03610
           XS(M,N)=ZERO                                                 XX 03620
           DO 236 L=1,NSX                                               XX 03630
 236       XS(M,N)=XS(M,N)+GH(M,L,NN)*XT(L,N)                           XX 03640
      IF(ABS(NH).GT.1) GO TO 260
C  NH = 0, 1, or -1 section.
      DO 240  M=1,NSX
      DO 240  N=1,M
        XR(M,N)=ZERO
        DO 242  I=1,NSX
 242    XR(M,N)=XR(M,N)+S3(I)*XS(M,I)*XS(N,I)
        XR(M,N)=XR(M,N)*F(M,N)
 240    XR(N,M)=XR(M,N)
      DO 244  N=1,NSX
      DO 244  M=1,NSX
 244    XS(M,N)=S1(N)/(TWO*S1(M))
      DO 246  MA=1,NSF
        DO 248  M=1,NSX                                                 XX 02390
        DO 248  N=1,NSX                                                 XX 02400
        IF(IA(M,MA)+IA(N,MA).EQ.2) THEN                                 XX 02420
           G(MA)=G(MA)+XR(M,N)                                          XX 02430
        ELSE IF(IA(M,MA).EQ.1.AND.IA(N,MA).EQ.0) THEN                   XX 02450
           G(MA)=G(MA)+XR(M,N)*XS(M,N)                                  XX 02460
        ELSE IF(IA(M,MA).EQ.0.AND.IA(N,MA).EQ.1) THEN                   XX 02480
           G(MA)=G(MA)+XR(M,N)*XS(N,M)                                  XX 02490
        END IF                                                          XX 02500
 248    CONTINUE                                                        XX 02510
 246  CONTINUE
      GO TO 210
C       NH=2 and NH=-2 section
 260  CONTINUE
      DO 262  N=1,NSX
      DO 262  M=1,NSX
 262    XT(M,N)=S1(N)/(TWO*S1(M))
      DO 264  MA=1,NSF
          DO 268  M=1,NSX                                               XX 02390
          DO 268  N=1,NSX                                               XX 02400
            D(N,M)=ZERO                                                 XX 02410
            IF(IA(M,MA)+IA(N,MA).EQ.2) THEN                             XX 02420
               D(N,M)=F(M,N)                                            XX 02430
            ELSE IF(IA(M,MA).EQ.1.AND.IA(N,MA).EQ.0) THEN               XX 02450
               D(N,M)=F(M,N)*XT(M,N)                                    XX 02460
            ELSE IF(IA(M,MA).EQ.0.AND.IA(N,MA).EQ.1) THEN               XX 02480
               D(N,M)=F(M,N)*XT(N,M)                                    XX 02490
            END IF                                                      XX 02500
 268        CONTINUE
        DO 270  J=1,NSX
        DO 270  M=1,NSX
        XX=ZERO
        DO 272  N=1,NSX
 272      XX=XX+D(N,M)*XS(N,J)
 270    XR(M,J)=XX
        DO 276  I=1,NSX
        DO 276  J=1,I
        XX=ZERO
        DO 274  M=1,NSX
 274      XX=XX+XS(M,I)*XR(M,J)
 276    TA(I,J,MA)=XX
 264  CONTINUE
C
      XX=WAVE*WAVE/TWO
      DO 278  MA=1,NSF
      DO 278  I=1,NSX
      JA(I,MA,NN)=TA(I,I,MA)*XX/W(I)
 278  G(MA)=G(MA)+S3(I)*TA(I,I,MA)
      IF(NH.EQ.-2) GO TO 289
C
      DO 280  M=1,NSX
      DO 280  N=1,M
        XR(M,N)=ZERO
        DO 282  I=1,NSX
 282    XR(M,N)=XR(M,N)+S3(I)*XS(M,I)*XS(N,I)
        XR(M,N)=XR(M,N)*F(M,N)
 280    XR(N,M)=XR(M,N)
      DO 284  M=1,NSX
      DO 284  N=1,NSX
        XS(M,N)=S1(N)/(FOUR*S1(M)**3)
 284    XT(M,N)=ONE/(FOUR*S1(M)*S1(N))
C
      DO 286  MA=1,NSF
      DO 286  MB=1,MA
        DO 288  M=1,NSX                                                 XX 02390
        DO 288  N=1,M                                                   XX 02400
        IR=IA(M,MA)+IA(N,MA)                                            XX 02410
        IS=IA(M,MB)+IA(N,MB)                                            XX 02410
        IT=IR*IS
        XX=ZERO
        IF(MA.EQ.MB) THEN
           IF(IA(M,MA).EQ.1.AND.IA(N,MA).EQ.0) THEN
             XX=-XR(M,N)*XS(M,N)                                        XX 02460
           ELSE IF(IA(M,MA).EQ.0.AND.IA(N,MA).EQ.1) THEN
             XX=-XR(M,N)*XS(N,M)                                        XX 02460
           END IF
        ELSE
           IF(IT.NE.0) THEN
             XX=XR(M,N)*XT(M,N)
           END IF
        END IF
        IF(M.NE.N) XX=TWO*XX
 288    HF(MA,MB)=HF(MA,MB)+XX
 286    CONTINUE
C
 289  DO 290  I=1,NSX
      XU(I+NSX)=WT(I,NN)/(TWO*S2(I))
 290  XU(I)=XU(I+NSX)*V(I,NN)/W(I)
      IF(NH.EQ.-2) GO TO 293
C  Note that forming XU overwrites XS
      DO 292  I=2,NSX
      DO 292  J=1,I-1
        IF(DABS(W(I)-W(J)).LT.VCUT) THEN
          XR(I,J)=ZERO
        ELSE
          XR(I,J)=TWO*(S3(I)-S3(J))/(S2(I)-S2(J))
        END IF
 292  CONTINUE
 293  DO 294  MA=1,NSF
      DO 294  MB=1,MA
        XX=ZERO
        XY=ZERO
        DO 296  I=1,NSX
        XZ=TA(I,I,MA)*TA(I,I,MB)
        XX=XX+XZ*XU(I)
 296    XY=XY+XZ*XU(I+NSX)
        IF(NH.GE.2) THEN
          HF(MA,MB)=HF(MA,MB)+XX
        ELSE IF (NH.EQ.-2) THEN
          HF(MA,MB)=HF(MA,MB)+XY
        END IF
        C(MA,MB)=C(MA,MB)+XY/TWO
 294  CONTINUE
      IF(NH.EQ.-2) GO TO 210
      DO 297  MA=1,NSF
      DO 297  MB=1,MA
        XX=ZERO
        DO 298  J=1,NSX
        DO 298  I=J+1,NSX
 298    XX=XX+TA(I,J,MA)*TA(I,J,MB)*XR(I,J)
        HF(MA,MB)=HF(MA,MB)+XX
 297  CONTINUE
 210  CONTINUE
C
C  Multiply by conversion factors to obtain final derivatives.
C
      XX=WAVE*WAVE
      DO 304  MA=1,NSF
 304  G(MA)=G(MA)*XX
C
      IF(NH.EQ.-1.AND.IFS.NE.0) THEN
         IF(IFS.EQ.1) THEN
            DO 308  MA=1,NSF
 308        HF(MA,IFD)=G(MA)
            IFS=-1
         ELSE IF (IFS.EQ.-1) THEN
            DO 310  MA=1,NSF
 310        HF(MA,IFD)=(HF(MA,IFD)-G(MA))/(TWO*FDISP)
            IFS=1
            FA(IFD)=FA(IFD)+FDISP
            IFD=IFD+1
         END IF
         IF(IFD.GT.NSF) THEN
            IFS=0
         END IF
         GO TO 200
      END IF
C
      IF(ABS(NH).GT.1) THEN
      DO 306  MA=1,NSF
      DO 306  MB=1,MA
      C(MA,MB)=C(MA,MB)*XX
      HF(MA,MB)=HF(MA,MB)*XX
      C(MB,MA)=C(MA,MB)
 306  HF(MB,MA)=HF(MA,MB)
      END IF
C
C  Scale factor update section                                          XX 02530
C
      IF(ABS(NH).GT.1.OR.NH.EQ.-1) GO TO 350
      IF(IOP.EQ.1) THEN                                                 XX 02540
        IF(NTER.EQ.1.AND.NH.EQ.0) THEN                                  XX 02550
          DO 314  I=1,NSF                                               XX 02560
  314     HF(I,I)=ONE/SUM1                                              XX 02570
        ELSE IF(NTER.GT.1) THEN                                         XX 02580
          DO 320  I=1,NSF                                               XX 02590
  320     W(I)=G(I)-G0(I)                                               XX 02600
          A0=ZERO                                                       XX 02610
          A1=ZERO                                                       XX 02620
          DO 324  I=1,NSF                                               XX 02630
          S2(I)=ZERO                                                    XX 02640
          DO 322  J=1,NSF                                               XX 02650
  322     S2(I)=S2(I)+HF(I,J)*W(J)                                      XX 02660
          A0=A0+S2(I)*W(I)                                              XX 02670
  324     A1=A1+DF(I)*W(I)                                              XX 02680
          DO 326  J=1,NSF                                               XX 02690
          DO 326  I=1,NSF                                               XX 02700
  326     HF(I,J)=HF(I,J)+XL*DF(I)*DF(J)/A1-S2(I)*S2(J)/A0              XX 02710
        END IF                                                          XX 02720
        A0=SUM                                                          XX 02730
        A1=ZERO                                                         XX 02740
        SNORM=ZERO                                                      XX 02750
        FNORM=ZERO                                                      XX 02750
        DO 330  I=1,NSF                                                 XX 02760
        S2(I)=G(I)*FA(I)/SUM                                            XX 02770
        SNORM=SNORM+S2(I)*S2(I)                                         XX 02780
        G0(I)=G(I)                                                      XX 02800
        F0(I)=FA(I)                                                     XX 02810
        DF(I)=ZERO                                                      XX 02820
        DO 332  J=1,NSF                                                 XX 02830
  332   DF(I)=DF(I)-HF(I,J)*G(J)                                        XX 02840
        A1=A1+DF(I)*G(I)                                                XX 02850
        XX=XL*DF(I)
        FNORM=FNORM+XX*XX
  330   FA(I)=FA(I)+XX                                                  XX 02860
        IOP=2                                                           XX 02870
        SNORM=DSQRT(SNORM)                                              XX 02880
        FNORM=DSQRT(FNORM)                                              XX 02880
        WRITE(IOUT,30)                                                  XX 02890
        WRITE(IOUT,31) NTER,SUM,A1,FNORM                                XX 02900
        WRITE(IOUT,32)                                                  XX 02910
        WRITE(IOUT,33) (S2(I),I=1,NSF)                                  XX 02920
        IF(SNORM.GT.GCONV.OR.DABS(ONE-SUM0/SUM).GT.SCONV) THEN          XX 02930
           SUM0=SUM
           GO TO 200
        END IF
      ELSE IF (IOP.EQ.2) THEN                                           XX 02940
        XX=ZERO                                                         XX 02950
        DO 334  I=1,NSF                                                 XX 02960
  334   XX=XX+DF(I)*G(I)                                                XX 02970
        A3=XL*(XX-A1)-TWO*(SUM-A0-XL*A1)                                XX 02980
        A3=A3/XL**3                                                     XX 02990
        A2=(XX-A1-THREE*A3*XL*XL)/(TWO*XL)                              XX 03000
        DISC=A2*A2-THREE*A1*A3                                          XX 03010
        XLN=(-A2+DSQRT(DABS(DISC)))/(THREE*A3)                          XX 03020
        IF(DISC.LT.ZERO.OR.XLN.LT.XLCUT) THEN                           XX 03030
          WRITE(IOUT,34) DISC,XLN                                       XX 03040
          GO TO 336                                                     XX 03050
        END IF                                                          XX 03060
        XL=XLN                                                          XX 03070
  336   WRITE(IOUT,35) XL                                               XX 03080
        DO 338  I=1,NSF                                                 XX 03090
  338   FA(I)=F0(I)+XL*DF(I)                                            XX 03100
        IOP=1                                                           XX 03110
        NTER=NTER+1                                                     XX 03120
        IF(NTER.LT.MAXIT) GO TO 200                                     XX 03130
      END IF                                                            XX 03140
      GO TO 380
C   NH = 2 and NH = -2 section
  350 CONTINUE
      WRITE(ICHECK,60) NTER
      DO 351  I=1,NSF
  351 WRITE(ICHECK,61) (HF(I,J),J=1,NSF)
      DO 352  J=1,NSF
      DO 352  I=1,NSF
      XR(I,J)=HF(I,J)
  352 XR(I,J+NSF)=ZERO
      DO 354  I=1,NSF
  354 XR(I,I+NSF)=ONE
C  Note that the above filling of XR overwrites XT.
C  Save HF before inversion for later use:
      II=0
      DO 353  I=1,NSF
      DO 353  J=1,I
      II=II+1
  353 XU(II)=HF(I,J)
      CALL FLIN(XR,NC,NSF,NSF,DET)
      DO 356  J=1,NSF
      DO 356  I=1,NSF
  356 HF(I,J)=XR(I,J+NSF)
      SNORM=ZERO
      FNORM=ZERO
      A1=ZERO
      DO 360  I=1,NSF
        S2(I)=G(I)*FA(I)/SUM                                            XX 02770
        SNORM=SNORM+S2(I)*S2(I)                                         XX 02780
        DF(I)=ZERO
        DO 358  J=1,NSF
  358   DF(I)=DF(I)-HF(I,J)*G(J)
        A1=A1+DF(I)*G(I)
        FA(I)=FA(I)+DF(I)
  360   FNORM=FNORM+DF(I)*DF(I)
        SNORM=DSQRT(SNORM)                                              XX 02880
        FNORM=DSQRT(FNORM)                                              XX 02880
        WRITE(IOUT,30)                                                  XX 02890
        WRITE(IOUT,31) NTER,SUM,A1,FNORM                                XX 02900
        WRITE(IOUT,32)                                                  XX 02910
        WRITE(IOUT,33) (S2(I),I=1,NSF)                                  XX 02920
        IF(SNORM.GT.GCONV2.OR.DABS(ONE-SUM0/SUM).GT.SCONV2) THEN        XX 02930
           IF(NTER.LT.MAXIT) THEN
             SUM0=SUM
             NTER=NTER+1
             IF(NH.EQ.-1) THEN
               IFS=1
               IFD=1
             END IF
             GO TO 200
           END IF
        END IF
  380 IF(NTER.GE.MAXIT) THEN                                            XX 03150
        IFLAG=3                                                         XX 03160
        RETURN                                                          XX 03170
      END IF                                                            XX 03180
      WRITE(IOUT,38)                                                    XX 03190
      WRITE(IOUT,39) (HF(I,I),I=1,NSF)                                  XX 03200
 400  CONTINUE                                                          XX 03220
C  Statistics section
      IF(NOPT.EQ.1.AND.ABS(NH).GE.2) THEN
        DO 502  MA=1,NSF
        DO 502  NA=1,NSF
           D(MA,NA)=ZERO
           DO 502  LA=1,NSF
 502       D(MA,NA)=D(MA,NA)+HF(MA,LA)*C(LA,NA)
        XX=ZERO
        XY=ZERO
        DO 504  MA=1,NSF
        XX=XX+D(MA,NA)
        DO 504  NA=1,NSF
 504    XY=XY+D(MA,NA)*D(MA,MA)
        DO 506  MA=1,NSF
        DO 506  NA=1,NSF
           C(MA,NA)=ZERO
           DO 506  LA=1,NSF
 506       C(MA,NA)=C(MA,NA)+D(MA,LA)*HF(LA,NA)
        XZ=FOUR*(XY-XX)
        STX=DSQRT(SUM/(NFT-XZ))
        DO 508  MA=1,NSF
 508    S2(MA)=DSQRT(C(MA,MA))*STX*TWO
        DO 510  MA=2,NSF
        DO 510  NA=1,MA-1
        XR(MA,NA)=C(MA,NA)/DSQRT(C(MA,MA)*C(NA,NA))
 510    XR(NA,MA)=XR(MA,NA)
        DO 512  MA=1,NSF
 512    XR(MA,MA)=ONE
        WRITE(IOUT,49)                                                  XX 03270
        WRITE(IOUT,50) (I,FA(I),S2(I),I=1,NSF)                          XX 03280
        WRITE(IOUT,51) STX,NFT-XZ,NFT-NSF                               XX 03270
        WRITE(IOUT,52)                                                  XX 03280
        CALL TABLE1(NC,NSF,XR)
        NVF=NSF*(NSF+1)/2
        CALL RSP(NC,NSF,NVF,XU,W,1,XT,S2,S3)
        WRITE(IOUT,53)
        CALL TABLE4(NC,NSF,NSF,W,XT)
      ELSE
        WRITE(IOUT,23)                                                  XX 03270
        WRITE(IOUT,24) (I,FA(I),I=1,NSF)                                XX 03280
      END IF
C  Scale force constants; Perform vibrational analysis for isotopomers  XX 03210
      DO 404  I=1,NSF                                                   XX 03230
      DO 404  K=1,NSC(I,0)                                              XX 03240
      J=NSC(I,K)                                                        XX 03250
 404  S1(J)=DSQRT(FA(I))                                                XX 03260
      DO 406  J=1,NSX                                                   XX 03290
      DO 406  I=1,NSX                                                   XX 03300
 406  F(I,J)=F(I,J)*S1(I)*S1(J)                                         XX 03310
      WRITE(IOUT,22)                                                    XX 03320
      CALL TABLE1(NSX,NSX,F)                                            XX 03330
C                                                                       XX 03340
      DO 408 M=1,NSX                                                    XX 03350
      DO 408 N=1,M                                                      XX 03360
       IF(DABS(F(M,N)).GT.F2CUT) WRITE(ICHECK,44) M,N,F(M,N)            XX 03370
 408  CONTINUE                                                          XX 03380
C                                                                       XX 03390
      DO 410  NN=1,NISO                                                 XX 03400
      DO 414  N=1,NSX                                                   XX 03410
      DO 414  M=1,NSX                                                   XX 03420
           XX=ZERO                                                      XX 03430
           DO 416  J=1,NSX                                              XX 03440
 416       XX=XX+F(M,J)*GH(J,N,NN)                                      XX 03450
 414  XT(M,N)=XX                                                        XX 03460
      II=0                                                              XX 03470
      DO 418  M=1,NSX                                                   XX 03480
      DO 418  N=1,M                                                     XX 03490
      II=II+1                                                           XX 03500
           XX=ZERO                                                      XX 03510
           DO 420  J=1,NSX                                              XX 03520
 420       XX=XX+GH(M,J,NN)*XT(J,N)                                     XX 03530
 418  XU(II)=XX                                                         XX 03540
      CALL RSP(NC,NSX,NV,XU,W,1,XT,XR(1,1),XR(1,2))                     XX 03550
      DO 430  I=1,NSX                                                   XX 03560
           XX=DSIGN(DSQRT(DABS(W(I))),W(I))                             XX 03570
 430       W(I)=XX*WAVE                                                 XX 03590
      DO 436  M=1,NSX                                                   XX 03600
      DO 436  N=1,NSX                                                   XX 03610
           XS(M,N)=ZERO                                                 XX 03620
           DO 440 L=1,NSX                                               XX 03630
 440       XS(M,N)=XS(M,N)+GH(M,L,NN)*XT(L,N)                           XX 03640
 436  CONTINUE                                                          XX 03650
      WRITE(IOUT,25) NN                                                 XX 03750
      CALL TABLE2(NC,NSX,NSX,W,XS)                                      XX 03760
      IF(NMODE.EQ.1) THEN
        DO 445  N=1,NSX                                                 XX 03660
             XX=ZERO                                                    XX 03670
             DO 446  M=1,NSX                                            XX 03680
             XR(M,N)=XS(M,N)*F(M,M)*XS(M,N)                             XX 03690
 446         XX=XX+DABS(XR(M,N))                                        XX 03700
             IF(XX.EQ.ZERO) GO TO 445                                   XX 03710
             DO 447  M=1,NSX                                            XX 03720
 447         XR(M,N)=ONEH*XR(M,N)/XX                                    XX 03730
 445    CONTINUE                                                        XX 03740
      ELSE
        DO 450  N=1,NSX
        DO 450  M=1,NSX
 450         XR(M,N)=XS(M,N)
C   Note that the following step overwrites the XT array in memory.
        DO 451  N=1,NSX
        DO 451  M=1,NSX
 451         XR(M,N+NSX)=ZERO
        DO 452  M=1,NSX
 452         XR(M,M+NSX)=ONE
        CALL FLIN(XR,NC,NSX,NSX,DET)
        DO 453  N=1,NSX
        DO 453  M=1,NSX
 453         XR(M,N)=XS(M,N)*XR(N,M+NSX)*ONEH
      END IF
      IF(NMODE.EQ.1) THEN
        WRITE(IOUT,26) NN                                               XX 03770
      ELSE
        WRITE(IOUT,45) NN                                               XX 03770
      END IF
      CALL TABLE3(NC,NSX,W,XR)                                          XX 03780
      DO 460  J=1,NSX                                                   XX 03790
      IA(5,J)=4                                                         XX 03800
      DO 460  I=1,4                                                     XX 03810
      IA(I,J)=0                                                         XX 03820
 460  XT(I,J)=ZERO                                                      XX 03830
      DO 465  J=1,NSX                                                   XX 03840
          DO 470  I=1,NSX                                               XX 03850
          IF(DABS(XR(I,J)).GT.DABS(XT(1,J))) THEN                       XX 03860
             XT(1,J)=XR(I,J)                                            XX 03870
             IA(1,J)=I                                                  XX 03880
          END IF                                                        XX 03890
 470      CONTINUE                                                      XX 03900
          DO 472  I=1,NSX                                               XX 03910
          ITEST=I.NE.IA(1,J)                                            XX 03920
          IF(DABS(XR(I,J)).GT.DABS(XT(2,J)).AND.ITEST) THEN             XX 03930
             XT(2,J)=XR(I,J)                                            XX 03940
             IA(2,J)=I                                                  XX 03950
          END IF                                                        XX 03960
 472      CONTINUE                                                      XX 03970
          DO 474  I=1,NSX                                               XX 03980
          ITEST=I.NE.IA(1,J).AND.I.NE.IA(2,J)                           XX 03990
          IF(DABS(XR(I,J)).GT.DABS(XT(3,J)).AND.ITEST) THEN             XX 04000
             XT(3,J)=XR(I,J)                                            XX 04010
             IA(3,J)=I                                                  XX 04020
          END IF                                                        XX 04030
 474      CONTINUE                                                      XX 04040
          DO 476  I=1,NSX                                               XX 04050
          ITEST=I.NE.IA(1,J).AND.I.NE.IA(2,J).AND.I.NE.IA(3,J)          XX 04060
          IF(DABS(XR(I,J)).GT.DABS(XT(4,J)).AND.ITEST) THEN             XX 04070
             XT(4,J)=XR(I,J)                                            XX 04080
             IA(4,J)=I                                                  XX 04090
          END IF                                                        XX 04100
 476      CONTINUE                                                      XX 04110
          DO 478  K=1,4                                                 XX 04120
             IF(DABS(XT(K,J)).LT.EDCUT) THEN                            XX 04130
             IA(5,J)=K-1                                                XX 04140
             GO TO 465                                                  XX 04150
             END IF                                                     XX 04160
 478      CONTINUE                                                      XX 04170
 465  CONTINUE                                                          XX 04180
      DO 480  J=1,NSX                                                   XX 04190
          DO 485  K=1,4                                                 XX 04200
          XX=XS(IA(K,J),J)                                              XX 04210
          IF(XX.LT.ZERO) IA(K,J)=-IA(K,J)                               XX 04220
 485      CONTINUE                                                      XX 04230
          IF(IA(1,J).LT.0) THEN                                         XX 04240
          DO 490  K=1,4                                                 XX 04250
 490      IA(K,J)=-IA(K,J)                                              XX 04260
          END IF
 480      CONTINUE                                                      XX 04280
      WRITE(IOUT,27) NN                                                 XX 04290
      IF(NOPT.EQ.1.AND.ABS(NH).GE.2) THEN
        DO  550  I=1,NSX
        WT(I,NN)=ZERO
            DO 552  MA=1,NSF
            DO 552  NA=1,NSF
 552        WT(I,NN)=WT(I,NN)+JA(I,MA,NN)*JA(I,NA,NN)*C(MA,NA)
 550    WT(I,NN)=DSQRT(WT(I,NN))*TWO*STX
        IF(NMODE.EQ.1) THEN
          WRITE(IOUT,56)                                                XX 04300
        ELSE
          WRITE(IOUT,57)                                                XX 04300
        END IF
        DO 560  J=1,NSX                                                 XX 04310
        IF(V(J,NN).NE.ZERO) THEN                                        XX 04320
        XX=V(J,NN)-W(J)
        XY=XX/V(J,NN)*ONEH
        WRITE(IOUT,58) J,V(J,NN),W(J),XX,XY,WT(J,NN),                   XX 04330
     $                  (IA(K,J),XT(K,J),K=1,IA(5,J))                   XX 04330
        ELSE                                                            XX 04340
        WRITE(IOUT,59) J,W(J),WT(J,NN),(IA(K,J),XT(K,J),K=1,IA(5,J))    XX 04350
        END IF                                                          XX 04360
 560    CONTINUE                                                        XX 04370
      ELSE
        IF(NMODE.EQ.1) THEN
          WRITE(IOUT,28)                                                XX 04300
        ELSE
          WRITE(IOUT,46)                                                XX 04300
        END IF
        DO 495  J=1,NSX                                                 XX 04310
        IF(V(J,NN).NE.ZERO) THEN                                        XX 04320
        XX=V(J,NN)-W(J)
        XY=XX/V(J,NN)*ONEH
      WRITE(IOUT,29) J,V(J,NN),W(J),XX,XY,(IA(K,J),XT(K,J),K=1,IA(5,J)) XX 04330
        ELSE                                                            XX 04340
          WRITE(IOUT,40) J,W(J),(IA(K,J),XT(K,J),K=1,IA(5,J))           XX 04350
        END IF                                                          XX 04360
 495  CONTINUE                                                          XX 04370
      END IF
C
 410  CONTINUE                                                          XX 04380
      RETURN                                                            XX 04390
      END                                                               XX 04400
