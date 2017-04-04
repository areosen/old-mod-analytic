      PROGRAM MODELLERING
C**********************************************************************
C***  Program that derives on DIR.F. Purpose is to model the        ***
C***  pressure and pressuregradient on a receiver surface lying     ***
C***  below an array of sources.                                    ***
C**********************************************************************
C**********************************************************************
      INTEGER      NX,NT,NMAX,NW
      PARAMETER    (NX=512,NT=1024,NW=NT/2+1,NMAX=10)

      INTEGER      NNN,NSIG
      PARAMETER    (NNN=600,NSIG=8)
C----------------------------------------------------------------------
C     Variabeldeklarasjon
C----------------------------------------------------------------------
      REAL         A
      INTEGER      NWST,I,N,II
      REAL         DT,DX,VEL,R0,ZS,ZR,T0,EPS,XS,DW,PI,WNYQ
      REAL         SORT(NT),SORS(NT)
      REAL         X(NMAX),Z(NMAX)
      REAL         DATA1(NT,NX),DATA2(NT,NX)
      REAL         DMAT1(NT,NX),DMAT2(NT,NX)
      CHARACTER*60 FILUT1,FILUT2,SIGNFIL
      COMPLEX      CA1(NT),CA2(NT),CA3(NT),CA4(NT)
      COMPLEX      SOURCE(NT),WARR(NW)
      COMPLEX      J
      INTEGER      FLAG
      REAL         ZZ(NX)

C-----Nye variabler.....

      REAL    EXTRA(NNN,NSIG)
      REAL    SIGDATA(NT,5)

C----------------------------------------------------------------------
C     Begynner hovedprogram.
C----------------------------------------------------------------------
      CALL LES_INN(DT,DX,NWST,VEL,R0,N,X,Z,ZR,T0,EPS,FILUT1,
     +     FILUT2,NMAX,NW,A,FLAG,SIGNFIL)

      J    = (0.,1.)
      PI   = 2.*ACOS(0.)
      WNYQ = PI/DT
      DW   = 2.*PI/(NT*DT)

C-----Tester om 'nullfrekvens' mindre enn Nyqvistfrekvens
      IF (NW.GT.NT/2+1) THEN
         PRINT *,'NW MUST BE LESS THAN NT/2+1'
         PRINT *,'NW    = ',NW
         PRINT *,'NT/2+1= ',NT/2+1
         GOTO 999
      ENDIF

C-----Leser inn signaturene
      OPEN(UNIT=11,FILE=SIGNFIL,FORM='UNFORMATTED')
      PRINT *,'Leser fil med signatur : ',SIGNFIL
      DO 101 I=1,1
         READ(11) (SIGDATA(IT,I),IT=1,256)
 101  CONTINUE

C-----Oppretter vektor med alle frekvenser
      DO 100 IW=1,NW
         WARR(IW) = (IW-1)*DW
 100  CONTINUE

C-----Dersom FLAG=1; leser inn streamerdyp, som ellers er konstant
      IF (FLAG .EQ. 0) THEN
         PRINT *,'Konstant streamer-dybde.'
      ELSE
         IF (FLAG .EQ. 1) THEN
            PRINT *,'Variabel streamer-dybde.'
            OPEN(UNIT=44,FILE='./fdm/sd.d',FORM='FORMATTED')
            PRINT *,'LESER FIL MED STREAMDYP : ./fdm/sd.d'
            DO 10 I=1,NX
               READ(44,*) ZZ(I)
 10         CONTINUE      
         ELSE
            PRINT *,'Feil innlest FLAG.'
            STOP
         ENDIF
      ENDIF      

C-----Looper over alle kildene, beregner trykk & norma.trykk.
      DO 200 I=1,N
         II=I
         CALL SOURFU(SOURCE,WARR,NW,T0,EPS,VEL,
     +        NNN,NT,NSIG,EXTRA,SIGDATA,CA1,DT,II,N)
         XS = X(I)
         ZS = Z(I)
C-----Beregner for gitt (x,z) og signatur => DATA1&2
         CALL BEREGN(DATA1,DATA2,NX,NT,NWST,A,
     +        DT,DX,VEL,R0,ZR,ZS,NW,
     +        SORT,SORS,CA3,CA4,XS,WARR,SOURCE,FLAG,ZZ)
C-----Summerer DATA1&2 for alle kildene => DMAT1&2
         CALL SUM_DATA(DMAT1,DMAT2,DATA1,DATA2,NX,NT)
 200  CONTINUE
      PRINT *,'Currently listing out on files.'
      CALL LIST_OUT(DMAT1,DMAT2,SORT,SORS,FILUT1,FILUT2,
     +     NT,NX,FSORS,FSORT)
C----------------------------------------------------------------------
C     Avslutter hovedprogram
C----------------------------------------------------------------------
 999  CONTINUE
      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     SUBROUTINES
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


C**********************************************************************
C     Subrutinen LES_INN
C**********************************************************************
      SUBROUTINE LES_INN(DT,DX,NWST,VEL,R0,N,X,Z,ZR,T0,EPS,FILUT1,
     +     FILUT2,NMAX,NW,A,FLAG,SIGNFIL)

C----------------------------------------------------------------------
C     Variabeldeklarasjon
C----------------------------------------------------------------------
      INTEGER      NWST,NW,N,NMAX
      REAL         DT,DX,ZR,R0,VEL,T0,EPS,A
      REAL         X(NMAX),Z(NMAX)
      CHARACTER*60 FILUT1,FILUT2,SIGNFIL
      INTEGER      FLAG
C----------------------------------------------------------------------
C     Lokale variable
C----------------------------------------------------------------------
      INTEGER      IX

      OPEN(UNIT=10,FILE='mod.job',FORM='FORMATTED')
      PRINT *,'Leser parameter fra fil : mod.job'
      READ (10,*) FLAG
      READ (10,*) DT,DX
      READ (10,*) A
      READ (10,*) NWST
      READ (10,*) VEL
      READ (10,*) R0
      READ (10,*) N
      DO 100 IX=1,N
         READ (10,*) X(IX),Z(IX)
 100  CONTINUE
      READ (10,*) ZR
      READ (10,*) T0,EPS
      READ (10,*) FILUT1
      READ (10,*) FILUT2
      READ (10,*) SIGNFIL
      RETURN
      END

C***************************************************************
C***  Subrutine som beregner kildefunksjonen,                ***
C***  dvs. den returnerer en matrise, SOURCE, som            ***
C***  inneholder kildefunksjonen for de signifikante         ***
C***  frekvensene.                                           ***
C***************************************************************
      SUBROUTINE SOURFU(SOURCE,WARR,NW,T0,SIGMA,VEL,
     +     NNN,NT,NSIG,EXTRA,SIGDATA,CA1,DT,II,N)

      INTEGER NW
      REAL    T0,SIGMA,VEL
      COMPLEX SOURCE(NW),WARR(NW)

C-----Nye variabler......

      INTEGER NNN,NSIG,NT,II,N
      REAL    EXTRA(NNN,NSIG)
      REAL    SIGDATA(NT,5)
      COMPLEX CA1(NT)
      INTEGER NSIGN
      REAL    DT
      INTEGER IT,I,K      

C-----Lokale variable i subrutinen

      INTEGER IW
      REAL    PI,FACT,TEMP
      COMPLEX J,W

C-----Definerer kontanter

      DATA PI / 3.141592654 /
      J = (0.,1.)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     GAUSSISK KIDESIGNATUR
C
c      DO 100 IW=2,NW
c         W          = WARR(IW)
c         SOURCE(IW) = CEXP(-W*W/(4.*SIGMA))*CEXP(J*W*T0)
c         SOURCE(IW) = J*W*SOURCE(IW)*SQRT(PI)/SIGMA
c         SOURCE(IW) = SOURCE(IW)
c 100  CONTINUE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     ENHETS-IMPULS KIDESIGNATUR
C
c      DO 100 IW=2,NW
c         W          = WARR(IW)
c         SOURCE(IW) = ( CEXP(J*W*0.2) - 1)/(J*W)
c 100  CONTINUE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      FACT=1.0
      IF (II .GT. 5) THEN
         IF (N .GT. 10) THEN
            PRINT *,'WRONG NUMBER OF SOURCES...'
            STOP
         ENDIF
         II=II-5
         print *,'Give the last sources half the amplitude'
         FACT=0.5
      ENDIF
C-----Tar Fouriertransformen av signaturen nr. II
      NSIGN=1
      DO 105 IT=1,NT
        CA1(IT)=SIGDATA(IT,II)
 105  CONTINUE
      CALL FOUR(CA1,NT,NT,DT,NSIGN)
      DO 106 IW=1,NW
         SOURCE(IW)=FACT*CA1(IW)
 106  CONTINUE
      RETURN
      END

C***************************************************************
C     Subrutine som beregner dataverdier
C***************************************************************
      SUBROUTINE BEREGN(DATA1,DATA2,NX,NT,NWST,A,
     +        DT,DX,VEL,R0,ZR,ZS,NW,
     +        SORT,SORS,CA3,CA4,XS,WARR,SORW,FLAG,ZZ)

      INTEGER NX,NT,NWST,NW
      REAL    ZZ(NX)
      INTEGER FLAG
      REAL    DT,DX,VEL,R0,ZS,ZR,XS,A
      REAL    SORT(NT),SORS(NT)
      REAL    DATA1(NT,NX),DATA2(NT,NX)
      COMPLEX WARR(NT),SORW(NT),CA3(NT),CA4(NT)
C---------------------------------------------------------------
C     Lokale variable i subrutinen
C---------------------------------------------------------------
      REAL    PI,X,RM,RP,DW
      INTEGER N,NSIGN
      REAL    WNYQ
      COMPLEX J,K,CARG1,CARG2,CJ0,CY0,CJ1,CY1
C---------------------------------------------------------------
C     Definerer kontanter
C---------------------------------------------------------------
      PI   = 2.*ACOS(0.)
      J    = (0.,1.)
      WNYQ = PI/DT
      DW   = 2.*WNYQ/NT

C---------------------------------------------------------------
C     THIS PART IS TO PRINT OUT THE SOURCE CARACTERISTICS !!!
C
c      DO 10 IW=1,NW
c         CA3(IW)  = SORW(IW)
c 10   CONTINUE
c      DO 100 IW=NT/2+2,NT
c         SORW(IW) = CONJG(SORW(NT-IW+2))
c         CA3(IW)  = SORW(IW)
c 100  CONTINUE
c      DO 20 IW=1,NT
c         SORS(IW) = SQRT(SORW(IW)*CONJG(SORW(IW)))
c 20   CONTINUE
c      NSIGN = -1
c      CALL FOUR(CA3,NT,NT,DT,NSIGN)
c      PRINT *,'Source sig.-data file: ./signature.d'
c      OPEN(UNIT=25,FILE='./signature.d',FORM='UNFORMATTED')
c         WRITE(25) (REAL(CA3(IT)),IT=1,NT)
c      CLOSE(25)
C--------------------------------------------------------------
      DO 1000 IX=1,NX
         X = ix*DX
ccc         x = 1.
ccc + (ix-1)*2*dx
C-----Beregner streamer-dyp avh. av verdi av FLAG
ccc         RM = SQRT((X - XS)**2+(ZR - ZS)**2)
ccc         RP = SQRT((X - XS)**2+(ZR + ZS)**2)
         IF (FLAG .EQ. 0) THEN
            RM = SQRT(( X - (XS+NX*DX/2))**2+(ZR - ZS)**2)
            RP = SQRT(( X - (XS+NX*DX/2.))**2+(ZR + ZS)**2)
         ELSE
            IF (FLAG .EQ. 1) THEN
               RM = SQRT((X-(XS+(NX/2)*DX))**2+(ZZ(IX)-ZS)**2)
               RP = SQRT((X-(XS+(NX/2)*DX))**2+(ZZ(IX)+ZS)**2)
            ENDIF            
         ENDIF      
         DO 400 IW= NWST,NW
            K     = WARR(IW)/VEL
            CARG1 = K*RM
            CALL BSSLJ(CARG1,0,CJ0)
            CALL BSSLY(CARG1,0,CY0)
            CARG1 = CJ0 + J*CY0
            
            CARG2 = K*RP
            CALL BSSLJ(CARG2,0,CJ0)
            CALL BSSLY(CARG2,0,CY0)
            CARG2 = CJ0 + J*CY0
            
            CA3(IW) = - CARG1 + R0*CARG2
            
            CA3(IW) = CA3(IW)*J*SORW(IW)/4.
 400     CONTINUE
         DO 500 IW=NWST,NW
            K     = WARR(IW)/VEL
            CARG1 = K*RM
            CALL BSSLJ(CARG1,1,CJ1)
            CALL BSSLY(CARG1,1,CY1)
            CARG1 = CJ1+J*CY1
            
            CARG2 = K*RP
            CALL BSSLJ(CARG2,1,CJ1)
            CALL BSSLY(CARG2,1,CY1)
            CARG2 = CJ1+J*CY1
            
            CA4(IW) = (ZR - ZS)*CARG1/RM - R0*(ZR + ZS)*CARG2/RP
            CA4(IW)=  K*CA4(IW)*J*SORW(IW)/4.
 500     CONTINUE
         DO 600 IW=1,NWST-1
            CA3(IW) = (0.,0.)
            CA4(IW) = (0.,0.)
 600     CONTINUE
         DO 700 IW=NW+1,NT-NW+1
            CA3(IW) = (0.,0.)
            CA4(IW) = (0.,0.)
 700     CONTINUE
         DO 800 IW=NT-NW+2,NT
            CA3(IW) = CONJG(CA3(NT-IW+2))
            CA4(IW) = CONJG(CA4(NT-IW+2))
 800     CONTINUE
         NSIGN = -1
         CALL FOUR(CA3,NT,NT,DT,NSIGN)
         CALL FOUR(CA4,NT,NT,DT,NSIGN)
C---------------------------------------------------------------
C     Overfoerer beregnvede verdier til matrisa DATA
C---------------------------------------------------------------
         DO 900 IT=1,NT
            DATA1(IT,IX) = REAL(CA3(IT))
            DATA2(IT,IX) = REAL(CA4(IT))
 900     CONTINUE
 1000 CONTINUE

      RETURN
      END

C***************************************************************
C     Subrutinen SUM_DATA
C***************************************************************
      SUBROUTINE SUM_DATA(DMAT1,DMAT2,DATA1,DATA2,NX,NT)

      REAL    DATA1(NT,NX),DATA2(NT,NX)
      REAL    DMAT1(NT,NX),DMAT2(NT,NX)
      INTEGER NX,NT
      
      INTEGER IT,IX
     
      DO 100 IT=1,NT
         DO 200 IX=1,NX
            DMAT1(IT,IX) = DMAT1(IT,IX) + DATA1(IT,IX)
            DMAT2(IT,IX) = DMAT2(IT,IX) + DATA2(IT,IX)
 200     CONTINUE
 100  CONTINUE

      RETURN
      END

C***************************************************************
C     Subrutinen LIST_OUT, som skriver ut paa filer
C***************************************************************
      SUBROUTINE LIST_OUT(DMAT1,DMAT2,SORT,SORS,FILUT1,
     +     FILUT2,NT,NX,FSORS,FSORT)

      REAL         DMAT1(NT,NX),DMAT2(NT,NX)
      INTEGER      NX,NT
      REAL         SORT(NT),SORS(NT)
      CHARACTER*60 FILUT1,FSORT,FSORS,FILUT2

      PRINT *,'LIST OVER PROGRAM OUTPUT:'
      PRINT *,'-------------------------------------'
      
      PRINT *,'Output press.-data file        : ',FILUT1
      OPEN(UNIT=15,FILE=FILUT1,FORM='UNFORMATTED')
      DO 103 IX=1,NX
         WRITE(15) (DMAT1(IT,IX),IT=1,NT)
 103  CONTINUE
      
      PRINT *,'Output norm.deriv.-data file   : ',FILUT2
      OPEN(UNIT=16,FILE=FILUT2,FORM='UNFORMATTED')
      DO 104 IX=1,NX
         WRITE(16) (DMAT2(IT,IX),IT=1,NT)
 104  CONTINUE

c      PRINT *,'Output source spectrum data file   : 
c     +        /konaos/Data/moddata/sors.d'
c      OPEN(UNIT=19,FILE='/homes/konaos/Data/
c     +        moddata/sors.d',FORM='UNFORMATTED')
c      WRITE(19) (SORS(IT),IT=1,NT)
      RETURN
      END
      











