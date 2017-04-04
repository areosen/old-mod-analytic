      PROGRAM AMODFS
C********************************************************************
C     This program models a strip reflector and gets the corresponding
C     field record. Free surface is included.
C********************************************************************
C     Calls subroutine SOURCE, returns source function values
C     in the common block FF.
C      
C     * V(NX,NZ) contains the velocity values in the coarse grid.
C     * P(NX,NZ)    ... pressure values at grid nodes
C     * PM1(NX,NZ)  ...
C     * PM1(NX,NZ)  ...
C     * TRACE(NT,NX)... Time series at some (?) nodes
C
C     Stability criterion : del t < (del x / sqrt(2).* maxvelocity)
C********************************************************************
      parameter(NT=5000)
      parameter(NX=1000)
      parameter(NZ=NX)

      dimension V(0:NX,0:NZ)
      dimension T(0:NT,0:NX),t1(0:nt,0:nx)
      dimension FF(200)
      dimension P(0:NX,0:NZ),PM1(0:NX,0:NZ),PM2(0:NX,0:NZ)

      common FF

      real    TRACE(0:NT,0:NX)
      integer P1,P2,P3,P4,P5,DTMS,DZMS
      real    W(0:NT)
C-------------------------------------
C     Set up the parameters
C-------------------------------------
      NX2=NX/2
      NX2P1=NX2+1
      NZ2=NZ/2
      NZ2P1=NZ2+1
      DELX=2.5
C-----This alogrithm work well only for  delz=delx
      DELZ=DELX
      DELT=.001
      NSX=NX/2+2
      NSZ=3
      NL1=NSZ+10
      NL2=NSZ+20
C-----These paramteters seem not to enter the program
C      NWL=NSX-20
C      DTMS=DELT*1000.
C      DZMS=DELZ*1000.
C      CLIP=4.

C--------------------------------------
C     Set up the structure of the model
C--------------------------------------
      PRINT *,'Set up the structure of the model'
C-----Layer one
      DO 20 J=0,NL1
         DO 30 K=0,NX
            V(K,J)=1500.
 30      CONTINUE
 20   CONTINUE
C-----Layer two
      DO 40 J=NL1+1,NL2
         DO 50 K=0,NX
      	    V(K,J)=1600.
 50      CONTINUE
C-----Possible strip reflector in layer two.
C         DO 60 K=NWL+1,NX
C            V(K,J)=1600.
C 60      CONTINUE
 40   CONTINUE
C-----Layer three
      DO 70 J=NL2+1,NZ
         DO 80 K=0,NX
            V(K,J)=1700.
 80      CONTINUE
 70   CONTINUE
C--------------------------------------
C     Put a gaussian wavelet on the data
C--------------------------------------
      DO 90 I=0,NT
         W(I)=0.
 90   CONTINUE
C--------------------------------------
C     Call subroutine. SOURCE
C--------------------------------------
      CALL SOURCE
      PRINT *,'Put a gaussian wavelet on the data'
      DO 91 I=1,200
         W(I-1)=FF(I)
 91   CONTINUE
C--------------------------------------
C     Now do the model (first initiate PM1/2)
C--------------------------------------
C-----Set up the first slice (PM2)
      DO 6 J=0,NZ
         DO 61 K=0,NX
            PM2(K,J)=0.
 61      CONTINUE
 6    CONTINUE

C-----Setup the second slice (PM1)
      DO 7 J=0,NZ
         DO 71 K=0,NX
            PM1(K,J)=0.
 71      CONTINUE
 7    CONTINUE
C--------------------------------------
C     Now ready for the general case
C--------------------------------------
C-----Begin loop over time
      PRINT *,'Begin loop over time'
      DO 100 I=0,nt
C--------------------------------------
C     Now calculate the field P
C--------------------------------------
         J=0
         DO 118 K=1,NX-1
            P(K,J)=0.
 118     CONTINUE
C-----Set endpoint pressure values to zero
         P(0,J)=0.
         P(NX,J)=0.
         DO 17 J=1,NZ-1
            DO 18 K=1,NX-1
               A=V(K,J)*DELT/DELX
               A2=A*A
               FACTOR=1.-2.*A2
               P(K,J)=2.*FACTOR*PM1(K,J)-PM2(K,J)
               SUM=PM1(K+1,J)+PM1(K-1,J)+PM1(K,J+1)+PM1(K,J-1)
               P(K,J)=P(K,J)+A2*SUM
 18         CONTINUE
C-----Set endpoint pressure values to zero
            P(0,J)=0.
            P(NX,J)=0.
 17      CONTINUE
C-----Set grid bottom pressure values to zero
         DO 31 K=0,NX
            P(K,NZ)=0.
 31      CONTINUE
         
C-----Add in a point source
         P(NSX,NSZ)=P(NSX,NSZ)+W(I)
         
C-----pick off data at z=nsz+3
         DO 41 K=0,NX
            T(I,K)=P(K,NSZ+3)
 41      CONTINUE

C-----pick off data at z=nsz+4
         DO 42 K=0,NX
            T1(I,K)=P(K,NSZ+4)
 42      CONTINUE         
C--------------------------------------
C     Shift data to  PM1 and PM2
C--------------------------------------
         DO 11 J=0,NZ
            DO 12 K=0,NX
               PM2(K,J)=PM1(K,J)
               PM1(K,J)=P(K,J)
 12         CONTINUE
 11      CONTINUE
CCC         DO 14 J=0,NZ
CCC            DO 15 K=0,NX
CCC               PM1(K,J)=P(K,J)
CCC 15         CONTINUE
CCC 14      CONTINUE
         P1=NT/5
         P2=2*NT/5
         P3=3*NT/5
         P4=4*NT/5
         P5=NT
         IF(I.EQ.P1.OR.I.EQ.P2.OR.I.EQ.P3.OR.I.EQ.P4) THEN
            NFLAG=0
            IF(I.NE.P1) THEN
               NFLAG=1
            ENDIF
            WMX=0.
            DO 221 K=0,NX
               DO 222 J=0,NZ
                  TRACE(J,K)=P(K,J)
 222           CONTINUE
 221        CONTINUE
            RMAX=1.
            NSAMP=NZ+1
         ENDIF
C-----End loop over time
 100  CONTINUE
C--------------------------------------
C     Convert and print out on files
C--------------------------------------      
      DO 44 K=0,NX
         DO 45 I=0,NT
            TRACE(I,K)=0.
 45      CONTINUE
 44   CONTINUE
C-----Get first data set into TRASE
      DO 21 K=0,NX
         DO 22 J=0,NT
            TRACE(J,K)=T(J,K)
 22      CONTINUE
 21   CONTINUE
      RMAX=1.
      NSAMP=NT
C-----Now TRACE contains the pressure data for R1
      PRINT *,'Output press.-data file        : ./pres1.d'
      OPEN(UNIT=15,FILE='./pres1.d',FORM='UNFORMATTED')
      DO 103 NTR=0,NX
         WRITE(15) (TRACE(I,NTR),I=1,NT)
 103  CONTINUE

      DO 544 K=0,NX
         DO 545 I=0,NT
            TRACE(I,K)=0.
 545     CONTINUE
 544  CONTINUE
C-----Get second data set into TRACE
      DO 521 K=0,NX
         DO 522 J=0,NT
            TRACE(J,K)=T1(J,K)
 522     CONTINUE
 521  CONTINUE
      NSAMP=NT
C-----Now TRACE contains the pressure data for R2
      PRINT *,'Output press.-data file        : ./pres2.d'
      OPEN(UNIT=16,FILE='./pres2.d',FORM='UNFORMATTED')
      DO 104 NTR=0,NX
         WRITE(16) (TRACE(I,NTR),I=1,NT)
 104  CONTINUE
C-----Terminate main program
      STOP
      END

C********************************************************************
C************** SUBROUTINES   ***************************************
C********************************************************************
      subroutine SOURCE
C********************************************************************
C
C     This program generates the source function - first derivative
C     of a gaussian.
C
C********************************************************************
C
C     Calculate the parameters needed to implement the first
C     derivative of a gaussian as a aource (cf. Alford et. al.
C     Geophysics, 1974).
C
C     SWIDTH = Negative peak to positive peak width of the pulse
C
C     SHIFT  = The number of SWITHS the zero crossover lies to
C	       the right of the time origin
C
C     SDELAY = The corresponding time delay
C
C     SLOPE  = Forces truncated gaussian derivative to zero
C              at times 0. and 2.*DELAY
C
C     SNORM  = Normalizes the source peaks to unit amplitude
C
C     RLOAD  = Is the desired amplitude of the source
C
C********************************************************************

C-----Local variable declaration
      implicit integer(I-N)
      dimension FF(200)
      common FF

      SWIDTH=.02
      RLOAD=1.

C-----Initialize the source function

      SNORM=1.
      SHIFT=2.
      SDELAY=2.*SWIDTH
      S2=SDELAY+SDELAY
      SLOPE=-EXP(-2.*4.)
      SNORM=SWIDTH*(.5*EXP(-.5)+SLOPE) 
      PERIOD=2.*SWIDTH                          
      CFREQ=1./PERIOD                           
      
      VMIN=6000.
      VMAX=9000.
     
      SQ2=SQRT(2.)
      DX=(VMIN*PERIOD)/12.
      DT=(VMIN*PERIOD)/(12.*VMAX*SQ2)

C-----WEIGHT(K) = The relative weight of the k'th load.
      WEIGHT=1.
      WEIGHT=WEIGHT*RLOAD/SNORM                                           

C-----Determine the maximum time step NFMAX for a nonzero source
      NFMAX=(2.*SDELAY)/DT+1.5                           
      DO 699 JJ=1,200 
         FF(JJ)=0.0
 699  CONTINUE
      T=0.0
      DO 1000 ITT=1,NFMAX
         T=T+DT
         P0=T-SDELAY
         P1=P0                                                        
         IF(P1+SDELAY .LT. 0.) GO TO 1000
         IF(P1-SDELAY .GT. 0.) GO TO 1000
         P2=P1/SWIDTH
         P3=2.*P2*P2
         FF(ITT)=P1*(EXP(-P3)+SLOPE)*WEIGHT         
 1000 CONTINUE
      PRINT *,'Output time source file : ./gau.d'
      OPEN(UNIT=18,FILE='./gau.d',FORM='UNFORMATTED')
      WRITE(18) (FF(ITT),ITT=1,200)
      CLOSE(18)
      RETURN
      END











































































