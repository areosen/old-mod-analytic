      program gausisk_kildesignatur

      integer nt,nw
      parameter (nt=1024,nw=nt/2+1)

      integer iw,it
      integer nsign
      real    warr(nw)
      real    sort(nt)
      complex ca(nw),ca2(nt)
      real    dw,pi,t0,eps,dt
      complex j
      complex w
      
C-----Setter parametre
      DT   = 0.001
      J    = (0.,1.)
      PI   = 2.*ACOS(0.)
      WNYQ = PI/DT
      DW   = 2.*PI/(NT*DT)
      EPS  = 10000.
      T0   = 0.06

C-----Oppretter vektor med alle frekvenser
      DO 100 IW=1,NW
         WARR(IW) = (IW-1)*DW
 100  CONTINUE

C-----Beregner kilden i frekvensdomenet
      DO 200 IW=2,NW
         W      = WARR(IW)
         CA(IW) = CEXP(-W*W/(4.*EPS))*CEXP(J*W*T0)
         CA(IW) = J*W*CA(IW)*SQRT(PI)/EPS
 200  CONTINUE

C-----Tar en invers fouriertransform
      DO 350 IT=1,NW
         CA2(IT)=CA(IT)
 350  CONTINUE
      DO 400 IW=NT-NW+2,NT
         CA2(IW) = CONJG(CA(NT-IW+2))
 400  CONTINUE
      NSIGN = -1
      CALL FOUR(CA2,NT,NT,DT,NSIGN)
      DO 500 IT=1,NT
         SORT(IT) = REAL(CA2(IT))
 500  CONTINUE
      
      DO 600 IT=1,NT
         PRINT *,SORT(IT)
 600  CONTINUE
      END




