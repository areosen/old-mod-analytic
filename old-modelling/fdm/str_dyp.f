      program streamer_depth
C**********************************************************************
C     Program som beregner streamer-dyp som funksjon av x.          ***
C     Tar utgangspunkt i artikkelen "Signature Determination by     ***
C     Inversion" av Landroe og Sollie (Geophys. 12, 1992)           ***
C                                                                   ***
C**********************************************************************
C     Bruker ligning:                                               ***
C                                                                   *** 
C           z = z0 + sum_{n=1}^{K} [ a_n cos(2 pi n x/L)            ***
C                                   + b_n sin(2 pi n x/L) ]         ***
C                                                                   ***
C     Her er L lengde av streamer, K er ant. Fourier koeffisienter, ***
C     z0 er gj.snittlig streamer-dyp.                               ***
C**********************************************************************
      integer   nx,K
      real      dx
      parameter (nx=256,dx=2.5,K=6)

      real      pi
      real      L
      real      a(K),b(K)
      real      x(nx),z(nx),z0

      a(1) = .01
      b(1) = .3
      a(2) = .2
      b(2) = .01
      a(3) = .3
      b(3) = .3
      a(4) = .05
      b(4) = .05
      a(5) = .03
      b(5) = .01
      a(6) = .09
      b(6) = .09

      pi = 2.*acos(0.)
      z0 = 15.
      L  = nx*dx

      do i=1,nx
         z(i) = z0
         x(i) = (i-1)*dx
         do j=1,K
            z(i) = z(i) + 0.5*(a(j)*cos(2.*pi*j*x(i)/L)
     +                  + b(j)*sin(2.*pi*j*x(i)/L))
         enddo
      enddo

      print *, Size(z)
      print *, MinVal(z), MinLoc(z)
      print *, MaxVal(z), MaxLoc(z)

      PRINT *,'Print to file in ascii-format'
      PRINT *,'Output streamer depth data file: ./test.d'
      OPEN(UNIT=1,FILE='./test.d',FORM='FORMATTED')
      do i=1,nx
         WRITE(1,*) i*dx,z(i)
      enddo
      CLOSE(1)
      end


