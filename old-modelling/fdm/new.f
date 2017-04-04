       PROGRAM TODIM
       
       INTEGER PDIMX,PDIMZ,PDIMT
       PARAMETER(PDIMX=512,PDIMZ=512,PDIMT=600)

       INTEGER I,J
       REAL    P(PDIMX,PDIMZ,PDIMT)

       do i=1,512
          do j=1,512
             print *,p(i,j,1)
          enddo
       enddo

       END
