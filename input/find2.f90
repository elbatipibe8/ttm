      Program find
        implicit none  
        REAL*8 da,db,a,b,eps,epsx
        INTEGER I,J,I0,J0,N

        OPEN(UNIT = 1,FILE = "find.out")
        da = 4.929689d0 ! water
        db = 4.084717d0 ! Au
        epsx = 5.0d0
        N = 2000

        DO I = 1,N
           DO J = 1,N
              a = I*da
              b = J*db
              eps = a - b
              IF(ABS(eps).LE.epsx) THEN
                 I0 = I
                 J0 = J
                 epsx = ABS(eps)
              ENDIF
           ENDDO
        ENDDO

        print *, epsx,I0,J0
      END Program find
