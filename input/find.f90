      Program find
        implicit none  
        REAL*8 da,db,a,b,eps
        INTEGER I,K,M

        OPEN(UNIt = 1,FILE = "find.out")
        da = 4.929689d0 ! water

!        da = 4.8415d0
        
        da = 4.9277d0


        db = 4.084717d0 ! Au

        a = da
        b = db
        K = 0
        M = 0

        DO I = 1,2000
           a = K*da
           b = M*db
           eps = a - b

           IF(eps.GE.db) THEN
              b = b + db
              M = M + 1
              eps = a - b
           ENDIF
           IF(eps.LE.-db) THEN
              a = a + da
              K = K + 1
              eps = a - b
           ENDIF
           
           WRITE(1,*) I,eps,K,M

           K = K + 1
           M = M + 1
        ENDDO

        CLOSE(UNIT = 1)
      END Program find
