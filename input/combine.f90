      Program Combine
        REAL*8 , DIMENSION (:,:), ALLOCATABLE :: XD,Q1D,XD2,Q1D2
        INTEGER, DIMENSION (:), ALLOCATABLE :: KTYPE,KHIST,KTYPE2,KHIST2
        INTEGER I,J,K,NP,NN1,NN2
        REAL*8 XCENTR,YCENTR,ZCENTR,XL,YL,ZL
        REAL*8 XCENTR2,YCENTR2,ZCENTR2,XL2,YL2,ZL2
        CHARACTER*5 chpid

        NP = 256 !Number of processors

        OPEN(UNIT = 980,FILE = 'visio.out')
        
        DO I=0,NP-1
              
           WRITE(chpid, FMT='(I5.5)') I
              
           OPEN(UNIT=969,FILE='./Au.'//trim(chpid))
           OPEN(UNIT=970,FILE='./H2O.'//trim(chpid))
                            
           READ(969,*) NN1,TIME,ZCRIT
           READ(969,*) XL,YL,ZL
           READ(969,*) XCENTR,YCENTR,ZCENTR

           READ(970,*) NN2,TIME2,ZCRIT2
           READ(970,*) XL2,YL2,ZL2
           READ(970,*) XCENTR2,YCENTR2,ZCENTR2
              
           ALLOCATE(XD(3,NN1),Q1D(3,NN1),KTYPE(NN1),KHIST(NN1))
           ALLOCATE(XD2(3,NN2),Q1D2(3,NN2),KTYPE2(NN2),KHIST2(NN2))
              
           READ(969,*) (KTYPE(J),XD(1,J),XD(2,J),XD(3,J),J=1,NN1)
           READ(969,*) (KHIST(J),Q1D(1,J),Q1D(2,J),Q1D(3,J),J=1,NN1)

           READ(970,*) (KTYPE2(J),XD2(1,J),XD2(2,J),XD2(3,J),J=1,NN2)
           READ(970,*) (KHIST2(J),Q1D2(1,J),Q1D2(2,J),Q1D2(3,J),J=1,NN2)
              
           CLOSE(UNIT = 969)
           CLOSE(UNIT = 970)

           OPEN(UNIT=971,FILE='./AuWT.'//trim(chpid))
              
           WRITE(971,*) NN1+NN2,0.0d0,ZCRIT
           WRITE(971,*) XL,YL,10000.0d0
           WRITE(971,*) XCENTR,YCENTR,3000.0d0
              
           WRITE(971,67) (KTYPE(J),XD(1,J),XD(2,J),XD(3,J),J=1,NN1)
           WRITE(971,67) (2,XD2(1,J),XD2(2,J),XD2(3,J) - ZCENTR2 + 0.5d0*ZL2 + 1.3032d0,J=1,NN2)
!           print *,- ZCENTR2 + 0.5d0*ZL2
           
           WRITE(971,68) (KHIST(J),Q1D(1,J),Q1D(2,J),Q1D(3,J),J=1,NN1)
           WRITE(971,68) (KHIST2(J),Q1D2(1,J),Q1D2(2,J),Q1D2(3,J),J=1,NN2)

           DO K = 1,NN1
              IF(ABS(XD(1,K)).LE.5) THEN
                 WRITE(980,69) XD(1,K),XD(2,K),XD(3,K),KTYPE(K)
              ENDIF
           ENDDO
           DO K = 1,NN2
              IF(ABS(XD2(1,K)).LE.5) THEN
                 WRITE(980,69) XD2(1,K),XD2(2,K),XD2(3,K) - ZCENTR2 + 0.5d0*ZL2 + 1.3032d0,2
              ENDIF
           ENDDO
           
67         FORMAT(I3,1x,F20.12,1x,F20.12,1x,F20.12)
68         FORMAT(I3,1x,F18.12,1x,F18.12,1x,F18.12)
69         FORMAT(F20.12,1x,F20.12,1x,F20.12,1x,I2)

           CLOSE(UNIT = 971)

           DEALLOCATE(XD,Q1D,KHIST,KTYPE)
           DEALLOCATE(XD2,Q1D2,KHIST2,KTYPE2)

           print *, "file AuWT ",chpid,"  prepared"
              
        ENDDO

        CLOSE(UNIT = 980)

      END Program Combine
