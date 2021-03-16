      Program Merge
        REAL*8 , DIMENSION (:,:), ALLOCATABLE :: XD
        REAL*8 , DIMENSION (:), ALLOCATABLE :: TlG,TeG,TP,RON,WALA,Q1
        INTEGER I,J,JJ,NN1,NO,NF,NP,ITIME,IT
        REAL*8 XCENTR,YCENTR,ZCENTR,DELTAL,DELTAR,YL
        CHARACTER*5 chpid,chtime

        ZCENTR = 200.0d0 ! [nm]
        YCENTR = 0.0d0
        XCENTR = 0.0d0
        YL     =  349.651775200000d0
        DELTAL  =  0.00d0   ! Left border in Angstrons
        DELTAR  =  0.5d0   ! Right border in Angstrons

        NI = 0 !Initial time
        NF = 200  !Final time

        NN1 = 9121536 ! a product of all the link cells

        NP = 256 !Number of processors

        Loop: DO ITIME = NI,NF

           IF(MOD(ITIME,25).NE.0) CYCLE Loop

           WRITE(chtime, FMT='(I3.3)') ITIME

!           OPEN(UNIT=1,FILE='./time'//trim(chtime)//'.d')
           OPEN(UNIT=2,FILE='./slice_ttm'//trim(chtime)//'.d   ')
              
           OPEN(UNIT=969,FILE='./ttm'//trim(chtime)//'.d')
              
           ALLOCATE(XD(3,NN1),TlG(NN1),TeG(NN1),TP(NN1),RON(NN1),WALA(NN1),Q1(NN1))

           READ(969,*) (XD(1,JJ),XD(2,JJ),XD(3,JJ), &
                TlG(JJ),TeG(JJ),TP(JJ),RON(JJ),WALA(JJ),Q1(JJ),JJ=1,NN1)
              
           print *,"reading",ITIME,"ps is over"
              
           CLOSE(UNIT = 969)
              
           DO J=1,NN1
              IF(ABS(XD(1,J)).LE.1.0d0.AND.XD(3,J).GE.-200.0d0) &
                    WRITE(2,166) XD(1,J),XD(2,J),XD(3,J), &
                         TlG(J),TeG(J),TP(J),RON(J),WALA(J),Q1(J)
!                 IF(MOD(J,1000000).EQ.0) print *,REAL(J)/REAL(NN1)*100," % done"
           ENDDO
              
166        FORMAT(F8.2,1x,F8.2,1x,F8.2,1x, &
                F6.0,1x,F7.0,1x,F7.3,1x, &
                F5.3,1x,F5.3,1x,F8.2)
              
           DEALLOCATE(XD,TlG,TeG,TP,RON,WALA,Q1)
              
           CLOSE(UNIT = 2)
           
           print *,"File ttm", chtime,"is sliced"
           print *," "
           
        ENDDO Loop

      END Program Merge
