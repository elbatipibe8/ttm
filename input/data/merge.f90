      Program Merge
        REAL*8 , DIMENSION (:,:), ALLOCATABLE :: XD,Q1D
        REAL*8 , DIMENSION (:), ALLOCATABLE :: POT,WALA,CSP
        INTEGER, DIMENSION (:), ALLOCATABLE :: KTYPE
        INTEGER I,J,JJ,NN1,NO,NF,NP,ITIME,IT
        REAL*8 XCENTR,YCENTR,ZCENTR,DELTAL,DELTAR
        CHARACTER*5 chpid,chtime

        ZCENTR = 200.0d0 ! [nm]
        YCENTR = 0.0d0
        XCENTR = 0.0d0
        DELTAL  =  0.00d0   ! Left border in Angstrons
        DELTAR  =  0.5d0   ! Right border in Angstrons

        NI = 1 !Initial time
        NF = 12 !Final time

        NP = 4 !Number of processors

        Loop: DO ITIME = NI,NF

!           IF(MOD(ITIME,5).NE.0) CYCLE Loop

!           IF(MOD(ITIME,5).NE.0) CYCLE Loop

           WRITE(chtime, FMT='(I3.3)') ITIME

!           OPEN(UNIT=1,FILE='./time'//trim(chtime)//'.d')
           OPEN(UNIT=2,FILE='./slice'//trim(chtime)//'.d')
           
           JJ = 0

           DO I=0,NP-1
              
              WRITE(chpid, FMT='(I5.5)') I
              
              OPEN(UNIT=969,FILE='./time'//trim(chtime)//'.'//trim(chpid))
              
              READ(969,*) NN1
              
              ALLOCATE(XD(3,NN1),Q1D(3,NN1),POT(NN1), &
                   WALA(NN1),CSP(NN1),KTYPE(NN1))

              READ(969,*) (XD(1,J),XD(2,J),XD(3,J), &
                   Q1D(1,J),Q1D(2,J),Q1D(3,J), &
                   POT(J),WALA(J),CSP(J),KTYPE(J),J=1,NN1)
              
!              WRITE(1,166)(XD(1,J),XD(2,J),XD(3,J),POT(J),WALA(J),CSP(J),J=1,NN1)!KHIST(J),J=1,NN1)
              
              CLOSE(UNIT = 969)
              
              DO J=1,NN1
                 IF(ABS(XD(1,J)).LE.0.5d0) THEN !)POT(J).GE.-3.2d0) THEN
                    IF(CSP(J).EQ.0.0d0) CSP(J) = 1.0d0
                    WRITE(2,166) XD(1,J),XD(2,J),XD(3,J), &
!                         Q1D(1,J),Q1D(2,J),Q1D(3,J), &
                         POT(J),WALA(J),CSP(J),KTYPE(J)
                    JJ = JJ + 1
                 ENDIF
              ENDDO
              
166           FORMAT(F9.3,1x,F9.3,1x,F9.3,1x, &
!                   F9.3,1x,F9.3,1x,F9.3,1x, &
                   F14.3,1x,F6.3,1x,F6.3,1x,I3)
             
              DEALLOCATE(XD,Q1D,POT,WALA,CSP,KTYPE)
              
              print *,I,"processing.....",ITIME,"ps", JJ, "atoms"

           ENDDO


!           CLOSE(UNIT = 1)
           CLOSE(UNIT = 2)
           
           print *,"File time", chtime,"is merged and sliced"
           print *,JJ,"atoms in the slice"
           print *," "
           
        ENDDO Loop

      END Program Merge
