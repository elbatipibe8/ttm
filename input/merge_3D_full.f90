      Program Merge
        REAL*8 , DIMENSION (:,:), ALLOCATABLE :: XD,Q1D
        REAL*8 , DIMENSION (:), ALLOCATABLE :: POT,WALA,CSP
        INTEGER, DIMENSION (:), ALLOCATABLE :: KTYPE
        INTEGER I,J,JJ,NN1,NO,NF,NP,ITIME,IT
        REAL*8 XCENTR,YCENTR,ZCENTR,DELTAL,DELTAR,SHIFT
        CHARACTER*5 chpid,chtime

        ZCENTR =  2000.5795287437300d0 ! [nm]
        YCENTR = 0.0d0
        XCENTR = 0.0d0
        DELTAL  =  19.6066416d0   ! Left border in Angstrons
        DELTAR  = 349.65177520000002   ! Right border in Angstrons
        SHIFT = 5.0d0


        NI =  0 !Initial time
        NF =  5 !Final time

        NP = 400 !Number of processors

        Loop: DO ITIME = NI,NF

           IF(MOD(ITIME,5).NE.0) CYCLE Loop

           WRITE(chtime, FMT='(I3.3)') ITIME

!           OPEN(UNIT=1,FILE='./time'//trim(chtime)//'.d')
           OPEN(UNIT=2,FILE='./slice'//trim(chtime)//'.d_3D')
           
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
              
              Loop2: DO J=1,NN1
!                 IF(XD(1,J).GE.5.0d0) XD(1,J) = XD(1,J) - 2.0d0*DELTAL
!                 IF(ABS(XD(1,J)).GE.1.40d0.OR.POT(J).GE.-3.2d0) THEN !)POT(J).GE.-3.2d0) THEN
!                 IF((XD(3,J).LE.-50.0d0).AND.(ABS(XD(2,J)).LE.60.0d0).AND.(ABS(XD(1,J)).LE.60.0d0)) Cycle Loop2
!                 IF(XD(3,J).GE.175.0d0.AND.XD(2,J).LE.130.0d0.AND.XD(1,J).LE.3.0d0) Cycle Loop2

!                 IF(XD(3,J).LE.-2.0d0.AND.ABS(XD(2,J)).GE.130.0d0 &
!                      .AND.XD(2,J).LE.248.0d0.AND.XD(1,J).LE.3.0d0) Cycle Loop2


!                    IF(XD(2,J).LE.0.0d0) XD(2,J) = XD(2,J) + DELTAR
!                    IF(XD(2,J)-0.5d0*DELTAR.GE.-100.0d0.AND.XD(2,J)-0.5d0*DELTAR.LE.0.d0 &
!                         .AND.XD(3,J).GE.-80.0d0)THEN
!                       IF(XD(1,J).LE.-DELTAL*0.20d0) XD(1,J) = XD(1,J)+DELTAL
!                       IF(ABS(XD(1,J)-0.3d0*DELTAL).LE.0.2d0) THEN
!                 IF(XD(3,J).GE.-125.0d0) THEN
!                 IF(XD(2,J).LE.0.0d0) THEN
                 IF(CSP(J).EQ.0.0d0) CSP(J) = 1.0d0
                    WRITE(2,166) XD(1,J),XD(2,J),XD(3,J), &
                         Q1D(1,J),Q1D(2,J),Q1D(3,J), &
                         POT(J),WALA(J),CSP(J),KTYPE(J)
                    JJ = JJ + 1
!                 ENDIF
!                    ENDIF
              ENDDO Loop2
              
166           FORMAT(F9.3,1x,F9.3,1x,F9.3,1x, &
                   F9.3,1x,F9.3,1x,F9.3,1x, &
                   F6.3,1x,F5.3,1x,F5.3,1x,I3)
             
              DEALLOCATE(XD,Q1D,POT,WALA,CSP,KTYPE)
              
              print *,I,"processing.....",ITIME,"ps", JJ, "atoms in"

           ENDDO


!           CLOSE(UNIT = 1)
           CLOSE(UNIT = 2)
           
           print *,"File time", chtime,"is merged and sliced"
           print *,JJ,"atoms in the slice"
           print *," "
           
        ENDDO Loop

      END Program Merge
