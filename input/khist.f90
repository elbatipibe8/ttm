      PROGRAM Khist_boundary
        IMPLICIT None
        INTEGER, PARAMETER :: LPMX = 125000000    !Maximum number of particles
        INTEGER, PARAMETER :: nproc = 4 

        REAL*8, DIMENSION(:,:), ALLOCATABLE :: RD !Coordinates, velocities
        INTEGER, DIMENSION(:), ALLOCATABLE :: KHIST !History of each atom
        INTEGER, DIMENSION(:), ALLOCATABLE :: KTYPE !Type of each atom

        REAL*8 CUBE(3)                   !Size of the computational cell
        REAL*8 CENTR(3)                  !Position of the center of C.C.

        REAL*8  time,Q1D,ZCRIT
        INTEGER NAN,mypid,I,J,KKK
        CHARACTER*5 chpid

        KKK = 0

        DO mypid = 0,nproc - 1

           ALLOCATE(RD(3,LPMX),KHIST(LPMX),KTYPE(LPMX))

           WRITE(chpid, FMT='(I5.5)') mypid

           OPEN(UNIT = 16, FILE = 'Al.khist.vel.'//trim(chpid))
           OPEN(UNIT = 17, FILE = 'Al.vel.'//trim(chpid))
           OPEN(UNIT = 18, FILE = 'test.vel.out')

!          Define the system boundary by KHIST

           READ(16,*) NAN,TIME,ZCRIT
           READ(16,*) CUBE(1:3)
           READ(16,*) CENTR(1:3)
           READ(16,*) (KTYPE(J),RD(1:3,J),J=1,NAN)
           READ(16,*) (KHIST(J),Q1D,Q1D,Q1D,J=1,NAN)
!           READ(16,*) (KHIST(J),Q1D,Q1D,Q1D,J=1,NAN)

           WRITE(17,*) NAN,TIME,ZCRIT
           WRITE(17,*) CUBE(1:3)
           WRITE(17,*) CENTR(1:3)
           WRITE(17,*) (KTYPE(J),RD(1:3,J),J=1,NAN)

           DO I=1,NAN


               GOTO 123
              
!            For the case of 100 nm thickness
!              IF(RD(3,I).LE.-999.0d0.AND.RD(3,I).GT.-1001.0d0) THEN
!                 KHIST(I) = 11
!                 GOTO 123
!              ENDIF
!              IF(RD(3,I).LE.-1001.0d0.AND.RD(3,I).GT.-1003.0d0) THEN
!                 KHIST(I) = 12
!                 GOTO 123
!              ENDIF
!              IF(RD(3,I).LE.-1003.0d0.AND.RD(3,I).GT.-1005.0d0) THEN
!                 KHIST(I) = 13
!                 GOTO 123
!              ENDIF


!             For the case of 150nm thickness
!              IF(RD(3,I).LE.-1497.0d0.AND.RD(3,I).GT.-1499.0d0) THEN
!                 KHIST(I) = 11
!                 GOTO 123
!              ENDIF
!              IF(RD(3,I).LE.-1499.0d0.AND.RD(3,I).GT.-1501.0d0) THEN
!                 KHIST(I) = 12
!                 GOTO 123
!              ENDIF
!              IF(RD(3,I).LE.-1501.0d0.AND.RD(3,I).GT.-1503.0d0) THEN
!                 KHIST(I) = 13
!                 GOTO 123
!              ENDIF


              GOTO 123

!             For the case of 200nm thickness
              IF(RD(3,I).LE.-524.0d0.AND.RD(3,I).GT.-525.0d0) THEN
                 KHIST(I) = 11
                 GOTO 123
              ENDIF
              IF(RD(3,I).LE.-525.0d0.AND.RD(3,I).GT.-526.0d0) THEN
                 KHIST(I) = 12
                 GOTO 123
              ENDIF
              IF(RD(3,I).LE.-526.0d0.AND.RD(3,I).GT.-528.0d0) THEN
                 KHIST(I) = 13
                 GOTO 123
              ENDIF


!              IF(RD(2,I).LE.-29994.0d0) KHIST(I) = 0
!              IF(RD(2,I).GE.22050.0d0.AND.RD(2,I).LT.22052.0d0) THEN
!                 KHIST(I) = 14
!                 GOTO 123
!              ENDIF
!              IF(RD(2,I).GE.22052.0d0.AND.RD(2,I).LT.22054.0d0) THEN
!                 KHIST(I) = 15
!                 GOTO 123
!              ENDIF
!              IF(RD(2,I).GE.22054.0d0.AND.RD(2,I).LT.22056.0d0) THEN
!                 KHIST(I) = 16
!                 GOTO 123
!              ENDIF
!              IF(RD(2,I).GE.29994.0d0) KHIST(I) = 0

123           CONTINUE

             IF(KHIST(I).NE.0) KKK = KKK + 1


          ENDDO

          WRITE(17,*) (KHIST(J),Q1D,Q1D,Q1D,J=1,NAN)

          DO J=1,NAN
             IF(DABS(RD(2,J)).LE.1.5d0) &
                  WRITE(18,111) RD(1:3,J),KTYPE(J),KHIST(J)
          ENDDO

          DEALLOCATE(RD,KHIST,KTYPE)
          
          CLOSE(UNIT = 16)
          CLOSE(UNIT = 17)
          
          print *,KKK,"Khist values were prescribed for ",mypid
          
       ENDDO

        print *,KKK,"Total number of KHIST"
        
111     FORMAT(3(1x,E14.7),2(1x,I2))
        

        CLOSE(UNIT = 18)

        STOP
      END PROGRAM KHIST_BOUNDARY
