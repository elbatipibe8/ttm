      PROGRAM cg_velocity
        IMPLICIT None
        INTEGER, PARAMETER :: LPMX  = 125000000    !Maximum number of particles
        INTEGER, PARAMETER :: nproc = 4
        REAL*8, PARAMETER  :: QTEM  = 294.0d0   !Distributed temperature
        REAL*8, PARAMETER :: ENUNIT = 1.036419026d-04 
        REAL*8, PARAMETER :: BK     = 8.617385d-05

        REAL*8, DIMENSION(:,:), ALLOCATABLE :: RD !Coordinates, velocities
        INTEGER, DIMENSION(:), ALLOCATABLE :: KHIST !History of each atom
        INTEGER, DIMENSION(:), ALLOCATABLE :: KTYPE !Type of each atom

        REAL*8 CUBE(3)                   !Size of the computational cell
        REAL*8 CENTR(3)                  !Position of the center of C.C.
        REAL*8 XMASS(3)

        REAL*8  time,RWW,RWW1,AWW,SK,P,VM,PR,T1,QIN,ZCRIT
        INTEGER NAN,mypid,I,J,NBOUND,NN
        CHARACTER*5 chpid

!        DATA XMASS/15.9994d0,47.867d0,87.62d0/
!        DATA XMASS/58.6934d0, 58.6934d0, 58.6934d0/
         DATA XMASS/26.981539d0,26.981539d0,26.981539d0/     ! Al
!         DATA XMASS/196.966569d0, 196.966569d0, 196.966569d0/ ! Au

        DO mypid = 0,nproc - 1

           ALLOCATE(RD(3,LPMX),KHIST(LPMX),KTYPE(LPMX))

           WRITE(chpid, FMT='(I5.5)') mypid

           OPEN(UNIT = 16, FILE = 'Al.vel.'//trim(chpid))
           OPEN(UNIT = 17, FILE = 'Al.'//trim(chpid))

!          Define the system boundary by KHIST

           READ(16,*) NAN,TIME,ZCRIT
           READ(16,*) CUBE(1:3)
           READ(16,*) CENTR(1:3)
           READ(16,*) (KTYPE(J),RD(1:3,J),J=1,NAN)
           READ(16,*) (KHIST(J),SK,SK,SK,J=1,NAN)

           WRITE(17,*) NAN,TIME,ZCRIT
           WRITE(17,*) CUBE(1:3)
           WRITE(17,*) CENTR(1:3)
           WRITE(17,*) (KTYPE(J),RD(1:3,J),J=1,NAN) 

           RD(1:3,1:NAN) = 0.0d0

           DO NN=1,3

              SK=0.0d0
              QIN=0.0d0
              NBOUND = 0

              DO J=1,NAN
                 IF(KHIST(J).NE.3) THEN
                    call random_number(rww)
                    RD(NN,J)=2.0d0*SQRT(QTEM*BK*RWW/ENUNIT/XMASS(1))
                    call random_number(rww1)
                    AWW=RWW1-0.5d0
                    RD(NN,J)=SIGN(RD(NN,J),AWW)
                    SK=SK+RWW
                 ELSE
                    NBOUND = NBOUND + 1
                 ENDIF
              ENDDO
                 
!             FINITE SIZE CORRECTION
              P=(NAN-NBOUND)/2.0d0/SK
              VM=0.0d0
              DO J=1,NAN
                 IF(KHIST(J).NE.3) THEN
                    RD(NN,J)=RD(NN,J)*P
                    VM=VM+RD(NN,J)*XMASS(1)
                 ENDIF
              ENDDO
                 
!             TOTAL MOMENTUM CORRECTION
              VM=VM/(NAN-NBOUND)
              DO J=1,NAN
                 IF(KHIST(J).NE.3) THEN
                    RD(NN,J)=RD(NN,J)-VM/XMASS(1)
                    QIN = QIN + 0.5d0*RD(NN,J)*RD(NN,J)*XMASS(1)
                 ENDIF
              ENDDO

!             TEMPERATURE CORRECTION
              T1=QIN*ENUNIT/BK/(NAN-NBOUND)
              PR=SQRT(QTEM/T1)
              QIN=0.0d0
              DO J=1,NAN
                 IF(KHIST(J).NE.3) THEN
                    RD(NN,J)=RD(NN,J)*PR
                    QIN = QIN + 0.5d0*RD(NN,J)*RD(NN,J)*XMASS(1)
                 ENDIF
              ENDDO
              T1=QIN*ENUNIT/BK/(NAN-NBOUND)
              print *, 'VEL--> T=',QTEM,',Corrected Tvel=' &
                   ,T1,'K in ',chpid,' unit'  

           ENDDO

           print *,chpid,"unit is completed"
           print *," " 

           WRITE(17,*) (KHIST(J),RD(1:3,J),J=1,NAN)

           DEALLOCATE(RD,KHIST,KTYPE)

           CLOSE(UNIT = 16)
           CLOSE(UNIT = 17)

        ENDDO

        print *,"Velocities were prescribed"

        STOP
      END PROGRAM CG_VELOCITY
