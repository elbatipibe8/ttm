      PROGRAM ZL
        IMPLICIT None
        INTEGER, PARAMETER :: LPMX = 125000000    !Maximum number of particles
        INTEGER, PARAMETER :: nproc = 84!72!804

        REAL*8, DIMENSION(:,:), ALLOCATABLE :: RD !Coordinates, velocities
        REAL*8, DIMENSION(:,:), ALLOCATABLE :: Q1D !Coordinates, velocities
        INTEGER, DIMENSION(:), ALLOCATABLE :: KHIST !History of each atom
        INTEGER, DIMENSION(:), ALLOCATABLE :: KTYPE !Type of each atom

        REAL*8 CUBE(3)                   !Size of the computational cell
        REAL*8 CENTR(3)                  !Position of the center of C.C.

        REAL*8  time,ZCRIT,ZD,DE,aa,bb,cc,dd,delta
        INTEGER NAN,mypid,I,J,KKK,K11,K12,K13
        CHARACTER*5 chpid

!        delta = 0.5d0*4.927595111d0
!        aa = -3.0d0*delta
!        bb = -2.0d0*delta + 0.0000d0
!        cc =    - delta - 0.00021d0
!        bb = cc - delta - 0.00057d0
!        aa = bb - delta - 0.00110d0

!        dd = cc + delta + 0.00740d0
        
        K11 = 0
        K12 = 0
        K13 = 0

        KKK = 0

        DO mypid = 0,nproc - 1

           ALLOCATE(RD(3,LPMX),Q1D(3,LPMX),KHIST(LPMX),KTYPE(LPMX))

           WRITE(chpid, FMT='(I5.5)') mypid

           OPEN(UNIT = 16, FILE = 'H2O.'//trim(chpid))
!           OPEN(UNIT = 17, FILE = 'H2O.hist.'//trim(chpid))
           OPEN(UNIT = 18, FILE = 'test.hist')

!          Define the system boundary by KHIST

           READ(16,*) NAN,TIME,ZCRIT
           READ(16,*) CUBE(1:3)
           READ(16,*) CENTR(1:3)
           READ(16,*) (KTYPE(J),RD(1:3,J),J=1,NAN)
           READ(16,*) (KHIST(J),Q1D(1:3,J),J=1,NAN)
!           READ(16,*) (KHIST(J),Q1D,Q1D,Q1D,J=1,NAN)

           CLOSE(UNIT = 16)
           CUBE(3) = 10000.0d0
           OPEN(UNIT = 17,FILE = 'H2O.'//trim(chpid))
           
           WRITE(17,*) NAN,TIME,ZCRIT
           WRITE(17,*) CUBE(1:3)
           WRITE(17,*) CENTR(1:3)
           WRITE(17,67) (KTYPE(J),RD(1:3,J),J=1,NAN)
!           WRITE(17,67  (KHIST(J),RD(1:3,J),J=1,NAN)

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
!              IF(RD(3,I).GT.-6.5d0) THEN
!                 KHIST(I) = 3
!                 GOTO 123
!              ENDIF

!              goto 123

              IF(RD(3,I).GT.aa.AND.RD(3,I).LE.bb) THEN
                 KHIST(I) = 11
                 K11 = K11 + 1
                 GOTO 123
              ENDIF
              IF(RD(3,I).GT.bb.AND.RD(3,I).LE.cc) THEN
                 KHIST(I) = 12
                 K12 = K12 + 1
                 GOTO 123
              ENDIF
              IF(RD(3,I).GT.cc) THEN
                 KHIST(I) = 13
                 K13 = K13 + 1
                 GOTO 123
              ENDIF
              GOTO 123
              
!             For the case of 150nm thickness
              IF(RD(3,I).LE.-1.0d0.AND.RD(3,I).GT.-2.0d0) THEN
                 KHIST(I) = 11
                 GOTO 123
              ENDIF
              IF(RD(3,I).LE.-3.0d0.AND.RD(3,I).GT.-4.0d0) THEN
                 KHIST(I) = 12
                 GOTO 123
              ENDIF
              IF(RD(3,I).LE.-

              IF(DABS(RD(1,J)).LE.1.5d0) &
                   WRITE(18,111) RD(1:3,J),KTYPE(J),KHIST(J)

!             IF(KHIST(I).NE.0) KKK = KKK + 1

!             IF(RD(3,I).LE.ZD.0d0) RD(3,I) = RD(3,I)-ZD
             

          ENDDO

          WRITE(17,68) (KHIST(J),Q1D(1:3,J),J=1,NAN)


          DEALLOCATE(RD,Q1D,KHIST,KTYPE)
          
          CLOSE(UNIT = 16)
!          CLOSE(UNIT = 17)
          
          print *,K11,K12,K13,"      Khist values were prescribed for ", mypid
          
       ENDDO

        print *,KKK,"Total number of KHIST"
        
111     FORMAT(3(1x,E14.7),2(1x,I2))
67      FORMAT(I3,1x,F20.12,1x,F20.12,1x,F20.12)
68      FORMAT(I3,1x,F18.12,1x,F18.12,1x,F18.12)
        

        CLOSE(UNIT = 18)

        STOP
      END PROGRAM ZL
