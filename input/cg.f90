      PROGRAM Au_Crystall
        IMPLICIT None

        REAL*8, DIMENSION(:,:), ALLOCATABLE :: RD !Coordinates, velocities
        INTEGER, DIMENSION(:), ALLOCATABLE :: KHIST !History of each atom
        INTEGER, DIMENSION(:), ALLOCATABLE :: KTYPE !Type of each atom
!        INTEGER, PARAMETER :: LPMX = 200000000
        
        REAL*8 CUBE(3)                   !Size of the computational cell
        REAL*8 CENTR(3)                  !Position of the center of C.C.
        REAL*8 dnode(3),rr(3)
        REAL*8 xlstart,ylstart,zlstart,ZCRIT
        REAL*8 xnmin,xnmax,ynmin,ynmax,znmin,znmax
 
        INTEGER, PARAMETER :: Ncell = 4  !Number ov atoms per unit cell

        REAL*8  BASIS(1:3,Ncell)        !Positions of atoms in unit cell
        REAL*8  SHIFT(3),TRANS(3),XINI(3)
        REAL*8  bx(Ncell),by(Ncell),bz(Ncell)
        REAL*8  alat,time,factor
        INTEGER NAN,NANG,iatom,i,j,k,jj,NKH
        INTEGER mynodex,mynodey,mynodez,mypid

        CHARACTER*5 chpid

        INTEGER LPMX
        INTEGER nn(3)

        INTEGER, PARAMETER :: nxnodes = 2
        INTEGER, PARAMETER :: nynodes = 2
        INTEGER, PARAMETER :: nznodes = 1
        INTEGER, PARAMETER :: nproc   = nxnodes*nynodes*nznodes

        DATA nn/20,20,493/
!        DATA nn/50,50,150/

        
        DATA bx/0, 0.5, 0.0, 0.5/
        DATA by/0, 0.5, 0.5, 0.0/
        DATA bz/0, 0.0, 0.5, 0.5/

        DATA SHIFT/0.25, 0.25, 0.25/

!        alat=3.53151d0   ! Ni at 300 K T_distr = 303
!        alat = 4.065d0   ! Au at T = 0
!        alat = 3.513d0 ! Ni at 0 K, ) GPa
!        alat = 3.523717d0
!        alat = 4.1606d0  ! Au at 1223K P=0, T_given = 1000 K
         alat = 4.056433 !  (Al)
       
        DO i=1,Ncell
           BASIS(1,i) = bx(i)
           BASIS(2,i) = by(i)
           BASIS(3,i) = bz(i)
        ENDDO

        print *,'building crystall with', nn(1:3)
        print *,'lattice constant (Ang):', alat
        print *,' ' 
        
!       Define the sizes and positions of the computational cell
        cube(1:3) = nn(1:3) * alat

        factor      = 10.0d0
        ZCRIT       = -0.25d0*alat + alat/SQRT(2.0d0)
        !- cube(3) + 0.25d0*alat - alat/SQRT(2.0d0) + 0.052863d0
        centr(1:2) = 0.0d0
!        centr(3)   = factor*cube(3)*0.5d0-cube(3)
!        cube(3)     = factor*cube(3)
        centr(3) = -0.5d0*cube(3)

        LPMX = INT(1.5d0*nn(1)*nn(2)*nn(3)*Ncell/nproc) !Number atoms/proc

        print *,"Averaged number of atoms per processors LPMX=",nn(1)*nn(2)*nn(3)*Ncell/nproc
        print *,"Correct the array boundaries for khist.f90 and vel.f90"
        print *," "  

        xlstart = centr(1) - 0.5d0*cube(1)
        ylstart = centr(2) - 0.5d0*cube(2)
        zlstart = centr(3) - 0.5d0*cube(3)

        XINI(1) = xlstart
        XINI(2) = ylstart
        XINI(3) = zlstart

        print *,XINI(1:3),"beginning of the crystall"

        dnode(1) = cube(1)/DBLE(nxnodes)
        dnode(2) = cube(2)/DBLE(nynodes)
        dnode(3) = cube(3)/DBLE(nznodes)

        WRITE(* ,910) cube(1:3)
        print *,"Real size is:",nn(1:3)*alat,"Angstons"
910     FORMAT (' The cube sizes in Ang are  ',3f12.6,/)
        
        NANG = 0
        NKH = 0

        DO mypid = 0,nproc - 1

           ALLOCATE(RD(3,LPMX),KHIST(LPMX),KTYPE(LPMX))

           WRITE(chpid, FMT='(I5.5)') mypid

           OPEN(UNIT = 16, FILE = 'Al.vel.'//trim(chpid))
           OPEN(UNIT = 17, FILE = 'test.vel.out')

!          Define the system geometry
           mynodex = MOD(mypid,nxnodes)
           mynodey = MOD((mypid/nxnodes),nynodes)
           mynodez = mypid/(nxnodes*nynodes)

           xnmin = xlstart + DBLE(mynodex)*dnode(1)
           xnmax = xlstart + DBLE(mynodex+1)*dnode(1)

           ynmin = ylstart + DBLE(mynodey)*dnode(2)
           ynmax = ylstart + DBLE(mynodey+1)*dnode(2)

           znmin = zlstart + DBLE(mynodez)*dnode(3)
           znmax = zlstart + DBLE(mynodez+1)*dnode(3)

!          Generating atoms per processor

           NAN = 0
           z_loop:   DO k=0,nn(3)-1

             IF((k+1)*alat+XINI(3).LT.znmin &
                  .OR.(k-1)*alat+XINI(3).GE.znmax) CYCLE z_loop
              
              y_loop:  DO j=0,nn(2)-1

                IF((j+1)*alat+XINI(2).LT.ynmin &
                     .OR.(j-1)*alat+XINI(2).GE.ynmax) CYCLE y_loop

                 x_loop: DO i=0,nn(1)-1

                   IF((i+1)*alat+XINI(1).LT.xnmin &
                        .OR.(i-1)*alat+XINI(1).GE.xnmax) CYCLE x_loop
                      
                       TRANS(1) = i*alat
                       TRANS(2) = j*alat
                       TRANS(3) = k*alat

                       DO iatom = 1, Ncell
                       
                       rr(1:3) = XINI(1:3) + &
                            (BASIS(1:3,iatom) + SHIFT(1:3))*alat + TRANS(1:3)

                       IF((rr(1).GE.xnmin.AND.rr(1).LT.xnmax) &
                            .AND.(rr(2).GE.ynmin.AND.rr(2).LT.ynmax) &
                            .AND.(rr(3).GE.znmin.AND.rr(3).LT.znmax)) THEN
                          
!                          IF(rr(3).GE.-506.0d0) THEN
                          
                          NAN = NAN + 1

                          rd(1:3,NAN) = rr(1:3)

!                          IF(iatom.EQ.1) THEN
                             KTYPE(NAN) = 1
!                          ELSEIF(iatom.eq.4) THEN
!                             KTYPE(NAN) = 1
!                          ELSE
!                             KTYPE(NAN) = 1
!                          ENDIF
                          
                          KHIST(NAN) = 0

!                          IF(rr(3).LE.-1523.0d0.AND.rr(3).GT.-1525.0d0) THEN
!                             KHIST(NAN) = 11
!                             NKH = NKH + 1
!                             GOTO 123
!                          ENDIF
!                          IF(rr(3).LE.-1525.0d0.AND.rr(3).GT.-1527.0d0) THEN
!                             KHIST(NAN) = 12
!                             NKH = NKH + 1
!                             GOTO 123
!                          ENDIF
!                          IF(rr(3).LE.-1527.0d0.AND.rr(3).GT.-1529.0d0) THEN
!                             KHIST(NAN) = 13
!                             NKH = NKH + 1
!                             GOTO 123
!                          ENDIF
!              IF(RD(2,I).LE.-29994.0d0) KHIST(I) = 0
!                          IF(rr(2).GE.22050.0d0.AND.rr(2).LT.22052.0d0) THEN
!                             KHIST(NAN) = 14
!                             NKH = NKH + 1
!                             GOTO 123
!                          ENDIF
!                          IF(rr(2).GE.22052.0d0.AND.rr(2).LT.22054.0d0) THEN
!                             KHIST(NAN) = 15
!                             NKH = NKH + 1
!                             GOTO 123
!                          ENDIF
!                          IF(rr(2).GE.22054.0d0.AND.rr(2).LT.22056.0d0) THEN
!                             KHIST(NAN) = 16
!                             NKH = NKH + 1
!                             GOTO 123
!                          ENDIF
!              IF(RD(2,I).GE.29994.0d0) KHIST(I) = 0

!                       ENDIF
                          
123                       CONTINUE

                       ENDIF

                    enddo
      
                 END DO x_loop
              END DO y_loop
           END DO z_loop
           
           TIME=0.0d0
!           ZCRIT = 0.0d0
!           CUBE(3) = 10000.0d0
!           CENTR(3)= 2000.0d0

           print *,mypid,"unit contains ",NAN,"particles",xnmin,znmin

!          Write output coordinate file
           REWIND 16
           WRITE(16,*) NAN,TIME,ZCRIT
           WRITE(16,*) CUBE(1:3)
           WRITE(16,*) CENTR(1),CENTR(2),CENTR(3)!-0.75d0*alat
           WRITE(16,*) (KTYPE(J),RD(1,J),RD(2,J),RD(3,J),J=1,NAN)
           WRITE(16,*) (KHIST(J),0.0d0,0.0d0,0.0d0,J=1,NAN)
           CLOSE(UNIT = 16)

           DO JJ=1,NAN
              IF(ABS(RD(2,JJ)).LE.1.5d0) &
                   WRITE(17,111) RD(1,JJ),RD(2,JJ),RD(3,JJ),KTYPE(JJ),KHIST(JJ)
           ENDDO
           
           NANG = NANG + NAN

           DEALLOCATE(RD,KHIST,KTYPE)

        ENDDO

        IF(NANG.NE.nn(1)*nn(2)*nn(3)*Ncell) print *,"The sampple is incomplete !!!!!!"

        print *,"The sample consist of",NANG,"atoms totally", NKH,"KHIST"
        
111     FORMAT(3(1x,E14.7),2(1x,I2))

        CLOSE(UNIT = 17)

        STOP
      END PROGRAM AU_CRYSTALL
