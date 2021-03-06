       SUBROUTINE FIRSTLIST(DENSM,DENSN)
  
        INCLUDE 'common.h'
        INCLUDE 'commonTTM.h'
        INCLUDE 'mpif.h'
        REAL*8 ARMG(6),ARM(6)
        INTEGER nodeix(26),nodeiy(26),nodeiz(26)
        INTEGER noresx(26),noresy(26),noresz(26)
        INTEGER status(MPI_STATUS_SIZE)
        INTEGER, DIMENSION(:), ALLOCATABLE :: LNONREF
        CHARACTER*5 chpid
        LOGICAL, PARAMETER :: NRBK = .FALSE.

!       Define coordinates of 26 relative neighbore for each processor (node)
        DATA noresx/1,0,-1,1,0,-1,1,0,-1,1,0,-1,1, &
             -1,1,0,-1,1, 0,-1,1,0,-1,1,0,-1/
        DATA noresy/1,1,1,0,0,0,-1,-1,-1,1,1,1,0, &
             0,-1,-1,-1,1,1,1,0,0,0,-1,-1,-1/
        DATA noresz/1,1,1,1,1,1,1,1,1,0,0,0,0, &
             0,0,0,0,-1,-1,-1,-1,-1,-1,-1,-1,-1/

        DATA nodeix/-1,0,1,-1,0,1,-1,0,1,-1,0,1,-1, &
             1,-1,0,1,-1,0,1,-1,0,1,-1,0,1/
        DATA nodeiy/-1,-1,-1,0,0,0,1,1,1,-1,-1,-1,0, &
             0,1,1,1,-1,-1,-1,0,0,0,1,1,1/
        DATA nodeiz/-1,-1,-1,-1,-1,-1,-1,-1,-1,0,0,0,0, &
             0,0,0,0,1,1,1,1,1,1,1,1,1/

        CALL MPI_BARRIER(MPI_COMM_WORLD,ipass0)
        IF(mypid.eq.0) THEN
           print *," "
           print *,"Start reading coordinates/velocities"
           print *,' '
!           call Flushout(6)
        ENDIF

        IF(mypid.eq.0) THEN
           CLOSE(UNIT = 15)
           CLOSE(UNIT = 16)
        ENDIF
	
        WRITE(chpid,FMT='(I5.5)') mypid

        IF(TIME.EQ.0.0d0) THEN
           OPEN(UNIT=16,FILE='./Al.'//trim(chpid))
        ELSE
          OPEN(UNIT=16,FILE='./Al.out.'//trim(chpid))
        ENDIF

        READ(16,*) NN1,TIME,ZCRIT
        READ(16,*) XL,YL,ZL
        READ(16,*) XCENTR,YCENTR,ZCENTR
        READ(16,*) (KTYPE(J),XD(1,J),XD(2,J),XD(3,J),J=1,NN1)
        READ(16,*) (KHIST(J),Q1D(1,J),Q1D(2,J),Q1D(3,J),J=1,NN1)

        IF(NRBK) THEN
           DO i=1,NN1
              IF(XD(3,i).GE.-100.0d0) XD(3,i) = &
                   (XD(3,i)+100.0d0)*(-0.02d0*XD(3,i)/100.0d0+0.98d0)-100.0d0
           ENDdO
        ENDIF

        IF(NN1.GT.LPMZ) THEN
           print *,"Number of particles in", mypid,"node"
           print *,NN1,"exceeded the maximum allowed value in"
           print *,"FirstList LPMZ=",LPMZ,"Increase RMD or ZSIZE in"
           print *,"parameters.h The porogram is halted now"
           STOP
        ENDIF
        
        CALL MPI_REDUCE(NN1,NAN,1,MPI_INTEGER,MPI_SUM, &
             0,MPI_COMM_WORLD,ierr)
        CALL MPI_BCAST(NAN,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        
        IF(NAN.GT.LPMX0) THEN
           IF(mypid.eq.0) THEN
              print *,"The number of particles NAN=",NAN," exceeds LPMX0=",LPMX0
              print *,"The program is halted"
           ENDIF
           CALL MPI_FINALIZE(info)
           STOP
        ENDIF

        CALL MPI_BARRIER(MPI_COMM_WORLD,ipass1)

        IF(LGOT.EQ.2.AND.TIME.NE.0.0d0) THEN
           IF(mypid.eq.0) THEN
              print *,"Read TTM input part from the previous run"
              print *," "
!              call Flushout(6)
           ENDIF

           READ(16,*) nlcx,nlcy,nlcz
           READ(16,*) XLREAL,YLREAL,ZLREAL
           READ(16,*) DENWAT,EOUTG,RODENO
           READ(16,*) vden,zdist,F_totalG

           DO iz = 1,nlcz
              DO iy = 1,nlcy
                 DO ix = 1,nlcx
                    j = ix + (nlcx + 2)*(iy + (nlcy + 2)*iz)
                    READ(16,*) Xsh(j),Ysh(j),Zsh(j),XED(j),YED(j),ZED(j), &
                         difx(j),dify(j),difz(j),T_l(j),T_e(j),IBULK(j)
                 ENDDO
              ENDDO
           ENDDO
           
        ENDIF

        CLOSE(UNIT = 16)

        CALL MPI_BARRIER(MPI_COMM_WORLD,ipass2)

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!       Define the real geometry of the system for correct pressure definition
        IF(.NOT.(LGOT.EQ.2.AND.TIME.NE.0.0d0)) THEN

!          If water is in calculations, we do it for water first
           IF(NTYPE.GE.2) THEN
              ARM(1) = XCENTR + XL
              ARM(2) = XCENTR - XL
              ARM(3) = YCENTR + YL
              ARM(4) = YCENTR - YL
              ARM(5) = ZCENTR + ZL
              ARM(6) = ZCENTR - ZL

              DO I=1,NN1
                 IF(KTYPE(I).EQ.KWT) THEN   ! If not water molecule
                    IF(ARM(1).GE.XD(1,I)) ARM(1) = XD(1,I)
                    IF(ARM(2).LE.XD(1,I)) ARM(2) = XD(1,I)
                    IF(ARM(3).GE.XD(2,I)) ARM(3) = XD(2,I)
                    IF(ARM(4).LE.XD(2,I)) ARM(4) = XD(2,I)
                    IF(ARM(5).GE.XD(3,I)) ARM(5) = XD(3,I)
                    IF(ARM(6).LE.XD(3,I)) ARM(6) = XD(3,I)
                 ENDIF
              ENDDO

              IF(mypid.ne.0) CALL MPI_SEND(ARM,6,MPI_DOUBLE_PRECISION, &
                   0,0,MPI_COMM_WORLD,ierror)

              IF(mypid.eq.0)  THEN
                 ARMG(1:6) = ARM(1:6)
                 DO nrec = 1,nnodes-1
                    CALL MPI_RECV(ARM,6,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,0, &
                         MPI_COMM_WORLD,status,ierror)
                    DO j=1,3
                       IF(ARMG(2*j-1).ge.ARM(2*j-1)) ARMG(2*j-1) = ARM(2*j-1)
                       IF(ARMG(2*j).le.ARM(2*j)) ARMG(2*j) = ARM(2*j)
                    ENDDO
                 ENDDO
              ENDIF

              CALL MPI_BCAST(ARMG,6,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

              IF(LIDZ.EQ.0) THEN
                 ZLREALW = ARMG(6) - ARMG(5)
              ELSE
                 ZLREALW = ZL
              ENDIF
              IF(LIDX.EQ.0) THEN
                 XLREALW = ARMG(2) - ARMG(1)
              ELSE
                 XLREALW = XL
              ENDIF
              IF(LIDY.EQ.0) THEN
                 YLREALW = ARMG(4) - ARMG(3)
              ELSE
                 YLREALW = YL
              ENDIF

!             Water density calculation
              VOLUMEW=XLREALW*YLREALW*ZLREALW*1.0d-24!-32 ! cm3
              TMASSW=0.0d0
              NWAT = 0
              DO J=1,NN1
                 IF(KTYPE(J).EQ.KWT) THEN
                    TMASSW=TMASSW+XMASS(KTYPE(J))       ! Total mass in amu
                    NWAT = NWAT + 1
                 ENDIF
              ENDDO

              CALL MPI_REDUCE(TMASSW,TMASSGW,1,MPI_DOUBLE_PRECISION, &
                   MPI_SUM,0,MPI_COMM_WORLD,ierr)
              CALL MPI_BCAST(TMASSGW,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
              CALL MPI_REDUCE(NWAT,NWATG,1,MPI_INTEGER, &
                   MPI_SUM,0,MPI_COMM_WORLD,ierr)
              CALL MPI_BCAST(NWATG,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
              
!             The following concwerns water only
              DENSMW=TMASSGW*AMUTOKG*1.0D+03/VOLUMEW  ! gr/cm3
              DENSNW=NWATG/VOLUMEW
              DENWAT = NWATG/(XLREALW*YLREALW*ZLREALW)
           ENDIF

!*************************************************************
!          Next, we do the same for metallic particles
           ARM(1) = XCENTR + XL
           ARM(2) = XCENTR - XL
           ARM(3) = YCENTR + YL
           ARM(4) = YCENTR - YL
           ARM(5) = ZCENTR + ZL
           ARM(6) = ZCENTR - ZL

           DO I=1,NN1
              IF(KTYPE(I).NE.KWT) THEN   ! If not water molecule
                 IF(ARM(1).GE.XD(1,I)) ARM(1) = XD(1,I)
                 IF(ARM(2).LE.XD(1,I)) ARM(2) = XD(1,I)
                 IF(ARM(3).GE.XD(2,I)) ARM(3) = XD(2,I)
                 IF(ARM(4).LE.XD(2,I)) ARM(4) = XD(2,I)
                 IF(ARM(5).GE.XD(3,I)) ARM(5) = XD(3,I)
                 IF(ARM(6).LE.XD(3,I)) ARM(6) = XD(3,I)
              ENDIF
           ENDDO

           IF(mypid.ne.0) CALL MPI_SEND(ARM,6,MPI_DOUBLE_PRECISION, &
        	    0,0,MPI_COMM_WORLD,ierror)	

           IF(mypid.eq.0)  THEN
              ARMG(1:6) = ARM(1:6)
              DO nrec = 1,nnodes-1
                 CALL MPI_RECV(ARM,6,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,0, &
                      MPI_COMM_WORLD,status,ierror)
                 DO j=1,3
                    IF(ARMG(2*j-1).ge.ARM(2*j-1)) ARMG(2*j-1) = ARM(2*j-1)
                    IF(ARMG(2*j).le.ARM(2*j)) ARMG(2*j) = ARM(2*j)
                 ENDDO
              ENDDO
           ENDIF

           CALL MPI_BCAST(ARMG,6,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

           IF(LGOT.EQ.2.AND.TIME.EQ.0.0d0) THEN
              XD(3,1:NN1) = XD(3,1:NN1) + DABS(ARMG(6))
              ZCRIT = ZCRIT + DABS(ARMG(6))
              ZCENTR = ZCENTR + DABS(ARMG(6))
           ENDIF
           
!           IF(LGOT.EQ.2.OR.KBOUND.EQ.2) THEN
!              IF(ARMG(5).GE.ZCRIT) ARMG(5) = ZCRIT
!              IF(ARMG(6).LE.ZCRIT) THEN
!                 IF(mypid.eq.0) print *,"an atom below ZCRIT"
!                 STOP
!              ENDIF
!           ENDIF

!          The following concerns only metal atoms
           IF(LIDZ.EQ.0) THEN
              ZLREAL = ARMG(6) - ARMG(5)
           ELSE
              ZLREAL = ZL
           ENDIF
           IF(LIDX.EQ.0) THEN
              XLREAL = ARMG(2) - ARMG(1)
           ELSE
              XLREAL = XL
           ENDIF
           IF(LIDY.EQ.0) THEN
              YLREAL = ARMG(4) - ARMG(3)
           ELSE
              YLREAL = YL
           ENDIF

!          Density calculation
           VOLUME=XLREAL*YLREAL*ZLREAL*1.0d-24!-32 ! cm3
           TMASS=0.0d0
           NMET = 0
           DO J=1,NN1
              IF(KTYPE(J).NE.KWT) THEN
                 TMASS=TMASS+XMASS(KTYPE(J))       ! Total mass in amu
                 NMET = NMET + 1
              ENDIF
           ENDDO
        
           CALL MPI_REDUCE(TMASS,TMASSG,1,MPI_DOUBLE_PRECISION, &
                MPI_SUM,0,MPI_COMM_WORLD,ierr) 
           CALL MPI_BCAST(TMASSG,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
           CALL MPI_REDUCE(NMET,NMETG,1,MPI_INTEGER, &
                MPI_SUM,0,MPI_COMM_WORLD,ierr)
           CALL MPI_BCAST(NMETG,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

!          The following concwerns metal only
           DENSM=TMASSG*AMUTOKG*1.0D+03/VOLUME  ! gr/cm3
           DENSN=NMETG/VOLUME                     ! molecules/cm3

           CALL MPI_BARRIER(MPI_COMM_WORLD,ipass3)

           IF(mypid.eq.0) THEN
              print *,"Initial setting complited"
              print *,"Number of metal atoms is:", NMETG
              print *,"Density of metal target is:", DENSM,"g/cm3"
              print *,"Density of metal particles is:", DENSN, "at/cm3"
              print *," "
              print *,"Number of water atoms is:", NWATG
              print *,"Density of water volume is:", DENSMW,"g/cm3"
              print *,"Density of water particles is:", DENSNW, "at/cm3"
              print *," "
              !           call Flushout(6)
           ENDIF
        ENDIF

        Q1D(1:3,1:NN1) = Q1D(1:3,1:NN1)*DELTA
        Q2(1:3*NN1) = 0.0d0
        Q3(1:3*NN1) = 0.0d0
        Q4(1:3*NN1) = 0.0d0
        Q5(1:3*NN1) = 0.0d0

!        IF(KBOUND.EQ.11) THEN
!           NCIRCC = 0
!           DO I=1,NN1
!              IF(KHIST(I).GT.10) NCIRCC = NCIRCC + 1
!           ENDDO
!              CALL MPI_REDUCE(NCIRCC,NCIRC,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
!              CALL MPI_BCAST(NCIRC,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!        ENDIF

        IF(KBOUND.EQ.2.OR.KBOUND.EQ.11) THEN
           DO KK = 11,13
              ALLOCATE(LNONREF(LPMZ))
              QQ = 0.0d0
              XDD = 0.0d0
              NCIR = 0
              DO I=1,NN1
                 IF(KHIST(I).EQ.KK) THEN
                       QQ = QQ + Q1D(3,I)
                       NCIR = NCIR + 1
                       XDD = XDD + XD(3,I)
                       LNONREF(NCIR) = I
                 ENDIF
              ENDDO
              
              CALL MPI_REDUCE(QQ,QQG,1,MPI_DOUBLE_PRECISION, &
                   MPI_SUM,0,MPI_COMM_WORLD,ierr)
              CALL MPI_REDUCE(XDD,XDDG,1,MPI_DOUBLE_PRECISION, &
                   MPI_SUM,0,MPI_COMM_WORLD,ierr)
              CALL MPI_REDUCE(NCIR,NCIRC,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
              
              CALL MPI_BARRIER(MPI_COMM_WORLD,ipass4)
              
              CALL MPI_BCAST(QQG,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
              CALL MPI_BCAST(XDDG,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
              CALL MPI_BCAST(NCIRC,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

              IF(mypid.eq.0) THEN
                 print *,"Velocity of NRB is set to",QQG/DELTA/DBLE(NCIRC),"for",KK,"KHIST"
                 print *,"Position of NRB is set to",XDDG/DBLE(NCIRC),"for",KK,"KHIST"
                 print *,"Number of particles in NRB with",KK,"KHIST is",NCIRC
!                 call Flushout(6)
              ENDIF
              
              DO II=1,NCIR
                 I = LNONREF(II)
                 Q1D(3,I) = QQG/DBLE(NCIRC)
!                 XD(3,I)  = XDDG/DBLE(NCIRC)
                 Q1D(1,I) = 0.0d0
                 Q1D(2,I) = 0.0d0
              ENDDO

              DEALLOCATE(LNONREF)

           ENDDO
        ENDIF

        ZNRB = XDDG/DBLE(NCIRC)+alat ! Initial position of NRB

        CALL MPI_BARRIER(MPI_COMM_WORLD,ipass5)

        IF(mypid.eq.0) THEN
           print *," "
           print *,"Assignment of the input data completed"
           print *," "
!           call Flushout(6)
        ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!   SETNODES   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       Calculate max and min x,y,z coordinates for
!       this node's region in true system coordinates

!       Define the system geography
        mynodex = MOD(mypid,nxnodes)
        mynodey = MOD((mypid/nxnodes),nynodes)
        mynodez = mypid/(nxnodes*nynodes)

        xdiff   = nxsign*XL/DBLE(nxnodes)
        xlstart = XCENTR - nxsign*XLHALF
        xnmin   = xlstart + DBLE(mynodex)*xdiff
        xnmax   = xnmin + xdiff

        ydiff   = nysign*YL/DBLE(nynodes)
        ylstart = YCENTR - nysign*YLHALF
        ynmin   = ylstart + DBLE(mynodey)*ydiff
        ynmax   = ynmin + ydiff

        zdiff   = nzsign*ZL/DBLE(nznodes)
        zlstart = ZCENTR - nzsign*ZLHALF
        znmin   = zlstart +  DBLE(mynodez)*zdiff
        znmax   = znmin + zdiff

        CALL MPI_BARRIER(MPI_COMM_WORLD,ipass2)

!       Define auxiliary system geometry in terms of
!       long, midle, and short coordinates
        nlong  = MAX(nznodes,nynodes,nxnodes)
        IF(nlong.eq.nznodes) THEN
           longnode = mynodez
           longarray = (/1,2,3,4,5,6,7,8,9,18,19,20,21,22,23,24,25,26/)
           nmidl = MAX(nxnodes,nynodes)
           IF(nmidl.eq.nxnodes) THEN
              midlnode = mynodex
              midlarray = (/10,12,13,14,15,17/)
              nshort = nynodes
              kornode = mynodey
              korarray = (/11,16/)
           ELSE
              midlnode = mynodey
              midlarray = (/10,11,12,15,16,17/)
              nshort = nxnodes
              kornode = mynodex
              korarray = (/13,14/)
           ENDIF
        ELSEIF(nlong.eq.nxnodes) THEN
           longnode = mynodex
           longarray = (/1,3,4,6,7,9,10,12,13,14,15,17,18,20,21,23,24,26/)
           nmidl = MAX(nynodes,nznodes)
           IF(nmidl.eq.nynodes) THEN
              midlnode = mynodey
              midlarray = (/ 2, 8,11,16,19,25/)
              nshort = nznodes
              kornode = mynodez
              korarray = (/ 5,22/)
           ELSE
              midlnode = mynodez
              midlarray = (/ 2, 5, 8,19,22,25/)
              nshort = nynodes
              kornode = mynodey
              korarray = (/11,16/)
           ENDIF
        ELSE
           longnode = mynodey
           longarray = (/1,2,3,7,8,9,10,11,12,15,16,17,18,19,20,24,25,26/)
           nmidl = MAX(nxnodes,nznodes)
           IF(nmidl.eq.nxnodes) THEN
              midlnode = mynodex
              midlarray = (/ 4, 6,13,14,21,23/)
              nshort = nznodes
              kornode = mynodez
              korarray = (/ 5,22/)
           ELSE
              midlnode = mynodez
              midlarray = (/ 4, 5, 6,21,22,23/)
              nshort = nxnodes
              kornode = mynodex
              korarray = (/13,14/)
           ENDIF
        ENDIF

!       We will use DXR(I) to pick a value of force or energy that
!       corresponds to a given interparticle distance from the table
        DO I=1,NPOTS
           DXR(I)=NT/RM(I)
           RLIST(I) = RM(I) + 0.5d0
           RLIST2(I) = RLIST(I)*RLIST(I)
        ENDDO

!       Step of Neighbore List renewal
        NLREN = 10

!       Define number of linkcells in each direction
        IF(LGOT.EQ.2.AND.TIME.NE.0.0d0) GOTO 555
        IF(LGOT.EQ.2) THEN
           rsize = RM(1) + 3.25d0
        ELSE
           rsize = RLIST(1)
        ENDIF
        nlcx = INT(XL/(nxnodes*rsize))
        nlcy = INT(YL/(nynodes*rsize))
        nlcz = INT(ZL/(nznodes*rsize))
555     CONTINUE

        REDCUT = MAX(XL/DBLE(nxnodes)/DBLE(nlcx), &
             YL/DBLE(nynodes)/DBLE(nlcy), &
             ZL/DBLE(nznodes)/DBLE(nlcz))

        IF(mypid.eq.0) print *,REDCUT,"size of a cell in Angstrons"

        IF(mypid.eq.0) THEN
           WRITE(17,*)"number of linked cells is based on cut-off",rsize
           WRITE(17,*)"Each node has:",nlcx,nlcy,nlcz
           WRITE(17,*)"linkcells in X, Y, and Z directions"
           WRITE(17,*)"The sample axes directions are: (positive/negative)"
           WRITE(17,*)"X =",nxsign,",Y =",nysign,", Z =",nzsign
           WRITE(17,*)"if you intended to have them differnt from that,"
           WRITE(17,*)"set them manually in SetInit.f90"
           WRITE(17,*)"      "
           call Flushout(17)
           print *,"number of linked cells is based on cut-off",rsize
           print *,"Each node has:",nlcx,nlcy,nlcz
           print *,"linkcells in X, Y, and Z directions"
           print *,"The sample axes directions are: (positive/negative)"
           print *,"X =",nxsign,",Y =",nysign,", Z =",nzsign
           print *,"if you intended to have them differnt from that,"
           print *,"set them manually in SetInit.f90"
           print *,"      "
!           call Flushout(6)
        ENDIF

!       Define some auxiliry parameters
        nlc = nlcx*nlcy*nlcz
        nlcx1 = nlcx+1
        nlcy1 = nlcy+1
        nlcz1 = nlcz+1
        nlcx2 = nlcx+2
        nlcy2 = nlcy+2
        nlcz2 = nlcz+2
        nlc2 = nlcx2*nlcy2*nlcz2
        nlcxg = nlcx * nxnodes
        nlcyg = nlcy * nynodes
        nlczg = nlcz * nznodes
        nlcg = nlcxg * nlcyg * nlczg

        IF(nlc.gt.mnlc) THEN
           IF(mypid.eq.0) THEN
              print *,"mnlc =",mnlc," is less then the expected number"
              print *,"of linked cells",nlcx*nlcy*nlcz
              print *,"Increse mnlc in parameters.h The program is halted now"
!              call Flushout(6)
           ENDIF
           CALL MPI_FINALIZE(info)
           STOP
        ENDIF

        IF(nlc2.gt.mnlc2) THEN
           IF(mypid.eq.0) THEN
              print *,"mnlc2 =",mnlc2," is less then the expected number"
              print *,"of linked cells",nlc2,"Increse mnlc2."
              print *,"The program is halted now"
!              call Flushout(6)
           ENDIF
           CALL MPI_FINALIZE(info)
           STOP
        ENDIF

        IF(nlcz.gt.mnlcz) THEN
           IF(mypid.eq.0) THEN
              print *,"mnlcz =",mnlcz," is less then the expected number"
              print *,"of linked cells",nlcz,"Increse mnlcz in parameters.h"
              print *,"The program is halted now"
!              call Flushout(6)
           ENDIF
           CALL MPI_FINALIZE(info)
           STOP
        ENDIF

        IF(nlcg.gt.mnlcg) THEN
           IF(mypid.eq.0) THEN
              print *,"mnlcg =",mnlcg," is less then the expected number"
              print *,"of linked cells",nlcg,"Increse mnlcg in parameters.h"
              print *,"The program is halted now"
!              call Flushout(6)
           ENDIF
           CALL MPI_FINALIZE(info)
           STOP
        ENDIF

!****************************************************************
!       Define some parameters for MPI
        IMIGRATE = IAHEAD
        ISHARE = ICOME
        ISWRITE = IWRITE
        ISWRITE2= IWRITE2

!       Want to reduce the amount of consumed memory but still
!       do not want to go beyond an array's boundary
!       Assume no more than ~% of particles can migrate accross
!       node's boundaries per one MD step
        MSGMIG = IMIGRATE * LPMZ !/ 1 !0 !0
        MSGSHA = ISHARE * LPMX
        MSGSWR = ISWRITE * LPMZ
        MSGSWR2 = IWRITE2* LPMZ
        MSGSWBS= ISWBS * nlcz

        IF(LGOT.EQ.2) THEN
           MSGTTM  = ITTM  * nlc2
           MSGTTM2 = ITTMS * nlc2
           MSGTTM3 = ITTM3 * nlc2
           MSGTTMD = ITTMD * nlc2
           MSGTTMDSL = ITTMDSL * nlc2
           MSGSWR  = IWRITE * LPMZ
           MSGSTTM = nlc * ITTMR + 3
           MSGSTTM2= nlc * ITTM2 + 3
           MSGTTMB = nlc2 * ITTMBULK
           MSGTTMB2= nlc2 * ITTMBULK2
        ENDIF

        IF(LIDX.EQ.1.AND.LIDY.EQ.1.AND.LIDZ.EQ.1) ALLPBC = .TRUE.

!       By uncommenting the following line, you will force the variable
!       ALLPBC to be true. That will save some calculation time, but you
!       have to be sure that the total number of particles NAN does not
!       change during all the calculation process. For example, you may
!       take ZL with some extra space to account for the system expansion.
!       Otherwise, you will accumulating errors anywhere you have
!       variable NAN involved.

!       ALLPBC = .TRUE.

!**************************************************************

!       Define IDs for all 26 neighbores of a processor
!       Will use these ID in SEND and RECEIVE routine
!       Also define the image of particles accordign to PBC
!       Loop through all the neighbors
        DO nb = 1,26
!          Prevent communication between nodes if no periodic
!          boundary in this direction is active
           ncom(nb) = 1

           nbnodex = mynodex + nodeix(nb)
           nodresx = mynodex + noresx(nb)
           ximage(nb) = 0.0d0

           IF(nbnodex.LT.0) THEN
              nbnodex = nxnodes - 1
              ximage(nb) = XL
              IF(LIDX.EQ.0) ncom(nb) = 0
           ELSEIF(nbnodex.gt.nxnodes-1) THEN
              nbnodex = 0
              ximage(nb) = -XL
              IF(LIDX.EQ.0) ncom(nb) = 0
           ENDIF

           IF(nodresx.LT.0) THEN
              nodresx = nxnodes - 1
           ELSEIF(nodresx.gt.nxnodes-1) THEN
              nodresx = 0
           ENDIF

           nbnodey = mynodey + nodeiy(nb)
           nodresy = mynodey + noresy(nb)
           yimage(nb) = 0.0d0

           IF(nbnodey.LT.0) THEN
              nbnodey = nynodes - 1
              yimage(nb) = YL
              IF(LIDY.EQ.0) ncom(nb) = 0
              ELSEIF(nbnodey.gt.nynodes-1) THEN
              nbnodey = 0
              yimage(nb) = -YL
              IF(LIDY.EQ.0) ncom(nb) = 0
           ENDIF

           IF(nodresy.LT.0) THEN
              nodresy = nynodes - 1
           ELSEIF(nodresy.gt.nynodes-1) THEN
              nodresy = 0
           ENDIF

           nbnodez = mynodez + nodeiz(nb)
           nodresz = mynodez + noresz(nb)
           zimage(nb) = 0.0d0

           IF(nbnodez.LT.0) THEN
              nbnodez = nznodes - 1
              zimage(nb) = ZL
              IF(LIDZ.EQ.0) ncom(nb) = 0
              ELSEIF(nbnodez.gt.nznodes-1) THEN
              nbnodez = 0
              zimage(nb) = -ZL
              IF(LIDZ.EQ.0) ncom(nb) = 0
           ENDIF

           IF(nodresz.LT.0) THEN
              nodresz = nznodes - 1
           ELSEIF(nodresz.gt.nznodes-1) THEN
              nodresz = 0
           ENDIF

!          Generate PID of neighboring node for sending part
           nbpid(nb) = nbnodex+nxnodes*(nbnodey+nynodes*nbnodez)

!          Generate PID of neighboring node for receiving part
           nbres(nb) = nodresx+nxnodes*(nodresy+nynodes*nodresz)

        ENDDO

!       Determine bounds of skin and core border region
        cell_list: DO nb = 1,26

           IF(nodeix(nb).EQ.0) THEN
              nblox = 1
              nbhix = nlcx
              kslox = 1
              kshix = nlcx
              ELSEIF(nodeix(nb).EQ.(-1)) THEN
              nblox = 1
              nbhix = 1
              kslox = 0
              kshix = 0
           ELSE
              nblox = nlcx
              nbhix = nlcx
              kslox = nlcx1
              kshix = nlcx1
           ENDIF

           IF(nodeiy(nb).EQ.0) THEN
              nbloy = 1
              nbhiy = nlcy
              ksloy = 1
              kshiy = nlcy
              ELSEIF(nodeiy(nb).EQ.(-1)) THEN
              nbloy = 1
              nbhiy = 1
              ksloy = 0
              kshiy = 0
           ELSE
              nbloy = nlcy
              nbhiy = nlcy
              ksloy = nlcy1
              kshiy = nlcy1
           ENDIF

           IF(nodeiz(nb).EQ.0) THEN
              nbloz = 1
              nbhiz = nlcz
              ksloz = 1
              kshiz = nlcz
              ELSEIF(nodeiz(nb).EQ.(-1)) THEN
              nbloz = 1
              nbhiz = 1
              ksloz = 0
              kshiz = 0
           ELSE
              nbloz = nlcz
              nbhiz = nlcz
              ksloz = nlcz1
              kshiz = nlcz1
           ENDIF

!          Generate list of skin cells (scell) and core border
!          cells (bcell) whose content must be sent to this neighbore
!          Number of such cells is stored in nscell/nbcell
!          Building list here simpllifies migrate and share subroutines

           next = 0

           DO iz = ksloz,kshiz
              DO iy = ksloy,kshiy
                 DO ix = kslox,kshix
                    ic = ix + nlcx2*(iy + nlcy2*iz)
                    next = next + 1
                    kscell(nb,next) = ic
                 ENDDO
              ENDDO
           ENDDO

           nscell(nb) = next
           next = 0

           DO iz = nbloz,nbhiz
              DO iy = nbloy,nbhiy
                 DO ix = nblox,nbhix
                    ic = ix + nlcx2*(iy + nlcy2*iz)
                    next = next + 1
                    kbcell(nb,next) = ic
                 ENDDO
              ENDDO
           ENDDO

           nbcell(nb) = next

        ENDDO cell_list

        CALL MPI_BARRIER(MPI_COMM_WORLD,ipass6)

        RETURN
      END SUBROUTINE FIRSTLIST
