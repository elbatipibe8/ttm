!     Write output files for further analysis for TTM
!     Leonid Zhigilei and Dmitri Ivanov, 2002
      SUBROUTINE SwriteTTM_Res()
        INCLUDE 'common.h'
        INCLUDE 'commonTTM.h'
        INCLUDE 'mpif.h'
        REAL*8, DIMENSION(:), ALLOCATABLE :: xmsgbuf
        REAL*8, DIMENSION(:), ALLOCATABLE :: TPSXG,RONG,FSOLG,Q1VG
        REAL*8, DIMENSION(:), ALLOCATABLE :: T_lG,T_eG,XEDT,YEDT,ZEDT
        INTEGER statusttm(MPI_STATUS_SIZE)
        CHARACTER*5 chtime, chpid

        IF(mypid.eq.0) CALL CPU_TIME(timestart0)
        
        Call Liquid()
        Call Symmetry()

        ITIME=ANINT(TIME)
        IF (ITIME.LT.1000) Then
           WRITE(chtime, FMT='(I3.3)') ITIME
        ELSEIF (ITIME.LT.10000) Then
           WRITE(chtime, FMT='(I4.4)') ITIME
        ELSEIF (ITIME.LT.100000) Then
           WRITE(chtime, FMT='(I5.5)') ITIME
        ELSE
           PRINT *, "Simulation is too long... Time=",Time," ps."
           STOP
        ENDIF

        WRITE(chpid,FMT='(I5.5)') mypid

        OPEN(UNIT=969,FILE='./'//DRNAME(1:LDR)//'/time'//trim(chtime)// &
             '.'//trim(chpid))

        CALL MPI_BARRIER(MPI_COMM_WORLD,ipass1)

!       Do not write everything to avoid too big output files
        REWIND 969
        WRITE(969,*) NN1
        WRITE(969,166)(0.1d0*XD(1,J),0.1d0*XD(2,J),0.1d0*XD(3,J), &
        100.0d0*Q1D(1,J)/DELTA,100.0d0*Q1D(2,J)/DELTA,100.0d0*Q1D(3,J)/DELTA, &
        POT(J)*ENUNIT,1.0d0-FD(1,J),FD(3,J),KTYPE(J),J=1,NN1)
166     FORMAT(F10.4,1x,F10.4,1x,F10.4,1x &
             ,F8.2,1x,F8.2,1x,F8.2,1x,F7.4,1x,F6.4,1x,F6.4,1x,I2)
        CLOSE(UNIT = 969)

        CALL MPI_BARRIER(MPI_COMM_WORLD,ipass2)

        IF(mypid.eq.0) THEN
           CALL CPU_TIME(time00)
           print *," "
           print *,"time.proc file written",time00-timestart0
           print *," "
!           call Flushout(6)
        ENDIF

!!!!    Write output coordinate file (restart)
        OPEN(UNIT=16,FILE='./Ni.out.'//trim(chpid))

        REWIND 16
        WRITE(16,*) NN1,TIME,ZCRIT
        WRITE(16,*) XL,YL,ZL
        WRITE(16,*) XCENTR,YCENTR,ZCENTR
        WRITE(16,67) (KTYPE(J),XD(1,J),XD(2,J),XD(3,J),J=1,NN1)
        WRITE(16,68) (KHIST(J),Q1D(1,J)/DELTA,Q1D(2,J)/DELTA,Q1D(3,J)/DELTA,J=1,NN1)
        WRITE(16,*) nlcx,nlcy,nlcz
        WRITE(16,*) XLREAL,YLREAL,ZLREAL
        WRITE(16,*) DENWAT,EOUTG,RODENO
        WRITE(16,*) vden,zdist,F_totalG

        DO iz = 1,nlcz
           DO iy = 1,nlcy
              DO ix = 1,nlcx
                 j = ix + (nlcx + 2)*(iy + (nlcy + 2)*iz)
                 WRITE(16,69) Xsh(j),Ysh(j),Zsh(j),XED(j),YED(j),ZED(j), &
                      difx(j),dify(j),difz(j),T_l(j),T_e(j),IBULK(j)
              ENDDO
           ENDDO
        ENDDO

        CLOSE(UNIT = 16)

67      FORMAT(I3,1x,F20.12,1x,F20.12,1x,F20.12)
68      FORMAT(I3,1x,F18.12,1x,F18.12,1x,F18.12)
69      FORMAT(6(F20.12,1x),3(E20.12,1x),2(F14.8,1x),1x,L2)

        CALL MPI_BARRIER(MPI_COMM_WORLD,ipass3)

        IF(mypid.eq.0) THEN
           CALL CPU_TIME(time00)
           print *," "
           print *,"Ni.out.proc file written", time00-timestart0
           print *," " 
!           call Flushout(6)
        ENDIF

!##########################################################################
!##########   Begin TTM writing now    ####################################

!       In case of water we need to collect some additional information here
        IF(NTYPE.GE.2) FLUID = .TRUE.

        CALL TTMTemper()
        CALL Share_TTM()

!       Calculation of pressure for the each link cell
        TPSX(0:nlc2-1) = 0.0d0
        Do iz = 1,nlcz
           Do iy = 1,nlcy
              Loop: Do ix = 1,nlcx
                 ip = ix + nlcx2*(iy + nlcy2*iz)
                 
                 IF(NVK(ip).EQ.0) CYCLE Loop

                 j = ltop(ip)
                 
31               IF(j.EQ.0) GOTO 32
                 
                 IF(.NOT.((KHIST(j).EQ.3.OR.KHIST(j).GT.10).AND.NOTEM)) &
                      TPSX(ip) = TPSX(ip) + &
                      (STEN(j,1,1)+STEN(j,2,2)+STEN(j,3,3))*ENUNIT/3.0d0
                 
                 j = link(j)
                 
                 GOTO 31
                 
32               CONTINUE
                 
              ENDDO Loop
           ENDDO
        ENDDO

        IF(FLUID) THEN
           NWT(0:nlc2-1) = 0
           DO ii = 1,ippl
              ip = idlc(ii)
              NWT(ip) = 1
              i = ltop(ip)

              IF(i.eq.0) GOTO 28

29            IF(KTYPE(I).EQ.KWT) TPSX(ip) = TPSX(ip) + &
                   (STEN(i,1,1)+STEN(i,2,2)+STEN(i,3,3))*ENUNIT/3.0d0

              i = link(i)
              IF(i.ne.0)GOTO 29

28            CONTINUE
           ENDDO
        ENDIF
        
        CALL Share_PR()
        CALL Smooth_PR()
        IF(FLUID) CALL Smooth_PRW()

        Do iz = 1,nlcz
           Do iy = 1,nlcy
              Loop2: Do ix = 1,nlcx
                 ip = ix + nlcx2*(iy + nlcy2*iz)

                 IF(.NOT.(ICOND(ip)))  TPSX(ip) = TPSX(ip) + RON(ip)*T_l(ip)*BK

              ENDDO Loop2
           ENDDO
        ENDDO
        
!       Convert Units [Pr.En.units/Ang^3] --> [J/m^3] == [Pa] => [GPa]
        TPSX(0:nlc2-1)=TPSX(0:nlc2-1)*EVTOJOU*1.0D-09/volcell(0:nlc2-1)

        IF(FLUID) THEN
           DO ii = 1,ippl
              ip = idlc(ii)
              RON(ip) = RON(ip)*RODENO/DENWAT
           ENDDO
        ENDIF

!       Void the aditional data determining until next writing event
        FLUID = .FALSE.
        
        ALLOCATE(xmsgbuf(MSGSTTM))

!       Preparation to write 3D TTM profile
        IF(mypid.ne.0) THEN
           k = 1
           xmsgbuf(k)     = DBLE(mynodex)
           xmsgbuf(k + 1) = DBLE(mynodey)
           xmsgbuf(k + 2) = DBLE(mynodez)
           k = k + 3
           Do iz = 1,nlcz
              Do iy = 1,nlcy
                 Do ix = 1,nlcx
                    ip = ix + nlcx2*(iy + nlcy2*iz)
                    
                    IF(ICOND(ip)) THEN
                       RDENS = 1.0d0
                    ELSE
                       RDENS = RON(ip)/RODENO
                    ENDIF
                    
                    xmsgbuf(k    ) = T_l(ip)
                    xmsgbuf(k + 1) = T_e(ip)
                    xmsgbuf(k + 2) = TPSX(ip)
                    xmsgbuf(k + 3) = RDENS
                    xmsgbuf(k + 4) = FSOL(ip)
                    xmsgbuf(k + 5) = Q1V(ip)*100.0d0 ! m/s
                    xmsgbuf(k + 6) = 0.1d0*XED(ip)   ! nm
                    xmsgbuf(k + 7) = 0.1d0*YED(ip)   ! nm
                    xmsgbuf(k + 8) = 0.1d0*ZED(ip)   ! nm
                    
                    k = k + ITTMR
                 ENDDO
              ENDDO
           ENDDO
           
           msglent = MSGSTTM

!          Send all the information to the processor with ID "0"
           CALL MPI_SEND(xmsgbuf,msglent,MPI_DOUBLE_PRECISION, &
                0,208+nnodes+mynodez+1,MPI_COMM_WORLD,ierr)

           DEALLOCATE(xmsgbuf)

        ENDIF

!       Processor with ID "0" recieves the information and writes it
!       into the output file
        IF(mypid.eq.0) THEN

           ALLOCATE(TPSXG(nlcg),RONG(nlcg),FSOLG(nlcg),T_lG(nlcg), &
                T_eG(nlcg),Q1VG(nlcg),XEDT(nlcg),YEDT(nlcg),ZEDT(nlcg))

           DO iz = 1,nlcz
              DO iy = 1,nlcy
                 DO ix = 1,nlcx
                    ic = ix + nlcxg*((iy-1)+nlcyg*(iz-1))
                    ip = ix + nlcx2*(iy + nlcy2*iz)

                    IF(ICOND(ip)) THEN
                       RDENS = 1.0d0
                    ELSE
                       RDENS = RON(ip)/RODENO
                    ENDIF

                    T_lG(ic)    = T_l(ip)
                    T_eG(ic)    = T_e(ip)
                    TPSXG(ic)   = TPSX(ip)
                    RONG(ic)    = RDENS
                    FSOLG(ic)   = FSOL(ip)
                    Q1VG(ic)    = Q1V(ip)*100.0d0 ! m/s
                    XEDT(ic)    = 0.1d0*XED(ip)   ! nm
                    YEDT(ic)    = 0.1d0*YED(ip)   ! nm
                    ZEDT(ic)    = 0.1d0*ZED(ip)   ! nm

                    k = k + ITTMR
                 ENDDO
              ENDDO
           ENDDO

           DO nz = 1,nznodes
              DO ny = 1,nynodes
                 KLS: DO nx = 1,nxnodes
                    if(nz.eq.1.and.ny.eq.1.and.nx.eq.1) CYCLE KLS
                    msglent = MSGSTTM
                    CALL MPI_RECV(xmsgbuf,msglent,MPI_DOUBLE_PRECISION, &
                         MPI_ANY_SOURCE,208+nnodes+nz,MPI_COMM_WORLD, &
                         statusttm,ierr)
                    ifromttm = statusttm(MPI_SOURCE)

!                   CALL CPU_TIME(time00)
!                   print *,time00-timestart0," SwriteTTM-Res received from",ifromttm
!                   call Flushout(6)

                    k = 1
                    idx = INT(xmsgbuf(k))
                    idy = INT(xmsgbuf(k + 1))
                    idz = INT(xmsgbuf(k + 2))
                    k = k + 3
                    DO iz = 1,nlcz
                       DO iy = 1,nlcy
                          DO ix = 1,nlcx
                             ic = (ix+idx*nlcx) + nlcxg*((iy+idy*nlcy-1) + &
                                  nlcyg*(iz+idz*nlcz-1))

                             T_lG(ic)    = xmsgbuf(k    )
                             T_eG(ic)    = xmsgbuf(k + 1)
                             TPSXG(ic)   = xmsgbuf(k + 2)
                             RONG(ic)    = xmsgbuf(k + 3)
                             FSOLG(ic)   = xmsgbuf(k + 4)
                             Q1VG(ic)    = xmsgbuf(k + 5)
                             XEDT(ic)    = xmsgbuf(k + 6)
                             YEDT(ic)    = xmsgbuf(k + 7)
                             ZEDT(ic)    = xmsgbuf(k + 8)

                             k = k + ITTMR
                          ENDDO
                       ENDDO
                    ENDDO
                 ENDDO KLS
              ENDDO
           ENDDO

           DEALLOCATE(xmsgbuf)

!!!          Search for an average temerature on the edge
!!           Tem_el = 0.0d0
!!           Tem_ph = 0.0d0
!!           kount = 0

!          3D TTTM profile
           IF(KONF.EQ.1) THEN
              OPEN(UNIT=989,FILE='./'//DRNAME(1:LDR)//'/ttmreal'//trim(chtime)//'.d')
              REWIND 989
           ELSE
              OPEN(UNIT=979,FILE='./'//DRNAME(1:LDR)//'/ttm'//trim(chtime)//'.d')
              REWIND 979
           ENDIF
!          To make the future analisys of the data more handy, we write
!          the evenly spaced part of the system into 979 unit
           DO iz = 1,nlczg
              DO iy = 1,nlcyg
                 DO ix = 1,nlcxg
                    ip = ix + nlcxg*((iy-1) + nlcyg*(iz-1))
                    IF(KONF.EQ.1) THEN
                       WRITE(989,515) XEDT(ip),YEDT(ip),ZEDT(ip), &
                            T_lG(ip),T_eG(ip),TPSXG(ip),RONG(ip),FSOLG(ip),Q1VG(ip)
!                       IF(iz.ge.izstart.and.iz.le.izend &
!                            .and.iy.ge.iystart.and.iy.le.iyend &
!                            .and.ix.ge.ixstart.and.ix.le.ixend) &
!                            WRITE(979,515) XEDT(ip),YEDT(ip),ZEDT(ip),T_lG(ip), &
!                            T_eG(ip),TPSXG(ip),RONG(ip),FSOLG(ip),Q1VG(ip)
!!                       IF((iy.eq.1.or.iy.eq.nlcyg) &
!!                           .and.10.0d0*ZEDT(ip).LE.-2.0d0*rstep &
!!                            .AND.10.0d0*ZEDT(ip).GE.(ZCRIT+2.0d0*rstep)) THEN
!!                          kount = kount + 1
!!                          Tem_el = Tem_el + T_eG(ip)
!!                          Tem_ph = Tem_ph + T_lG(ip)
!!                       ENDIF
                    ELSE
                       WRITE(979,515) XEDT(ip),YEDT(ip),ZEDT(ip),T_lG(ip), &
                            T_eG(ip),TPSXG(ip),RONG(ip),FSOLG(ip),Q1VG(ip)
                    ENDIF
                 ENDDO
              ENDDO
           ENDDO
515        FORMAT(3(F9.2,1x),F10.0,1x,F10.0,1x,F7.3,1x,2(F6.3,1x),F8.2)
518        FORMAT(3(F9.2,1x),F10.0,1x,F10.0,1x,F5.3)

           IF(KONF.EQ.1) THEN
              CLOSE(UNIT = 989)
           ELSE
              CLOSE(UNIT = 79)
           ENDIF
           
           DEALLOCATE(TPSXG,RONG,FSOLG,Q1VG,T_eG,T_lG,XEDT,YEDT,ZEDT)

           CALL CPU_TIME(time00)
           print *," "
           print *,"ttm.d and ttmreal.d files written",time00-timestart0
           print *," "
!           call Flushout(6)

!!           IF(KONF.EQ.1) THEN
!!              IF(kount.ne.0) THEN
!!                 print *,"Temperature of phonons on the edge ",Tem_ph/kount,"K"
!!                 print *,"Temperature of electrons on the edge ",Tem_el/kount,"K"
!!                 print *,"Averaging is done over ",kount,"cells"
!!                 print *," "
!!!                 call Flushout(6)
!!              ELSE
!!                 print *,"Temperature on the edge of continuum part is not defined"
!!                 print *,"Number of averaged cells ",kount
!!                 print *," "
!!!                 call Flushout(6)
!!              ENDIF
!!           ENDIF
           
        ENDIF
        
        CALL MPI_BARRIER(MPI_COMM_WORLD,ipass4)
        
        RETURN
      END SUBROUTINE SwriteTTM_Res
