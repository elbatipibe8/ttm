!     Write output files for further analysis
!     For pair potentials and EAM
!     Leonid Zhigilei, 2002
      SUBROUTINE Swrite_Res()
        INCLUDE 'common.h' 
        INCLUDE 'mpif.h'
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
!       Available information: J,KHISTG(J),XDG(1,J),XDG(2,J),XDG(3,J),
!       Q1DG(1,J),Q1DG(2,J),Q1DG(3,J),POTG(J),
!       QING(J),ENG(J),TEMPAG(J),TAtPressG(J),IDGG(J),WALG(J),WALAG(J)
!       If you add anything, correct the format of writing as well
        REWIND 969
        WRITE(969,*) NN1
        WRITE(969,166)(0.1d0*XD(1,J),0.1d0*XD(2,J),0.1d0*XD(3,J), &
        100.0d0*Q1D(1,J)/DELTA,100.0d0*Q1D(2,J)/DELTA,100.0d0*Q1D(3,J)/DELTA, &
        POT(J)*ENUNIT,1.0d0-FD(1,J),FD(3,J),KHIST(J),J=1,NN1)
166     FORMAT(F10.4,1x,F10.4,1x,F10.4,1x &
             ,F8.2,1x,F8.2,1x,F8.2,1x,F7.4,1x,F6.4,1x,F6.4,1x,I2)
        CLOSE(UNIT = 969)
        
        CALL MPI_BARRIER(MPI_COMM_WORLD,ipass2)

        IF(mypid.eq.0) THEN
           CALL CPU_TIME(time00)  
           print *,"time.proc file written", time00-timestart0
           print *," "
!           call Flushout(6)
        ENDIF

        OPEN(UNIT=16,FILE='./Al.'//trim(chpid))

        TIME = 0.0d0
         ZL = 3000.0d0
!        ZL = 10000.0d0
!        ZCENTR = -4500.0d0
!        XL = 257.195265530d0
!        YL = XL
!        ZL = ZL

        
!!!!    Write output coordinate file (restart)
        REWIND 16
        WRITE(16,*) NN1,TIME,ZCRIT
        WRITE(16,*) XL,YL,ZL
        WRITE(16,*) XCENTR,YCENTR,ZCENTR
        WRITE(16,67) (KTYPE(J),XD(1,J),XD(2,J),XD(3,J),J=1,NN1)
        WRITE(16,68) (KHIST(J),Q1D(1,J)/DELTA,Q1D(2,J)/DELTA,Q1D(3,J)/DELTA,J=1,NN1)
        CLOSE(UNIT = 16)
67      FORMAT(I3,1x,F20.12,1x,F20.12,1x,F20.12)
68      FORMAT(I3,1x,F18.12,1x,F18.12,1x,F18.12)

        CALL MPI_BARRIER(MPI_COMM_WORLD,ipass3)

        IF(mypid.eq.0) THEN
           CALL CPU_TIME(time00)
           print *,"Ni.proc file written", time00-timestart0
           print *," "
!           call Flushout(6)
        ENDIF

        RETURN
      END SUBROUTINE Swrite_Res
