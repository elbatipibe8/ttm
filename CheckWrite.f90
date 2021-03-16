      Subroutine CheckWrite()
        INCLUDE 'common.h'
        INCLUDE 'mpif.h'

        IF(mypid.eq.0) THEN
           CLOSE(UNIT = 14)
           OPEN(UNIT = 14,FILE = 'md.input')
           REWIND 14
           READ(14,100) NSTEP,MATER,NEPRT,NWRITE,NTYPE,KFLAG,LFLAG, &
                KEYBS,KEYDEP,NDEP,LIDX,LIDY,LIDZ,KBOUND,LGOT,NDIM
           READ(14,200) QTEM,QPRESS,AddEn,DELTA
        ENDIF
100     FORMAT(I10)
200     FORMAT(D15.5)

        CALL MPI_BCAST(NSTEP,NINTEGER,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        CALL MPI_BCAST(TIME,NREAL,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

        CALL MPI_BARRIER(MPI_COMM_WORLD,ipass)

        IF(mypid.eq.0) THEN
           print *,"Indicator of NRB is set to",KBOUND,"NSTEP is set to",NSTEP
           print *," "
!           call Flushout(6)
        ENDIF

        RETURN

      END Subroutine CheckWrite
