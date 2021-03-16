      SUBROUTINE BCAST_POS()
        INCLUDE 'common.h'
        INCLUDE 'commonTTM.h'
        INCLUDE 'mpif.h'

!       Processor with ID "0" broadcasts all the input variables except 
!       the velocity related ones to all the processors
        CALL MPI_BCAST(NSTEP,NINTEGER,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        CALL MPI_BCAST(TIME,NREAL,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

!       For input/output
        CALL MPI_BCAST(LDR,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        CALL MPI_BCAST(DRNAME,LDR,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)

!       Nonreflective boundary conditions
        IF(KBOUND.EQ.2) CALL MPI_BCAST(FSTART,NNONREFL, &
             MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

!       Variables that have to be passed due to TTM
        IF(LGOT.EQ.2) THEN
           CALL MPI_BCAST(NTST,NTTMIN,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
           CALL MPI_BCAST(WIDTH,NTTMDD,MPI_DOUBLE_PRECISION,0, &
                MPI_COMM_WORLD,ierr)
        ENDIF

        CALL MPI_BARRIER(MPI_COMM_WORLD,ipass)

        RETURN
      END SUBROUTINE BCAST_POS
