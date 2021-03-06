!     Performes initial setting for TTM model
      SUBROUTINE SetTTM()
        INCLUDE 'common.h'
        INCLUDE 'commonTTM.h'
        INCLUDE 'mpif.h'
        CHARACTER*5 chpid

        IF(mypid.eq.0) THEN
           IF(TIME.EQ.0) THEN
              E_pereG   = 0.0d0
              F_totalG  = 0.0d0
           ENDIF
           OPEN (UNIT = 100,File='TTMenergy.out')  ! Energies/T vs time
        ENDIF
        
        IF(TIME.EQ.0.0d0) RODENO = DBLE(NMETG)*xdiff*ydiff*zdiff/ &
             (XLREAL*YLREAL*ZLREAL)/DBLE(nlcx*nlcy*nlcz)

        IF(TIME.EQ.0.0d0.AND.NTYPE.GE.2) DENWAT = DBLE(NWATG)*xdiff*ydiff*zdiff/ &
             (XLREALW*YLREALW*ZLREALW)/DBLE(nlcx*nlcy*nlcz)
        
        FLUID = .FALSE.

        CALL TTMInit()

        IF(TIME.EQ.0.0d0) THEN
           T_e(0:nlc2-1) = 0.0d0
           DO iz = 1,nlcz
              DO iy = 1,nlcy
                 DO ix = 1,nlcx
                    ic = ix + nlcx2*(iy + nlcy2*iz)
                    IF(NVK(ic).EQ.1) T_e(ic) = Teminit
                    IF(IBULK(ic).AND.RON(ic).LT.RMINR*RODENO) T_l(ic)=Teminit
                 ENDDO
              ENDDO
           ENDDO
        ENDIF

!       Drawing initial atomic density: neigbores/atom
        dmcount = 0
        dncount = 0
        DO iz = 1,nlcz
           DO iy = 1,nlcy
              Tem: DO ix = 1,nlcx

                 ip = ix + nlcx2*(iy + nlcy2*iz)

                 IF(ICOND(ip).OR.NVK(ip).EQ.0) CYCLE Tem

                 i = ltop(ip)
                 IF(i.eq.0) GOTO 24

23               dmcount = dmcount + DBLE(NNG(i))
                 dncount = dncount + 1.0d0

                 i = link(i)
                 IF(i.ne.0) GOTO 23

24               CONTINUE

              ENDDO Tem
           ENDDO
        ENDDO

        CALL MPI_REDUCE(dmcount,dmcountg,1,MPI_DOUBLE_PRECISION, &
             MPI_SUM,0,MPI_COMM_WORLD,ierr)
        CALL MPI_REDUCE(dncount,dncountg,1,MPI_DOUBLE_PRECISION, &
             MPI_SUM,0,MPI_COMM_WORLD,ierr)
        vden = dmcountg/dncountg

        CALL MPI_BCAST(vden,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        IF(mypid.eq.0.and.dncountg.EQ.0) THEN
           print *,"There is a problem with initial density"
           print *,"The program will be halted"
           CALL MPI_FINALIZE(info)
           STOP
        ENDIF
        
        CALL TTMTemper()
        CALL Share_TTM()

        CALL MPI_BARRIER(MPI_COMM_WORLD,ipass)

        EGG(1) = 0.5d0*A_e*TeG(1)*TeG(1)
        DO k = 2,NGT
           EGG(k) = EGG(k-1) + 0.5d0*(CeN(k)+CeN(k-1))*(TeG(k)-TeG(k-1))
        ENDDO

        CALL Solid()

        CALL Coeff()
        
        RETURN  
      END SUBROUTINE SetTTM
