      SUBROUTINE Collect()
        INCLUDE 'common.h'
        INCLUDE 'commonTTM.h'
        INCLUDE 'mpif.h'

        CALL Temper(QINT,POTT,TEMPTR)
        CALL Pressure(TPRESS)

!       Coolect the following information for the future analysis into
!       the processor with ID "0"
        CALL MPI_REDUCE(QINT,QINTG,1,MPI_DOUBLE_PRECISION, &
             MPI_SUM,0,MPI_COMM_WORLD,ierr)
        CALL MPI_REDUCE(POTT,POTTG,1,MPI_DOUBLE_PRECISION, &
             MPI_SUM,0,MPI_COMM_WORLD,ierr)
        CALL MPI_REDUCE(TEMPTR,TEMPTRG,1,MPI_DOUBLE_PRECISION, &
             MPI_SUM,0,MPI_COMM_WORLD,ierr)
        CALL MPI_REDUCE(TPRESS,TPRESSG,1,MPI_DOUBLE_PRECISION, &
             MPI_SUM,0,MPI_COMM_WORLD,ierr)
        CALL MPI_REDUCE(Etotal_e,Etotal_e0,1,MPI_DOUBLE_PRECISION, &
             MPI_SUM,0,MPI_COMM_WORLD,ierr)
        CALL MPI_REDUCE(ETTM,ETTMG,1,MPI_DOUBLE_PRECISION, &
             MPI_SUM,0,MPI_COMM_WORLD,ierr)
        CALL MPI_REDUCE(E_pere,E_pere0,1,MPI_DOUBLE_PRECISION, &
             MPI_SUM,0,MPI_COMM_WORLD,ierr)
        CALL MPI_REDUCE(ZZLIQ,ZZLIQG,1,MPI_DOUBLE_PRECISION, &
             MPI_SUM,0,MPI_COMM_WORLD,ierr)
        CALL MPI_REDUCE(EOUT,EOUT0,1,MPI_DOUBLE_PRECISION, &
             MPI_SUM,0,MPI_COMM_WORLD,ierr)
        CALL MPI_REDUCE(F_total,F_total0,1,MPI_DOUBLE_PRECISION, &
             MPI_SUM,0,MPI_COMM_WORLD,ierr)

        E_pere   = 0.0d0
        Etotal_e = 0.0d0
        F_total  = 0.0d0

        IF(mypid.eq.0) THEN
           WRITE(17,9) ISTEP,QINTG+POTTG,QINTG,POTTG,TEMPTRG,TPRESSG*1.0d-09 
           CALL Flushout(17)

!          EOUTG: Energy of the pressure wave transferred through the NRB, eV
!          Etotal_eG: Energy of electrons, eV
!          EcrG: Thermal energy of continuum part, eV
!          E_pereG: Energy transferred from electrons to lattice, eV
!          E_meltingG: Energy of melting, eV
!          EMOVG: Kinetic energy of elastic vibr.(coolective motion),eV
!          ETHERG: Kinetic energy of thermal motion of atoms, eV 

           E_pereG = E_pereG + E_pere0*XJOUTOEV

           IF(LGOT.EQ.2) THEN
              F_totalG = F_totalG + F_total0/(XL*YL*1.0d-20) !J/m2
           ELSE
              F_totalG = 0.0d0
           ENDIF
           Etotal_eG = Etotal_e0*XJOUTOEV
           EOUTG = EOUTG + EOUT0
           EcrG = ETTMG*XJOUTOEV
           E_meltingG = ZZLIQG*FUSION
           
           WRITE(100,91) TIME,(QINTG+POTTG+Etotal_eG+EcrG+EOUTG), &
                QINTG+POTTG,QINTG,POTTG,Etotal_eG,EOUTG,ZL, &!EcrG, &
                F_totalG,TEMPTRG,TPRESSG,(Etotal_eG+E_pereG), &
                (QINTG+POTTG+EcrG+EOUTG-E_pereG),E_pereG,E_meltingG, &
                Etotal_eG+ENTG-P_ZER+EcrG+EOUTG
           CALL Flushout(100)
           print *,"Total Energy =",QINTG+POTTG+Etotal_eG+EcrG+EOUTG, "Total Fluence =",F_totalG/10.0d0,"mJ/cm2"
           print *,"Total number of particles =",NAN
           print *,"Temperature =",TEMPTRG,"K ; Pressure = ",TPRESSG*1.0d-09,"GPa"
           print *,"      "
!           call Flushout(6)

           IF (ISTEP.GE.INT(0.8d0*NSTEP)) THEN
              A_TEM=(A_TEM*DBLE(NTTPR)+TEMPTRG)/DBLE(NTTPR+1)
              A_Press=(A_Press*DBLE(NTTPR)+TPressG)/DBLE(NTTPR+1)
              A_XL = (A_XL*DBLE(NTTPR)+XL)/DBLE(NTTPR+1)
              NTTPR = NTTPR + 1
           ENDIF
        ENDIF
9       FORMAT(I7,3(1X,E14.7),1X,F6.1,1X,F7.4)
91      FORMAT(16(1x,D20.12))

        RETURN
      END SUBROUTINE Collect
