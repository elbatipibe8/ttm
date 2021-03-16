!     Implementing Berendsen barostat for keeping constant P
!     At each time step coordinates are scaled by the factor Dzeta.
!     TPRESS is external pressure in GPa
!     Berendsen et al., J. Chem. Phys., 81 (1984), p. 3684.
      SUBROUTINE BERE_P (TPRESS)
        INCLUDE 'common.h'
        INCLUDE 'mpif.h'
        REAL*8, PARAMETER :: pbeta=0.5d0

        CALL MPI_REDUCE(TPRESS,TPRESSG,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
             MPI_COMM_WORLD,ierr)
        CALL MPI_BCAST(TPRESSG,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        
        DZETA=1.0d0-pbeta*DELTA*(QPRESS-TPRESSG*1.0d-09)

        ZL=ZL*DZETA
        YL=YL*DZETA
        XL=XL*DZETA
        xlstart=XCENTR-XL/2.0d0
        ylstart=YCENTR-YL/2.0d0
        zlstart=ZCENTR-ZL/2.0d0
        XD(1,1:NN1)=(XD(1,1:NN1)-XCENTR)*DZETA+XCENTR
        XD(2,1:NN1)=(XD(2,1:NN1)-YCENTR)*DZETA+YCENTR
        XD(3,1:NN1)=(XD(3,1:NN1)-ZCENTR)*DZETA+ZCENTR
        ximage(1:26)=ximage(1:26)*DZETA 
        yimage(1:26)=yimage(1:26)*DZETA
        zimage(1:26)=zimage(1:26)*DZETA
        
        RETURN
      END SUBROUTINE BERE_P

!     Implementing van Gunstern-Berendsen thermostat for keeping constant T
!     in a layer of atoms with KHIST(J)=2.  At each time step velocities are 
!     scaled by the factor Vscale.
!     Berendsen, J. Chem. Phys., 81 (1984), p. 3684.
!     Although this method does not reproduce canonical ensemble, it is 
!     widely used, and usually gives the same results as more rigorous methods
!     such as GLE implemented below. It reproduces the correct average energy,
!     but the distribution is wrong. Therefore, averages are usually correct,
!     but fluctuations are not. 
!     Leonid Zhigilei, 2002
      SUBROUTINE BERE_T(TEMPTR)
        INCLUDE 'common.h'
        INCLUDE 'mpif.h'
        
        CALL MPI_REDUCE(TEMPTR,TEMPTRG,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
             MPI_COMM_WORLD,ierr)
        CALL MPI_BCAST(TEMPTRG,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        
!     TAU is the time constant that defines the strength of coupling to 
!     the thermal bath.  It is usually choosen from the range of 0.1 - 2 ps.
        TAU=0.5d0  
        Vscale=SQRT(1.0d0+(QTEM/TEMPTRG-1.0d0)*DELTA/TAU)
        
        DO I=1,NN1
           Q1D(1,I)=Q1D(1,I)*Vscale
           Q1D(2,I)=Q1D(2,I)*Vscale
           Q1D(3,I)=Q1D(3,I)*Vscale
        ENDDO

        RETURN
      END SUBROUTINE BERE_T
