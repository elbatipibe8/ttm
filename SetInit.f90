!     THIS SUBROUTINE SETS THE INITIAL PARAMETERS
!     Perform this with all involved processors
!     Leonid Zhigilei, 2000
      SUBROUTINE SetInit()
        INCLUDE 'common.h'
        INCLUDE 'commonTTM.h' 
        
!       Number of types of particle pairs
        Npots=Ntype*(Ntype+1)/2  

!        IF(LFLAG.EQ.1) THEN
!           IF(.NOT.(LIDX.EQ.1.AND.LIDY.EQ.1.AND.LIDZ.EQ.1))THEN
!              print *,"Set all boundary to periodic"
!              print *,"The progaram is halted now"
!              call Flushout(6)
!              CALL MPI_FINALIZE(info)
!              STOP
!           ENDIF
!        ENDIF

!       Let's define a few variables
        XLHALF=XL*0.5d0
        YLHALF=YL*0.5d0
        ZLHALF=ZL*0.5d0
        
!       Define the sign of sample position with respect to the origin
        IF(XCENTR.NE.0.0d0) THEN
           nxsign  = XCENTR/ABS(XCENTR)
        ELSE
           nxsign = 1
        ENDIF
        IF(YCENTR.NE.0.0d0) THEN
           nysign  = YCENTR/ABS(YCENTR)
        ELSE
           nysign = 1
        ENDIF
        IF(ZCENTR.NE.0.0d0) THEN
           nzsign  = ZCENTR/ABS(ZCENTR)
        ELSE
           nzsign = 1
        ENDIF
!       If this authomatical position definition does not give you what
!       you wanted, you may set theses axes manually by uncommentig the
!       following statements
        nxsign = 1
        nysign = 1
        nzsign = 1

        IF(mypid.eq.0) THEN
           print *,"nxsign =",nxsign
           print *,"nysign =",nysign
           print *,"nzsign =",nzsign
           print *,"  "
        ENDIF
        
!       Cutoff for liquid/solid distinguishing
        CUTLIQ = 0.05d0
        CUTLIQ2= 0.11d0

!       If have a specific layer in which no temperature should be defined
        IF((KBOUND.EQ.2.AND.LIDZ.EQ.0).OR.KBOUND.EQ.1.OR.LFLAG.EQ.2)THEN
           NOTEM = .TRUE.
        ELSE
           NOTEM = .FALSE.
        ENDIF

!       Initial value for energy loss due to nonreflective boundary conditions
        IF(TIME.EQ.0.0d0) EOUTG = 0.0d0

!       Let's assume that we have up to 3 types of atoms
!       interacting via up to 6 pair potentials
        IJINDEX(1,1)=1
        IJINDEX(1,2)=3
        IJINDEX(1,3)=6
        IJINDEX(2,1)=3
        IJINDEX(2,2)=2
        IJINDEX(2,3)=5
        IJINDEX(3,1)=6
        IJINDEX(3,2)=5
        IJINDEX(3,3)=4
        
!       G1 is used in Nordsieck Integration method
        DO J=1,NTYPE
           G1(J)=0.5d0*DELTA*DELTA/XMASS(J)
        ENDDO

!       Nordsieck predictor-corrector coefficients
        CC1=3.d0/16.d0
        C2=251.d0/360.d0
        C3=11.d0/18.d0
        C4=1.d0/6.d0
        C5=1.d0/60.d0
        
        RETURN
      END SUBROUTINE SetInit
