!     Calculation of Forces and Energies with Johnson's EAM potential
      SUBROUTINE F_EAM_br_test()
        INCLUDE 'common.h'
        INCLUDE 'commonTTM.h'
        INCLUDE 'mpif.h'
        Real*8  RD,RR,PENER,P,DZ,DY,DX,Rhoi,Rhoj,Fi,dFi,dFj,DDenI,DDenJ
        Integer  IJ,KF,RKFi,RKFj,KTYPEI,KTYPEJ,I,K,J,N,M
        INTEGER, DIMENSION(:), ALLOCATABLE :: LNONREF
        REAL*8 FDD(NN1)
        REAL*8 FVZ(9),FVZG(9)

        ALLOCATE(LNONREF(LPMZ))

        CALL EAM_Density()

!       To calculate forces need also information about density in the skin
        CALL SHARE_EAM()

        FDD(1:NN1) = 0.0d0
        NRL = 0 !Renulling of the number of particles in the NRBC
        FVZ(1:9) = 0.0d0 ! Renulling the accumulating force due to NRB
        outer_loop: DO I=1,NN1
           IF (NNG(I).gt.MAXNNB) THEN
              print *, ' TOO MANY NEIGHBOURS: ',NNG(I),' > ', maxnnb
!              call Flushout(6)
              STOP
           ENDIF
           ktypeI=KTYPE(I)
!          (i) EAM Local Electron Density (calculated in CalcDensity())
           Rhoi=RhoDen(I)*NT/RhoM(ktypei)
           RKFi=Rhoi		! RKFi is used below
           Rhoi=Rhoi-RKFi		! Rhoi is used below
!          Callulate the distance to the centre if this particle is in NRB
!           IF(KHIST(I).EQ.10) RI = XD(2,I) - YCENTR

           inner_loop: DO K=1,NNG(I)
              J=NNNG(k,i)
              ktypeJ=KTYPE(J)
              DX=XD(1,J)-XD(1,I)
              DY=XD(2,J)-XD(2,I)
              DZ=XD(3,J)-XD(3,I)
              
              IJ = ijindex(ktypeI, ktypeJ)
              RR=DX*DX+DY*DY+DZ*DZ
              RD=SQRT(RR)
              
              IF (RD.GE.RM(IJ)) CYCLE inner_loop
              
              IF (RD.LE.RM(IJ)/NT) THEN
                 Print *, 'Atoms ',I,' and ',J,' are too close.',XD(1:3,I),XD(1:3,J), &
                      'Program will stop.',NN1,mypid
!                 call Flushout(6)
                 STOP
              ENDIF

!             (j) EAM Local Electron Density (calculated in CalcDensity())
              Rhoj=RhoDen(j)*NT/RhoM(ktypeJ)
              RKFj=Rhoj
              Rhoj=Rhoj-RKFj
              IF(RKFj.GE.1) Then
                 dFj=DEFT(ktypeJ,RKFj)+(DEFT(ktypeJ,RKFj+1)- &
                      DEFT(ktypeJ,RKFj))*Rhoj
              ELSE 
                 dFj=DEFT(ktypeJ,1)+(DEFT(ktypeJ,1)-DEFT(ktypeJ,2))*(1-Rhoj)
              ENDIF  ! This takes care for low el. density, e.g. for dimers

!             (i) EAM Local Electron Density (calculated in CalcDensity())
              IF(RKFi.GE.1) Then
                 dFi=DEFT(ktypeI,RKFi)+(DEFT(ktypeI,RKFi+1)- &
                      DEFT(ktypeI,RKFi))*Rhoi
              ELSE 
                 dFi=DEFT(ktypeI,1)+(DEFT(ktypeI,1)-DEFT(ktypeI,2))*(1-Rhoi)
              ENDIF  ! This takes care for low el. density, e.g. for dimers

              P=RD*DXR(IJ)
              KF=P
              P=P-KF

!             Pair Potential energy (give half the potential to this particle)
              PENER=0.5d0*(UT(IJ,KF)+(UT(IJ,KF+1)-UT(IJ,KF))*P)
              POT(I)=POT(I)+PENER
!             End of potential energy

              DDenI=DDFT(ktypeI,KF)+(DDFT(ktypeI,KF+1)-DDFT(ktypeI,KF))*P
              DDenJ=DDFT(ktypeJ,KF)+(DDFT(ktypeJ,KF+1)-DDFT(ktypeJ,KF))*P
              DDenI=DDenI/RD
              DDenJ=DDenJ/RD
              
              P=FT(IJ,KF)+(FT(IJ,KF+1)-FT(IJ,KF))*P
              P=P/RD
              
!             For calculation of static portion of stress tensor (Force*rij)
!             The sign is oposite (multiply by -1 if stresses are of interest)
!             (we are not taking into account the terminating forces here)
              STEN(I,1,1)=STEN(I,1,1)+(P-(dFi*DDenJ+dFj*DDenI))*DX*DX
              STEN(I,2,2)=STEN(I,2,2)+(P-(dFi*DDenJ+dFj*DDenI))*DY*DY
              STEN(I,3,3)=STEN(I,3,3)+(P-(dFi*DDenJ+dFj*DDenI))*DZ*DZ

!             Force Calculation (in neg direction)
!             Non-reflecting boundary - no interaction from ABOVE
              IF(KHIST(I).GT.10) THEN
                 DJI = XD(3,J)-XD(3,I)
                 IF(DJI.LT.0.0d0.AND.KHIST(I).NE.KHIST(J)) FDD(I) = &
                      FDD(I) + (dFi*DDenJ+dFj*DDenI)*DZ - P*DZ
              ENDIF
              
              FD(1,I) = FD(1,I) + (dFi*DDenJ+dFj*DDenI)*DX - P*DX
              FD(2,I) = FD(2,I) + (dFi*DDenJ+dFj*DDenI)*DY - P*DY
              FD(3,I) = FD(3,I) + (dFi*DDenJ+dFj*DDenI)*DZ - P*DZ              
           End do inner_loop

!          Potential from EAM Embedding Function (with Rho as input)
           IF(RKFi.GE.1) Then
              Fi=EFT(ktypeI,RKFi)+(EFT(ktypeI,RKFi+1)-EFT(ktypeI,RKFi))*Rhoi
           ELSE 
              Fi=EFT(ktypeI,1)+(EFT(ktypeI,1)-EFT(ktypeI,2))*(1-Rhoi)
              IF (Rhoi.EQ.0.0) Fi=0.0d0
           ENDIF  ! This takes care for low el. density, e.g. for dimers

           POT(I)=POT(I)+Fi

!          We shal now complete a search for the nonreflective layer particles,
!          NRL is set to 0 in the beginning of this subroutine.
!          Accumulate the radial forces and velocities here on RIC projection
           IF(KHIST(I).GT.10) THEN
              KK = KHIST(I)-10

              FVZ(KK) = FVZ(KK) + FD(3,I)
              
              FVZ(KK+3) = FVZ(KK+3) + Q1D(3,I)/DELTA

              FVZ(KK+6) = FVZ(KK+6) + FDD(I)
              
              NRL = NRL + 1
              LNONREF(NRL) = I
           ENDIF
           
        End do outer_loop

!       There is no reasonable way to redict in which node the part of NR
!       boundaries happens to be so that we prefer usual MPI_REDUCE over
!       blocking MPI_SEND/MPI_RECEIVE
        CALL MPI_REDUCE(FVZ,FVZG,9,MPI_DOUBLE_PRECISION, &
             MPI_SUM,0,MPI_COMM_WORLD,ierr)
        CALL MPI_BCAST(FVZG,9,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        CALL MPI_BARRIER(MPI_COMM_WORLD,ipass)

        FVZG(1:9) = FVZG(1:9)/DBLE(NCIRC)

!       Plot V3 versusu V2 to see the force you need to apply versus velocity
        IF(mypid.eq.0) THEN
           WRITE(411,412)FVZG(1:9),ISTEP
412        FORMAT(9(E14.7,1x),I6)
           CALL Flushout(411)
        ENDIF

        DO II=1,NRL
           I = LNONREF(II)
!           RIC = SQRT((XD(1,I) - XCENTR)*(XD(1,I) - XCENTR) + &
!                (XD(2,I) - YCENTR)*(XD(2,I) - YCENTR))
           KK = KHIST(I)-10
           FD(3,I) = FVZG(KK)
           Q1D(3,I) = FVZG(KK+3)*DELTA
        ENDDO

        DEALLOCATE(LNONREF)

        RETURN
      END SUBROUTINE F_EAM_br_test
