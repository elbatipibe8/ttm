!     Calculation of Forces and Energies with Johnson's EAM potential
      SUBROUTINE F_EAM_br()
        INCLUDE 'common.h'
        INCLUDE 'commonTTM.h'
        INCLUDE 'mpif.h'
        Real *8  RD,RR,PENER,P,DZ,DY,DX,Rhoi,Rhoj,Fi,dFi,dFj,DDenI,DDenJ
        Integer  IJ,KF,RKFi,RKFj,KTYPEI,KTYPEJ,I,K,J,N,M
        INTEGER, DIMENSION(:), ALLOCATABLE :: LNONREF
        DIMENSION FVZ(2),FVZG(2)

        ALLOCATE(LNONREF(LPMZ))

        CALL EAM_Density()

!       To calculate forces need also information about density in the skin
        CALL SHARE_EAM()

        NRL = 0 !Renulling of the number of particles in the NRBC
        FVZ(1:2) = 0.0d0 ! Renulling the accumulating force due to NRB
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
           IF(KHIST(I).EQ.10) THEN
              RI2 = (XD(1,I) - XCENTR)*(XD(1,I) - XCENTR) + &
                   (XD(2,I) - YCENTR)*(XD(2,I) - YCENTR)
              RIC = SQRT(RI2)
           ENDIF

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
                 Print *, 'Atoms ',I,' and ',J,' are too close.', &
                      'Program will stop.'
!                call Flushout(6)
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
              PENER=.5*(UT(IJ,KF)+(UT(IJ,KF+1)-UT(IJ,KF))*P)
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
              IF(KHIST(I).EQ.10) THEN
                 RJ2 = (XD(1,J) - XCENTR)*(XD(1,J) - XCENTR) + &
                      (XD(2,J) - YCENTR)*(XD(2,J) - YCENTR)
                 DJI2 = RJ2 - RI2
                 IF(DJI2.LE.0.0d0) THEN
                    FD(1,I) = FD(1,I) + (dFi*DDenJ+dFj*DDenI)*DX - P*DX
                    FD(2,I) = FD(2,I) + (dFi*DDenJ+dFj*DDenI)*DY - P*DY
                 ENDIF
              ELSE
                 FD(1,I) = FD(1,I) + (dFi*DDenJ+dFj*DDenI)*DX - P*DX
                 FD(2,I) = FD(2,I) + (dFi*DDenJ+dFj*DDenI)*DY - P*DY
              ENDIF

!             Add F component as it is
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
           IF(KHIST(I).EQ.10) THEN
              FFr=(FD(1,I)*(XD(1,I)-XCENTR)+ FD(2,I)*(XD(2,I)-YCENTR))/RIC
              FVZ(1) = FVZ(1) + FFr
              QQr=(Q1D(1,I)*(XD(1,I)-XCENTR)+Q1D(2,I)*(XD(2,I)-YCENTR))/RIC
              FVZ(2) = FVZ(2) + QQr/DELTA
              NRL = NRL + 1
              LNONREF(NRL) = I
           ENDIF

        End do outer_loop

!       There is no reasonable way to predict in which node the part of NR
!       boundaries happens to be so that we prefer usual MPI_REDUCE over
!       blocking MPI_SEND/MPI_RECEIVE
        CALL MPI_REDUCE(FVZ,FVZG,2,MPI_DOUBLE_PRECISION, &
             MPI_SUM,0,MPI_COMM_WORLD,ierr)
        CALL MPI_BCAST(FVZG,2,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        
        FVZG(1:2) = FVZG(1:2)/DBLE(NCIRC)

        IF(mypid.eq.0) THEN
           WRITE(411,412) FVZG(1) + FSTART - GIMPED*FVZG(2), &
                FVZG(2),FVZG(1),FSTART - GIMPED*FVZG(2)
412        FORMAT(4(E14.7,1x))
           call Flushout(411)
        ENDIF

        FVZG(1) = FVZG(1) + FSTART - GIMPED*FVZG(2)

        EOUT = 0.0d0
!       Adding/substracting forces due to the NRB
        DO II=1,NRL
           I = LNONREF(II)
           RIC = SQRT((XD(1,I) - XCENTR)*(XD(1,I) - XCENTR) + &
                (XD(2,I) - YCENTR)*(XD(2,I) - YCENTR))
           FD(1,I) = FVZG(1)*(XD(1,I) - XCENTR)/RIC
           FD(2,I) = FVZG(1)*(XD(2,I) - YCENTR)/RIC
!          The final energy eout should be divided to the number of
!          monolayers in the  boundary zone, but the definition of a
!          monolayer is questionable here. Based on the initial sample
!          construction We may approximately chose equal to 2
           EOUT = EOUT + GIMPED*FVZG(2)*FVZG(2)*ENUNIT*DELTA/2.0d0
        ENDDO
        
        DEALLOCATE(LNONREF)

        RETURN
      END SUBROUTINE F_EAM_br
