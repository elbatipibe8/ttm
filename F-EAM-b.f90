!     Calculation of Forces and Energies with Johnson's EAM potential
      SUBROUTINE F_EAM_b()
        INCLUDE 'common.h'
        INCLUDE 'commonTTM.h'
        INCLUDE 'mpif.h'
        Real *8  RD,RR,PENER,P,DZ,DY,DX,Rhoi,Rhoj,Fi,dFi,dFj,DDenI,DDenJ
        Integer  IJ,KF,RKFi,RKFj,KTYPEI,KTYPEJ,I,K,J,N,M
        DIMENSION FVZ(2),FVZA(2),FVZG(2),NSEQ(nxnodes*nynodes)
        DIMENSION LNONREF(LPMZ),LNONREFA(LPMZ)
        INTEGER status1(MPI_STATUS_SIZE),status2(MPI_STATUS_SIZE)

        CALL EAM_Density()

!       To calculate forces need also information about density in the skin
        CALL SHARE_EAM()

        NRL = 0 ! Renulling of the number of particles in the NRBC in this node
        NSH = 0 !Renulling of the number of particles in the piston in the node
        outer_loop: DO I=1,NN1
           IF (NNG(I).gt.MAXNNB) THEN
              print *, ' TOO MANY NEIGHBOURS: ',NNG(I),' > ', maxnnb
              call flush(6)
              STOP
           ENDIF
           ktypeI=KTYPE(I)
!          (i) EAM Local Electron Density (calculated in CalcDensity())
           Rhoi=RhoDen(I)*NT/RhoM(ktypei)
           RKFi=Rhoi		! RKFi is used below
           Rhoi=Rhoi-RKFi		! Rhoi is used below
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
                 call flush(6)
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
              FD(1,I)=FD(1,I) + (dFi*DDenJ+dFj*DDenI)*DX - P*DX
              FD(2,I)=FD(2,I) + (dFi*DDenJ+dFj*DDenI)*DY - P*DY

!             Non-reflecting boundary - no interaction from ABOVE
              IF(KHIST(I).LE.10) THEN
                 FD(3,I) = FD(3,I) + (dFi*DDenJ+dFj*DDenI)*DZ - P*DZ
              ELSE
                 IF(nrbsign*DZ.LT.0.0d0) FD(3,I) = &
                      FD(3,I) + (dFi*DDenJ+dFj*DDenI)*DZ - P*DZ   
              ENDIF
              
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
           IF(KHIST(I).GT.10) THEN
              NRL = NRL + 1
              LNONREF(NRL) = I
           ENDIF
           
           IF(LFLAG.EQ.2) THEN
!          Since calculation of shocks includes TTM, and TTM is used for 
!          metals only, and metals are treated with EAM, one can be sure
!          that if shock is the case, the following statement will be passed
!          and simplify the shock forces calculation. We shal now complite a
!          search for the piston particles, NSH is set to 0 in the begining    
!          of this subroutine
              IF(KHIST(I).GT.10) THEN
                 NSH = NSH + 1
                 LSHOCK(NSH) = I
              ENDIF
           ENDIF

        End do outer_loop

        EOUT = 0.0d0

!       Adding Forces due to the Non-reflecting Boundary Condition
        DO KK = 1,Nmono
           KKK = KK + 10
           INONREF = 0
           FVZ(1:2) = 0.0d0
           NRLIN = 0
           DO II=1,NRL
              I = LNONREF(II)
              IF (KHIST(I).EQ.KKK) THEN
                 FVZ(1) = FVZ(1) + FD(3,I)
                 FVZ(2) = FVZ(2) + Q1D(3,I)/DELTA
                 NRLIN = NRLIN + 1
                 LNONREFA(NRLIN) = I
                 IF(INONREF.EQ.0) INONREF = 1
              ENDIF
           ENDDO
           
           IF(mypid.NE.0.AND.INONREF.EQ.1) CALL MPI_SEND(FVZ, &
                2,MPI_DOUBLE_PRECISION,0,0,MPI_COMM_WORLD,ierror)
           
           IF(mypid.eq.0)  THEN   
              FVZG(1:2) = 0.0d0
              IF(INONREF.EQ.1) FVZG(1:2) = FVZG(1:2) + FVZ(1:2)
              DO nrec = 1,(nxnodes*nynodes - INONREF)
                 CALL MPI_RECV(FVZA,2,MPI_DOUBLE_PRECISION, &
                      MPI_ANY_SOURCE,0,MPI_COMM_WORLD,status1,ierror)
                 NSEQ(nrec) = status1(MPI_SOURCE)
                 FVZG(1:2) = FVZG(1:2) + FVZA(1:2)
              ENDDO
              
              DO nrec = 1,(nxnodes*nynodes - INONREF)
                 CALL MPI_SEND(FVZG,2,MPI_DOUBLE_PRECISION, &
                      NSEQ(nrec),0,MPI_COMM_WORLD,ierror)
              ENDDO
              IF(INONREF.EQ.1) FVZ(1:2) = FVZG(1:2)
           ENDIF
           
           IF(mypid.NE.0.AND.INONREF.EQ.1) CALL MPI_RECV(FVZ, &
                2,MPI_DOUBLE_PRECISION,0,0,MPI_COMM_WORLD,status2,ierror)
           
           FVZ(1:2) = FVZ(1:2)/Ponrefl(KK)
           FVZ(1) = FVZ(1) + FSTART - GIMPED*FVZ(2)

           DO II = 1,NRLIN
              I = LNONREFA(II)
              IF(KHIST(I).EQ.KKK) THEN
                 FD(3,I) = FVZ(1)
                 EOUT = EOUT + GIMPED*FVZ(2)*FVZ(2)*ENUNIT*DELTA/DBLE(Nmono)
              ENDIF
           ENDDO
        ENDDO

        RETURN
      END SUBROUTINE F_EAM_b
