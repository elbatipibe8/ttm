!     Calculation of Forces and Energies with Johnson's EAM potential
      SUBROUTINE F_EAM_br()
        INCLUDE 'common.h'
        INCLUDE 'commonTTM.h'
        INCLUDE 'mpif.h'
        Real *8  RD,RR,PENER,P,DZ,DY,DX,Rhoi,Rhoj,Fi,dFi,dFj,DDenI,DDenJ
        Integer  IJ,KF,RKFi,RKFj,KTYPEI,KTYPEJ,I,K,J,N,M
        INTEGER, DIMENSION(:), ALLOCATABLE :: LNONREF
        DIMENSION FVZ(6),FVZG(6),FVZA(6)
        DIMENSION NSEQ(3*nxnodes*nynodes)
        INTEGER status1(MPI_STATUS_SIZE),status2(MPI_STATUS_SIZE)

        ALLOCATE(LNONREF(LPMZ))

        CALL EAM_Density()

!       To calculate forces need also information about density in the skin
        CALL SHARE_EAM()

        NRL = 0 !Renulling of the number of particles in the NRBC
        FVZ(1:6) = 0.0d0 ! Renulling the accumulating force due to NRB
        INON = 0
        outer_loop: DO I=1,NN1

           IF(KTYPE(I).EQ.KWT) CYCLE outer_loop

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
!           IF(KHIST(I).GT.10) RI = (XD(2,I) - YCENTR)

           inner_loop: DO K=1,NNG(I)
              J=NNNG(k,i)

              IF(KTYPE(J).EQ.KWT) CYCLE inner_loop
              
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
              IF(KHIST(I).GT.10) THEN
!                 RJ = XD(2,J) - YCENTR
                 DJI = XD(3,J)-XD(3,I)!DABS(RJ) - DABS(RI)
                 IF(DJI.GT.0.0d0.AND.KHIST(I).NE.KHIST(J)) FD(3,I) = &
                      FD(3,I) + (dFi*DDenJ+dFj*DDenI)*DZ - P*DZ
              ELSE
                 FD(3,I) = FD(3,I) + (dFi*DDenJ+dFj*DDenI)*DZ - P*DZ
              ENDIF

!             Add F component X and Z as it is
              FD(1,I) = FD(1,I) + (dFi*DDenJ+dFj*DDenI)*DX - P*DX
              FD(2,I) = FD(2,I) + (dFi*DDenJ+dFj*DDenI)*DY - P*DY
           End do inner_loop

!          Potential from EAM Embedding Function (with Rho as input)
           IF(RKFi.GE.1) Then
              Fi=EFT(ktypeI,RKFi)+(EFT(ktypeI,RKFi+1)-EFT(ktypeI,RKFi))*Rhoi
           ELSE
              Fi=EFT(ktypeI,1)+(EFT(ktypeI,1)-EFT(ktypeI,2))*(1-Rhoi)
              IF (Rhoi.EQ.0.0) Fi=0.0d0
           ENDIF  ! This takes care for low el. density, e.g. for dimers

           POT(I)=POT(I)+Fi

           IF(KHIST(I).GT.10) THEN
              KK = KHIST(I)-10
              FVZ(KK) = FVZ(KK) + FD(3,I)
              FVZ(KK+3) = FVZ(KK+3) + Q1D(3,I)/DELTA
              NRL = NRL + 1
              LNONREF(NRL) = I
           ENDIF

        End do outer_loop

        IF(NRL.NE.0) INON = 1

        CALL MPI_REDUCE(INON,INONG,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)

        IF(mypid.NE.0.AND.INON.EQ.1) CALL MPI_SEND(FVZ, &
             6,MPI_DOUBLE_PRECISION,0,0,MPI_COMM_WORLD,ierror)

        IF(mypid.eq.0)  THEN
           FVZG(1:6) = 0.0d0
           FVZG(1:6) = FVZG(1:6) + FVZ(1:6)*INON
           DO nrec = 1,INONG-INON
              CALL MPI_RECV(FVZA,6,MPI_DOUBLE_PRECISION, &
                   MPI_ANY_SOURCE,0,MPI_COMM_WORLD,status1,ierror)
              NSEQ(nrec) = status1(MPI_SOURCE)
              FVZG(1:6) = FVZG(1:6) + FVZA(1:6)
           ENDDO

           FVZG(1:6) = FVZG(1:6)/DBLE(NCIRC)

           DO nrec = 1,INONG-INON
              CALL MPI_SEND(FVZG,6,MPI_DOUBLE_PRECISION, &
                   NSEQ(nrec),0,MPI_COMM_WORLD,ierror)
           ENDDO
        ENDIF


        IF(mypid.NE.0.AND.INON.EQ.1) CALL MPI_RECV(FVZG, &
             6,MPI_DOUBLE_PRECISION,0,0,MPI_COMM_WORLD,status2,ierror)

        IF(mypid.eq.0) THEN
           WRITE(411,412) FVZG(1) - FSTART - GIMPED*FVZG(4), &
                FVZG(4),FVZG(1),-FSTART - GIMPED*FVZG(4),ISTEP
412        FORMAT(4(E14.7,1x),I6)
           call Flushout(411)
        ENDIF

        EOUT = 0.0d0
!       Adding/substracting forces due to the NRB
        DO II=1,NRL
           I = LNONREF(II)
           KK = KHIST(I) - 10
           FD(3,I) = FVZG(KK) - FSTART - GIMPED*FVZG(KK+3)
!          The final energy eout should be divided to the number of
!          monolayers in the  boundary zone
           EOUT = EOUT + GIMPED*FVZG(KK+3)*FVZG(KK+3)*ENUNIT*DELTA/3.0d0
        ENDDO

        DEALLOCATE(LNONREF)

        RETURN
      END SUBROUTINE F_EAM_br
