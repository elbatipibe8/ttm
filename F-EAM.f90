!     Calculation of Forces and Energies with Johnson's EAM potential
      SUBROUTINE F_EAM()
        INCLUDE 'common.h'
        INCLUDE 'commonTTM.h'
        Real *8  RD,RR,PENER,P,DZ,DY,DX,Rhoi,Rhoj,Fi,dFi,dFj,DDenI,DDenJ
        Integer  IJ,KF,RKFi,RKFj,KTYPEI,KTYPEJ,I,K,J,N,M
        
        CALL EAM_Density()
!       To calculate forces need also information about density in the skin
        CALL SHARE_EAM()

        outer_loop: DO I=1,NN1

!           IF(KTYPE(I).EQ.KWT) CYCLE outer_loop
           
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
           inner_loop: DO K=1,NNG(I)
              J=NNNG(k,i)

!              IF(KTYPE(J).EQ.KWT) CYCLE inner_loop
              
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

!             Force Calculation (in neg direction)
              FD(1,I)=FD(1,I)+(dFi*DDenJ+dFj*DDenI)*DX
              FD(2,I)=FD(2,I)+(dFi*DDenJ+dFj*DDenI)*DY
              FD(3,I)=FD(3,I)+(dFi*DDenJ+dFj*DDenI)*DZ
              
              P=FT(IJ,KF)+(FT(IJ,KF+1)-FT(IJ,KF))*P
              P=P/RD
              
!             For calculation of static portion of stress tensor (Force*rij)
!             The sign is oposite ( multiply by -1 if stresses are of interest)
              STEN(I,1,1)=STEN(I,1,1)+(P-(dFi*DDenJ+dFj*DDenI))*DX*DX
              STEN(I,2,2)=STEN(I,2,2)+(P-(dFi*DDenJ+dFj*DDenI))*DY*DY
              STEN(I,3,3)=STEN(I,3,3)+(P-(dFi*DDenJ+dFj*DDenI))*DZ*DZ
            
              DX=P*DX                          ! x component of force
              FD(1,I)=FD(1,I)-DX
              DY=P*DY                          ! y component of force
              FD(2,I)=FD(2,I)-DY
              DZ=P*DZ                          ! z component of force
              FD(3,I)=FD(3,I)-DZ
              
           End do inner_loop

!          Potential from EAM Embedding Function (with Rho as input)

           IF(RKFi.GE.1) Then
              Fi=EFT(ktypeI,RKFi)+(EFT(ktypeI,RKFi+1)-EFT(ktypeI,RKFi))*Rhoi
           ELSE 
              Fi=EFT(ktypeI,1)+(EFT(ktypeI,1)-EFT(ktypeI,2))*(1-Rhoi)
              IF (Rhoi.EQ.0.0) Fi=0.0d0
           ENDIF  ! This takes care for low el. density, e.g. for dimers

           POT(I)=POT(I)+Fi
           
        End do outer_loop
        
        RETURN
      END SUBROUTINE F_EAM
