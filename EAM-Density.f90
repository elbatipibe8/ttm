!     Local Electron Density calculation (inefficient)
!     Electron density is in eV/A, do not convert to pr. units
      SUBROUTINE EAM_Density()
        INCLUDE 'common.h'
        Real *8  RD,RR,P,DZ,DY,DX
        Integer  IJ,KF,KTYPEI,KTYPEJ,I,K,J
        
        RhoDen(1:NN1)=0.0d0
        
        outer_loop: DO I=1,NN1

!           IF(KTYPE(I).EQ.KWT) CYCLE outer_loop
           
           IF (NNG(I).gt.MAXNNB) THEN
              print *, ' TOO MANY NEIGHBOURS: ',NNG(I),' > ', maxnnb
!              call Flushout(6)
              STOP
           ENDIF
           ktypei=KTYPE(I)
           inner_loop: DO K=1,NNG(I)
              J=NNNG(k,i)

!              IF(KTYPE(J).EQ.KWT) CYCLE inner_loop
              
              ktypeJ=KTYPE(J)
              DX=XD(1,J)-XD(1,I)
              DY=XD(2,J)-XD(2,I)
              DZ=XD(3,J)-XD(3,I)
              
              IJ = ijindex(KTYPEI, KTYPEJ)
              RR=DX*DX+DY*DY+DZ*DZ
              RD=DSQRT(RR)
              IF (RD.GE.RM(IJ)) CYCLE inner_loop
              
              IF (RD.LE.RM(IJ)/NT) THEN
                 Print *, 'Atoms ',I,' and ',J,' are too close.', XD(1:3,I),XD(1:3,J),&
                      'Program will stop.',xnmin,xnmax,ynmin,ynmax,znmin,znmax,"borders",NN1,mypid
!                 call Flushout(6)
                 STOP
              ENDIF
              
              P=RD*DXR(IJ)
              KF=P
              P=P-KF
              
!             Local Electron Density for atom I 
              RhoDen(I)=RhoDen(I)+DFT(ktypeJ,KF)+ &
                   (DFT(ktypeJ,KF+1)-DFT(ktypeJ,KF))*P

           End do inner_loop

!          Check and see if our Rho is too large
           IF(RhoDen(I).gt.RhoM(ktypeI)) THEN
              Print *, RhoDen(I), '>', RhoM(ktypeI)
              Print *, 'Electron Density is too high for molecule ',I,'.', &
                   ' Program will stop.',NNG(I),'neighbores',mypid, &
                   'process',RD,'distans in Ang'
 !             call Flushout(6)
              STOP
           END IF

        End do outer_loop
        
        RETURN
      END SUBROUTINE EAM_Density

