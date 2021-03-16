      SUBROUTINE F_pair()
        INCLUDE 'common.h'
        
        outer_loop: DO I=1,NN1
           IF (NNG(I).gt.MAXNNB) THEN
              print *, ' TOO MANY NEIGHBOURS: ',NNG(I),' > ', maxnnb
!              call Flushout(6)
              STOP
           ENDIF
           ktypei=KTYPE(I)
           inner_loop: DO K=1,NNG(I)
              J=NNNG(k,i)
              
              IF(KTYPE(I).EQ.1.AND.KTYPE(J).EQ.1.AND.KEYBS.EQ.4.) CYCLE inner_loop
              
              DX=XD(1,J)-XD(1,I)
              DY=XD(2,J)-XD(2,I)
              DZ=XD(3,J)-XD(3,I)

              IJ = ijindex(KTYPEI,KTYPE(J))
              RR=DX*DX+DY*DY+DZ*DZ
              RD=SQRT(RR)
              IF (RD.GE.RM(IJ)) CYCLE inner_loop

!             If your run crashes, you can try to uncomment this part
              IF (RD.LE.RM(IJ)/NT) THEN
                 Print *, 'Molecules ',I,' and ',J,' are too close.', &
                      'Program will stop.'
!                 call Flushout(6)
                 STOP
              ENDIF

              P=RD*DXR(IJ)
              KF=P
              P=P-KF
              
!             Potential energy (give half the potential to each particle)
              PENER=.5*(UT(IJ,KF)+(UT(IJ,KF+1)-UT(IJ,KF))*P)
              POT(I)=POT(I)+PENER
!              POT(J)=POT(J)+PENER
!             End of potential energy

              P=FT(IJ,KF)+(FT(IJ,KF+1)-FT(IJ,KF))*P
              P=P/RD

!             Static Portion of Stress Tensor
              STEN(I,1,1)=STEN(I,1,1)+P*DX*DX
              STEN(I,2,2)=STEN(I,2,2)+P*DY*DY
              STEN(I,3,3)=STEN(I,3,3)+P*DZ*DZ
              
!              STEN(J,1,1)=STEN(J,1,1)+P*DX*DX
!              STEN(J,2,2)=STEN(J,2,2)+P*DY*DY
!              STEN(J,3,3)=STEN(J,3,3)+P*DZ*DZ
              
              DX=P*DX                          ! x component of force
              FD(1,I)=FD(1,I)-DX
!              FD(1,J)=FD(1,J)+DX
              DZ=P*DZ                          ! z component of force
              FD(3,I)=FD(3,I)-DZ
!              FD(3,J)=FD(3,J)+DZ
              DY=P*DY                          ! y component of force
              FD(2,I)=FD(2,I)-DY
!              FD(2,J)=FD(2,J)+DY
           End do inner_loop
        End do outer_loop
        
        RETURN
      END SUBROUTINE F_pair
