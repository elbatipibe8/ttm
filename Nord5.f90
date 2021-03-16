!     Nordsieck Integrator (fifth-order predictor-corrector algorithm)
!     [Nordsieck A. Math.Comput.,v.16,p.22,1962]
!     [Melker A.I. et. al. A.F.Ioffe Phys.Tech.Inst. Preprint 661,1980]
!     [Allen & Tildesley, p.340]
!     [Lecture notes for MSE 524]
!     Kinetic energy is calculated here as well.
!     Leonid Zhigilei, 2000
      SUBROUTINE NORD5()           
        INCLUDE 'common.h'
        INCLUDE 'commonTTM.h'
        INCLUDE 'mpif.h'

        trans_loop: DO I=1,NN13 
           
           I3 = (I-0.5)/3.+1.
           M3 = MOD(I,3)
           IF(M3.EQ.0) M3=3

!          Rigid boundary
           IF(KBOUND.EQ.1.AND.KHIST(I3).EQ.3) CYCLE trans_loop

!          Nonreflective Boundary
           IF((KBOUND.EQ.2.OR.KBOUND.EQ.11).AND.KHIST(I3).GT.10.AND.M3.NE.3) CYCLE trans_loop

           P1=Q2(I)                                                  
           P2=Q3(I)               
           P3=P2+P2+P2
           P4=Q4(I)
           P5=P4+P4
           P6=P5+P5
           P7=Q5(I)
           P8=(P7+P7+P7)+(P7+P7)                 
           P9=P8+P8                            
           XD(M3,I3)=P7+P4+P2+P1+Q1D(M3,I3)+XD(M3,I3)       
           Q1D(M3,I3)=P8+P6+P3+P1+P1+Q1D(M3,I3)      
           Q2(I)=P9+P6+P5+P3+P1          
           Q3(I)=P9+P6+P2        
           Q4(I)=P8+P4
        END DO trans_loop

        IF(MOD(ISTEP,NLREN).EQ.0) Then
           CALL LINKLIST()
           CALL MIGRATE()
        ENDIF

        NN13 = NN1*3
        CALL SHARE()

        IF(MOD(ISTEP,NLREN).EQ.0) CALL NbList_MPI()

!       TTM and laser pulse
        IF(LGOT.EQ.2) Then
           CALL TTMTemper()
           CALL SHARE_TTM()
           CALL TTMDiffuse()
        ENDIF

!       To calculate forces for EAM need also information about density
!       in the skin, performed in F_EAM() and shared in Share-EAM.f90.
        CALL Forces()
        
        QIN(1:NN1)=0.0d0

        trans2_loop: DO I=1,NN13

           I3 = (I-0.5)/3.+1.
           M3 = MOD(I,3)
           IF(M3.EQ.0) M3 = 3

!          Rigid boundary
           IF(KBOUND.EQ.1.AND.KHIST(I3).EQ.3) CYCLE trans2_loop

!          Nonreflective Boundary
           IF((KBOUND.EQ.2.OR.KBOUND.EQ.11).AND.KHIST(I3).GT.10.AND.M3.NE.3) CYCLE trans2_loop

           ktypei = KTYPE(I3)
           GG=G1(ktypei)
           P=GG*FD(M3,I3)-Q2(I)                                
           XD(M3,I3)=XD(M3,I3)+CC1*P                               
           Q1D(M3,I3)=Q1D(M3,I3)+C2*P                           
           Q2(I)=Q2(I)+P                            
           Q3(I)=Q3(I)+C3*P                       
           Q4(I)=Q4(I)+C4*P 
           Q5(I)=Q5(I)+C5*P     
           QIN(I3)=QIN(I3)+Q1D(M3,I3)*Q1D(M3,I3)*0.25/GG
        END DO trans2_loop

        RETURN 
      END SUBROUTINE NORD5

