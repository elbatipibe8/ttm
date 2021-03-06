!     Making the Energy and Force tables
!     for Jonhson's new EAM potential
!     Dewey Murdick and Leonid Zhigilei, 2000, 2002
!     Number of pair potential defined here should be 
!     Npots=Ntype*(Ntype+1)/2,  Ntype is number of particle types

      SUBROUTINE EF1_EAM(KTEF,KE1,KE2)
        INCLUDE 'common.h'
!       KTEF - type of the potential
!       KE1, KE2 - elements from the list Element(N)
!       Units: XT in Angs, UT in eV, FT in eV/Angs, mass in amu

        INTEGER I
	
	! the next parameters are taken from "Molecular dynamics simulation of femtosecond ablation and spallation with different interatomic potentials" by V.V. Zhakhovskiia
!	real*8, parameter :: aa1 = 2.7471581015136728d0
!        real*8, parameter :: aa2 = 5.3593750000000000d0
!        real*8, parameter :: aa3 = 3.2500000000000000d0

!        real*8, parameter :: bb1 = 8.2311259601633768d0
!        real*8, parameter :: bb2 =-382.38931538388255d0
!        real*8, parameter :: bb3 = 16.250071667347235d0
!        real*8, parameter :: bb4 = 1.4586663896542300d0

!        real*8, parameter :: cc1 = 3.0697898737897571d0
!        real*8, parameter :: cc2 = 20.750105835621902d0

        real*8 :: rc = 6.875   ! cutoff in angstroem
	real*8 xx, xc, rhomax
	! /these parameters are taken from "Molecular dynamics simulation of femtosecond ablation and spallation with different interatomic potentials" by V.V. Zhakhovskiia

        REAL *8 XT, DX, Rho, DRho, Rhoin, Rhoout
        REAL *8 re, fe, Rhoe, Alpha, Beta, A, B, Kappa, Lambda
        REAL *8 Fi0,Fi1,Fi2,Fi3,Fm0,Fm1,Fm2,Fm3,niu,Fn
        
        CHARACTER*2, DIMENSION(8), PARAMETER:: Element =   &
             (/'Cu','Au','Ni','Fe','Zr','Al','ST','HO'/)
        
!        REAL*8, PARAMETER:: m=20.0d0, n=20.0d0
!        REAL*8, DIMENSION(5), PARAMETER:: Tre=   &
!             (/ 2.556162d0,  2.885034d0,  2.488746d0,  2.481987d0,  3.199978d0/)
!        REAL*8, DIMENSION(5), PARAMETER:: Tfe=   &
!             (/ 1.554485d0,  1.529021d0,  2.007018d0,  1.885957d0,  2.230909d0/)
!        REAL*8, DIMENSION(5), PARAMETER:: TRhoe=   &
!             (/22.150141d0, 21.319637d0, 27.984706d0, 20.041463d0, 30.879991d0/)
!        REAL*8, DIMENSION(5), PARAMETER:: TAlpha=   &
!             (/ 7.669911d0,  8.086176d0,  8.029633d0,  9.818270d0,  8.559190d0/)
!        REAL*8, DIMENSION(5), PARAMETER:: TBeta=   &
!             (/ 4.090619d0,  4.312627d0,  4.282471d0,  5.236411d0,  4.564902d0/)
!        REAL*8, DIMENSION(5), PARAMETER:: TA=   &
!             (/ 0.327584d0,  0.230728d0,  0.439664d0,  0.392811d0,  0.424667d0/)
!        REAL*8, DIMENSION(5), PARAMETER:: TB=   &
!             (/ 0.468735d0,  0.336695d0,  0.632771d0,  0.646243d0,  0.640054d0/)
!        REAL*8, DIMENSION(5), PARAMETER:: TKappa=   &
!             (/ 0.431307d0,  0.420755d0,  0.413436d0,  0.170306d0,  0.500000d0/)
!        REAL*8, DIMENSION(5), PARAMETER:: TLambda=   &
!             (/ 0.862140d0,  0.841511d0,  0.826873d0,  0.340613d0,  1.000000d0/)
!        REAL*8, DIMENSION(5), PARAMETER:: TFi0=   &
!             (/-2.176490d0, -2.930281d0, -2.693996d0, -2.534992d0, -4.485793d0/)
!        REAL*8, DIMENSION(5), PARAMETER:: TFi1=   &
!             (/-0.140035d0, -0.554034d0, -0.066073d0, -0.059605d0, -0.293129d0/)
!        REAL*8, DIMENSION(5), PARAMETER:: TFi2=   &
!\             (/ 0.285621d0,  1.489437d0,  0.170482d0,  0.193065d0,  0.990148d0/)
!        REAL*8, DIMENSION(5), PARAMETER:: TFi3=   &
!             (/-1.750834d0, -0.886809d0, -2.457442d0, -2.282322d0, -3.202516d0/)
!        REAL*8, DIMENSION(5), PARAMETER:: TFm0=   &
!             (/-2.190000d0, -2.980000d0, -2.700000d0, -2.540000d0, -4.510000d0/)
!        REAL*8, DIMENSION(5), PARAMETER:: TFm1=   &
!             (/ 0.000000d0,  0.000000d0,  0.000000d0,  0.000000d0,  0.000000d0/)
!        REAL*8, DIMENSION(5), PARAMETER:: TFm2=   &
!             (/ 0.702991d0,  2.283863d0,  0.282257d0,  0.200269d0,  0.928602d0/)
!        REAL*8, DIMENSION(5), PARAMETER:: TFm3=   &
!             (/ 0.683705d0,  0.494127d0,  0.102879d0, -0.148770d0, -0.981870d0/)
!        REAL*8, DIMENSION(5), PARAMETER:: Tniu=   &
!             (/ 0.921150d0,  1.286960d0,  0.509860d0,  0.391750d0,  0.597133d0/)
!        REAL*8, DIMENSION(5), PARAMETER:: TFn=   &
!             (/-2.191675d0, -2.981365d0, -2.700493d0, -2.539945d0, -4.509025d0/)
        
        REAL*8, DIMENSION(8), PARAMETER:: Tmass=   &
             (/63.6d0, 196.966569d0,  58.6934d0,  55.9d0, 91.2d0, &
             26.981539d0, 183.49d0, 18.02d0/)
        
!        IF(mypid.eq.0) write (17,*) 'Creating tables for ', Element(KE1),'-', &
!             Element(KE2),' EAM potential' 
        
        IF(KE1.EQ.KE2) Then   
!       going from the type of the potential to the type of atom
           IF(KTEF.EQ.1) KF=1  
           IF(KTEF.EQ.2) KF=2
           IF(KTEF.EQ.4) KF=3
           IF(KTEF.EQ.3.OR.KTEF.GT.4) Then
              If(mypid.eq.0) Print *,"The values of KTEF=",KTEF,", KE1=",KE1,", and KE2=",KE2, &
                   " are uncompatible. Program will stop."
              STOP
           ENDIF
!           re=Tre(KE1)
!           fe=Tfe(KE1)
!           Rhoe=TRhoe(KE1)
!           Alpha=TAlpha(KE1)
!           Beta=TBeta(KE1)
!           A=TA(KE1)
!           B=TB(KE1)
!           Kappa=TKappa(KE1)
!           Lambda=TLambda(KE1)
!           Fi0=TFi0(KE1)
!           Fi1=TFi1(KE1)
!           Fi2=TFi2(KE1)
!           Fi3=TFi3(KE1)
!           Fm0=TFm0(KE1)
!           Fm1=TFm1(KE1)
!           Fm2=TFm2(KE1)
!           Fm3=TFm3(KE1)
!           niu=Tniu(KE1)
!           Fn=TFn(KE1)
!           Rhoin=0.85d0*Rhoe
!           Rhoout=1.15d0*Rhoe
           XMass(KF)=Tmass(KE1)      ! Mass of the particle in pr.unit [aem]
           RM(KTEF) = rc		       ! Cutoff dist. for 2-body potential [A]
           alat = 4.056433d0
        ELSE
           If(mypid.eq.0)Print *,"Multicomponent EAM is not implemented yet. Program will stop."
           STOP
        ENDIF
        
        RhoM(KF) = 8.0d0

!        XT = 0.d0
        DX = RM(KTEF)/NT
!        rcc = 0.1d0*rc    ! changing units for easier calculation of energy: angstroems -> nm

!        DO I=1,NT
!           XT=XT+DX
!           XTT=XT*0.1d0  ! changing units for easier calculation of energy: an!gstroems -> nm
!           DFT(KF,I)  =  cc1*(XTT*XTT-rcc*rcc)**2/(1.d0+(cc2*XTT*XTT)**3)
!           if (DFT(KF,I).GT.RhoM(KF)) RhoM(KF) = DFT(KF,I)
!        enddo

!        RhoM(KF) = RhoM(KF)*4.5d0 ! Empirically found form summation of densities


!W      Uncomment these lines if you want to plot the potential
!        IF(mypid.eq.0) OPEN (UNIT = 111,FILE='EAM.out')

        DRho=RhoM(KF)/NT
!        XT=0.0d0
!        Rho=0.0d0
!	xc = aa1*rcc*rcc

        OPEN(UNIT = 111,FILE = 'allum.data')

        READ(111,*)(XT,UT(KTEF,I),FT(KTEF,I),DFT(KF,I),DDFT(KF,I), &
             Rho,EFT(KF,I),DEFT(KF,I),I=1,NT)
        
!        DO I=1,NT
!           XT=XT+DX
!            XTT=XT*0.1d0  ! changing units for easier calculation of energy: angstroems -> nm
!           Rho=Rho+DRho
!	   xx = aa1 *XTT*XTT
	   
!!          Pair Potential
!           UT(KTEF,I) = (  (1.d0/xx - aa2) * (xx - xc)**10 * ((xx - xc)**18 + aa3)  ) 
!           FT(KTEF,I) = -( -1.d0/(xx*xx) * (xx - xc)**10 * ((xx - xc)**18 + aa3) + &
!                (1.d0/xx - aa2) * 10.d0*(xx - xc)**9 * ((xx - xc)**18 + aa3) + &
!                (1.d0/xx - aa2) * (xx - xc)**10 * 18.d0*(xx - xc)**17  ) &
!                *2.d0* XTT*aa1
		    
!!          Electron density function
!           DFT(KF,I)  =  cc1*(XTT*XTT-rcc*rcc)**2/(1.d0+(cc2*XTT*XTT)**3)
!           DDFT(KF,I) = cc1*2.d0*(XTT*XTT-rcc*rcc)*2.d0*XTT/(1.d0+(cc2*XTT*XTT)**3) &
!			   -  cc1*(XTT*XTT-rcc*rcc)**2/(1.d0+(cc2*XTT*XTT)**3)**2 &
!                           *3.d0*(cc2*XTT*XTT)**2 * cc2*2.d0*XTT

!!          Embedding energy function
!           EFT(KF,I)  =  (bb1*Rho*(bb2+(bb3+Rho)**2)/(1.d0+bb4*Rho))
!           DEFT(KF,I) = (bb1       *(bb2+(bb3+Rho)**2)/(1.d0+bb4*Rho)  &
!                + bb1*Rho*2.d0*(bb3+Rho)/(1.d0+bb4*Rho)  &
!                - bb1*Rho*(bb2+(bb3+Rho)**2)/(1.d0+bb4*Rho)**2 * bb4)
    

!W         Uncomment this line if want to plot the potential
!           if (mypid.eq.0)WRITE (111,102) XT,UT(KTEF,I),FT(KTEF,I),EFT(KF,I),DEFT(KF,I),&
!                Rho,DFT(KF,I),DDFT(KF,I)
!102        Format (8(1x,E14.7))

!          Transfer to program units
!          Electron density and its derivative is in eV/A, do not convert to pr.un.
        DO I = 1,NT
           UT(KTEF,I) = UT(KTEF,I)/ENUNIT!* &
!                (0.010382059800664451827242524916944d0)  ! KJ/mol to eV/atom
           FT(KTEF,I) = FT(KTEF,I)/ENUNIT!* &
!                (0.0010382059800664451827242524916944d0) !  kJ/(mol*nm) to eV/(atom*angstroem)
           EFT(KF,I) = EFT(KF,I)/ENUNIT!* &
!                (0.010382059800664451827242524916944d0)  ! KJ/mol to eV/atom
           DEFT(KF,I) = DEFT(KF,I)/ENUNIT!* &
!                (0.0010382059800664451827242524916944d0) !  kJ/(mol*nm) to eV/(atom*angstroem)
        END DO
	
!W      Uncomment this line if want to plot the potential
        CLOSE (UNIT = 111)
        RETURN
      END SUBROUTINE EF1_EAM
