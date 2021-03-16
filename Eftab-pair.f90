!     Making the Energy and Force tables
!     Leonid Zhigilei, 2000
!     Number of pair potentials defined here should be
!     Npots=Ntype*(Ntype+1)/2,  Ntype is number of particle types
!     Lennard - Jones is defined here so far.

      SUBROUTINE EF1LJ(KTEF,KE1,KE2)
        INCLUDE 'common.h' 
!       KTEF - type of the potential
!       KE1, KE2 - elements from the list Element(N)
!       LJeps is in eV, LJsig is in A, LJmass is in amu
!       [D.V.Matyushov and R. Schmid, J.Chem.Phys. 104, 8627 (1996)] - Ar
!       [T.Halicioglu and G.M.Pound, Phys.Stat.Sol.A 30, 619 (1975)] - metals
!       Cross-species parameters are fixed using the Lorentz-Berthelot rules
!       [see, e.g. J. S. Rowlinson, Liquid and liquid mixtures 
!       (Butterworth Scientific, London, 1982]

        REAL*8 LJCeps, LJCsig

        CHARACTER*3, DIMENSION(20), PARAMETER:: Element =   &
             (/'Ar ','Al ','Ca ','Au ','Pb ','Ni ','Pd ','Pt ','Ag ','Cu ', &
             'Cr ','Fe ','Li ','Mo ','W  ','Na ','K  ','Si ','Ge ','H2O'/)

        REAL*8, DIMENSION(20), PARAMETER:: LJeps=   &
             (/0.0103d0, 0.3922d0, 0.2152d0, 0.4415d0, 0.2364d0, 0.5197d0, &
             0.4267d0, 0.6817d0, 0.3448d0, 0.4094d0, 0.5020d0, 0.5264d0, &
             0.2054d0, 0.8379d0, 1.0678d0, 0.1379d0, 0.1144d0, 3.3900d0, &
             2.7400d0, 0.026d0/)!0.00589692d0/)!0.020d0/)!0.00606634d0/)

        REAL*8, DIMENSION(20), PARAMETER:: LJsig=   &
             (/3.405d0,  2.62d0,   3.60d0,   2.637d0,  3.197d0,  2.282d0, &
             2.52d0,   2.542d0,  2.644d0,  2.338d0,  2.336d0,  2.321d0, &
             2.839d0,  2.551d0,  2.562d0,  3.475d0,  4.285d0,  2.0936d0, &
             2.1827d0, 2.90d0/)!3.10456955d0/)!2.75d0/)!3.041066914d0/)!2.865672727d0/) !This one worked for matching the lattise but gave too high pressur
        
        REAL*8, DIMENSION(20), PARAMETER:: mass=   &
             (/40.0d0,  27.0d0,   40.1d0,   197.0d0,  207.2d0,  58.7d0, &
             106.4d0,   195.1d0,  107.9d0,  63.6d0,  52.0d0,  55.9d0, &
             6.94d0,  95.9d0,  183.85d0,  23.0d0,  39.1d0,  28.09d0, &
             72.59d0, 18.02d0/)

!        Write (17,*) 'Creating tables for ', Element(KE1),'-', &
!             Element(KE2),' potential' 

        IF(KE1.EQ.KE2) Then
           LJCeps=LJeps(KE1)
           LJCsig=LJsig(KE1)
           XMass(KTEF) = mass(KE1)   ! Mass of the particle in pr.unit [aem]
!           IF(KTEF.EQ.4) XMass(3) = mass(KE1)
           RM(KTEF) = 2.5*LJCsig     ! Cutoff distance for interactin potential [A]
        ELSE
           LJCeps=SQRT(LJeps(KE1)*LJeps(KE2))
           LJCsig=(LJsig(KE1)+LJsig(KE2))/2.0d0
           RM(KTEF) = 2.5*LJCsig 
        ENDIF

!       Uncomment this line if want to plot the potential
        OPEN (UNIT = 111,FILE='LJ-'//Element(KE1)//'-'//Element(KE2)//'.pot')
        sig6 = LJCsig**6
        DX=RM(KTEF)/NT
        XT=0.0d0
        DO I=1,NT
           XT=XT+DX
           XT6=XT**6
           sr6 = sig6/XT6                                    ! (sigma/r)**6
           UT(KTEF,I) = 4.0d0 * LJCeps * sr6 * (sr6 - 1.0d0)     ! potential
           FT(KTEF,I) = LJCeps * 48.0d0/XT * sr6 * (sr6 - 0.5d0) ! force
!          It is a common practice to tabulate force/r rather than force
!          But I am not doing this here to make code more transparent
!          Uncomment this line if want to plot the potential
           WRITE (111,*) XT,UT(KTEF,I),FT(KTEF,I)
!          Transfer to program units
           UT(KTEF,I) = UT(KTEF,I)/ENUNIT
           FT(KTEF,I) = FT(KTEF,I)/ENUNIT
        END DO
!W      Uncomment this line if want to plot the potential
        CLOSE (UNIT = 111)
        RETURN
      END SUBROUTINE EF1LJ











