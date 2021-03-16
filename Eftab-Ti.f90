      SUBROUTINE EFTAB_Ti()
        INCLUDE 'common.h' 
!       KE1, KE2 - elements from the list Element(N)
!       Units: XT in Angs, UT in eV, FT in eV/Angs, mass in amu

        INTEGER I
        REAL *8 XT, DX
        
        CHARACTER*2, DIMENSION(3), PARAMETER:: Element = (/'OX','Ti','Sr'/)
        
        REAL*8, DIMENSION(3), PARAMETER:: Zet=   &
             (/ -1.4d0,     2.2d0,    2.0d0    /)
        REAL*8, DIMENSION(3), PARAMETER:: anm=   &
             (/  1.85d0,    1.385d0,  1.846d0  /)
        REAL*8, DIMENSION(3), PARAMETER:: bnm=   &
             (/  0.18d0,    0.1d0,    0.13d0   /)
        REAL*8, DIMENSION(3), PARAMETER:: cpar=   &
             (/ 8.55372d0,  0.0d0,    5.90614d0/)
        REAL*8, DIMENSION(3), PARAMETER:: TXmass=  &
             (/  15.9994d0, 47.867d0, 87.62d0  /)
        REAL*8, PARAMETER:: dpar  =  1.123692581 ! ev/at
        
        REAL*8, PARAMETER:: beta  =  2.0d0

        REAL*8, PARAMETER:: rstar = 1.6d0

        REAL*8, PARAMETER:: fzero = 0.043364394d0  !N

        IF(mypid.eq.0) Write (17,*) 'Creating potencial table'

        
        alat=SQRT(2.0d0)*(anm(1)+anm(2)+ &
             anm(3))/3.0d0  ! FCC lattice parameter (used in Liquid.f90) [A]

        NIJ = 0

        DO I=1,NTYPE

           XMASS(I) = TXmass(I)

           DO J = I,NTYPE

              NIJ = NIJ + 1

              RM(NIJ) = 15.0d0  ! Cut-off distance in An

!W         Uncomment these lines if you want to plot the potential
              IF(mypid.eq.0) OPEN (UNIT = 111, FILE = &
                   './'//trim(Element(I))//'-'//trim(Element(J))//'.out')
              
              DX=RM(NIJ)/NT

              XT=0.0d0

              IJINDEX(I,J) = NIJ
              IJINDEX(J,I) = NIJ

              DO II=1,NT
                 XT=XT+DX

                 UT(NIJ,II) = Zet(I)*Zet(J)/XT + &
                      fzero*(bnm(I) + &
                      bnm(J))*EXP((anm(I)+anm(J)-XT)/(bnm(I)+bnm(J))) +&
                      cpar(I)*cpar(J)/XT**6.0d0
                      
                 FT(NIJ,II) = Zet(I)*Zet(J)/XT/XT + &
                      fzero*EXP((anm(I)+anm(J)-XT)/(bnm(I)+bnm(J))) + &
                      6.0d0*cpar(I)*cpar(J)/XT**7.0d0

                 IF(NIJ.EQ.2) THEN

                    UT(NIJ,II) = UT(NIJ,II) + &
                         dpar*(EXP(-2.0d0*beta*(XT-rstar)) - &
                         2.0d0*EXP(-beta*(XT-rstar)))


                    FT(NIJ,II) = FT(NIJ,II) + &
                         dpar*beta*(2.0d0*EXP(-2.0d0*beta*(XT-rstar)) - &
                         2.0d0*EXP(-beta*(XT-rstar)))

                 ENDIF

!W               Uncomment this line if want to plot the potential
                 IF(mypid.eq.0) WRITE (111,102) XT,UT(NIJ,II),FT(NIJ,II)
102              Format (3(1x,D10.4))

!                Transfer to program units
!                Electron density and its derivative is in eV/A, do not
!                convert it into pr.un.

                 UT(NIJ,II) = UT(NIJ,II)/ENUNIT
                 FT(NIJ,II) = FT(NIJ,II)/ENUNIT
              END DO

!W            Uncomment this line if want to plot the potential
              CLOSE (UNIT = 111)

           ENDDO
        ENDDO

        
        RETURN
      END SUBROUTINE EFTAB_TI
