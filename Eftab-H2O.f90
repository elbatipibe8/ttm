      SUBROUTINE EFTAB_H2O()
        INCLUDE 'common.h' 
!       Units: XT in Angs, UT in eV, FT in eV/Angs, mass in amu
        INTEGER I
        REAL*8 XT, DX
        CHARACTER*2, DIMENSION(2), PARAMETER:: Element = (/'H','O'/) 
   
        REAL*8, DIMENSION(2), PARAMETER:: anm=   &
             (/  0.0d0,  1388.7730d0   /)    ! eV
        REAL*8, DIMENSION(2), PARAMETER:: bnm=   &
             (/  0.0d0,     2.7690d0   /)    ! A-1 
        REAL*8, DIMENSION(2), PARAMETER:: cpar=   &  
             (/ 0.0d0,     175.000d0   /)    ! eV A+6
        REAL*8, DIMENSION(2), PARAMETER:: TXmass=  & 
             (/  28.085d0,  15.9994d0  /)     ! at. u.
        REAL*8, DIMENSION(2), PARAMETER:: Eps=  & 
             (/ 0.0d0,     0.101411d0  /)    ! kJ mol-1
        REAL*8, DIMENSION(2), PARAMETER:: Sig=  & 
             (/ 0.0d0,     1.779239d0  /)            ! A 

	fconst= 14.39964391d0            !A ev-1
        fzero = 4.1868d+03/1.60217646d-19/6.02214179d+23 ! ev/A
!       Parameter for bound. cond. in Ewald Sum, EPSR = 1 is for vacuum
        EPSR = 1.0d0
        ftest = 1.0d+03/1.60217646d-19/6.02214179d+23 ! ev 

	Zet = (/ 2.4d0 , -1.2d0, 0.0d0 /)         ! ev

        IF(mypid.eq.0) Write (17,*) 'Creating potencial table'

        alat = 7.389285d0 ! [A]
        NIJ = 0

        DO I=1,NTYPE
           XMASS(I) = TXmass(I)
           DO J = I,NTYPE
              NIJ = NIJ + 1

              IF(I.NE.J) THEN
                 EpsIJ = 0.2989065d0   ! kJ mol-1
                 SigIJ = 1.313635d0    ! A
                 aij   = 18003.7572d0  ! eV
                 bij   = 4.87318d0     ! A-1
                 cij   = 133.5381d0    ! eV A+6
              ELSE
                 EpsIJ = Eps(I)
                 SigIJ = Sig(I)
                 aij   = anm(I)
                 bij   = bnm(I)
                 cij   = cpar(I)
              ENDIF

              RM(NIJ) = 15.0d0 ! Cut-off in [A]

!W            Uncomment these lines if you want to plot the potential
              OPEN (UNIT = 111, FILE = &
                   './'//trim(Element(I))//'-'//trim(Element(J))//'.out')
              
              DX=RM(NIJ)/NT
              XT = 0.0d0

              IJINDEX(I,J) = NIJ
              IJINDEX(J,I) = NIJ

              DO II=1,NT
                 XT=XT+DX
                 ERFC(II) = ERFCC(falfa*XT)
                 UT(NIJ,II) = aij*EXP(-XT*bij) - cij/XT**6 + &
                      ftest*4.0d0*EpsIJ*((SigIJ/XT)**30.0d0-(SigIJ/XT)**6.0d0)
                      
                 UTK(NIJ,II) = fconst*Zet(I)*Zet(J)*ERFC(II)/XT
                 
                 FT(NIJ,II)=aij*bij*EXP(-XT*bij)- 6.0d0*cij/XT**7.0d0 + &
              ftest*24.0d0*EpsIJ*(5.0d0*(SigIJ/XT)**30.0d0-(SigIJ/XT)**6.0d0)/XT

                 FTK(NIJ,II) = fconst*Zet(I)*Zet(J)*(ERFC(II) + &
                      2.0d0*falfa*EXP(-falfa*falfa*XT*XT)*XT/SQRT(PI))/XT/XT

!W               Uncomment this line if want to plot the potential
                 WRITE (111,102) XT,UT(NIJ,II)+UTK(NIJ,II),UT(NIJ,II), &
		    UTK(NIJ,II),FT(NIJ,II)+FTK(NIJ,II),FT(NIJ,II),FTK(NIJ,II)
102              Format (7(1x,E14.7))

                 UT(NIJ,II) = UT(NIJ,II) + UTK(NIJ,II)
                 FT(NIJ,II) = FT(NIJ,II) + FTK(NIJ,II)

!                Transfer to program units
                 UT(NIJ,II) = UT(NIJ,II)/ENUNIT
                 FT(NIJ,II) = FT(NIJ,II)/ENUNIT
              END DO

!W            Uncomment this line if want to plot the potential
              CLOSE (UNIT = 111)

           ENDDO
        ENDDO
        
!	CALL Set_Ewald()

        RETURN
      END SUBROUTINE EFTAB_GLASS
