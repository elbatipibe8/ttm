!     Calculation of the total pressure (April 2002)
      SUBROUTINE Pressure(TPRESS)
        INCLUDE 'common.h'
        REAL*8, DIMENSION(:), ALLOCATABLE :: AtPress,VAtPress,TAtPress

        ALLOCATE(AtPress(LPMZ),VAtPress(LPMZ),TAtPress(LPMZ))

        SPress=0.0d0
        VPress=0.0d0
!       Calculating virial portion first
!       Min Zhou claims that only this part makes sense (Cauchy stress)
        outer_loop: DO I=1,NN1
!          Calculate Trace to obtain local atomic pressure
!          The sign is oposite (multiply by -1 if stresses are of interest)
           IF(KTYPE(I).NE.KWT) THEN
              AtPress(i)=(STEN(I,1,1)+STEN(I,2,2)+STEN(I,3,3))/3.0d0
              SPress=SPress+AtPress(I)
           ENDIF
        End do outer_loop

!       Adding the velocity portion
        VAtPress(1:NN1)=0.0d0
!       This is the kinetic portion of the atomic pressure [in program units]
        i_loop: DO i=1,NN1
           IF(KTYPE(I).NE.KWT) THEN
              m_loop: DO M=1,3
                 VAtPress(i)=VAtPress(i)+ &
                      XMASS(KTYPE(I))*Q1D(M,I)*Q1D(M,I)/DELTA/DELTA/3.0d0
              END DO m_loop
              VPress=VPress+VAtPress(I)
           ENDIF
        END DO i_loop

!       Convert Units [Pr.En.units/Ang^3] --> [J/m^3] == [Pa]
        SPress=SPress*ENUNIT*EVTOJOU*1.0d+30/(XLREAL*YLREAL*ZLREAL)
        VPress=VPress*ENUNIT*EVTOJOU*1.0d+30/(XLREAL*YLREAL*ZLREAL)
        TPress=SPress+VPress

        DEALLOCATE(AtPress,VAtPress,TAtPress)

        RETURN
      END SUBROUTINE Pressure
