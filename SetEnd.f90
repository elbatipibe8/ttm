!     THIS SUBROUTINE CLEANS UP AFTER THE MD LOOP
!     Leonid Zhigilei, 2000
      SUBROUTINE SetEnd()
        INCLUDE 'common.h' 

!       Since we multiplied velocities by timestep DELTA
!       at the beginning of the run (in SetInit.f90) we have
!       to divide by DELTA here to get real velocities
        DO I=1,NN13
           Q1(I)=Q1(I)/DELTA
        ENDDO

!!!        DO I=1,NN1
!!!           R1(I)=R1(I)/DELTA
!!!        ENDDO

        RETURN
      END SUBROUTINE SetEnd
