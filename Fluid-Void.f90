!     For the analisys of the results we give the microscopic water variable
!     only for the writing the output data. Following the writing procedure,
!     all these variable must be nullified with Fluid-Void.f90 procedure
      SUBROUTINE Fluid_Void()
        INCLUDE 'common.h'
        INCLUDE 'commonTTM.h'

        DO ii = 1,ippl
           ip = idlc(ippl)
           RON(ip) = 0.0d0
           T_l(ip) = 0.0d0
           TPSX(ip)= 0.0d0
           Q1V(ip) = 0.0d0
        ENDDO
        
      END SUBROUTINE Fluid_Void
