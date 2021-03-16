
!12.02.2019
!
!CDF (complex dielectric function) could be taken from
     !mod_CDF=1 (paper Vial PRB 2005 Improved analytical fit of gold dispersion, Ep=[1.24–2.48 eV]),
     !mod_CDF=2 (paper H.S. Sehmi PRB 2017,fit for the paper S. Babar and J. H. Weaver, 
     !Appl. Opt. 54, 477 (2015) with 4 Lorentz terms, Ep=[0.1–6.0 eV] 
     
     !mod_CDF=3 (paper H.S. Sehmi PRB 2017 with a fit to the paper P. B. Johnson and R. W. Christy, 
     !Phys. Rev. B 6, 4370 (1972) with 3 Lorentz terms and epsilon_infinity <> 1),
     !Ep=[0.64–6.6 eV]
     
     !mod_CDF=4 (for Au, paper R.L. Olmon PRB 2012,Ep=[0.04974–4.133 eV], it means lambda=[300 - 24930 nm]

!*******************************************************************************************************************
!Function Source is a fully analytical description of a TTM source with SPP at normal incidence

      !Source term at normal incidence is described by a Gaussian profile for temporal and space distribution.
      !Poynting vector + SPP      
      !"-div(S)" 
      subroutine Source_SPP(x,y,z,time,beta,delta,rp,source_res_spp_div) 
      use Module_mod
      use Module_coeff
      
      implicit none
      real*8 x,y,z,time,source_res_spp_div
      real*8 Part0,Part1,Part2,Part3,Part4
      real*8 tetta,sigma,beta,delta,rp
      real*8 f6
      
	  real*8, parameter :: pi = 2.d0*dacos(0.d0)
      
      !*********************************************************************
      
      sigma = 4.0d0*dlog(2.0d0)
      !t_laser = 200.0d-15  !sec
      !tstart=3.0d0*t_laser !sec
      !lambda_las=800.00d-09 !m
      
      Part0=dsqrt(sigma/pi)*dexp(-sigma*(time-tstart)*(time-tstart)/(t_laser*t_laser))/t_laser
      Part1=dexp(-(x*x+y*y)/(rp*rp))
      Part2=koeff11*dexp(2.0d0*koeff14*z)
      
      if (x.ge.(0.0d0)) then
      
        f6=koeff5*x+koeff13*z+koeff15 
        Part3=beta*(2.0d0/koeff10)*dexp(koeff14*z)*dexp(-x*koeff6+z*koeff3)*(f13*dcos(f6)+f14*dsin(f6))
        Part4=2.00d0*beta*beta*((koeff6*f1+koeff3*f9)/(koeff12))*dexp(2.0d0*(-x*koeff6+z*koeff3)) 
        
      elseif (x.lt.(0.0d0)) then
      
        f6=-koeff5*x+koeff13*z+koeff15 
        Part3=beta*(2.0d0/koeff10)*dexp(koeff14*z)*dexp(x*koeff6+z*koeff3)*(f13*dcos(f6)+f14*dsin(f6))
        Part4=2.00d0*beta*beta*((koeff6*f1+koeff3*f9)/(koeff12))*dexp(2.0d0*(x*koeff6+z*koeff3)) 

      endif  
      
       source_res_spp_div=Fluence*Part0*Part1*(Part2-Part3+Part4)
      
      return 
      end
!*******************************************************************************************************************

!***********************************************************************************************************
!Function CDF of Au as a function of laser frequency.
      
	  complex*16 function CDF_m(omega2)
	    use Module_mod
	    use Module_CDF
	    
	    implicit none
        real*8 omega2,eps,omega_d,gamma_d
        real*8 sigma_d,omega_11,omega_12,omega_21,omega_22,omega_31,omega_32,omega_41,omega_42
        real*8 sigma_11,sigma_12,sigma_21,sigma_22,sigma_31,sigma_32,sigma_41,sigma_42
        real*8 e1_res,e2_res
        complex*16 term0,term1,term2,term3,term4
        
        
        complex*16, parameter :: ci=(0.0d0,1.0d0),c1=(1.d0,0.0d0)
        real*8, parameter :: pi = 2.d0*dacos(0.d0)
        real*8, parameter :: hcross=6.582119514d-16 !eV*sec; Planck constant over 2 pi
        real*8, parameter :: c=299792458.0d0 !m/s


            if (mod_CDF.eq.1) then
                eps = 9.0685
                gamma_d = 18.36d0*2.0d0*pi*1.0d12 ! Hz=1/sec
                omega_d = 2155.6d0*2.0d0*pi*1.0d12 ! Hz=1/sec
                
                CDF_m=c1*eps-((c1*omega_d*omega_d)/(c1*omega2*omega2+ci*omega2*gamma_d))
                
            elseif (mod_CDF.eq.2) then
                eps = 0.83409
                gamma_d = 0.02334d0/hcross !1/sec
                sigma_d = 3134.5d0/hcross !1/sec
                omega_d = dsqrt(gamma_d*sigma_d) !1/sec
                
                omega_11 = 2.6905d0/hcross !1/sec
                omega_12 = -0.16645d0/hcross !1/sec
                sigma_11 = -0.01743d0/hcross !1/sec
                sigma_12 = 0.3059d0/hcross !1/sec
                
                omega_21 = 2.8772d0/hcross !1/sec
                omega_22 = -0.44473d0/hcross !1/sec
                sigma_21 = 1.0349d0/hcross !1/sec
                sigma_22 = 1.2919d0/hcross !1/sec
                
                omega_31 = 3.7911d0/hcross !1/sec
                omega_32 = -0.81981d0/hcross !1/sec
                sigma_31 = 1.2274d0/hcross !1/sec
                sigma_32 = 2.5605d0/hcross !1/sec
                
                omega_41 = 4.8532d0/hcross !1/sec
                omega_42 = -13.891d0/hcross !1/sec
                sigma_41 = 9.85d0/hcross !1/sec
                sigma_42 = 37.614d0/hcross !1/sec
                
                term0 = ((c1*gamma_d*sigma_d)/(c1*omega2*omega2+ci*omega2*gamma_d))
                term1 = ((ci*(c1*sigma_11+ci*sigma_12))/(c1*omega2-(c1*omega_11+ci*omega_12))) + &
                        ((ci*(c1*sigma_11-ci*sigma_12))/(c1*omega2+(c1*omega_11-ci*omega_12)))
                term2 = ((ci*(c1*sigma_21+ci*sigma_22))/(c1*omega2-(c1*omega_21+ci*omega_22))) + &
                        ((ci*(c1*sigma_21-ci*sigma_22))/(c1*omega2+(c1*omega_21-ci*omega_22)))
                term3 = ((ci*(c1*sigma_31+ci*sigma_32))/(c1*omega2-(c1*omega_31+ci*omega_32))) + &
                        ((ci*(c1*sigma_31-ci*sigma_32))/(c1*omega2+(c1*omega_31-ci*omega_32)))
                term4 = ((ci*(c1*sigma_41+ci*sigma_42))/(c1*omega2-(c1*omega_41+ci*omega_42))) + &
                        ((ci*(c1*sigma_41-ci*sigma_42))/(c1*omega2+(c1*omega_41-ci*omega_42)))
                
                CDF_m=c1*eps-term0+term1+term2+term3+term4
                
            elseif (mod_CDF.eq.3) then
                eps = -10.534d0
                gamma_d = 0.07373d0/hcross !1/sec
                sigma_d = 997.41d0/hcross !1/sec
                omega_d = dsqrt(gamma_d*sigma_d) !1/sec
                
                omega_11 = 2.5997d0/hcross !1/sec
                omega_12 = -0.43057d0/hcross !1/sec
                sigma_11 = 1.4835d0/hcross !1/sec
                sigma_12 = 0.88382d0/hcross !1/sec
                
                omega_21 = 3.7429d0/hcross !1/sec
                omega_22 = -1.2267d0/hcross !1/sec
                sigma_21 = 1.1372d0/hcross !1/sec
                sigma_22 = 3.8223d0/hcross !1/sec
                
                omega_31 = 7.3145d0/hcross !1/sec
                omega_32 = -21.843d0/hcross !1/sec
                sigma_31 = 225.27d0/hcross !1/sec
                sigma_32 = -193.27d0/hcross !1/sec
                
                term0 = ((c1*gamma_d*sigma_d)/(c1*omega2*omega2+ci*omega2*gamma_d))
                term1 = ((ci*(c1*sigma_11+ci*sigma_12))/(c1*omega2-(c1*omega_11+ci*omega_12))) + &
                        ((ci*(c1*sigma_11-ci*sigma_12))/(c1*omega2+(c1*omega_11-ci*omega_12)))
                term2 = ((ci*(c1*sigma_21+ci*sigma_22))/(c1*omega2-(c1*omega_21+ci*omega_22))) + &
                        ((ci*(c1*sigma_21-ci*sigma_22))/(c1*omega2+(c1*omega_21-ci*omega_22)))
                term3 = ((ci*(c1*sigma_31+ci*sigma_32))/(c1*omega2-(c1*omega_31+ci*omega_32))) + &
                        ((ci*(c1*sigma_31-ci*sigma_32))/(c1*omega2+(c1*omega_31-ci*omega_32)))
                
                CDF_m=c1*eps-term0+term1+term2+term3
            
            elseif (mod_CDF.eq.4) then
            
                if (omega2.ge.omega_CDF_tab(1).and.omega2.le.omega_CDF_tab(N_CDF)) then
                    call spline(omega2,omega_CDF_tab,e1_tab,b_e1_koef,c_e1_koef,d_e1_koef,N_CDF,e1_res)
                    call spline(omega2,omega_CDF_tab,e2_tab,b_e2_koef,c_e2_koef,d_e2_koef,N_CDF,e2_res)
                
                    CDF_m=c1*e1_res+ci*e2_res
                else
                   write(*,*)"You are out of the experimental data range"
                   stop
                endif  
            endif

          
	  end function CDF_m
!***********************************************************************************************************

!***********************************************************************************************************
!Function n1=Re(n) is a real part of a complex refractive index for Au as a function of laser frequency
      
	  real*8 function n1(omega3)
	    use Module_mod
	    
	    implicit none
        real*8 omega3,ReCDF,ImCDF
        
        complex*16 CDF_m
        external CDF_m
        
        ReCDF=dreal(CDF_m(omega3))
        ImCDF=dimag(CDF_m(omega3))
        
        n1=dsqrt((ReCDF+dsqrt(ReCDF*ReCDF+ImCDF*ImCDF))/(2.0d0))
 
	  end function n1
!***********************************************************************************************************

!***********************************************************************************************************
!Function k1=Im(n) is an imaginary part of a complex refractive index for Au as a function of laser frequency
      
	  real*8 function k1(omega4)
	    use Module_mod
	    
	    implicit none
        real*8 omega4,ReCDF,ImCDF
        
        complex*16 CDF_m
        external CDF_m
        
        ReCDF=dreal(CDF_m(omega4))
        ImCDF=dimag(CDF_m(omega4))
        
        k1=dsqrt((-ReCDF+dsqrt(ReCDF*ReCDF+ImCDF*ImCDF))/(2.0d0))
 
	  end function k1
!***********************************************************************************************************

!***********************************************************************************************************
!Function kx of the SPP for Au as a function of laser frequency
      
	  complex*16 function kx(omega3)
	    use Module_mod
	    
	    implicit none
        real*8 omega3
        
        complex*16 CDF_m
        external CDF_m
        
        complex*16, parameter :: ci=(0.0d0,1.0d0),c1=(1.d0,0.0d0)
        real*8, parameter :: c = 299792458.0d0 !m/sec
        
        !CDF_d=1 (air)
        kx=c1*(omega3/c)*cdsqrt((CDF_m(omega3))/(1.0d0+CDF_m(omega3)))
 
	  end function kx
!***********************************************************************************************************

!***********************************************************************************************************
!Function kz_m of the SPP for Au as a function of laser frequency
      
	  complex*16 function kz_m(omega4)
	    use Module_mod
	    
	    implicit none
        real*8 omega4
        
        complex*16 CDF_m
        external CDF_m
        
        complex*16, parameter :: ci=(0.0d0,1.0d0),c1=(1.d0,0.0d0)
        real*8, parameter :: c = 299792458.0d0 !m/sec       
        
        !CDF_d=1 (air)
        kz_m=c1*(omega4/c)*cdsqrt((-CDF_m(omega4)*CDF_m(omega4))/(1.0d0+CDF_m(omega4)))
 
	  end function kz_m
!***********************************************************************************************************

!***********************************************************************************************************
subroutine spline_koef(x,y,b,c,d,n)
!======================================================================
!  Calculate the coefficients b(i), c(i), and d(i), i=1,2,...,n
!  for cubic spline interpolation
!  s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
!  for  x(i) <= x <= x(i+1)
!  Alex G: January 2010
!----------------------------------------------------------------------
!  input..
!  x = the arrays of data abscissas (in strictly increasing order)
!  y = the arrays of data ordinates
!  n = size of the arrays xi() and yi() (n>=2)
!  output..
!  b, c, d  = arrays of spline coefficients
!  comments ...
!  spline.f90 program is based on fortran version of program spline.f
!  the accompanying function fspline can be used for interpolation
!======================================================================
    implicit none
    integer n
    real*8,dimension(1:n)::x,y,b,c,d
    integer i,j,gap
    real*8 h

    gap = n-1
    ! check input
    if ( n < 2 ) return
    if ( n < 3 ) then
      b(1) = (y(2)-y(1))/(x(2)-x(1))   ! linear interpolation
      c(1) = 0.0d0
      d(1) = 0.0d0
      b(2) = b(1)
      c(2) = 0.0d0
      d(2) = 0.0d0
      return
    end if
    !
    ! step 1: preparation
    !
    d(1) = x(2) - x(1)
    c(2) = (y(2) - y(1))/d(1)
    do i = 2, gap
      d(i) = x(i+1) - x(i)
      b(i) = 2.0d0*(d(i-1) + d(i))
      c(i+1) = (y(i+1) - y(i))/d(i)
      c(i) = c(i+1) - c(i)
    end do
    !
    ! step 2: end conditions 
    !
    b(1) = -d(1)
    b(n) = -d(n-1)
    c(1) = 0.0d0
    c(n) = 0.0d0
    if(n /= 3) then
      c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))
      c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))
      c(1) = c(1)*d(1)**2/(x(4)-x(1))
      c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))
    end if
    !
    ! step 3: forward elimination 
    !
    do i = 2,n
      h = d(i-1)/b(i-1)
      b(i) = b(i) - h*d(i-1)
      c(i) = c(i) - h*c(i-1)
    end do
    !
    ! step 4: back substitution
    !
    c(n) = c(n)/b(n)
    do j = 1, gap
      i = n-j
      c(i) = (c(i) - d(i)*c(i+1))/b(i)
    end do
    !
    ! step 5: compute spline coefficients
    !
    b(n) = (y(n) - y(gap))/d(gap) + d(gap)*(c(gap) + 2.0*c(n))
    do i = 1, gap
      b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2.0*c(i))
      d(i) = (c(i+1) - c(i))/d(i)
      c(i) = 3.*c(i)
    end do
    c(n) = 3.0d0*c(n)
    d(n) = d(n-1)

return
end subroutine spline_koef
!***********************************************************************************************************

!***********************************************************************************************************
  subroutine spline(u,x,y,b,c,d,n,res)
!======================================================================
! subroutine spline evaluates the cubic spline interpolation at point u
! ispline = y(i)+b(i)*(u-x(i))+c(i)*(u-x(i))**2+d(i)*(u-x(i))**3
! where  x(i) <= u <= x(i+1)
!----------------------------------------------------------------------
! input:
! u       = the abscissa at which the spline is to be evaluated
! x, y    = the arrays of given data points
! b, c, d = arrays of spline coefficients computed by spline
! n       = the number of data points
! output:
! res = interpolated value at point u
!=======================================================================
    implicit none
    integer n
    real*8,dimension(1:n)::x,y,b,c,d
    integer i,j,k
    real*8 u,dx,res

    ! if u is ouside the x() interval take a boundary value (left or right)
    if(u <= x(1)) then
      res = y(1)
      return
    end if
    if(u >= x(n)) then
      res = y(n)
      return
    end if

    !*
    !  binary search for for i, such that x(i) <= u <= x(i+1)
    !*
    i = 1
    j = n+1
    do while (j > i+1)
      k = (i+j)/2
      if(u < x(k)) then
        j=k
        else
        i=k
       end if
    end do
    !*
    !  evaluate spline interpolation
    !*
    dx = u - x(i)
    res = y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))
    
    end subroutine spline
!***********************************************************************************************************
