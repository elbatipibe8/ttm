
!12.02.2019
!Add to the main code

module Module_mod
implicit none
  
  integer mod_CDF,mod_opt,mod_fun,mod_m,mod_d,mod_source
  real*8 Fluence,t_laser,tstart,alfa
  
end module Module_mod

module Module_coeff
        implicit none
          
           real*8 f1,f2,f7,f8,f9,f10,f11,f12,f13,f14
           real*8 phi1,phi2,phi7,phi8,phi9,phi10,phi11,phi12,phi13,phi14,phi18,phi19,phi20,phi21,phi22,phi23
           real*8 koeff0,koeff1,koeff2,koeff3,koeff4,koeff5,koeff6
           real*8 koeff7,koeff8,koeff9,koeff10,koeff11,koeff12,koeff13,koeff14,koeff15,koeff16
          
end module Module_coeff

module Module_CDF
implicit none
  
  integer N_CDF
  real*8,dimension(:),allocatable :: En_CDF_tab,omega_CDF_tab,e1_tab,e2_tab
  real*8,dimension(:),allocatable :: b_e1_koef,c_e1_koef,d_e1_koef,b_e2_koef,c_e2_koef,d_e2_koef
  
end module Module_CDF

      Program The_name_of_your_code
      use Module_mod
      use Module_coeff
      use Module_CDF
      
      real*8 omega_fix,beta_fix,delta_fix
      
      !add the part below after the refinition of variables
       
        N_CDF = 448 !Number of elements in the file with the experimental data of Olmon PRB 2012
        allocate(En_CDF_tab(1:N_CDF),stat=err)
	    allocate(omega_CDF_tab(1:N_CDF),stat=err)
	    allocate(e1_tab(1:N_CDF),stat=err)
	    allocate(e2_tab(1:N_CDF),stat=err)
	    allocate(b_e1_koef(1:N_CDF),stat=err)
	    allocate(c_e1_koef(1:N_CDF),stat=err)
	    allocate(d_e1_koef(1:N_CDF),stat=err)
	    allocate(b_e2_koef(1:N_CDF),stat=err)
	    allocate(c_e2_koef(1:N_CDF),stat=err)
	    allocate(d_e2_koef(1:N_CDF),stat=err)
        
        !Laser parameters
	    Fluence = 40000.0d0 !J/m2 ! Fluence = 40000.0d0 J/m2 for lambda=800nm and Fluence = 500.0d0 J/m2 for lambda=300nm
	    t_laser = 200.0d-15 ! [sec]
	    tstart=3.0d0*t_laser
	    rp_fix=40000.00d-09 !m, radius of the pulse 
	    
	    !======================================== 
         !Fix paremeters  
          omega_fix=2.0d0*pi*c/(300.00d-09) ! 1/sec
          beta_fix=0.00357585d0 !non-dimensional, 
          delta_fix=0.0d0 ! degree
          !at lambda=800 beta=0.00357585 at CDF_4 and delta=0, at lambda=300 beta=0.122802619 at CDF_4 and delta=0 
         !======================================== 
	    
	    mod_CDF=4
        !!mod_CDF=1 (paper Vial PRB 2005 Improved analytical fit of gold dispersion, Ep=[1.24–2.48 eV]),
        !mod_CDF=2 (paper H.S. Sehmi PRB 2017,fit for the paper S. Babar and J. H. Weaver, 
        !Appl. Opt. 54, 477 (2015) with 4 Lorentz terms, Ep=[0.1–6.0 eV] 
     
        !mod_CDF=3 (paper H.S. Sehmi PRB 2017 with a fit to the paper P. B. Johnson and R. W. Christy, 
        !Phys. Rev. B 6, 4370 (1972) with 3 Lorentz terms and epsilon_infinity <> 1),
        !Ep=[0.64–6.6 eV]
     
        !mod_CDF=4 (for Au, paper R.L. Olmon PRB 2012,Ep=[0.04974–4.133 eV], it means lambda=[300 - 24930 nm]
        
        
        !********************************************************************************************************
        !Reading the experimental data of Olmon PRB 2012 and calculation spline coefficients
        open(unit=100,file='CDF_Au_EV_en_e1_e2_Olmon.dat',status='old')
         do i=1,N_CDF
                read(100,*) En_CDF_tab(i),e1_tab(i),e2_tab(i) ! ev
                omega_CDF_tab(i)=En_CDF_tab(i)/hcross ! [1/sec]
         enddo
        close(100)
      
        !Calculation of spline coefficients
         call spline_koef(omega_CDF_tab,e1_tab,b_e1_koef,c_e1_koef,d_e1_koef,N_CDF)
         call spline_koef(omega_CDF_tab,e2_tab,b_e2_koef,c_e2_koef,d_e2_koef,N_CDF)
        !********************************************************************************************************
        
        
        !====================================================================================================================
         !There are coefficients for the source term
          
          koeff0=omega_fix/c
          koeff1=n1(omega_fix) 
          koeff2=k1(omega_fix)
          koeff3=dreal(kz_m(omega_fix)) 
          koeff4=dimag(kz_m(omega_fix)) 
          koeff5=dreal(kx(omega_fix))
          koeff6=dimag(kx(omega_fix)) 
          koeff7=dreal(CDF_m(omega_fix)) 
          koeff8=dimag(CDF_m(omega_fix)) 
          koeff9=koeff1*koeff1+koeff2*koeff2 ! |n|^2
          koeff10=(1.0d0+koeff1)*(1.0d0+koeff1)+koeff2*koeff2 
          koeff11=(8.0d0*koeff0*koeff1*koeff2)/(koeff10) 
          koeff12=koeff0*(koeff7*koeff7+koeff8*koeff8) ! k0*|eps_m|^2
          koeff13=koeff4+koeff0*koeff1
          koeff14=koeff0*koeff2
          koeff15=delta_fix*pi/180.0d0  

          f1=koeff5*koeff7+koeff6*koeff8 
          f2=koeff5*koeff8-koeff6*koeff7 
          f7=f1*(koeff1+koeff9)-f2*koeff2 
          f8=f2*(koeff1+koeff9)+f1*koeff2 
          f9=-koeff7*koeff4+koeff8*koeff3 
          f10=koeff7*koeff3+koeff8*koeff4 
          f11=1.0d0+koeff1+(f9*(koeff1+koeff9)+f10*koeff2)/(koeff12) 
          f12=-koeff2+(-f10*(koeff1+koeff9)+f9*koeff2)/(koeff12) 
          f13=((koeff6*f7-koeff5*f8)/(koeff12))+(koeff3+koeff0*koeff2)*f11+(koeff4+koeff0*koeff1)*f12 
          f14=((koeff6*f8+koeff5*f7)/(koeff12))+(koeff3+koeff0*koeff2)*f12-(koeff4+koeff0*koeff1)*f11 
          
         !====================================================================================================================