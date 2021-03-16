!     Leonid Zhigilei, December 2001
!     Programm solves 1D diffusion equation by the finite difference 
!     method (the explicit forward-time difference formulation).
      SUBROUTINE TTMDiffuse()
        INCLUDE 'common.h'
        INCLUDE 'commonTTM.h'
        REAL*8, DIMENSION(:), ALLOCATABLE :: T_e_new, Ediff
        REAL*8, DIMENSION(:), ALLOCATABLE :: RIL
        INTEGER ixcell(3),iycell(3),izcell(3)

        DATA ixcell/-1, 0, 0/
        DATA iycell/ 0,-1, 0/
        DATA izcell/ 0, 0,-1/
        
        ALLOCATE(T_e_new(0:nlc2-1),Ediff(0:nlc2-1))
        ALLOCATE(RIL(1:nlcz1))

        Time_curr=TIME*1.0E-12
        Etrans(0:nlc2-1)  = 0.0d0
        T_e_new(0:nlc2-1) = T_e(0:nlc2-1)
        
        Time_loop: DO

           Ediff(0:nlc2-1) = 0.0d0

           Do iz = 1,nlcz1
              DO iy = 1,nlcy1
                 Outer: DO ix = 1,nlcx1
            
                    ic = ix + nlcx2*(iy + nlcy2*iz)

                    IF(NVK(ic).eq.0) CYCLE Outer

                    Conduct: DO k = 1,3
                       jx = ix + ixcell(k)
                       jy = iy + iycell(k)
                       jz = iz + izcell(k)

                       jc = jx + nlcx2*(jy + nlcy2*jz)

                       IF(NVK(jc).NE.0) THEN
                          e_cond = Cond(ic,jc)
                       ELSE
                          CYCLE Conduct
                       ENDIF

                       IF(k.eq.1) THEN
                          fmult = 2.0d0*dift/(difx(ic)+difx(jc))
                          fcoic = difx(ic)
                          fcojc = difx(jc)
                       ELSEIF(k.eq.2) THEN
                          fmult = 2.0d0*dift/(dify(ic)+dify(jc))
                          fcoic = dify(ic)
                          fcojc = dify(jc)
                       ELSE
                          fmult = 2.0d0*dift/(difz(ic)+difz(jc))
                          fcoic = difz(ic)
                          fcojc = difz(jc)
                       ENDIF

                       Ediff(ic) = Ediff(ic) + fmult*e_cond*(T_e(jc)-T_e(ic))/fcoic
                       Ediff(jc) = Ediff(jc) + fmult*e_cond*(T_e(ic)-T_e(jc))/fcojc
                    ENDDO Conduct

                 ENDDO Outer
              ENDDO
           ENDDO

           DO ix = 1,nlcx
              DO iy = 1,nlcy
                 RIL(1:nlcz1) = 0.0d0
                 IL0 = 0
                 Inter: DO iz = nlcz,1,-1

                    ic = ix + nlcx2*(iy + nlcy2*iz)

                    IF(NVK(ic).eq.0) THEN
                       T_e_new(ic)=0.0d0
                       CYCLE Inter
                    ENDIF
                    
!                   Distribute energy in MD/bulk connection in proportiions
                    IF(ICOND(ic)) THEN
                       RDENS = 1.0d0
                    ELSE
                       RDENS = RON(ic)/RODENO
                    ENDIF
                    
                    IF(MATER.EQ.2.or.MATER.EQ.3.OR.mater.eq.6) THEN
                       CeF  = (T_e(ic)-TeG(1))*DGR + 1.0d0
                       L    = CeF
                       IF(L.GE.1.AND.L.LT.NGT) THEN
                          PP   = CeF - L
                          G    = GGN(L) + (GGN(L+1)-GGN(L))*PP
                          C_e  = CeN(L) + (CeN(L+1)-CeN(L))*PP
!                          RLT  = RTe(L) + (RTe(L+1)-RTe(L))*PP
                       ELSEIF(L.GE.NGT) THEN
                          G   = GGN(NGT) 
                          C_e = CeN(NGT)
!                          RLT = RTe(NGT)
                       ELSE
                          G   = GGN(1)
                          C_e = CeN(1) !A_e*T_e(ic)
!                          RLT = RTe(1)
                       ENDIF
                    ELSE
                       G   = G !GGN(1)
                       C_e = A_e*T_e(ic)
!                       RLT = 0.95d0
                    ENDIF

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!                   ***** The followin g makes the reflectivity constant ****
                    RLT = 0.869d0 ! for Ni at 800 nm

                    Ecoupl=dift*G*(T_e(ic)-T_l(ic))*RDENS
                    Etrans(ic)=Etrans(ic) - Ecoupl
                    E_pere = E_pere + Ecoupl*volcell(ic)

                    IF(Time_curr.LE.2.0d0*tstart) THEN
!                       IF(ZED(ic).GE.ZNRB.AND.mynodez.eq.1) THEN
                          IF(IL0.EQ.0) THEN
                             RIL(iz+1) = Source_t(Time_curr+0.5d0*dift, &
                                  XED(ic)*1.0d-10,YED(ic)*1.0d-10,RLT)
                             IL0 = 1
                          ENDIF

                          ball = 0.0d0!RELAX(T_e(ic),T_l(ic))
                         
                          alfa = 1.0d0/(ball+OPTPEN)
                          Rstar = RIL(iz+1) + 0.5d0*RDENS*difz(ic)*Source_z(RIL(iz+1),alfa)
                          RIL(iz) = RIL(iz+1) +     RDENS*difz(ic)*Source_z(Rstar,alfa)
                          Rlaser = 0.5d0*(RIL(iz+1) + RIL(iz))
                          Esource = alfa*Rlaser*dift* &
                               RDENS/(1.0d0-(1-KONF)*EXP(-ZLREAL*1.0d-10*alfa))
!                       ELSE
!                          ball = 0.0d0
!                          Esource=dift*Source(ZED(ic)*1.0d-10,&
!                               Time_curr + dift*0.5D0,RLT,ball,zpov(ic)*1.0d-10, &
!                               XED(ic)*1.0d-10,YED(ic)*1.0d-10)* &
!                               RDENS/(1.0d0-(1-KONF)*EXP(-ZLREAL*1.0d-10*alfa))



!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!                       ball = Relax(T_e(ic),T_l(ic))
!                       CALL Source_SPP_8(XED(ic)*1.0d-10,YED(ic)*1.0d-10,ZED(ic)*1.0d-10,&
!                            Time_curr + dift*0.5D0,beta_fix,ball,source_res)
                       
!                       CALL Source_SPP_6(XED(ic)*1.0d-10,YED(ic)*1.0d-10,ZED(ic)*1.0d-10,&
!                            Time_curr + dift*0.5D0,beta_fix,source_res)
                       
!                       Esource=dift*source_res*RDENS ! if beta = 0, we have the same 
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                                              
!                       ENDIF
                    ELSE
                       Esource = 0.0d0
                    ENDIF
                    
                    Eadd_e = Ediff(ic) - Ecoupl + Esource  ! J/m^3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                   Converting Energy to temperature increase based on Ce
                    IF(MATER.EQ.2.or.MATER.EQ.3.OR.mater.eq.6) THEN
                       T_e_new(ic) = T_e(ic) + Eadd_e/C_e/RDENS 
                    ELSE
                       T_e_new(ic)=SQRT(T_e(ic)*T_e(ic) + &
                            2.0D0*Eadd_e/(RDENS*A_e))
                    ENDIF

                    F_total=F_total+Esource*volcell(ic)  ! J

                 ENDDO Inter
              ENDDO
           ENDDO
           
           T_e(0:nlc2-1)=T_e_new(0:nlc2-1)

!          We need to refresh the skin layer with the information
!          about electronic temperature each TTM time-step
           Call Share_Te()

           Time_curr=Time_curr+dift
           
           IF (Time_curr.GT.((TIME+DELTA)*1.0E-12-0.1*dift)) EXIT Time_loop
           
        End do Time_loop


        Etotal_e = 0.0d0
        ETTM = 0.0d0
        DO iz = 1,nlcz
           DO iy = 1,nlcy
              Inside: DO ix = 1,nlcx
                 ic = ix + nlcx2*(iy + nlcy2*iz)

                 IF(NVK(ic).EQ.0) CYCLE Inside

                 IF(ICOND(ic)) THEN
                    RDENS = 1.0d0
                    T_l(ic) = T_l(ic) - Etrans(ic)/(attm + bttm*T_l(ic))
                    ETTM = ETTM + (attm*T_l(ic) + 0.5d0*bttm*T_l(ic) * &
                         T_l(ic))*volcell(ic)*(1.0d0-RON(ic)/RODENO) ! J
                 ELSE
                    RDENS = RON(ic)/RODENO
                 ENDIF

                 IF(MATER.EQ.2.or.MATER.EQ.3.OR.mater.eq.6) THEN
                    CeF  = (T_e(ic)-TeG(1))*DGR + 1.0d0
                    L    = CeF
                    IF(L.GE.1.AND.L.LT.NGT) THEN
                       PP   = CeF - L
                       Eadd = EGG(L) + (EGG(L+1)-EGG(L))*PP
                    ELSEIF(L.GE.NGT) THEN
                       Eadd = EGG(NGT) + CeN(NGT)*(T_e(ic)-TeG(NGT))
                    ELSE
                       Eadd = 0.5d0*A_e*T_e(ic)*T_e(ic)
                    ENDIF
                 ELSE
                    Eadd = 0.5d0*A_e*T_e(ic)*T_e(ic)
                 ENDIF

                 Etotal_e = Etotal_e + Eadd*volcell(ic)*RDENS

                 Etrans(ic) = Etrans(ic)/DELTA 

              ENDDO Inside
           ENDDO
        ENDDO
        
        DEALLOCATE(T_e_new,Ediff)
        DEALLOCATE(RIL)

        RETURN
      END SUBROUTINE TTMDiffuse

      REAL*8 Function Source_t(time_las,xpos,ypos,Refl)
        INCLUDE 'common.h'
        INCLUDE 'commonTTM.h'
!        diar = (xpos-XCENTR)*(xpos-XCENTR) + (ypos-YCENTR)*(ypos-YCENTR)
        Source_t = SourceL*(1.0d0-Refl)*exp(-(time_las- &
             tstart)*(time_las-tstart)/sigmaL2)!*exp(-diar/HDI/HDI)
        
!        argum = PI/(YL*1.0d-10)*ypos
!        Source_t = 2.0d0*Source_tt*COS(argum)*COS(argum)
        return
      END Function Source_t

      REAL*8 Function Source_z(St_z,alfa)
        INCLUDE 'common.h'
        INCLUDE 'commonTTM.h'
        Source_z = -alfa*St_z
        return
      END Function Source_z

      !This is a source WITHOUT sigma=4*ln(2) parameter in the Gauss function over spase
      !The same Gauss fanction over time for the laser and for the SPP.
      subroutine Source_SPP(x,y,z,time_las,beta,source_res)
       INCLUDE 'common.h'
        INCLUDE 'commonTTM.h'
!        use Module_mod
!        use Module_coeff

!        implicit none
        real*8 x,y,z,time,source_res
        real*8 Part0,Part1,Part2,Part3,Part4
        real*8 tetta,sigma,beta,rp
        real*8 f6

!        real*8, parameter :: pi = 2.d0*dacos(0.d0)

      !*********************************************************************

        sigma = 4.0d0*dlog(2.0d0)
        t_laser = DURAT*1.0d-12  !sec
        rp = HDI
      !tstart=3.0d0*t_laser !sec
      !lambda_las=800.00d-09 !m

        Part0=dsqrt(sigma/pi)*dexp(-sigma*(time_las-tstart)*(time_las-tstart)/(t_laser*t_laser))/t_laser
        Part1=dexp(-(x*x+y*y)/(rp*rp))

        Part2=koeff11_1*koeff11_3*dexp(koeff11_3*z)
!        Part2=koeff11*dexp(2.0d0*koeff14*z)

        if (x.ge.(0.0d0)) then

           f6=koeff5*x+koeff13*z+koeff15
           Part3=beta*(2.0d0/koeff10)*dexp(koeff14*z)*dexp(-x*koeff6+z*koeff3)*(f13*dcos(f6)+f14*dsin(f6))
           Part4=2.00d0*beta*beta*((koeff6*f1+koeff3*f9)/(koeff12))*dexp(2.0d0*(-x*koeff6+z*koeff3))

        elseif (x.lt.(0.0d0)) then

           f6=-koeff5*x+koeff13*z+koeff15
           Part3=beta*(2.0d0/koeff10)*dexp(koeff14*z)*dexp(x*koeff6+z*koeff3)*(f13*dcos(f6)+f14*dsin(f6))
           Part4=2.00d0*beta*beta*((koeff6*f1+koeff3*f9)/(koeff12))*dexp(2.0d0*(x*koeff6+z*koeff3))
           
        endif
        source_res=Fluence*Part0*Part1*(Part2-Part3+Part4)
        
        return
      end subroutine Source_SPP

      !This is a source WITH sigma=4*ln(2) parameter in the Gauss function over spase.
      !In this case rp - is a radius of the source at half maximum
      !The same Gauss fanction over time for the laser and for the SPP.
      subroutine Source_SPP_2(x,y,z,time_las,beta,source_res)
        INCLUDE 'common.h'
        INCLUDE 'commonTTM.h'
!        use Module_mod
!        use Module_coeff

!        implicit none
        real*8 x,y,z,time,source_res
        real*8 Part0,Part1,Part2,Part3,Part4
        real*8 tetta,sigma,beta,rp
        real*8 f6

!        real*8, parameter :: pi = 2.d0*dacos(0.d0)

      !*********************************************************************

        sigma = 4.0d0*dlog(2.0d0)
        t_laser = DURAT*1.0d-12  !sec
        rp = HDI
      !tstart=3.0d0*t_laser !sec
      !lambda_las=800.00d-09 !m

        Part0=dsqrt(sigma/pi)*dexp(-sigma*(time_las-tstart)*(time_las-tstart)/(t_laser*t_laser))/t_laser
        Part1=0.25d0*sigma*dexp(-sigma*(x*x+y*y)/(4.0d0*rp*rp)) !It is a difference here in comparison with the Source_SPP
!        Part2=koeff11*dexp(2.0d0*koeff14*z)
        Part2=koeff11_1*koeff11_3*dexp(koeff11_3*z)
        
        if (x.ge.(0.0d0)) then

           f6=koeff5*x+koeff13*z+koeff15
           Part3=beta*(2.0d0/koeff10)*dexp(koeff14*z)*dexp(-x*koeff6+z*koeff3)*(f13*dcos(f6)+f14*dsin(f6))
           Part4=2.00d0*beta*beta*((koeff6*f1+koeff3*f9)/(koeff12))*dexp(2.0d0*(-x*koeff6+z*koeff3))
           
        elseif (x.lt.(0.0d0)) then

           f6=-koeff5*x+koeff13*z+koeff15
           Part3=beta*(2.0d0/koeff10)*dexp(koeff14*z)*dexp(x*koeff6+z*koeff3)*(f13*dcos(f6)+f14*dsin(f6))
           Part4=2.00d0*beta*beta*((koeff6*f1+koeff3*f9)/(koeff12))*dexp(2.0d0*(x*koeff6+z*koeff3))
           
        endif
        source_res=Fluence*Part0*Part1*(Part2-Part3+Part4)
        
        return
      end subroutine Source_SPP_2
      
      !This is a source WITHOUT sigma=4*ln(2) parameter in the Gauss function over spase
      !DIFFERENT Gauss fanctions over time for the laser and for the SPP.
      subroutine Source_SPP_3(x,y,z,time_las,beta,source_res)
       INCLUDE 'common.h'
        INCLUDE 'commonTTM.h'
!        use Module_mod
!        use Module_coeff

!        implicit none
        real*8 x,y,z,time,source_res
        real*8 Part01,Part02,Part03,Part11,Part12,Part13,Part14,Part15,Part2,Part3,Part4
        real*8 tetta,sigma,beta,rp
        real*8 f6

!        real*8, parameter :: pi = 2.d0*dacos(0.d0)

      !*********************************************************************

        sigma = 4.0d0*dlog(2.0d0)
        t_laser = DURAT*1.0d-12  !sec
        rp = HDI
      !tstart=3.0d0*t_laser !sec
      !lambda_las=800.00d-09 !m

        Part01=dsqrt(sigma/pi)*dexp(-sigma*(time_las-tstart)*(time_las-tstart)/(t_laser*t_laser))/t_laser
        Part11=dexp(-(x*x+y*y)/(rp*rp))
        Part14=dexp(-(x*x+y*y)/(2.0d0*rp*rp))*dexp(-(y*y)/(2.0d0*rp*rp))
        Part15=dexp(-(y*y)/(rp*rp))

        Part2=koeff11_1*koeff11_3*dexp(koeff11_3*z)
!        Part2=koeff11*dexp(2.0d0*koeff14*z)

        if (x.ge.(0.0d0)) then !&
           
           Part02=dsqrt(sigma/pi)*dexp(-sigma*(time_las-tstart-x/vgr)* &
           (time_las-tstart-x/vgr)/(t_laser*t_laser))/t_laser
           
           Part03=dsqrt(sigma/pi)*dexp(-sigma*(time_las-tstart)* &
           (time_las-tstart)/(2.0d0*t_laser*t_laser))*dexp(-sigma*(time_las-tstart-x/vgr)* &
           (time_las-tstart-x/vgr)/(2.0d0*t_laser*t_laser))/t_laser
           
           f6=koeff5*x+koeff13*z+koeff15
           Part3=beta*(2.0d0/koeff10)*dexp(koeff14*z)*dexp(-x*koeff6+z*koeff3)*(f13*dcos(f6)+f14*dsin(f6))
           Part4=2.00d0*beta*beta*((koeff6*f1+koeff3*f9)/(koeff12))*dexp(2.0d0*(-x*koeff6+z*koeff3))

        elseif (x.lt.(0.0d0)) then
           
           Part02=dsqrt(sigma/pi)*dexp(-sigma*(time_las-tstart+x/vgr)* &
           (time_las-tstart+x/vgr)/(t_laser*t_laser))/t_laser
           
           Part03=dsqrt(sigma/pi)*dexp(-sigma*(time_las-tstart)* &
           (time_las-tstart)/(2.0d0*t_laser*t_laser))*dexp(-sigma*(time_las-tstart+x/vgr)* &
           (time_las-tstart+x/vgr)/(2.0d0*t_laser*t_laser))/t_laser
           
           f6=-koeff5*x+koeff13*z+koeff15
           Part3=beta*(2.0d0/koeff10)*dexp(koeff14*z)*dexp(x*koeff6+z*koeff3)*(f13*dcos(f6)+f14*dsin(f6))
           Part4=2.00d0*beta*beta*((koeff6*f1+koeff3*f9)/(koeff12))*dexp(2.0d0*(x*koeff6+z*koeff3))
           
        endif
        source_res=Fluence*(Part01*Part11*Part2-Part03*Part14*Part3+Part02*Part15*Part4)
        
        return
      end subroutine Source_SPP_3

            !This is a source WITH sigma=4*ln(2) parameter in the Gauss function over spase
      !DIFFERENT Gauss fanctions over time for the laser and for the SPP.
      subroutine Source_SPP_4(x,y,z,time_las,beta,source_res)
       INCLUDE 'common.h'
        INCLUDE 'commonTTM.h'
!        use Module_mod
!        use Module_coeff

!        implicit none
        real*8 x,y,z,time,source_res
        real*8 Part01,Part02,Part03,Part11,Part12,Part13,Part14,Part15,Part2,Part3,Part4
        real*8 tetta,sigma,beta,rp
        real*8 f6

!        real*8, parameter :: pi = 2.d0*dacos(0.d0)

      !*********************************************************************

        sigma = 4.0d0*dlog(2.0d0)
        t_laser = DURAT*1.0d-12  !sec
        rp = HDI
      !tstart=3.0d0*t_laser !sec
      !lambda_las=800.00d-09 !m

        Part01=dsqrt(sigma/pi)*dexp(-sigma*(time_las-tstart)*(time_las-tstart)/(t_laser*t_laser))/t_laser
        Part11=0.25d0*sigma*dexp(-sigma*(x*x+y*y)/(4.0d0*rp*rp))
        Part14=0.25d0*sigma*dexp(-sigma*(x*x+y*y)/(8.0d0*rp*rp))*dexp(-sigma*(y*y)/(8.0d0*rp*rp))
        Part15=0.25d0*sigma*dexp(-sigma*(y*y)/(4.0d0*rp*rp))
        
!        Part2=koeff11*dexp(2.0d0*koeff14*z)
        Part2=koeff11_1*koeff11_3*dexp(koeff11_3*z)
        
        if (x.ge.(0.0d0)) then !&
           
           Part02=dsqrt(sigma/pi)*dexp(-sigma*(time_las-tstart-x/vgr)* &
           (time_las-tstart-x/vgr)/(t_laser*t_laser))/t_laser
           
           Part03=dsqrt(sigma/pi)*dexp(-sigma*(time_las-tstart)* &
           (time_las-tstart)/(2.0d0*t_laser*t_laser))*dexp(-sigma*(time_las-tstart-x/vgr)* &
           (time_las-tstart-x/vgr)/(2.0d0*t_laser*t_laser))/t_laser
           
           f6=koeff5*x+koeff13*z+koeff15
           Part3=beta*(2.0d0/koeff10)*dexp(koeff14*z)*dexp(-x*koeff6+z*koeff3)*(f13*dcos(f6)+f14*dsin(f6))
           Part4=2.00d0*beta*beta*((koeff6*f1+koeff3*f9)/(koeff12))*dexp(2.0d0*(-x*koeff6+z*koeff3))

        elseif (x.lt.(0.0d0)) then
           
           Part02=dsqrt(sigma/pi)*dexp(-sigma*(time_las-tstart+x/vgr)* &
           (time_las-tstart+x/vgr)/(t_laser*t_laser))/t_laser
           
           Part03=dsqrt(sigma/pi)*dexp(-sigma*(time_las-tstart)* &
           (time_las-tstart)/(2.0d0*t_laser*t_laser))*dexp(-sigma*(time_las-tstart+x/vgr)* &
           (time_las-tstart+x/vgr)/(2.0d0*t_laser*t_laser))/t_laser
           
           f6=-koeff5*x+koeff13*z+koeff15
           Part3=beta*(2.0d0/koeff10)*dexp(koeff14*z)*dexp(x*koeff6+z*koeff3)*(f13*dcos(f6)+f14*dsin(f6))
           Part4=2.00d0*beta*beta*((koeff6*f1+koeff3*f9)/(koeff12))*dexp(2.0d0*(x*koeff6+z*koeff3))
           
        endif
        source_res=Fluence*(Part01*Part11*Part2-Part03*Part14*Part3+Part02*Part15*Part4)
        
        return
      end subroutine Source_SPP_4

      !This is a source WITH 1/e**2 decay in the Gauss function over spase
      !DIFFERENT Gauss fanctions over time for the laser and for the SPP.
      !This is a source with a ballistic range of the 1 model
      subroutine Source_SPP_5(x,y,z,time_las,beta,source_res)
       INCLUDE 'common.h'
        INCLUDE 'commonTTM.h'
!        use Module_mod
!        use Module_coeff

!        implicit none
        real*8 x,y,z,time,source_res
        real*8 Part01,Part02,Part03,Part11,Part12,Part13,Part14,Part15,Part2,Part3,Part4
        real*8 tetta,sigma,beta,rp
        real*8 f6

!        real*8, parameter :: pi = 2.d0*dacos(0.d0)

      !*********************************************************************

        sigma = 4.0d0*dlog(2.0d0)
        t_laser = DURAT*1.0d-12  !sec
        rp = HDI
      !tstart=3.0d0*t_laser !sec
      !lambda_las=800.00d-09 !m

        Part01=dsqrt(sigma/pi)*dexp(-sigma*(time_las-tstart)*(time_las-tstart)/(t_laser*t_laser))/t_laser
        Part11=dexp(-2.0d0*(x*x+y*y)/(rp*rp))
        Part14=dexp(-(x*x+y*y)/(rp*rp))*dexp(-(y*y)/(rp*rp))
        Part15=dexp(-2.0d0*(y*y)/(rp*rp))
        
        Part2=koeff11_1*koeff11_3*dexp(koeff11_3*z)

        if (x.ge.(0.0d0)) then !&
           
           Part02=dsqrt(sigma/pi)*dexp(-sigma*(time_las-tstart-x/vgr)* &
           (time_las-tstart-x/vgr)/(t_laser*t_laser))/t_laser
           
           Part03=dsqrt(sigma/pi)*dexp(-sigma*(time_las-tstart)* &
           (time_las-tstart)/(2.0d0*t_laser*t_laser))*dexp(-sigma*(time_las-tstart-x/vgr)* &
           (time_las-tstart-x/vgr)/(2.0d0*t_laser*t_laser))/t_laser
           
           f6=koeff5*x+koeff13*z+koeff15
           Part3=beta*(2.0d0/koeff10)*dexp(koeff14*z)*dexp(-x*koeff6+z*koeff3)*(f13*dcos(f6)+f14*dsin(f6))
           Part4=2.00d0*beta*beta*((koeff6*f1+koeff3*f9)/(koeff12))*dexp(2.0d0*(-x*koeff6+z*koeff3))

        elseif (x.lt.(0.0d0)) then
           
           Part02=dsqrt(sigma/pi)*dexp(-sigma*(time_las-tstart+x/vgr)* &
           (time_las-tstart+x/vgr)/(t_laser*t_laser))/t_laser
           
           Part03=dsqrt(sigma/pi)*dexp(-sigma*(time_las-tstart)* &
           (time_las-tstart)/(2.0d0*t_laser*t_laser))*dexp(-sigma*(time_las-tstart+x/vgr)* &
           (time_las-tstart+x/vgr)/(2.0d0*t_laser*t_laser))/t_laser
           
           f6=-koeff5*x+koeff13*z+koeff15
           Part3=beta*(2.0d0/koeff10)*dexp(koeff14*z)*dexp(x*koeff6+z*koeff3)*(f13*dcos(f6)+f14*dsin(f6))
           Part4=2.00d0*beta*beta*((koeff6*f1+koeff3*f9)/(koeff12))*dexp(2.0d0*(x*koeff6+z*koeff3))
           
        endif
        source_res=Fluence*(Part01*Part11*Part2-Part03*Part14*Part3+Part02*Part15*Part4)

!        IF(time_las.GE.0.2d-12) THEN
!           IF(myppid.eq.0) print *, source_res,x,y,z,beta
!           CALL MPI_FINALIZE(info)
!           STOP
!        ENDIF
        
        return
      end subroutine Source_SPP_5

     !This is a source WITH the step function over spase,
      !I suppose that the laser spot size is big and
      !we are inside, where the step function gives 1
      !DIFFERENT Gauss fanctions over time for the laser and for the SPP.
      !This is a source with a ballistic range of the 1 model
      subroutine Source_SPP_6(x,y,z,time_las,beta,source_res)
       INCLUDE 'common.h'
        INCLUDE 'commonTTM.h'
!        use Module_mod
!        use Module_coeff

!        implicit none
        real*8 x,y,z,time,source_res
        real*8 Part01,Part02,Part03,Part11,Part12,Part13,Part14,Part15,Part2,Part3,Part4
        real*8 tetta,sigma,beta,rp
        real*8 f6

!        real*8, parameter :: pi = 2.d0*dacos(0.d0)

      !*********************************************************************

        sigma = 4.0d0*dlog(2.0d0)
        t_laser = DURAT*1.0d-12  !sec
        rp = HDI
      !tstart=3.0d0*t_laser !sec
      !lambda_las=800.00d-09 !m

        Part01=dsqrt(sigma/pi)*dexp(-sigma*(time_las-tstart)*(time_las-tstart)/(t_laser*t_laser))/t_laser
        !Part11=dexp(-2.0d0*(x*x+y*y)/(rp*rp))
        !Part14=dexp(-(x*x+y*y)/(rp*rp))*dexp(-(y*y)/(rp*rp))
        !Part15=dexp(-2.0d0*(y*y)/(rp*rp))
        
        Part2=koeff11_1*koeff11_3*dexp(koeff11_3*z)

        if (x.ge.(0.0d0)) then !&
           
           Part02=dsqrt(sigma/pi)*dexp(-sigma*(time_las-tstart-x/vgr)* &
           (time_las-tstart-x/vgr)/(t_laser*t_laser))/t_laser
           
           Part03=dsqrt(sigma/pi)*dexp(-sigma*(time_las-tstart)* &
           (time_las-tstart)/(2.0d0*t_laser*t_laser))*dexp(-sigma*(time_las-tstart-x/vgr)* &
           (time_las-tstart-x/vgr)/(2.0d0*t_laser*t_laser))/t_laser
           
           f6=koeff5*x+koeff13*z+koeff15
           Part3=beta*(2.0d0/koeff10)*dexp(koeff14*z)*dexp(-x*koeff6+z*koeff3)*(f13*dcos(f6)+f14*dsin(f6))
           Part4=2.00d0*beta*beta*((koeff6*f1+koeff3*f9)/(koeff12))*dexp(2.0d0*(-x*koeff6+z*koeff3))

        elseif (x.lt.(0.0d0)) then
           
           Part02=dsqrt(sigma/pi)*dexp(-sigma*(time_las-tstart+x/vgr)* &
           (time_las-tstart+x/vgr)/(t_laser*t_laser))/t_laser
           
           Part03=dsqrt(sigma/pi)*dexp(-sigma*(time_las-tstart)* &
           (time_las-tstart)/(2.0d0*t_laser*t_laser))*dexp(-sigma*(time_las-tstart+x/vgr)* &
           (time_las-tstart+x/vgr)/(2.0d0*t_laser*t_laser))/t_laser
           
           f6=-koeff5*x+koeff13*z+koeff15
           Part3=beta*(2.0d0/koeff10)*dexp(koeff14*z)*dexp(x*koeff6+z*koeff3)*(f13*dcos(f6)+f14*dsin(f6))
           Part4=2.00d0*beta*beta*((koeff6*f1+koeff3*f9)/(koeff12))*dexp(2.0d0*(x*koeff6+z*koeff3))
           
        endif
        source_res=Fluence*(Part01*Part2-Part03*Part3+Part02*Part4)
        
        !source_res=Fluence*(Part01*Part11*Part2-Part03*Part14*Part3+Part02*Part15*Part4)

!        IF(time_las.GE.0.2d-12) THEN
!           IF(myppid.eq.0) print *, source_res,x,y,z,beta
!           CALL MPI_FINALIZE(info)
!           STOP
!        ENDIF
        
        return
      end subroutine Source_SPP_6
      
      !This is a source WITH 1/e**2 decay in the Gauss function over spase
      !DIFFERENT Gauss fanctions over time for the laser and for the SPP.
      !This is a source with a ballistic range of the 1 model as Relax(Tel,Tph)
      subroutine Source_SPP_7(x,y,z,time_las,beta,ball_range,source_res)
       INCLUDE 'common.h'
        INCLUDE 'commonTTM.h'
!        use Module_mod
!        use Module_coeff

!        implicit none
        real*8 x,y,z,time,source_res
        real*8 Part01,Part02,Part03,Part11,Part12,Part13,Part14,Part15,Part2,Part3,Part4
        real*8 tetta,sigma,beta,rp
        real*8 f6,koeff_ball,ball_range

!        real*8, parameter :: pi = 2.d0*dacos(0.d0)

      !*********************************************************************

        sigma = 4.0d0*dlog(2.0d0)
        t_laser = DURAT*1.0d-12  !sec
        rp = HDI
      !tstart=3.0d0*t_laser !sec
      !lambda_las=800.00d-09 !m

        Part01=dsqrt(sigma/pi)*dexp(-sigma*(time_las-tstart)*(time_las-tstart)/(t_laser*t_laser))/t_laser
        Part11=dexp(-2.0d0*(x*x+y*y)/(rp*rp))
        Part14=dexp(-(x*x+y*y)/(rp*rp))*dexp(-(y*y)/(rp*rp))
        Part15=dexp(-2.0d0*(y*y)/(rp*rp))
        
        koeff_ball=koeff11_2/(1.000d0+koeff11_2*ball_range)
        Part2=koeff11_1*koeff_ball*dexp(koeff_ball*z)

        if (x.ge.(0.0d0)) then !&
           
           Part02=dsqrt(sigma/pi)*dexp(-sigma*(time_las-tstart-x/vgr)* &
           (time_las-tstart-x/vgr)/(t_laser*t_laser))/t_laser
           
           Part03=dsqrt(sigma/pi)*dexp(-sigma*(time_las-tstart)* &
           (time_las-tstart)/(2.0d0*t_laser*t_laser))*dexp(-sigma*(time_las-tstart-x/vgr)* &
           (time_las-tstart-x/vgr)/(2.0d0*t_laser*t_laser))/t_laser
           
           f6=koeff5*x+koeff13*z+koeff15
           Part3=beta*(2.0d0/koeff10)*dexp(koeff14*z)*dexp(-x*koeff6+z*koeff3)*(f13*dcos(f6)+f14*dsin(f6))
           Part4=2.00d0*beta*beta*((koeff6*f1+koeff3*f9)/(koeff12))*dexp(2.0d0*(-x*koeff6+z*koeff3))

        elseif (x.lt.(0.0d0)) then
           
           Part02=dsqrt(sigma/pi)*dexp(-sigma*(time_las-tstart+x/vgr)* &
           (time_las-tstart+x/vgr)/(t_laser*t_laser))/t_laser
           
           Part03=dsqrt(sigma/pi)*dexp(-sigma*(time_las-tstart)* &
           (time_las-tstart)/(2.0d0*t_laser*t_laser))*dexp(-sigma*(time_las-tstart+x/vgr)* &
           (time_las-tstart+x/vgr)/(2.0d0*t_laser*t_laser))/t_laser
           
           f6=-koeff5*x+koeff13*z+koeff15
           Part3=beta*(2.0d0/koeff10)*dexp(koeff14*z)*dexp(x*koeff6+z*koeff3)*(f13*dcos(f6)+f14*dsin(f6))
           Part4=2.00d0*beta*beta*((koeff6*f1+koeff3*f9)/(koeff12))*dexp(2.0d0*(x*koeff6+z*koeff3))
           
        endif
        source_res=Fluence*(Part01*Part11*Part2-Part03*Part14*Part3+Part02*Part15*Part4)

!        IF(time_las.GE.0.2d-12) THEN
!           IF(myppid.eq.0) print *, source_res,x,y,z,beta
!           CALL MPI_FINALIZE(info)
!           STOP
!        ENDIF
        
        return
      end subroutine Source_SPP_7
      
      !This is a source WITH the step function over spase,
      !I suppose that the laser spot size is big and
      !we are inside, where the step function gives 1
      !DIFFERENT Gauss fanctions over time for the laser and for the SPP.
      !This is a source with a ballistic range of the 1 model as Relax(Tel,Tph)
      subroutine Source_SPP_8(x,y,z,time_las,beta,ball_range,source_res)
       INCLUDE 'common.h'
        INCLUDE 'commonTTM.h'
!        use Module_mod
!        use Module_coeff

!        implicit none
        real*8 x,y,z,time,source_res
        real*8 Part01,Part02,Part03,Part11,Part12,Part13,Part14,Part15,Part2,Part3,Part4
        real*8 tetta,sigma,beta,rp
        real*8 f6,koeff_ball,ball_range

!        real*8, parameter :: pi = 2.d0*dacos(0.d0)

      !*********************************************************************

        sigma = 4.0d0*dlog(2.0d0)
        t_laser = DURAT*1.0d-12  !sec
        rp = HDI
      !tstart=3.0d0*t_laser !sec
      !lambda_las=800.00d-09 !m

        Part01=dsqrt(sigma/pi)*dexp(-sigma*(time_las-tstart)*(time_las-tstart)/(t_laser*t_laser))/t_laser
!        Part11=dexp(-2.0d0*(x*x+y*y)/(rp*rp))
!        Part14=dexp(-(x*x+y*y)/(rp*rp))*dexp(-(y*y)/(rp*rp))
!        Part15=dexp(-2.0d0*(y*y)/(rp*rp))
        
        koeff_ball=koeff11_2/(1.000d0+koeff11_2*ball_range)
        Part2=koeff11_1*koeff_ball*dexp(koeff_ball*z)

        if (x.ge.(0.0d0)) then !&
           
           Part02=dsqrt(sigma/pi)*dexp(-sigma*(time_las-tstart-x/vgr)* &
           (time_las-tstart-x/vgr)/(t_laser*t_laser))/t_laser
           
           Part03=dsqrt(sigma/pi)*dexp(-sigma*(time_las-tstart)* &
           (time_las-tstart)/(2.0d0*t_laser*t_laser))*dexp(-sigma*(time_las-tstart-x/vgr)* &
           (time_las-tstart-x/vgr)/(2.0d0*t_laser*t_laser))/t_laser
           
           f6=koeff5*x+koeff13*z+koeff15
           Part3=beta*(2.0d0/koeff10)*dexp(koeff14*z)*dexp(-x*koeff6+z*koeff3)*(f13*dcos(f6)+f14*dsin(f6))
           Part4=2.00d0*beta*beta*((koeff6*f1+koeff3*f9)/(koeff12))*dexp(2.0d0*(-x*koeff6+z*koeff3))

        elseif (x.lt.(0.0d0)) then
           
           Part02=dsqrt(sigma/pi)*dexp(-sigma*(time_las-tstart+x/vgr)* &
           (time_las-tstart+x/vgr)/(t_laser*t_laser))/t_laser
           
           Part03=dsqrt(sigma/pi)*dexp(-sigma*(time_las-tstart)* &
           (time_las-tstart)/(2.0d0*t_laser*t_laser))*dexp(-sigma*(time_las-tstart+x/vgr)* &
           (time_las-tstart+x/vgr)/(2.0d0*t_laser*t_laser))/t_laser
           
           f6=-koeff5*x+koeff13*z+koeff15
           Part3=beta*(2.0d0/koeff10)*dexp(koeff14*z)*dexp(x*koeff6+z*koeff3)*(f13*dcos(f6)+f14*dsin(f6))
           Part4=2.00d0*beta*beta*((koeff6*f1+koeff3*f9)/(koeff12))*dexp(2.0d0*(x*koeff6+z*koeff3))
           
        endif
        source_res=Fluence*(Part01*Part2-Part03*Part3+Part02*Part4)
        
        !source_res=Fluence*(Part01*Part11*Part2-Part03*Part14*Part3+Part02*Part15*Part4)

!        IF(time_las.GE.0.2d-12) THEN
!           IF(myppid.eq.0) print *, source_res,x,y,z,beta
!           CALL MPI_FINALIZE(info)
!           STOP
!        ENDIF
        
        return
      end subroutine Source_SPP_8
      
      FUNCTION Source(zpos,timel,Reflect,brange,zpol,xpos,ypos)
!       Source term is described by the Lambert-Beer law for spatial 
!       distribution and Gaussian profile for temporal distribution.      
        INCLUDE 'common.h'
        INCLUDE 'commonTTM.h'
        alfas   = 1.0d0/(OPTPEN+brange)
        diar = (xpos-XCENTR)*(xpos-XCENTR) + (ypos-YCENTR)*(ypos-YCENTR)
        SourceLL = alfas*SourceL*(1.0d0-Reflect)*exp((zpos-zpol)*alfas)
        Source=SourceLL*exp(-(timel-tstart)*(timel-tstart)/sigmaL2)*exp(-diar/HDI/HDI)
!*****        argum = PI/(YL*1.0d-10)*ypos
!*****        Source = SourceLLLL*COS(argum)*COS(argum)
        RETURN
      END FUNCTION Source


      FUNCTION Relax(Tel,Tph)
        INCLUDE 'common.h'
        INCLUDE 'commonTTM.h'
        xsik = 11749500.0d0
        Relax = xsik/(Ac_e*Tel*Tel+Bc_e*Tph)
        RETURN
      END FUNCTION Relax


      FUNCTION Cond(K,M)
!       Conductivity of electron system can be approximated as (1):K0*T_e/T_l
!       (2):1/3*Vfer**2*A_e*T_e/(B*T_l+A*T_e**2) (3) Anisimov approx
!       Take average value of conductivity from left and right nodes
        INCLUDE 'common.h'
        INCLUDE 'commonTTM.h'
        
        Cond = 0.0d0
           
           Select case (KON)
           
           case (1)   
              
              ek_k=ABS(1.0d0 - FSOL(K))*ek_liq + FSOL(K)*ek
              ek_m=ABS(1.0d0 - FSOL(M))*ek_liq + FSOL(M)*ek

              IF(NPR(K).EQ.0.AND.NPR(M).EQ.1) THEN
                 Cond = 0.5D0*(ek_k*T_e(K)/T_l(K) + &
                      ek_m*T_e(M)/T_l(M))*RON(M)/RODENO
              ELSEIF(NPR(K).EQ.1.AND.NPR(M).EQ.0) THEN
                 Cond = 0.5D0*(ek_k*T_e(K)/T_l(K) +  &
                      ek_m*T_e(M)/T_l(M))*RON(K)/RODENO
              ELSEIF(NPR(K).NE.0.AND.NPR(M).NE.0) THEN
                 Cond = 0.5D0*(ek_k*RON(K)*T_e(K)/T_l(K) + &
                      ek_m*RON(M)*T_e(M)/T_l(M))/RODENO
              ENDIF
              
              IF(ICOND(K).OR.ICOND(M)) Cond = &
                   0.5d0*ek*(T_e(K)/T_l(K) + T_e(M)/T_l(M))

           case (2)

              IF(NPR(K).EQ.0.AND.NPR(M).EQ.1) THEN
                 Cond = Vfer*Vfer*A_e*(T_e(K)/(Bc_e*T_l(K) + &
                      Ac_e*T_e(K)*T_e(K))+T_e(M)/(Bc_e*T_l(M) + &
                      Ac_e*T_e(M)*T_e(M)))*RON(M)/RODENO/6.0d0
              ELSEIF(NPR(K).EQ.1.AND.NPR(M).EQ.0) THEN
                 Cond = Vfer*Vfer*A_e*(T_e(K)/(Bc_e*T_l(K) + &
                      Ac_e*T_e(K)*T_e(K))+T_e(M)/(Bc_e*T_l(M) + &
                      Ac_e*T_e(M)*T_e(M)))*RON(K)/RODENO/6.0d0
              ELSEIF(NPR(K).NE.0.AND.NPR(M).NE.0) THEN
                 Cond = Vfer*Vfer*A_e*(RON(K)*T_e(K)/(Bc_e*T_l(K) + &
                      Ac_e*T_e(K)*T_e(K))+RON(M)*T_e(M)/(Bc_e*T_l(M) + &
                      Ac_e*T_e(M)*T_e(M)))/RODENO/6.0d0
              ENDIF

              IF(ICOND(K).OR.ICOND(M)) Cond = &
                   Vfer*Vfer*A_e*(T_e(K)/(Bc_e*T_l(K) + Ac_e*T_e(K)*T_e(K)) + &
                   T_e(M)/(Bc_e*T_l(M)+Ac_e*T_e(M)*T_e(M)))/6.0d0

           case (3)

              te_k=BK*T_e(K)/Efer
              te_m=BK*T_e(M)/Efer
              tl_k=BK*T_l(K)/Efer
              tl_m=BK*T_l(M)/Efer
              cond_k=xappa*(te_k*te_k+0.16d0)**(5.0d0/4.0d0)*(te_k*te_k+ &
                   0.44d0)*te_k/(SQRT(te_k*te_k+0.092d0)*(te_k*te_k+xnum*tl_k))
              cond_m=xappa*(te_m*te_m+0.16d0)**(5.0d0/4.0d0)*(te_m*te_m+ &
                   0.44d0)*te_m/(SQRT(te_m*te_m+0.092d0)*(te_m*te_m+xnum*tl_m))


              IF(NPR(K).EQ.0.AND.NPR(M).EQ.1) THEN
                 Cond = 0.5d0*(cond_k + cond_m)*RON(M)/RODENO
                 ELSEIF(NPR(K).EQ.1.AND.NPR(M).EQ.0) THEN
                 Cond = 0.5d0*(cond_k + cond_m)*RON(K)/RODENO
              ELSEIF(NPR(K).NE.0.AND.NPR(M).NE.0) THEN
                 Cond = 0.5d0*(RON(K)*cond_k + RON(M)*cond_m)/RODENO
              ENDIF
              
              IF(ICOND(K).OR.ICOND(M)) Cond = 0.5d0*(cond_k + cond_m)

           case (4)

!              a04 = 9.294d0  !!!
!              b04 = 0.03d0!!!
!              b14 =-0.2688d0!!!
!              b24 = 0.9722d0!!!
              
              IF(NPR(K).EQ.0.AND.NPR(M).EQ.1) THEN
                 r_k  = rden(M)
                 te_k = t6*T_e(M)/r_k
                 conds_k = conds(M)
                 condl_k = condl(M)
              ELSEIF(NPR(K).EQ.1.AND.NPR(M).EQ.0) THEN
                 r_k  = rden(K)
                 te_k = t6*T_e(K)/r_k
                 conds_k = conds(K)
                 condl_k = condl(K)
              ELSEIF(NPR(K).EQ.1.AND.NPR(M).EQ.1) THEN
                 r_k  = 0.5d0*(rden(K) + rden(M))
                 IF(ICOND(K)) THEN
                    RDENSK = RODENO
                 ELSE
                    RDENSK = RON(K)
                 ENDIF
                 IF(ICOND(M)) THEN
                    RDENSM = RODENO
                 ELSE
                    RDENSM = RON(M)
                 ENDIF
                 te_k = t6*(RDENSK*T_e(K) + RDENSM*T_e(M))/(RDENSK+RDENSM)/r_k
                 conds_k = (RDENSK*conds(K)+RDENSM*conds(M))/(RDENSK+RDENSM)
                 condl_k = (RDENSK*condl(K)+RDENSM*condl(M))/(RDENSK+RDENSM)
              ELSE
                 GOTO 765
              ENDIF
              
!              te_k = t6*T_e(K)/r_k
!              te_m = t6*T_e(M)/r_m
              
!              c_0_k  =  131.0d0*te_k* (1.0d0 + 3.07d0 *te_k *te_k)/(1.0d0 + 1.08d0 *te_k**(2.07d0))
!              c_0_m  =  131.0d0*te_m* (1.0d0 + 3.07d0 *te_m *te_m)/(1.0d0 + 1.08d0 *te_m**(2.07d0))
              
              conde_k = (1.0d-4*r_k**(-4.0d0/3.0d0)*9.294d0*te_k)/ &
                   (1.0d0 + 0.03d0*SQRT(te_k) - 0.2688d0*te_k + 0.9722d0*te_k*te_k)
              
!              conde_m = (1.0d-4*r_m**(-4.0d0/3.0d0)*a04*te_m)/ &
!                   (1.0d0 + b04*SQRT(te_m) + b14*te_m + b24*te_m*te_m)
              
!              Cond_k = c_0_k/(c_0_k*conde_k + conds(K) + condl(K))
!              Cond_m = c_0_m/(c_0_m*conde_m + conds(M) + condl(M))

              Cond = 1.0d0/(conde_k + conds_k + condl_k)
!              Cond_m = 1.0d0/(conde_m + conds(M) + condl(M))
              
!              Cond = 0.5d0*(Cond_k + Cond_m)

765           CONTINUE
           end select

        RETURN
      END FUNCTION Cond
