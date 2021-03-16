!     Defines lattice T at each node and 
!     the center of mass velocity of each node.
      SUBROUTINE TTMTemper()
        INCLUDE 'common.h'
        INCLUDE 'commonTTM.h'
        INTEGER ixcell(6),iycell(6),izcell(6)
        
        DATA ixcell/-1, 1, 0, 0, 0, 0/
        DATA iycell/ 0, 0,-1, 1, 0, 0/
        DATA izcell/ 0, 0, 0, 0,-1, 1/

        Q1V(0:nlc2-1) = 0.0d0 ! Velocity of the center of mass       
        RON(0:nlc2-1) = 0.0d0 ! Density of the cell
        NPR(0:nlc2-1) = 0
!        T_l(0:nlc2-1) = 0.0d0
        
        DO iz = 1,nlcz
           DO iy = 1,nlcy
              DO ix = 1,nlcx                 
                 ip = ix + nlcx2*(iy + nlcy2*iz)

                 mcount = 0
                 XMASSN = 0.0d0
                 
                 i = ltop(ip)
                 IF(i.eq.0) GOTO 21

22               IF(.NOT.((KHIST(i).EQ.3.OR.KHIST(i).GE.10) &
                      .AND.NOTEM.OR.KTYPE(I).EQ.KWT))THEN
                    mcount = mcount + 1
                    IF(.NOT.(IMMOV)) THEN
                       Q1V(ip)=Q1V(ip)+Q1D(3,i)*XMASS(KTYPE(I))/DELTA
                       XMASSN=XMASSN+XMASS(KTYPE(I))
                    ENDIF
                 ENDIF

                 IF(.NOT.(IBULK(ip)).AND.KHIST(i).GE.10.AND.KONF.EQ.1 &
                      .AND.ZED(ip).LE.2.0d0*rstep.AND.ZED(ip) &
                         .GE.(ZCRIT+0.25d0*rstep)) IBULK(ip) = .TRUE.
                 
                 i = link(i)
                 IF(i.ne.0)GOTO 22

                 RON(ip) = DBLE(mcount)
                 
21               CONTINUE 
                 
!                Define the velocity of center of mass for activated cells
                 IF(.NOT.(IMMOV)) THEN
                    IF(RON(ip).GE.RMINL*RODENO) THEN
                       Q1V(ip)= Q1V(ip)/XMASSN
                    ELSE
                       Q1V(ip) = 0.0d0
                    ENDIF
                 ENDIF
                 
                 IF(IBULK(ip).AND.RON(ip).LT.RMINR*RODENO) Q1V(ip) = 0.0d0
                 
                 IF(RON(ip).GE.RMINH*RODENO.OR.IBULK(ip)) THEN
                    NPR(ip) = 1
                 ELSE
                    NPR(ip) = 0
                 ENDIF

                 IF(IBULK(ip).AND.RON(ip).LT.RMINR*RODENO) THEN
                    ICOND(ip) = .TRUE.
                 ELSE
                    ICOND(ip) = .FALSE.
                 ENDIF

              ENDDO
           ENDDO
        ENDDO

!       Additional information about liquid when writing out data
        IF(FLUID) THEN
        idlc(0:nlc2-1) = 0
        ippl = 0
           DO iz = 1,nlcz
              DO iy = 1,nlcy
                 Fluid1: DO ix = 1,nlcx
                    ip = ix + nlcx2*(iy + nlcy2*iz)
                    
                    IF(NVK(ip).EQ.1.OR.IBULK(ip)) CYCLE Fluid1

                    mcount = 0
                    XMASSW = 0.0d0
                    Q1W = 0.0d0

                    i = ltop(ip)
                    IF(i.eq.0) GOTO 25

26                  IF(KTYPE(I).EQ.KWT)THEN
                       mcount = mcount + 1
                       Q1W = Q1W + Q1D(3,i)*XMASS(KTYPE(I))/DELTA
                       XMASSW=XMASSW+XMASS(KTYPE(I))
                    ENDIF

                    i = link(i)
                    IF(i.ne.0)GOTO 26

                    IF(mcount.GE.RMINH*DENWAT.AND.RON(ip).LT.RMINL*RODENO) THEN
                       ippl = ippl + 1
                       idlc(ippl) = ip
                       Q1V(ip)= Q1W/XMASSW
                       RON(ip) = DBLE(mcount)
                    ENDIF

25                  CONTINUE
                 ENDDO Fluid1
              ENDDO
           ENDDO           
        ENDIF
        
!        IF(.NOT.(IMMOV)) THEN
           CALL SHARE_Q1V()
           CALL Smooth_Q1V()
!        ENDIF
        
        IF(FLUID) THEN
           DO ii=1,ippl
              ip=idlc(ii)
              RON(ip) = 0.0d0
           ENDDO
        ENDIF

!       Temperature calculation
        QINK(0:nlc2-1)=0.0d0
        
        DO iz = 1,nlcz
           DO iy = 1,nlcy
              Tem: DO ix = 1,nlcx
                 
                 ip = ix + nlcx2*(iy + nlcy2*iz)

                 IF(ICOND(ip)) CYCLE Tem
                    
                 IF(.NOT.(IBULK(ip)).AND.RON(ip).LT.RMINL*RODENO) THEN
                    T_e(ip) = 0.0d0
                    T_l(ip) = 0.0d0
                    NVK(ip) = 0
                    CYCLE Tem
                 ENDIF

                 i = ltop(ip)
                 IF(i.eq.0) GOTO 24
                 
23               IF(.NOT.((KHIST(i).EQ.3.OR.KHIST(i).GE.10) &
                      .AND.NOTEM.OR.KTYPE(I).EQ.KWT)) THEN
                    QINK(ip)=QINK(ip)+ &
                         0.5*XMASS(KTYPE(I))*Q1D(1,i)*Q1D(1,i)/(DELTA*DELTA)+ &
                         0.5*XMASS(KTYPE(I))*Q1D(2,i)*Q1D(2,i)/(DELTA*DELTA)+ &
                         0.5*XMASS(KTYPE(I))*(Q1D(3,i)/DELTA - &
                         Q1V(ip))*(Q1D(3,i)/DELTA-Q1V(ip))
                 ENDIF
                    
                 i = link(i)
                 IF(i.ne.0)GOTO 23

24               CONTINUE
                 
                 T_l(ip) = QINK(ip)*ENUNIT*2.0d0/BK/DBLE(NDIM)/RON(ip)

!!!!!!!!!!!!!!!!!!!!!!!!!!!! "Shold never happen, but just in case )) !!!!!!!!!!!!!!!!!!
                 IF(NVK(ip).eq.1.and.T_e(ip).eq.0.0d0) T_e(ip) = T_l(ip)

!                If there is a new cell as a result of mass transfer, it has to
!                be given some temperature. The temperature value is an average
!                temperature over its 6 neighbores or equal to the temperature
!                of the lattice in this cell if there are no neighbores around.
                 IF(NVK(ip).EQ.0) THEN

                    xcount = 0.0d0

                    DO k = 1,6

                       jx = ix + ixcell(k)
                       jy = iy + iycell(k)
                       jz = iz + izcell(k)                     
                       
                       jc = jx + nlcx2*(jy + nlcy2*jz)
                       
                       IF(NVK(jc).EQ.1)THEN
                          IF(ICOND(jc)) THEN
                             RDENS = RODENO
                          ELSE
                             RDENS = RON(jc)
                          ENDIF
                          xcount = xcount + RDENS
                          T_e(ip) = T_e(ip) + T_e(jc)*RDENS
                       ENDIF
                    ENDDO
                    
                    IF(xcount.EQ.0) THEN
                       T_e(ip) = T_l(ip)
                    ELSE
                       T_e(ip) = T_e(ip)/xcount
                    ENDIF
                    NVK(ip) = 1
                 ENDIF

               ENDDO Tem
           ENDDO
        ENDDO

        IF(FLUID) THEN
           DO ii = 1,ippl
              ip = idlc(ii)
              i = ltop(ip)
              mcount = 0
              
              IF(i.eq.0) GOTO 27


28            IF(KTYPE(I).EQ.KWT) THEN
                 mcount = mcount + 1
                 QINK(ip)=QINK(ip)+ &
                      0.5*XMASS(KTYPE(I))*Q1D(1,i)*Q1D(1,i)/(DELTA*DELTA)+ &
                      0.5*XMASS(KTYPE(I))*Q1D(2,i)*Q1D(2,i)/(DELTA*DELTA)+ &
                      0.5*XMASS(KTYPE(I))*(Q1D(3,i)/DELTA - &
                      Q1V(ip))*(Q1D(3,i)/DELTA-Q1V(ip))
              ENDIF
                 
              i = link(i)
              IF(i.ne.0)GOTO 28

              RON(ip) = DBLE(mcount)
              T_l(ip) = QINK(ip)*ENUNIT*2.0d0/BK/DBLE(NDIM)/RON(ip)
              
27            CONTINUE
           ENDDO
        ENDIF


        zpov(0:nlc2-1) = 0.0d0
        IF(Time_curr.LE.2.0d0*tstart) THEN
           DO ix = 1,nlcx
              DO iy = 1,nlcy
                 zloop: DO iz = nlcz,1,-1
                    ic = ix + nlcx2*(iy + nlcy2*iz)
                    IF(NVK(ic).EQ.1.AND.iz.lt.nlcz.and.mynodez.eq.1) THEN
                       zpov(ic) = ZED(ic)
                       EXIT zloop
                    ENDIF
                 ENDDO zloop
              ENDDO
           ENDDO
        ENDIF
        
!       Drawing current relative density per atom
        IF(KON.EQ.4) THEN
           rden(0:nlc2-1) = 0.0d0
 
           DO iz = 1,nlcz
              DO iy = 1,nlcy
                 Tem2: DO ix = 1,nlcx
                    
                    ip = ix + nlcx2*(iy + nlcy2*iz)
                    
                    IF(NVK(ip).EQ.0) CYCLE Tem2
                    
                    IF(ICOND(ip)) THEN
                       rden(ip) = 1.0d0
                       CYCLE Tem2
                    ENDIF
                       
                    dmcount = 0.0d0
                    dncount = 0.0d0
                    
                    i = ltop(ip)
                    IF(i.eq.0) GOTO 30
                    
29                  dmcount = dmcount + DBLE(NNG(i))
                    dncount = dncount + 1.0d0
                    
                    i = link(i)
                    IF(i.ne.0) GOTO 29
                    
30                  CONTINUE
                    
                    rden(ip) = dmcount/dncount/vden ! will share it in Share_Tl
                    
                 ENDDO Tem2
              ENDDO
           ENDDO
           
           CALL Share_TlD()
           CALL Smooth_TlD()

           Efer4 = 9.2d0!!!!!
           r_rt = 19.3d0/19.5d0!
           T_rt = 293.0d0!
           t6 = 6.0d0*BK/Efer4!
           te_rt = t6*T_rt/r_rt!

           a4 =3.920d0!
           b4 = 1.95d0!
           alpha4 = 2.0d0*a4+1.0d0!
           beta4 = a4 + 1.0d0!
           c_ab = (a4-b4)/(b4+1.0d0)!
           c_0_rt=   131.0d0*te_rt*(1.0d0 + 3.07d0*te_rt*te_rt)/(1.0d0 + 1.08d0*te_rt**(2.07d0))
           y_rt =  (1.0d0 + c_ab)*r_rt**alpha4/(1.0d0 + c_ab*r_rt**beta4)
           
           rden(0:nlc2-1) = r_rt*rden(0:nlc2-1)
           condl(0:nlc2-1) = 0.0d0
           conds(0:nlc2-1) = 0.0d0


           DO iz = 1,nlcz
              DO iy = 1,nlcy
                 Tem3: DO ix = 1,nlcx

                    ip = ix + nlcx2*(iy + nlcy2*iz)

                    IF(NVK(ip).EQ.0) CYCLE Tem3

                    te_ip = t6*T_e(ip)/rden(ip)!
                    c_0_ip = 131.0d0*te_ip* (1.0d0 + 3.07d0 *te_ip *te_ip)/ &
                         (1.0d0 + 1.08d0 *te_ip**(2.07d0))

                    IF(FSOL(ip).GE.0.9d0) THEN
                       y_s  =  (1.0d0 + c_ab)* rden(ip)**alpha4/(1.0d0 +  c_ab*rden(ip)**beta4)
                       conds(ip) =318.0d0*rden(ip)**(4.0d0/3.0d0)*(y_s/y_rt)*(T_rt/T_l(ip))*c_0_ip/c_0_rt
                       conds(ip) = 1.0d0/conds(ip)
                       condl(ip) = 0.0d0
                    ELSEIF(FSOL(ip).LE.0.1d0) THEN
                       t3_s = T_l(ip)*1.0d-3
                       xl_s = 0.887179d0 - 0.0328321d0*(t3_s - 1.337d0) - &
                            0.0030982d0*(t3_s - 1.337d0)**2.0d0 - &
                            0.000164884d0*(t3_s - 1.337d0)**3.0d0
                       rl_s = 148.5d0 + 119.3d0*t3_s*15.337d0/(14.0d0+t3_s)
                       condl(ip) = c_0_ip*(3254.0d0/rl_s)*rden(ip)*(rden(ip)/xl_s)**2.0d0
                       condl(ip) = 1.0d0/condl(ip)
                       conds(ip) = 0.0d0
                    ELSE
                       y_s  =  (1.0d0 + c_ab)* rden(ip)**alpha4/(1.0d0 +  c_ab*rden(ip)**beta4)
                       conds(ip) =318.0d0*rden(ip)**(4.0d0/3.0d0)*(y_s/y_rt)*(T_rt/T_l(ip))*c_0_ip/c_0_rt
                       t3_s = T_l(ip)*1.0d-3
                       xl_s = 0.887179d0 - 0.0328321d0*(t3_s - 1.337d0) - &
                            0.0030982d0*(t3_s - 1.337d0)**2.0d0 - &
                            0.000164884d0*(t3_s - 1.337d0)**3.0d0
                       rl_s = 148.5d0 + 119.3d0*t3_s*15.337d0/(14.0d0+t3_s)
                       condl(ip) = c_0_ip*(3254.0d0/rl_s)*rden(ip)*(rden(ip)/xl_s)**2.0d0
                       conds(ip) = FSOL(ip)/conds(ip)
                       condl(ip) = (1.0d0-FSOL(ip))/condl(ip)
                    ENDIF
          
                 ENDDO Tem3
              ENDDO
           ENDDO

           CALL Share_DSL()
           
        ELSE
           Call Share_Tl()
           Call Smooth_Tl()
        ENDIF
           
       RETURN
      END SUBROUTINE TTMTemper
