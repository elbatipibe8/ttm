!     Program brakes down film for nodes

!     Setting initial data for chosen material
      Subroutine TTMInit()
!     Parameters of material
!     ek is thermal conductivity [W m-1 K-1]
!     Cp_e=A_e*T_e
!     G is in [W K-1 m-3]
!     Vfer is Fermi velosity for the given material [m s-1]
!     Ac_e is specific constant for conductivity [s-1 K-2]
!     Bc_e is specific constant for conductivity [s-1 K-1]
!     H_fus is fusion Enthalpy in [J m-3]
!     T_melt is melting temperature [k]
!     H_abl is vaporization enthalpy in [J m-3]
!     T_abl is vaporization temperature [K]     

        INCLUDE 'common.h'
        INCLUDE 'commonTTM.h' 
        INCLUDE 'mpif.h'
!        INTEGER status(MPI_STATUS_SIZE)

        Select case (MATER)
        case (1)   ! The following parameters are for Cu
           A_e=96.6D0    ! In electron heat capasity
           ek=401.0D0    ! In electron conductivity
           ek_liq=163.0d0! [Mills et al.,Int. Mater. Rev. 41, #6, p.209 (1996)]
           G=2.0D+17     ! Coupling constant
           Vfer=1.57D+06 ! Fermi velocity
           Ac_e=1.0D+07  ! In electron conductivity
           Bc_e=1.08D+11 ! In electron conductivity
           xappa=377.0d0 ! In electron conductivity
           xnum=0.139d0  ! In electron conductivity
           Efer=7.03d0   ! Fermi energy
           gren=0.92     ! Guenizer coefficient, not used currently
           P_ZER=-0.00d0*NAN    ! Potential Energy at O K for EAM Cu
           FUSION=0.00d0        ! Enthalpy of Fusion in eV/atom for EAM Cu
           ballis=0.0d0  ! ballistic energy range
           OPTPEN=14.9d0 ! optical penetration depth, [nm]
!          Thermal energy approximation constants (experimental) 
!          for bulk calculations
           attm=3.5d+06  ! J/m3/K
           bttm=0.0d0    ! J/m3/K2
           TTMBULK = 10000.0d0 ! Min lengh [A] of the ttm part for bulk, approx
           bulkin = 1.2d0! Power of TTM part increase in bulk section
!          Cp_l=3.5D+06, H_fus=1.788D+9, H_abl=4.202D+10 ! experimental
!          T_melt=1357.00D0, T_abl=2840.00D0             ! Not used here

        case (2)   ! The following parameters are for Au [Wellershoff et al.]
           A_e=67.6D0    ! In electron heat capasity
           ek=318.0D0    ! In electron conductivity
           ek_liq=105.0d0! [Mills et al,Int. Mater. Rev. 41, #6, p.209 (1996)]
           G=2.1D+16     ! Coupling constant
!          Reads tabulated data on G_Au calculated by Zhibin Lin with DOS
           IF(mypid.eq.0) THEN
              OPEN(UNIT = 567,FILE = 'G_Ce_R_Au.data')
              READ(567,*) NGT
              READ(567,*) (TeG(i),GGN(i),CeN(i),RTe(i),i=1,NGT)
              CLOSE(UNIT = 567)
           ENDIF
           CALL MPI_BCAST(NGT,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr0)
           CALL MPI_BCAST(TeG,NGT,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr1)
           CALL MPI_BCAST(GGN,NGT,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr2)
           CALL MPI_BCAST(CeN,NGT,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr3)
           CALL MPI_BCAST(RTe,NGT,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr3)
           CALL MPI_BARRIER(MPI_COMM_WORLD,ibar)
           DGR = DBLE(NGT-1)/(TeG(NGT)-TeG(1))
!          Use tabulated G_Au in TTMDiffuse
           Vfer=1.3D+06  ! Fermi velocity
           Ac_e=1.2D+07  ! In electron conductivity
           Bc_e=1.23D+11 ! In electron conductivity
           xappa=353.0d0 ! In electron conductivity
           xnum=0.16d0   ! In electron conductivity
           Efer=5.53d0   ! Fermi energy
           gren=1.6d0    ! Guenizer coefficient, not used currently
           P_ZER=-3.81653d0*NAN ! Potential Energy at O K for EAM Au
           FUSION=0.08259d0 !! Enthalpy of Fusion in eV/atom for EAM Au (0.087 Johns)
           ballis=75.0d0 !Integral value !105.0d0! ballistic range [nm]
           Refl = 0.949d0! Refelctivity at 300 K
           OPTPEN=13.2d0 ! optical penetration dept at 800 nm, [nm]
!          Thermal energy approximation constants (experimental) 
!          for bulk calculations
           attm=2.327679322d+06  ! J/m3/K 
           bttm=0.0d0    ! J/m3/K2
           TTMBULK = 500000.0d0 ! Min lengh [A] of the ttm part for bulk, approx
           bulkin = 1.06d0! Power of TTM part increase in bulk section
!          Cp_l=2.5D+06, H_fus=1.235D+9, H_abl=3.416D+10  ! experimental
!          T_melt=1337.330D0, T_abl=3153.0D0              ! Not used here

        case (3)   !The following parameters are for Ni

           IF(mypid.eq.0) THEN
              OPEN(UNIT = 567,FILE = 'G_Ce_R_Ni.data')
              READ(567,*) NGT
              READ(567,*) (TeG(i),GGN(i),CeN(i),RTe(i),i=1,NGT)
              CLOSE(UNIT = 567)
           ENDIF
           CALL MPI_BCAST(NGT,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr0)
           CALL MPI_BCAST(TeG,NGT,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr1)
           CALL MPI_BCAST(GGN,NGT,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr2)
           CALL MPI_BCAST(CeN,NGT,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr3)
           CALL MPI_BCAST(RTe,NGT,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr3)
           CALL MPI_BARRIER(MPI_COMM_WORLD,ibar)
           DGR = DBLE(NGT-1)/(TeG(NGT)-TeG(1))
           
           A_e=1065.0D0   ! In electron heat capasity
           ek=91.0D0      ! In electron conductivity
           ek_liq=60.0d0  ! [Mills et al.,Int. Mater. Rev. 41, #6,p.209 (1996)]
           G=3.6D+17      ! Coupling constant
           Vfer=2.03D+06  ! Fermi velocity
           Ac_e=1.4D+07   ! In electron conductivity, estimated
           Bc_e=162.24D+11! In electron conductivity, estimated   
           xappa=100.0d0  ! In electron conductivity, estimated
           xnum=0.1d0     ! In electron conductivity, estimated
           Efer=11.72d0   ! Fermi energy
           gren=2.1d0     ! Guenizer coefficient, not used currently
           P_ZER=-4.447276284d0*NAN  ! Potential Energy at O K for EAM Ni
           FUSION=0.1485d0           ! Enthalpy of Fusion in eV/atom for EAM Ni
           ballis=0.0d0   ! ballistic energy range
           OPTPEN=13.5d0  ! optical penetration depth, [nm]
!          Heat capasity approximation constants for model EAM Ni material
!          for bulk calculations. The constant were obtained at equlibrium
!          conditions. Basically, one can distinguish two cases for heat
!          capasity approximation 1) at constant volume, leading to Cv
!          capasity and 2) at constant pressure, leading to Cp capasity.
!          In a real simulations however, the processes that determine
!          thermal energy flow and thus heat capasity are rather complex
!          and can even be result in efective compression near NRB an
!          therefore decresing of heat capasity. It is espesialy gets
!          noticable at high fluences when rise of lattice temperature at
!          the back end of the sample becomes significant resulting in
           A = -1.29065838619645860000E+002
           B =  2.87922030743130270000E-003
           C =  2.24292472502200760000E-006
!          discontinuity in lattice temperature determined form MD and TTM.
!          To avoid such dicontinuity, one may increase the lengh of MD part.
!          You can chose one or the order type of approximation depending on
!          the expected contitions.
!          Cv
           attm = 3.801446251247d+06 ! J/m3/K
           bttm = 1.199079665326d+02 ! J/m3/K2
!          Cp
!          attm = 5.191446805401d+06 ! J/m3/K
!          bttm = 3.523223716278d+03 ! J/m3/K2
           TTMBULK = 2400.0d0     ! Min lengh [A] of the ttm part for bulk
           bulkin = 1.05d0! Power of TTM part increase in bulk section
!          Cp_l=4.1D+06, H_fus=1.658D+9, H_abl=3.727D+10  ! experiment
!          T_melt=1728.15D0, T_abl=3186.15D0              ! Not used here

        case (6)   !The following parameters are for Al
           IF(mypid.eq.0) THEN
              OPEN(UNIT = 567,FILE = 'G_Ce_R_Al.data')
              READ(567,*) NGT
              READ(567,*) (TeG(i),GGN(i),CeN(i),RTe(i),i=1,NGT)
              CLOSE(UNIT = 567)
           ENDIF
           CALL MPI_BCAST(NGT,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr0)
           CALL MPI_BCAST(TeG,NGT,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr1)
           CALL MPI_BCAST(GGN,NGT,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr2)
           CALL MPI_BCAST(CeN,NGT,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr3)
           CALL MPI_BCAST(RTe,NGT,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr3)
           CALL MPI_BARRIER(MPI_COMM_WORLD,ibar)
           DGR = DBLE(NGT-1)/(TeG(NGT)-TeG(1))
           
           A_e=92.0D0    ! In electron heat capasity
           ek=238.0d0    ! In electron conductivity
           ek_liq=91.0d0 ! [Mills et al.,Int. Mater. Rev. 41, #6, p.209 (1996)]
           G=4.9D+17     ! Coupling constant
           Vfer=2.028D+06! Fermi velocity
           Ac_e=1.390D+06  ! In electron conductivity
           Bc_e=5.924D+11  ! In electron conductivity
           xappa=100.0d0 ! In electron conductivity, estimated
           xnum=0.1d0    ! In electron conductivity, estimated
           Efer=11.7d0   ! Fermi energy
           gren=2.1d0    ! Guenizer coefficient, not used currently, estimated
           P_ZER=-3.50d0*NAN! ???????????   Potential Energy at O K for EAM Al
           FUSION=0.015d0   !????????? Enthalpy of Fusion in eV/atom for EAM Al
           ballis=16.0d0  ! ballistic energy range
           OPTPEN=8.0d0 ! optical penetration depth, [nm] not precise
!          Thermal energy approximation constants (experimental) 
!          for bulk calculations
           attm=2.44d+06 ! J/m3/K
           bttm=0.0d0    ! J/m3/K2
           TTMBULK = 300000.0d0 ! Min lengh [A] of the ttm part for bulk
           bulkin = 1.18d0! Power of TTM part increase in bulk section
!          Cp_l=2.42D+06, H_fus=1.068D+9, H_abl=2.94D+10  ! experimental
!          T_melt=933.25D0, T_abl=2740.0D0                ! Not used here

        end select

!       This will get rid off the effect of electrons if needed
        IF (LCOND.EQ.0) G=0.0d0

        CALL MPI_BARRIER(MPI_COMM_WORLD,ipass2)

!       Regularly we start modeling of laser-solid intercation at room 
!       temperature Teminit = 300.0d0. At lower temperatuture, however,
!       a smaller TTM time-step can follow from the stability criteria.
!       Find the current temperature of the system for TIME = 0.0d0 or
!       refer to QTEM as initial temperature if TIME is not 0.0d0. QTEM
!       parameter can be set in md.input input file.

!       Defining preliminary step of discretization
        adifx = DABS(xdiff)*1.0D-10/DBLE(nlcx)               ! [m]
        adify = DABS(ydiff)*1.0D-10/DBLE(nlcy)               ! [m]
        adifz = DABS(zdiff)*1.0D-10/DBLE(nlcz)               ! [m]

        IF(TIME.EQ.0.0d0) THEN

           difst2 = MIN(adifx*adifx,adify*adify,adifz*adifz)

           TEMPTR = 0.0d0
           tempe: DO I=1,NN1
              IF(KHIST(I).GT.10) CYCLE tempe
              TEMPTR = TEMPTR + &
                   0.5*XMASS(1)*Q1D(1,i)*Q1D(1,i)/DELTA/DELTA + &
                   0.5*XMASS(1)*Q1D(2,i)*Q1D(2,i)/DELTA/DELTA + &
                   0.5*XMASS(1)*Q1D(3,i)*Q1D(3,i)/DELTA/DELTA
           ENDDO tempe

           TEMPTR = 2.0d0*TEMPTR*ENUNIT/BK/NDIM/DBLE(NAN-3*NCIRC)

           CALL MPI_REDUCE(TEMPTR,Teminit,1,MPI_DOUBLE_PRECISION, &
                MPI_SUM,0,MPI_COMM_WORLD,ierr)
           CALL MPI_BCAST(Teminit,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

!          Use Neumann stability criteria in 3-D for TTM time-step definition
           Select case (KON)
           case (1)
              dift0 = difst2*A_e*Teminit/ek/DBLE(NDIM)/4.0d0
           case (2)
              dift0 = 1.5d0*difst2*Bc_e*Teminit/(Vfer*Vfer)/DBLE(NDIM)/2.0d0
           case (3)              
              IF(MATER.EQ.2.OR.MATER.EQ.3.OR.MATER.EQ.6) THEN
                 CeF = (Teminit-TeG(1))*DGR + 1.0d0
                 L = CeF
                 IF(L.GE.1.AND.L.LT.NGT) THEN
                    PP  = CeF - L
                    C_e  = CeN(L) + (CeN(L+1)-CeN(L))*PP
                 ELSEIF(L.GE.NGT) THEN
                    C_e = CeN(NGT)
                 ELSE
                    C_e = CeN(1)
                 ENDIF
              ELSE
                 C_e = A_e*Teminit
              ENDIF
              tetae=BK*Teminit/Efer
              difco=xappa*(tetae*tetae+0.16d0)*(5.0d0/4.0d0)*(tetae*tetae+ &
                   0.44d0)*tetae/(SQRT(tetae*tetae+ &
                   0.092)*(tetae*tetae+xnum*tetae))
              dift0 = difst2*C_e/difco/3.0d0 !DBLE(NDIM)/2.0d0
           end select

           dift = DELTA*1.0d-12/DBLE((INT(DELTA*1.0d-12/dift0)+1))
           dift = DELTA*1.0d-12/100.0d0

!           IF(mypid.eq.0) print *, dift,"time TTM"

!          Becouse of the temperature drops due to fluctuations may want to
!          assure for the TTM time step small enough to avoid instabilities
!          dift = dift/2.0d0    ! [s] for TTM part

!          Defining the step of discretization
           difx(0:nlc2-1) = adifx    ! [m]
           dify(0:nlc2-1) = adify    ! [m]
           difz(0:nlc2-1) = adifz    ! [m]

           CALL MPI_BARRIER(MPI_COMM_WORLD,ipass31)
        
!          To simplify information exchage in TTM part, make the mapping size of
!          all cells equal but use actual size in diffusion equation
!          Xsh, Ysh, Zsh in MPI; XED, YED, ZED in TTM equation
           DO iz = 1,nlcz
              DO iy = 1,nlcy
                 DO ix = 1,nlcx
                    ip = ix + nlcx2*(iy + nlcy2*iz)
                    
                    Xsh(ip) = xnmin + (DBLE(ix) - 0.5d0) * adifx * 1.0d+10
                    Ysh(ip) = ynmin + (DBLE(iy) - 0.5d0) * adify * 1.0d+10
                    Zsh(ip) = znmin + (DBLE(iz) - 0.5d0) * adifz * 1.0d+10
                    XED(ip) = xnmin + (DBLE(ix) - 0.5d0) * difx(ip) * 1.0d+10
                    YED(ip) = ynmin + (DBLE(iy) - 0.5d0) * dify(ip) * 1.0d+10
                    ZED(ip) = znmin + (DBLE(iz) - 0.5d0) * difz(ip) * 1.0d+10

                 ENDDO
              ENDDO
           ENDDO

        ENDIF

        IF(TIME.NE.0.0d0) dift = DELTA*1.0d-12/100.0d0
        if(mypid.eq.0) print *,dift,"TTMstep now"

        CALL MPI_BARRIER(MPI_COMM_WORLD,ipass3)

!       Bulk margines to account for NRB fluctuations, found empirically
!       4 link cell layers is usually enough for bulk margines.
        marbulk = 4 ! Used for bulk simulations

!       Perform bulk settings
        IF(KONF.EQ.1) THEN
           IF(TIME.EQ.0.0d0) CALL SetBulk()
           CALL Share_Bulk()
        ENDIF

!       The scaling of linkcells position in geometrical progressino does help
!       for reasonable memory consumption but is still does not help to avoid
!       very big output files. Therefore, we cut significat part of bulk out
!       of output files that do not bear significat temperature rise. To do
!       this remember initial and final position of a link cell for output
!       file formation. You may also decrease/increase the writing part by
!       chosing an appropriate mshift variable
        mshift   = 0
        mzoff    = 1
        iystart = 1!INT((YCENTR-ydist-ylstart)*nysign*nlcyg/YL+1)-mshift
        iyend   = nlcyg!INT((YCENTR+ydist-ylstart)*nysign*nlcyg/YL+1)+mshift
        ixstart = 1
        ixend = nlcxg
        izstart = 1
        izend = nlczg

        IF(mypid.eq.0) THEN
           WRITE(17,*) "Recommended sizes of a rectangular box for triangulation:"
           WRITE(17,*) "X =",ixend-ixstart+1,"Y =",iyend-iystart+1,"Z =",izend
           WRITE(17,*) " "
           print *,"Recommended sizes of a rectangular box for triangulation:"
           print *,"X =",ixend-ixstart+1,"Y =",iyend-iystart+1,"Z =",izend
           print *," " 
           print *,"Initial position of NRB below the surface: ", 0.1d0*ZNRB," nm"
           print *," "
!           call Flushout(6)
        ENDIF

        DO iz = 0,nlcz1
           DO iy = 0,nlcy1
              DO ix = 0,nlcx1
                 ip = ix + nlcx2*(iy + nlcy2*iz)
                 volcell(ip) = difx(ip)*dify(ip)*difz(ip) ! [m^3]
              ENDDO
           ENDDO
        ENDDO

        CALL MPI_BARRIER(MPI_COMM_WORLD,ipass4)

        RMINH=0.2d0 ! relative density for switching on conductivity
        RMINL=0.1d0 ! relative density for switching on node/coupling
        RMINR=0.9d0 ! relative density for switching between MD/TTM in bulk

!       Will deduct the velocity of center of mass for ablation processes only
        IMMOV = (LGOT.NE.2).OR.(ENRIN.EQ.0.0d0.OR.DURAT.EQ.0.0d0)

        F_total  = 0.0d0    ! Absorbed fluence   
        E_pere   = 0.0d0    ! Energy transfered from electrons
        Etotal_e = 0.0d0    ! Total energy of the lectrons
        
        IF(TIME.EQ.0d0) THEN
           T_e(0:nlc2-1) = 0.0d0  ! Electron temperature
           T_l(0:nlc2-1) = 0.0d0  ! Lattice temperture
           Etrans(0:nlc2-1)=0.0d0 ! Energy trnasfered to the lattice from e-ns
           NVK(0:nlc2-1) = 0      ! Initial value of linkcell presence
        ENDIF

        IF(TIME.EQ.0.0d0) IBULK(0:nlc2-1) = .false. ! Bulk/MD marker
        RON(0:nlc2-1)   =  0.0D0  ! Current number of atomes in the cell

        rstep = adifz*1.0d+10

        IF(TIME.EQ.0.0d0) THEN
           fcoef = 2.0d0 ! Coefficient of margins inside for bulk definition
           rdist = zdist - fcoef*marbulk*rstep
        ENDIF

        IF(TIME.EQ.0) THEN
           ZSUR = ZCENTR - ZL
        ELSE
           ZSUR = 0.0d0
        ENDIF
!       Assign density value for each linkcell
        DO iz = 1,nlcz
           DO iy = 1,nlcy
              LED: DO ix = 1,nlcx
                 ip = ix + nlcx2*(iy + nlcy2*iz)
                 mcount = 0
                 
                 i = ltop(ip)
                 IF(i.eq.0) GOTO 221
222              IF(.NOT.((KHIST(i).EQ.3.OR.KHIST(i).GE.10) &
                      .AND.NOTEM)) mcount = mcount + 1
                 i = link(i)
                 IF(i.ne.0)GOTO 222
                 
221              CONTINUE

                 RON(ip) = DBLE(mcount)

!                Account for energy deposition in TTMDIffuse atrating at Z = 0.0d0
                 IF(TIME.EQ.0.0d0.AND.RON(ip).GE.RMINL*RODENO &
                      .AND.ZED(ip).GE.ZSUR) ZSUR = ZED(ip)

!                Definition of link cells belonging to bulk
                 IF(TIME.EQ.0.0d0.AND.KONF.EQ.1) THEN
                    rtest = ZCENTR-ZED(ip)
                    IF(rtest.ge.rdist) IBULK(ip) = .true.

!                    IF(ZED(ip).GT.2.0d0*rstep.OR.ZED(ip) &
!                         .LT.(ZCRIT+0.25d0*rstep)) IBULK(ip) = .false.

                 ENDIF

                 IF(IBULK(ip).OR.RON(ip).GE.RMINL*RODENO) NVK(ip) = 1

              ENDDO LED
           ENDDO
        ENDDO

        CALL MPI_BARRIER(MPI_COMM_WORLD,ipass5)

        IF(TIME.EQ.0.0d0) THEN
           IF(mypid.ne.0) CALL MPI_SEND(ZSUR,1,MPI_DOUBLE_PRECISION, &
                0,0,MPI_COMM_WORLD,ierror)

           IF(mypid.eq.0)  THEN
              ZSURG = ZSUR
              DO nrec = 1,nnodes-1
                 CALL MPI_RECV(ZSUR,1,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,0, &
                      MPI_COMM_WORLD,status,ierror)
                 IF(ZSUR.GE.ZSURG) ZSURG = ZSUR
              ENDDO
           ENDIF
           CALL MPI_BCAST(ZSURG,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
           ZSURG = ZSURG*1.0d-10 ! [M] Shift in Source to the surface
        ENDIF

        IF(KONF.EQ.1) CALL Share_Bulk2()

!       Define some source factors here for TTM
        IF(ballis.GT.0.25d0*WIDTH) THEN
           IF(mypid.eq.0) print *,"The ballistic range is greater than quarter of pulse width."
           IF(mypid.eq.0) print *,"The run is set to test mode. Ballistic range is discarded."
           ballis = 0.001d0*adify
        ENDIF
        OPTPEN  = OPTPEN*1.0D-09        ! [M]
        ballis  = ballis*1.0d-09        ! [M]
!        alfas   = 1.0d0/(OPTPEN+ballis)
        tstart  = 2.5D0*DURAT*1.0D-12   ! [s]
        HDI     = WIDTH*0.5d-09        ! half width [m]
        sigmaL2 = DURAT*DURAT*1.0E-24/(4.0D0*Log(2.0D0))
        SourceL = ENRIN*SQRT(4.0d0*Log(2.0d0)/PI)/(XL*YL*DURAT*1.0d-32)
        fluence = ENRIN/(XL*YL*1.0d-20) ! incident fluence J/m2
        
        IF(mypid.eq.0) THEN
           print *,"MD step = ",DELTA*1.0d+03,"fs"
           print *,"TTM step = ",dift*1.0d+15,"fs"
           print *," "
           print *,"Average density per cell =",RODENO,"atoms"
           print *," "
           print *,"The incident fluence is:", fluence/10.0d0, " mJ/cm2"
           print *," "
!           print *,"The everaged fluence over the MD cell is",(YL/XL)*(WIDTH*WIDTH/(XL*YL/100.0d0))* &
!                (1.0d0 - EXP(-(XL*XL/100.0d0)/(WIDTH*WIDTH)))*fluence/10.0d0," mJ/cm2"
           print *,"The everaged fluence over the MD cell is",fluence/10.0d0," mJ/cm2"
           print *," "
!           call Flushout(6)
        ENDIF

        IF (adifz.GE.OPTPEN) Then
           print *,"dx=,dy=,dz=",adifx,adify,adifz,"is larger than Lp=", &
                OPTPEN,".The program is  halted." 
!           call Flushout(6)
           CALL MPI_FINALIZE(info)
           STOP
        ENDIF


        CALL MPI_BARRIER(MPI_COMM_WORLD,ipass6)

        RETURN
      END Subroutine TTMInit
