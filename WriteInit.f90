!     THIS SUBROUTINE WRITES INITIAL INFORMATIOIN FOR THE RUN
!     Leonid Zhigilei, 2000
      SUBROUTINE WriteInit(DENSM,DENSN)
        INCLUDE 'common.h' 
        INCLUDE 'commonTTM.h' 

        WRITE(17,15) NSTEP,DELTA,KFLAG,LFLAG,KEYBS,KEYDEP,NDEP,LIDX,LIDY, &
             LIDZ,KBOUND,NDIM,QTEM,QPRESS,AddEn,NEPRT,NWRITE,NT
15      FORMAT( 6X,'  -----> MD RUN:   '/ &
             'NSTEP=',I7,' (Number of steps);'/ &
             'DELTA=',F7.5,' (Integration time step in [psec]);',/, &  	
             'KFLAG=',I1,' (1-Quench,2-Vel,3-Heating,4-GLEQ,5-Bere_T,6-AddEn);',/, &
             'LFLAG=',I1,' ((1-Berendsen Const. pressure, 2-shock wave);',/, &
             'KEYBS=',I1,' (0-p.p,1-BS,2-rigid sph,3-BS+Polymer,4-EAM,5-SW);',/, &
             'KEYDEP=',I1,' (Cluster deposition: 0-No, 1-Yes);',/, &
             'NDEP=',I5,' (Step of cluster deposition);',/, &  
             'LIDX=',I1,' (1-periodicity in X direct.,0-free boundaries);'/ &
             'LIDY=',I1,' (1-periodicity in Y direct.,0-free boundaries);'/ &
             'LIDZ=',I1,' (1-Z PBC,0-Z free or boundary layer);'/ &
             'KBOUND=',I2,' (0-Free,1-rigid,2-non-reflecting);'/ &
             'NDIM=',I2,' (3 - 3D siulation,2- 2D -not implemented for BS);'/ &
             'QTEM=',F6.1,' (Temperature used by VEL, HEATING or GLEQ,[K]);'/ &
             'QPRESS=',F6.1,' (Pressure / shock wave (set LGOT to 2) in GPa);'/ &
             'AddEn=',F6.1,' (Energy added by AddEnergy (KFLAG=6) [eV]);'/ &
             'NEPRT=',I4,' (Step of printing output information);'/ &
             'NWRITE=',I5,' (Step of writing output information);'/ &
             'NT=',I4,' (Number of points in the energy & force tables);'/)
             
             WRITE(17,17) XL,YL,ZL,DENSM,DENSN,NAN,NTYPE
17      FORMAT( 6X,'  -----> MATERIAL:   '/ &
             'XL=',D11.5,' (X size of the computational cell [A]);'/ &
             'YL=',D11.5,' (Y size of the computational cell [A]);'/ &
             'ZL=',D11.5,' (Z size of the computational cell [A]);'/ &
             'DENSM=',D11.5,' (Density of the material [gm/cm3]);'/ &
             'DENSN=',D11.5,' (Density of the material [molec./cm3]);'/ &
             'NAN=',I8,' (Number of particles in the computational cell);'/ &
             'NTYPE=',I2,'(Number of particle types);'// &
             6X,'  -----> TYPES OF PARTICLES:   ')
        
        WRITE(17,18) (I,XMASS(I),I=1,NTYPE)
18      FORMAT('Type ',I2,':  XMASS =',D11.5,' [amu].'/)

        IF(LGOT.EQ.1) THEN
           WRITE(17,16) LTYPE,WINT,WTIME,WL,WABS,WUSED,WFLUE,WFLUX,WAREA,VHV
16         FORMAT( 6X,'  -----> LASER:   '/ &
                'LTYPE=',I1,' (1-Laser en. goes to c.o.m. motion,2-vibr.en.);'/ &
                'WINT=',D11.5,' (Intensity(Irradiance)-Power p.un.ar.[W/cm2]);'/ &
                'WTIME=',D11.5,' (Pulse Duration [psec]);'/ &
                'WL=',D11.5,' (Laser Wawelength [nm]);'/ &
                'WABS=',D11.5,' (Reciprocal Abs. coef.(penetr. depth) [nm]);'/ &
                'WUSED=',f4.0,' (En.,transfered to vibr. & c.o.m. motion [%]);'/ &
                'WFLUE=',D11.5,' (Fluence-Pulse en. per unit area [J/cm2]);'/ &
                'WFLUX=',D11.5,' (Flux-Number of phot. per unit area [cm-2]);'/ &
                'WAREA=',f4.0,' (Irradeated area) [%];'/ &
                'VHV=',D11.5,' (Photon Energy) [eV];'/)
        ENDIF

        IF(LGOT.EQ.2) THEN
           WRITE(17,14) MATER,NTST,LCOND,KON,KONF,WIDTH,DURAT,ENRIN,OPTPEN
14         FORMAT( 6X,'  -----> TTM-LASER:   '/ &
                'MATER=',I1,' (Type of material [1 - Cu, 2 - Au, 3 - Ni]);'/ &
                'NTST=',I4,' (Step of printing data for T/P contour plots);'/ &
                'LCOND=',I4,' (e- therm conductin is off/0, on/1);'/ &
                'KON=',I4,' (Type of thermal conductivity approximation);'/ &
                'KONF=',I4,' (0-film, 1-bulk target);'/ &
                'WIDTH=',D11.5,' (Width of the laser strip [nm]);'/ &
                'DURAT=',D11.5,' (Laser pulse duration [psec]);'/ &
                'ENRIN=',D11.5,' (Laser energy J);'/ &
                'OPTPEN=',D11.5,' (Optic penetration depth [nm])'/)
        ENDIF
        
        RETURN
      END SUBROUTINE WriteInit
