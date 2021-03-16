!----- The following parameters are aiming to decrease the memory consumption
!----- The MPI parallelization is realized in such a way that one needs eiher
!----- one or even number of processors in each direction. Because of some
!----- empty space in the case_ of free boundary conditions, the actual density
!----- of the system can be higherand should be accounted in by setting ZSIZE.
!----- Because of density fluctuations, sometimes, trhe number of particles per
!----- node can be significantly higher than its averaged value and that should
!----- be accounted in by setting RMD. In the case_ of TTM-MD, a bigger size of
!----- the link cell leads to a grater number of particles per node with skin
!----- and should be accounted in by setting RTTM. Recommended value for_ RTTM
!----- is 2.0 when dealing with TTM-MD of relatively small size (like 8 linc
!------cells in one of the directions), 1.5 when delaing with TTM-MD/BS of big
!------size and 1.2 when dealing with plain MD.Setting RMD to 1.1 is enough to
!------account for_ most of the situations. So, the variables that one hss to
!------change more less frequently are nxnodes, nynodes, nznodes, LPMX0. If
!------memory is running out, set RTTM to 1.2 when dealing with plain MD and
!------ZSIZE to 1 when dealing with periodic boundary conditions (PBC) in
!------all diretions. Also pay attention to the parameter "mdev" which is
!------responsible for_ finding the maximum number of link cells depending
!------on the cutoff distnace for_ the used potential. MAXNNB is responsible
!------for_ the size of biggest array in the code NNNG(MAXNNB,LMPZ). Therefor,
!------you may want to have MAXNNB as small sa it posible. However,the boundary
!------of such array could be exceeded if_ teh value of MAXNNB is too small.
!------Different system reacts differently on such an oocasion. Particularly
!------piranha does not even care if_ those boundaries wre exceeded, it simply
!------takes the closest free array and keep filling that one. As a result you
!------may get simingly right but erroneouse modeling whothout even knowing
!------about that. To prevent it, you may estimate in advance the maximum
!------possible MAXNNB but in the meantime keep it small enough to gain in the
!------memory consumption. For metals and silicon, MAXNNB = 100 usually works
!------well. For BS model MAXNNB = 150 was found sutiable, but probably may be
!------decreased down to 100 as well. One should chek it on.

      Parameter (nxnodes = 2)            ! number of nodes in x direction
      Parameter (nynodes = 2)            ! number of nodes in y direction
      Parameter (nznodes = 1)            ! number of nodes in z direction
      Parameter (nnodes=nxnodes*nynodes*nznodes)   ! total number of nodes

      Parameter (LPMX0 = 788800)       !Maximum # of particles in zeros node
      Parameter (ZSIZE = 1.0d0)         !Ratio LZ to the actual size of system
      Parameter (RTTM  = 2.5d0)         !Sets the array boundary for_ skin
      Parameter (RMD =   2.2d0)         !Implies the max density fluctuation

      Parameter (LPMX  = INT(RTTM*RMD*ZSIZE*LPMX0/nnodes)) !At/node with skin
      Parameter (LPMX3 = LPMX*3)
      Parameter (LPMZ  = INT(RMD*ZSIZE*LPMX0/nnodes))  !At/node without skin
      Parameter (LPMZ3= LPMZ*3)         !

      Parameter (KWT = 2)           ! TYPE marker of water molecule  

!------------------------Other MD parameters------------------------------
      Parameter (LPhv = 400)           !Maximum number of photons
      Parameter (KTMX = 3)               !Maximum number of particle types
      Parameter (KPMX = KTMX*(KTMX+1)/2) !Maximum number of potential types
      Parameter (MAXNNB = 150)            !Max. number of neighbour pairs
      Parameter (NT = 8192)             !Points in the energy & force tables
      Parameter (NPTTM = 2000)           ! Points of Ce and G of electrons
      Parameter (NSYM  = 15)             ! Number of neighbores in CSP
      Parameter (N_CDF = 448)           !Number of points in the experimental data of Olmon PRB 2012
  
!-----LinkCells parameters------------------------------------------------
      Parameter (mdev =   1)            ! 2 for_ SW; 10 for_ EAM/TTM and BS
      Parameter (mnlc  =  LPMZ/mdev) 
      Parameter (mnlc2 =  LPMX/mdev)
      Parameter (mnlcg =  2*LPMX0*ZSIZE/mdev)
      Parameter (mnlcz =  2*mnlc/mdev)

!---- Fundamental constants ----------------------------------------------
      Parameter (PI = 3.1415926535897932D+00) !
      Parameter (EVTOJOU = 1.60219D-19)       !J/eV
      Parameter (AMUTOKG = 1.6605402D-27)     !kg/amu
      Parameter (BK = 8.617385D-05)           !Boltzman constant, eV/K
      Parameter (XJOUTOEV = 1.d0/EVTOJOU)     !eV/J
      Parameter (CSPEED = 2.99792458D+08)     !speed of light in m/s
      Parameter (HPLANCK = 6.6260755D-34)     !plancks constant in J*s
      Parameter (hcross = 6.582119514d-16)    !eV*sec; Planck constant over 2*pi
      Parameter (wave = 343.0d-09)        ! wavelenght of the laser pulse
      Parameter (omega_fix = 2.0d0*pi*CSPEED/wave) !Laser beam frequency 1/sec 
      Parameter (beta_fix = 0.11651d0)   ! 11651d0)     ! cooupling of the plasmon to the surface 
      Parameter (delta_fix = 0.0d0)           ! Incident angle of the laser pulse

      Parameter (mod_CDF = 4)
	!mod_CDF=1 (paper Vial PRB 2005 Improved analytical fit of gold dispersion, Ep=[1.24 - 2.48 eV])
	!mod_CDF=2 (paper H.S. Sehmi PRB 2017,fit for the paper S. Babar and J. H. Weaver, &
		    !Appl. Opt. 54, 477 (2015) with 4 Lorentz terms, Ep = [0.1 - 6.0 eV],
	!mod_CDF=3 (paper H.S. Sehmi PRB 2017 with a fit to the paper P. B. Johnson and R. W. Christy, &
		    !Phys. Rev. B 6, 4370 (1972) with 3 Lorentz terms and epsilon_infinity Ep = [0.64 - 6.6 eV]
	!mod_CDF=4 (paper R.L. Olmon PRB 2012 , Ep = [0.04974 - 4.133 eV], it means lambda = [300 - 24930 nm]

						       
!---- Transfer to program units ------------------------------------------
      Parameter (ENUNIT = AMUTOKG*1.d4*XJOUTOEV)  !eV/pr.u.
