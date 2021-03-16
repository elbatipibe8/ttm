!     *************************************************************
!     ** Molecular Dynamic Code with Electronic Heat Conduction  **
!     ** based on TTM treatment of electrons                     **
!     *************************************************************
!     ** Version 1.0 by Leonid V. Zhigilei, Fall 2000            **
!     ** Supports multi-component systems,                       **
!     ** periodic, free, and rigid boundary conditions,          **
!     ** Breathing Spheres added 25 November 2000                **
!     ** Cleaned up and SMP added - August 2001                  **
!     ** LJ for multiple materials added - November 2001         **
!     ** Masses and cutoff distances are defined in Eftab        **
!     Copyright (c) 2000-2004 Computational Materials Group, UVa **
!     Leonid Zhigilei, http://www.faculty.virginia.edu/CompMat/  **
!     *************************************************************
       PROGRAM MD
         INCLUDE 'common.h' 
         INCLUDE 'commonTTM.h'

!        Include the following internal fortran module whenever MPI is involved

         INCLUDE 'mpif.h'
         INTEGER HDONE,HDONEI,EDONE
         CHARACTER*8 openf
!**********************************************************************
!    UNITS:
!    LENGTH      - 1A
!    TIME        - 1psec
!    MASS        - 1amu = 1.66057D-27 Kg
!    TEMPERATURE - K
!    ENERGY      - 1.0364381D-04 eV = 1.66057D-23 J
!
!**********************************************************************
! All input/output files that you use should be listed in file md.rc **
! UNIT 7  | PARAMETERS FOR BOUNDARY CONDITION             |READ
! UNIT 8  | CLUSTER FOR CLUSTER DEPOSITION                |READ
! UNIT 10 | NAMES OF INPUT/OUTPUT FILES                   |READ
! UNIT 11 | PARAMETERS FOR TTM MODEL                      |READ
! UNIT 12 | EXTERNAL IMPACT PARAMETERS                    |READ
! UNIT 13 | SEEDS FOR THE RUNDOM NUMBERS GENERATOR        |READ & WRITE
! UNIT 14 | INPUT DATA (THE MATERIAL AND RUN PARAMETERS)  |READ
! UNIT 15 | INPUT COORDINATES AND VELOCITIES              |READ
! UNIT 16 | OUTPUT COORDINATES AND VELOCITIES             |       WRITE
! UNIT 17 | OUTPUT FILE (ENERGY PER ATOM,TEMP.,TIME, etc.)|       WRITE
! UNIT 18 | FILE FOR MAKING SNAPSHOTS (COORDINATES)       |       WRITE
! UNIT 111| FORCE & ENERGY TABLES (ONLY FOR TESTING)      |       WRITE
!**********************************************************************
!  NAN is total number of particles
!  Ntype is the number of particle types
!  Npots=Ntype*(Ntype+1)/2 is number of types of particle pairs
!  For example, for Ntype=3, Npots=6 and for molecules I of type Ktype(I)
!  and J of type Ktype(J) IJindex is defined as
!
!           KTYPE(J) = 1:Ntype
!
!          | 1  2  3
!      K  ----------
!      T  1| 1  3  6
!      Y   |           IJINDEX(KTYPE(I),KTYPE(J)) = 1:Npots
!      P  2| 3  2  5
!      E   |
!     (I) 3| 6  5  4
!
!**********************************************************************
!  HISTORY VARIABLE DEFINES PROPERTIES OF THE PARTICLE OTHER THAN ITS TYPE
!
!     KHIST(J)=1 - full dynamics
!     KHIST(J)=2 - temperature control
!     KHIST(J)=3 - rigid, non-reflecting boundary
!     KHIST(J)>10 - non-refl. bound. or shock with EAM (by monolayers)
!
!**********************************************************************
!  KFLAG DEFINES TEMPERATURE CONTROL
!  KFLAG= 1 -Quench, 2-Vel, 3-Heeating, 4-GLEQ, 5-Bere_T, 6-AddEnergy 
!**********************************************************************
!  KEYBS = 0 -pair potential, 1 -breathing spheres, 2 -rigid spheres,
!          3 -breathing spheres with polymer molecules, 4 -EAM
!**********************************************************************
!  LGOT  = 0 -Impact/laser is switched off, 1-on, 2-TTM laser
!**********************************************************************
!        Initialize MPI

         CALL MPI_INIT(ierr)
         CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nprocess,ierr)
         CALL MPI_COMM_RANK(MPI_COMM_WORLD,mypid,ierr)
!         CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nprocess,ierr)

!        Initial values of MD step, used for auxiliary settings
         ISTEP = 0

!        Checking the number of processors requested and used in compilation
         IF(nprocess.ne.nnodes)THEN
            IF(mypid.eq.o) THEN
               print *,"The job was compiled for ",nnodes,"processors"
               print *,"but was submitted for",nprocess
               print *,"Check the geometry. The job is halted."
            ENDIF
            CALL MPI_FINALIZE(info)
            STOP
         ENDIF
         print *, "Process",mypid,"of",nprocess,"is active"
!         call Flushout(6)

!        Open and read the namelist datafile
         IF(mypid.eq.0) THEN
            openf = "md.rc" 

            CALL OpenFiles(openf)
            
!           Processor with ID "0" reads all the information
            CALL ReadFiles()     ! reading input files

            IF(LGOT.NE.2) THEN
               OPEN (UNIT = 100,FILE='energy.out')
               OPEN (UNIT = 234,FILE='liquid.out')
            ENDIF

            IF(KEYBS.EQ.4.AND.(KBOUND.EQ.2.OR.KBOUND.EQ.11)) &
                 OPEN(UNIT = 411,FILE = 'NRB.out')
         ENDIF

!        Processor with ID "0" broadcasts all the input variables except 
!        the velocity related ones to all the processors
         CALL BCAST_POS()
         
!        Create energies & forces tables
!        Pair potentials and forces are tabulated
         IF (KEYBS.EQ.4) Then
            CALL EF1_EAM(1,MATER,MATER)  ! Cu(1),Au(2),Ni(3),Fe(4),Zr(5),Al(6)
            IF(NTYPE.GE.2) THEN
               CALL EF1LJ(2,20,20)
               CALL EF1LJ(3,4,20)
            ENDIF
         ELSEIF (KEYBS.EQ.0) THEN
            CALL EF1LJ(1,20,20)
         ELSE
            Print *, "KEYBS=",KEYBS," is not defined. Program will stop."
            Stop
         ENDIF
         
!         CALL EFTAB_TI()
         
         CALL SetInit()    ! defining some initial parameters

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        Assign each particle to related processor

         IF(mypid.eq.0) CALL CPU_TIME(t_ini)
         
         CALL FirstList(DENSM,DENSN)   
         
         IF(mypid.eq.0) THEN
            CALL CPU_TIME(t_fi)
            print *,">>>--->>>--->>>   Inicialization Time:", t_fi-t_ini
            print *," "
!            call Flushout(6)
         ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!        Perform LinkList procedure
         CALL LINKLIST()

!        In this program we have a lot of communications involved by mean
!        of MPI. Each processor communicates with its 26 neigbors and for
!        each sending/receive pair we cpecify both the destination/source
!        and tag for the message. Additionally, to avoid any misdestination
!        in MPI a blocking version of MPI was implemented. Each meassgae has
!        its tag that we specify in such a way so that there is no overlapping
!        between them would be possible. For example, in subroutine Migrate
!        the tag of messages changes from 1 to 26. In Share, however, that tag
!        runs from 27 to 52 and so on. For very short information exchange
!        we use "0" tag, like in nonreflective boundary implementation or
!        in shocks. Every time we do it, we must be sure that those messages
!        will not cross over. In the case of nonreflective boundaries for BS
!        the additional usage of MPI can be seen in NordBS. Despite the fact
!        that those messages do not time-consuming at all, they are very
!        intense and, therefore, we apply different tags "53" and "54" that
!        should not hitch up any other message that could be possible involved
!        into BS calculations.

         CALL MIGRATE()
         
!        Take positions of particles from the skin layer to calculate NbList
         CALL SHARE()

!        13 neighboring link cells are used to construct full NbList 
         CALL NbList_MPI()
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         ZCRIT = ZCRIT + 0.001d0*80.0d0
         
         GOTO 456
         
         IF(mypid.eq.0) open(unit = 300, file = 'testenergy.out')
         xdelta = 0.001d0
         CALL FORCES()
         CALL collect()
         CALL MPI_BCAST(POTTG,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
         POTT0 = POTTG
         Loop_outer: DO J = 1,1000
            ZCRIT = ZCRIT + xdelta
            CALL FORCES()
            CALL COLLECT()
            CALL MPI_BCAST(POTTG,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
               IF(POTTG.GT.POTT0) THEN
                  print *, POTTG,POTT0
                  exit Loop_outer
                  CALL MPI_FINALIZE(info)
                  STOP
               ELSE
                  POTT0 = POTTG
                  IF(mypid.eq.0) THEN
                     WRITE(300,300) J,POTTG,J*delta
                     print *, J, J*delta,POTTG
300                  FORMAT(I6,1x,E20.10,1x,E20.10)
                  ENDIF
               ENDIF
            ENDDO Loop_outer
            CALL MPI_FINALIZE(info)
            STOP
456         CONTINUE
                     
!        Initialization of the laser pulse and the TTM model
!        Sample was shifted in ReadFiles so that the surface is at z = 0
!        Before starting laser time should be set to zero
         IF(LGOT.EQ.2) CALL SetTTM()

         IF(mypid.eq.0) THEN
            CALL WriteInit(DENSM,DENSN)  ! writing initial output before MD
            CALL Flushout(17)
         ENDIF

         NN13 = NN1*3            ! for 1-D array of Q
         A_TEM = 0.0d0           ! temperature averaging
         A_Press =0.0d0          ! pressure averaging
         A_XL = 0.0d0
         NTTPR = 0               ! step of T/P averaging
        
         IF(mypid.eq.0) CALL CPU_TIME(timestart)
         nntime = 0
         ttt = 0.0d0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        Beginning of MD loop
         main_loop: DO WHILE (ISTEP.LT.NSTEP)
            ISTEP = ISTEP + 1
            TIME=TIME+DELTA
            
!           For synchronization and timimg
            CALL MPI_BARRIER(MPI_COMM_WORLD,ipass)
            IF(mypid.eq.0) THEN
               print *,"Time = ",TIME,"ps"
               CALL CPU_TIME(timet)
               print *,"last time step",timet-time00," sec"
               CALL CPU_TIME(time00)
               print *,time00-timestart," sec Total time"
!               call Flushout(6)
            ENDIF
            
            CALL Nord5()


            IF(.NOT.ALLPBC.AND.MOD(ISTEP,NEPRT).EQ.0) THEN
               CALL MPI_REDUCE(NN1,NAN,1,MPI_INTEGER,MPI_SUM, &
                    0,MPI_COMM_WORLD,ierr)
               CALL MPI_BCAST(NAN,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            ENDIF
            
            IF(LFLAG.EQ.1) THEN
               CALL Pressure(TPRESS)
               CALL BERE_P(TPRESS)
            ENDIF
            IF(KFLAG.EQ.5) THEN
               CALL Temper(QINT,POTT,TEMPTR)
               CALL BERE_T(TEMPTR)
            ENDIF

!           Printing output information
            IF (MOD(ISTEP,NEPRT).EQ.0.OR.ISTEP.EQ.1) CALL Collect()

!           Writing melting front curve for TTM
            IF(LGOT.EQ.2.AND.MOD(ISTEP,NTST).EQ.0) CALL Solid()

!           Writing data for future analysis
            IF(MOD(ISTEP,NWRITE).EQ.0.OR.ISTEP.EQ.NSTEP.OR.ISTEP.EQ.1) Then

               CALL CheckWrite()

               IF (LGOT.EQ.2) THEN
                  IF(MOD(ISTEP,2000).EQ.0.OR.ISTEP.EQ.NSTEP.AND.NSTEP.NE.1) THEN
                     CALL SwriteTTM_Res()
                     NFSTEP = ISTEP
                  ELSE
!                     IF(ISTEP.NE.1) CALL SwriteTTM()
                     CALL SwriteTTM()
                  ENDIF
               ELSE
                  IF(MOD(ISTEP,2000).EQ.0.OR.ISTEP.EQ.NSTEP.AND.NSTEP.NE.1) THEN
                     NFSTEP = ISTEP
                     CALL Swrite_Res()
                  ELSE
!                     IF(ISTEP.NE.1) CALL Swrite()
	  	       CALL Swrite()
                  ENDIF
               ENDIF

            ENDIF

         END DO main_loop
!        End of MD loop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         IF(mypid.eq.0) THEN

            CALL WriteEnd()

            CALL CPU_TIME(timeend)

            print *,"Average Temperature=",A_TEM
            print *,"Average Pressure=",A_Press
            print *,"Average XL =", A_XL
            print *,"Total CPU time =",timeend-timestart
            print *," "
!            call Flushout(6)

         ENDIF
         
         CALL MPI_BARRIER(MPI_COMM_WORLD,imhe)
         
!        Finalize all MPI stuff
         CALL MPI_FINALIZE(info)

         STOP
       END PROGRAM MD
