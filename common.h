      IMPLICIT REAL*8(A-H,O-Z)
      character*14 drname
      LOGICAL ALLPBC,NOTEM
      INCLUDE 'parameters.h'

!     Coordinates, velocities, and forces
      COMMON/XXX/ XD(3,LPMX),Q1D(3,LPMZ),FD(3,LPMX)

!     Energies per atom (total, potential, and kinetic)
!     VW - kinetic energy at the previous step (used in quenching)
      COMMON/ENY/ POT(LPMX),QIN(LPMZ)

!     Force and Energy tables for_ pair potentials
      COMMON/TTT/ FT(KPMX,NT),UT(KPMX,NT)

!     Pressure and stress
      COMMON/PRST/ STEN(LPMX,3,3)

!     Index that defines type of the pair potential for_ two atoms
      COMMON/IND/ IJINDEX(KTMX,KTMX)

!     Cutoff distances for_ pair potentials and neighbor lists
      COMMON/CUT/ DXR(KPMX),RM(KPMX),RList(KPMX),RList2(KPMX),Rskin

!     Parameters for_ EAM
      COMMON/EAM/ RhoM(KTMX),RhoDen(LPMX)
      COMMON/EAM2/ DFT(KTMX,NT),DDFT(KTMX,NT),EFT(KTMX,NT),DEFT(KTMX,NT)

!     Neighbor lists arrays
      COMMON/NBR/ NNG(LPMZ),NNNG(MAXNNB,LPMZ)

!     Sizes of the computational cell
      COMMON/CELL/ XLHALF,YLHALF,ZLHALF

!     Characteristics of the atoms
      COMMON/KTP/ KHIST(LPMX),KTYPE(LPMX),NMETG,NWATG

!     Order parameter for_ identification of liquid
      COMMON/LIQ/ CUTLIQ,CUTLIQ2,alat,ExpCoef

!     Masses of the atoms
      COMMON/MSS/ XMASS(KTMX),G1(KTMX)

!     Higher time derivatives of the coordinates and parameters for
!     Nordsieck integrator
      COMMON/NORD/ Q2(LPMZ3),Q3(LPMZ3),Q4(LPMZ3),Q5(LPMZ3),CC1,C2,C3,C4,C5

      COMMON/NNN/ NN1,NNF,NN13,NCIRCr,NCIRCl,NCIRC,NPOTS,ISTEP,NLREN,NOTEM
      COMMON/LSR1/ NhvCr,NhvTot
      COMMON/LSR2/ NNGOT(LPhv)
      COMMON/NBC/ EOUT,EOUTG,A_TEM,A_Press,A_XL,POTTG
      COMMON/NMNLAY/ NTTPR,Nmono,nrbsign,nbssign

      COMMON/NRBC/ Ponrefl(5),Ponbs,ZNRB

      PARAMETER(IAHEAD   = 20)
      PARAMETER(ICOME    =  5)
      PARAMETER(IWRITE   =  3)
      PARAMETER(IWRITE2  =  4)
      PARAMETER(IBS      =  6)
      PARAMETER(ISWBS    =  3)

      COMMON/MSGL/ MSGMIG,MSGSHA,MSGSWR,MSGSWR2,MSGSWP,MSGSWBS

      COMMON/MESP/ IMIGRATE,ISHARE,ISWRITE,ISWRITE2,ALLPBC

!############################################################################
!     Input parameters for_ MPI
!     if_ you want to broadcast these variables, add them here,
!     otherwise chose another place for_ variables adding 
      PARAMETER (NINTEGER = 17 )
      COMMON/NNM/ NSTEP,MATER,NEPRT,NWRITE,NTYPE,KFLAG,LFLAG, &
      KEYBS,KEYDEP,NDEP,LIDX,LIDY,LIDZ,KBOUND,LGOT,NDIM,NAN

      PARAMETER (NREAL = 12 )
      COMMON/XNN/ TIME,QTEM,QPRESS,DELTA,AddEn,XL,YL,ZL, &
      XCENTR,YCENTR,ZCENTR,ZCRIT

      PARAMETER (NNONREFL = 2)
      COMMON/FNN/ FSTART,GIMPED

      PARAMETER (LASBI = 1)
      COMMON/LAS1/ LTYPE

      PARAMETER (LASBR = 10)
      COMMON/LAS2/ WFLUE,WTIME,WL,WABS,WAREA,WUSED,CHFR,WINT,VHV,WFLUX

      PARAMETER (NINTARR2 = LPMX0)

      PARAMETER (NPOSARR = 3*LPMX0)

      PARAMETER (NVELARR = 3*LPMX0)

      PARAMETER (NLCELL = 3)
      COMMON/CELLBCAST/ nlcx,nlcy,nlcz

      PARAMETER (NINDIM = 3 )
      COMMON/DIMBCAST/ XLREAL,YLREAL,ZLREAL,XLREALW,YLREALW,ZLREALW
!############################################################################

!     Variables related to MPI version of the CODE
      COMMON/STNODES/ xlstart,xdiff,xnmin,xnmax,ylstart,ydiff,ynmin,ynmax, &
      zlstart,zdiff,znmin,znmax

      COMMON/NSIGVAR/ mypid,mynodex,mynodey,mynodez,nxsign,nysign,nzsign, &
      nlong,nmidl,nshort,longnode,midlnode,kornode

      COMMON/NBNODES/ kscell(26,mnlc2),kbcell(26,mnlc),nscell(26),nbcell(26), &
      nbpid(26),nbres(26),ncom(26),imstart(0:nnodes-1,1:26),longarray(18), &
      midlarray(6),korarray(2)

      COMMON/XNODESET/ ximage(26),yimage(26),zimage(26)

      COMMON/LINKCELL/ ltop(0:mnlc2-1),link(LPMX)

      COMMON/CELLLIST/ nlcx1,nlcy1,nlcz1,nlcx2,nlcy2,nlcz2, &
      nlc,nlc2,nlcxg,nlcyg,nlczg,nlcg,marbulk,ICASE

!     Variables used to define names of input/output files
      COMMON/dr1/ LDR
      COMMON/dr2/ DRNAME

!     Seed for_ random number generator
      COMMON/RNDN/ ISEED(64), U1, U2


