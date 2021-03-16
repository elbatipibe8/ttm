      LOGICAL IMMOV,FLUID,IBULK,ICOND
      REAL*8 koeff0,koeff1,koeff2,koeff3,koeff4,koeff5,koeff6,koeff7, &
      koeff8,koeff9,koeff10,koeff11,koeff12,koeff13,koeff14,koeff15,koeff16, &
      koeff11_1,koeff11_2,koeff11_3

!      COMPLEX*16 CDF_m,kx,kz_m

!      COMMON/COMPL/CDF_m,kx,kz_m
	
      COMMON/NPRE/ FLUID,ippl,NPR(0:mnlc2),NVK(0:mnlc2),IBULK(0:mnlc2), &
      ICOND(0:mnlc2),idlc(0:mnlc2),NWT(0:mnlc2)

      COMMON/XYZDEP/ XED(0:mnlc2),YED(0:mnlc2),ZED(0:mnlc2), &
      Xsh(0:mnlc2),Ysh(0:mnlc2),Zsh(0:mnlc2),zpov(0:mnlc2)

      PARAMETER(ITTM      = 8)
      PARAMETER(ITTMS     = 4)
      PARAMETER(ITTMD     = 5)
      PARAMETER(ITTMDSL   = 6)
      PARAMETER(ITTMR     = 9)
      PARAMETER(ITTM2     = 9)
      PARAMETER(ITTM3     = 6)
      PARAMETER(ITTMBULK  = 6)
      PARAMETER(ITTMBULK2 = 4)

      COMMON/AAA/ tstart,sigmaL2,sigmaXY,SourceL,alfas,RODENO,DENWAT, &
      Etotal_e,F_total,ZZLIQ,ETTM,attm,bttm,F_totalG,HDI,ZSURG, &
      E_pere,E_pere0,F_total0,Etotal_e0,Etotal_eG,E_pereG,DGR,Fluence

      COMMON/DIF/ T_e(0:mnlc2),T_l(0:mnlc2),QINK(0:mnlc2), &
      Etrans(0:mnlc2),Q1V(0:mnlc2),RON(0:mnlc2),FSOL(0:mnlc2), &
      TPSX(0:mnlc2)

      COMMON/TGC/ TeG(NPTTM),GGN(NPTTM),CeN(NPTTM),RTe(NPTTM),EGG(NPTTM)

      COMMON/MSGTT/ MSGTTM,MSGTTM2,MSGTTM3,MSGTTMD,MSGSW2T,MSGSW3T,MSGSTTM, &
      MSGTTMB,MSGTTMB2,MSGSTTM2,MSGTTMDSL,NGT

      COMMON/NSPEC/ nodespec,NFIRST,NMIDLE,NLAST,NSH,Nmonosh,nsign,IMMOV, &
      ixstart,ixend,iystart,iyend,izstart,izend

      COMMON/TTMCONST/ A_e,ek,ek_liq,G,RMINH,RMINL,RMINR,TTMBULK,Teminit, &
      bulkin,Ac_e,Bc_e,Vfer,xappa,xnum,Efer,P_ZER,FUSION,OPTPEN,ydist,ballis

      COMMON/DIFFST/ dift,difx(0:mnlc2),dify(0:mnlc2),rden(0:mnlc2), &
      difz(0:mnlc2),volcell(0:mnlc2),adifx,adify,adifz,rstep,zdist, &
      vden,conds(0:mnlc2),condl(0:mnlc2),t6

      COMMON/MODU/ mod_opt,mod_fun,mod_m,mod_d,mod_source

      COMMON/ffs/ f1,f2,f7,f8,f9,f10,f11,f12,f13,f14
      COMMON/phiphi/ phi1,phi2,phi7,phi8,phi9,phi10,phi11,phi12,phi13,phi14, &
      phi18,phi19,phi20,phi21,phi22,phi23

      COMMON/xkoeff/ koeff0,koeff1,koeff2,koeff3,koeff4,koeff5,koeff6,koeff7, &
      koeff8,koeff9,koeff10,koeff11,koeff12,koeff13,koeff14,koeff15,koeff16, &
      koeff11_1,koeff11_2,koeff11_3,vgr

      COMMON/Ecoeff/ En_CDF_tab(N_CDF),omega_CDF_tab(N_CDF),e1_tab(N_CDF), &
      e2_tab(N_CDF),b_e1_koef(N_CDF),c_e1_koef(N_CDF),d_e1_koef(N_CDF), &
      b_e2_koef(N_CDF),c_e2_koef(N_CDF),d_e2_koef(N_CDF)
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Input parameters
      PARAMETER (NTTMIN = 4)
      COMMON/NINTTM/ NTST,LCOND,KON,KONF
      PARAMETER (NTTMDD = 3)
      COMMON/DDD/ WIDTH,DURAT,ENRIN
