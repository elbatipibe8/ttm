!     THIS SUBROUTINE READS MOST OF THE INPUT FILES
!     Leonid Zhigilei, 2000
      SUBROUTINE ReadFiles()
        INCLUDE 'common.h' 
        INCLUDE 'commonTTM.h'
        INTEGER, DIMENSION (:), ALLOCATABLE :: PUT
 
!       READ INPUT DATA
        REWIND 14
        READ(14,100) NSTEP,MATER,NEPRT,NWRITE,NTYPE,KFLAG,LFLAG, &
             KEYBS,KEYDEP,NDEP,LIDX,LIDY,LIDZ,KBOUND,LGOT,NDIM
        READ(14,200) QTEM,QPRESS,AddEn,DELTA
        
        IF(LGOT.EQ.2) THEN 
           REWIND 11
           READ(11,100) NTST,LCOND,KON,KONF
           READ(11,200) WIDTH,DURAT,ENRIN
        ENDIF
        
        REWIND 15
        READ(15,*) NN1,TIME,ZCRIT
        READ(15,*) XL,YL,ZL
        READ(15,*) XCENTR,YCENTR,ZCENTR
        
!       READ PARAMETERS FOR BOUNDARY CONDITION
        IF(KBOUND.EQ.2) THEN 
           REWIND 7
           READ(7,200) FSTART,GIMPED
        ENDIF
        
!       READ RANDOM NUMBERS FOR VEL SUBROUTINE
!        CALL Random_seed(SIZE = KSEED)
!        REWIND 13
!        READ(13,*)(ISEED(i),i=1,KSEED)    ! for random_number(rww1)
!        READ(13,*) U1,U2                  ! for rn()
!        ALLOCATE(PUT(KSEED))
!        CALL Random_seed(PUT = ISEED)
!        DEALLOCATE(PUT)
        
100     FORMAT (I10)
200     FORMAT (D15.5)

        RETURN
      END SUBROUTINE ReadFiles
