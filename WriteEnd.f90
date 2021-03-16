!     THIS SUBROUTINE WRITES DATA AT THE END OF THE RUN
!     Leonid Zhigilei, 2000
      SUBROUTINE WriteEnd()
        INCLUDE 'common.h'
        INTEGER, DIMENSION (:), ALLOCATABLE :: GET

!       Because
        CLOSE(UNIT = 16)
        
!       Write random numbers
        REWIND 13
        CALL Random_seed(SIZE = KSEED)
        print *,"KSEED is chosen for md.random",KSEED
        print *," "
!        call Flushout(6)
        ALLOCATE(GET(KSEED))
        CALL Random_seed(GET = ISEED)
        print *,"Random seed is written in md.random",ISEED
        print *," "
!        call Flushout(6)
        WRITE(13,*) (ISEED(i), i=1,KSEED)  ! for random_number(rww1)
        WRITE(13,*) U1,U2                  ! for rn()
        DEALLOCATE(GET)
        
!       Close the files
        CLOSE (UNIT =  7)
        CLOSE (UNIT =  8)
        CLOSE (UNIT =  9)
        CLOSE (UNIT = 10)
        CLOSE (UNIT = 11)
        CLOSE (UNIT = 12)
        CLOSE (UNIT = 13)
        CLOSE (UNIT = 14)
        CLOSE (UNIT = 15)
!        CLOSE (UNIT = 16)
        CLOSE (UNIT = 17)
        CLOSE (UNIT = 18)
        CLOSE (UNIT = 19)
        CLOSE(UNIT = 100)
        CLOSE(UNIT = 234)
        
        IF(KEYBS.EQ.4.AND.(KBOUND.EQ.2.OR.KBOUND.EQ.11)) &
             CLOSE(UNIT = 411)    ! NRB.out
        
        RETURN
      END SUBROUTINE WriteEnd
