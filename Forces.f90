!     Force calculation
!     Leonid Zhigilei, 2001
      SUBROUTINE Forces()
        INCLUDE 'common.h'
        INCLUDE 'mpif.h'

        FD(1:3,1:NNF)=0.0d0
        POT(1:NNF)=0.0d0
        STEN(1:NNF,1:3,1:3)=0.0d0

        IF(KEYBS.EQ.4) THEN
           IF(KBOUND.EQ.2) THEN
              CALL F_EAM_br()
!              CALL F_pair()
           ELSEIF(KBOUND.EQ.11) THEN
              CALL F_EAM_br_test()
           ELSE
              CALL F_EAM()
!              CALL F_pair()
           ENDIF
        ELSEIF(KEYBS.EQ.6) THEN
!           CALL F_pair()
           IF (mypid.EQ.0) print *,"Not implemented"
           CALL MPI_FINALIZE(info)
           STOP
        ELSEIF(KEYBS.EQ.0) THEN
           CALL F_pair()
        ELSE
           IF(mypid.eq.0) THEN
              print *,"The program implemets EAM and pair potentials only"
              print *,"for with KEYBS = 4"
!              call Flushout(6)
           ENDIF
           CALL MPI_FINALIZE(info)
           STOP
        ENDIF

!       Additional forces from TTM
        IF (LGOT.EQ.2) CALL F_TTM()
        
!       Additional forces from matter-substrate interactions, used here
!       specifically for the project with Chichkov
        IF(NDEP.EQ.1) CALL F_SUB()

!       STEN is the static portion of stress tensor x atomic volume
        STEN(1:NN1,1:3,1:3)=STEN(1:NN1,1:3,1:3)*0.5d0

        RETURN
      END SUBROUTINE Forces
