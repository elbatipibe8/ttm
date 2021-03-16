!     Writing melting front curve for TTM
      SUBROUTINE Solid()
        Include 'common.h'
        Include 'commonTTM.h'
        REAL*8, DIMENSION(:), ALLOCATABLE :: WS

        ALLOCATE(WS(0:nlc2-1))

!        Call Liquid()  ! calculation of the local order parameter
        Call Symmetry() ! Calculation of the central symmetry parameter

        FSOL(0:nlc2-1)= 0.0d0
        WS(0:nlc2-1)  = 0.0d0
        ZZLIQ=0.0d0

!       Preparation for the averaging inside the processor yet
        DO iz = 1,nlcz
           DO iy = 1,nlcy
              Inter: DO ix = 1,nlcx
                 
                 ic = ix + nlcx2*(iy + nlcy2*iz)
                 
                 IF(.NOT.(IBULK(ic)).AND.RON(ic).EQ.0.0d0) CYCLE Inter

                 i = ltop(ic)
                 
                 IF(i.eq.0) GOTO 26

!!                Determination of the solid part according to LOP
!25               IF(.NOT.((KHIST(i).EQ.3.OR.KHIST(i).GT.10).AND.NOTEM) &
!                      .AND.FD(1,i).GT.CUTLIQ) WS(ic)=WS(ic) + 1.0d0

!                Determination of the solid part according to CSP
25               IF(.NOT.((KHIST(i).EQ.3.OR.KHIST(i).GT.10).AND.NOTEM) &
                      .AND.FD(3,i).LT.CUTLIQ2) WS(ic)=WS(ic) + 1.0d0


                 i = link(i)
                 IF(i.ne.0)GOTO 25
                 
26               CONTINUE

                 IF(ICOND(ic))THEN
                    FSOL(ic) = 1.0d0
                 ELSE
                    FSOL(ic) = WS(ic)/RON(ic)
                    ZZLIQ    = ZZLIQ + RON(ic) - WS(ic)
                 ENDIF

              ENDDO Inter
           ENDDO
        ENDDO
        
        DEALLOCATE(WS)

        RETURN
      END SUBROUTINE Solid
