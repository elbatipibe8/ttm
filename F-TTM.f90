      SUBROUTINE F_TTM()
        INCLUDE 'common.h'
        INCLUDE 'commonTTM.h'

        DO iz = 1,nlcz
           DO iy = 1,nlcy
              Inter: DO ix = 1,nlcx
                 ic = ix + nlcx2*(iy + nlcy2*iz)
                 
                 IF(NVK(ic).EQ.0.OR.ICOND(ic)) CYCLE Inter

                 QINK(ic)=0.5d0*T_l(ic)*BK*DBLE(NDIM)*RON(ic)*EVTOJOU/volcell(ic) ![J/m^3]
                 
                 COEFX=0.5d0*Etrans(ic)/QINK(ic)
                 
                 j = ltop(ic)
                 
111              CONTINUE
                 
                 IF(KHIST(J).LT.10.AND.KTYPE(J).NE.KWT) THEN
                    FD(1,J)=FD(1,J)-COEFX*XMASS(KTYPE(J))*Q1D(1,J)/DELTA
                    FD(2,J)=FD(2,J)-COEFX*XMASS(KTYPE(J))*Q1D(2,J)/DELTA
                    FD(3,J)=FD(3,J)-COEFX*XMASS(KTYPE(J))*(Q1D(3,J)/DELTA- &
                         Q1V(ic))
                 ENDIF

                 j = link(j)
                 
                 IF(j.ne.0) GOTO 111

              ENDDO Inter
           ENDDO
        ENDDO
        
        RETURN
      END SUBROUTINE F_TTM
