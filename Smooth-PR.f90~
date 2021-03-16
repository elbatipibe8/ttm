      Subroutine Smooth_PR()
        INCLUDE 'common.h'
        INCLUDE 'commonTTM.h'
        
        REAL*8 P_a(0:nlc2-1),R_a(0:nlc2-1),dcount(0:nlc2-1)
        INTEGER lcellx(13),lcelly(13),lcellz(13)
        
        DATA lcellx/-1,0,1,-1,0,1,-1,0,1,-1,0,1,-1/
        DATA lcelly/-1,-1,-1,0,0,0,1,1,1,-1,-1,-1,0/
        DATA lcellz/-1,-1,-1,-1,-1,-1,-1,-1,-1,0,0,0,0/

        P_a(0:nlc2-1)   = 0.0d0
        R_a(0:nlc2-1)   = 0.0d0
        dcount(0:nlc2-1)= 0.0d0

        
        DO iz = 1,nlcz1
           DO iy = 0,nlcy1
              Inter: DO ix = 0,nlcx1

                 ic = ix + nlcx2*(iy + nlcy2*iz)

                 IF(NVK(ic).eq.0) CYCLE Inter

                 IF(ix.eq.0.or.ix.eq.nlcx1.or.iy.eq.0 &
                      .or.iy.eq.nlcy1.or.iz.eq.nlcz1) GOTO 333

                 IF(.NOT.ICOND(ic)) THEN
                    P_a(ic)   = P_a(ic) + TPSX(ic)*RON(ic)
                    R_a(ic)   = R_a(ic) + RON(ic)
                    dcount(ic)= dcount(ic) + 1.0d0
                 ENDIF

333              Neigh: DO nb = 1,13

                    jx = ix + lcellx(nb)
                    jy = iy + lcelly(nb)
                    jz = iz + lcellz(nb)

                    jc = jx + nlcx2*(jy + nlcy2*jz)

                    IF(NVK(jc).eq.0) CYCLE Neigh

                    IF(jx.ge.0.and.jx.le.nlcx1 &
                         .and.jy.ge.0.and.jy.le.nlcy1) THEN

                       IF(.NOT.ICOND(jc)) THEN
                          P_a(ic) = P_a(ic) + TPSX(jc)*RON(jc)
                          R_a(ic) = R_a(ic) + RON(jc)
                          dcount(ic)= dcount(ic) + 1.0d0
                          P_a(jc) = P_a(jc) + TPSX(ic)*RON(ic)
                          R_a(jc) = R_a(jc) + RON(ic)
                          dcount(jc)= dcount(jc) + 1.0d0
                       ENDIF

                    ENDIF

                 ENDDO Neigh

              ENDDO Inter
           ENDDO
        ENDDO

        DO iz = 1,nlcz
           DO iy = 1,nlcy
              DO ix = 1,nlcx

                 ic = ix + nlcx2*(iy + nlcy2*iz)

                 IF(.NOT.ICOND(ic).AND.NVK(ic).EQ.1) THEN
                    TPSX(ic)=P_a(ic)/R_a(ic)
                    RON(ic) = R_a(ic)/dcount(ic)
                 ENDIF
                    
              ENDDO
           ENDDO
        ENDDO

      END Subroutine Smooth_PR
      
