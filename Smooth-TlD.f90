      Subroutine Smooth_TlD()
        INCLUDE 'common.h'
        INCLUDE 'commonTTM.h'

        REAL*8 T_a(0:nlc2-1),R_a(0:nlc2-1),d_a(0:nlc2-1)
        INTEGER lcellx(13),lcelly(13),lcellz(13)

        DATA lcellx/-1,0,1,-1,0,1,-1,0,1,-1,0,1,-1/
        DATA lcelly/-1,-1,-1,0,0,0,1,1,1,-1,-1,-1,0/
        DATA lcellz/-1,-1,-1,-1,-1,-1,-1,-1,-1,0,0,0,0/

        T_a(0:nlc2-1) = 0.0d0
        d_a(0:nlc2-1) = 0.0d0
        R_a(0:nlc2-1) = 0.0d0

        DO iz = 1,nlcz1
           DO iy = 0,nlcy1
              Inter: DO ix = 0,nlcx1

                 ic = ix + nlcx2*(iy + nlcy2*iz)

                 IF(NVK(ic).eq.0) CYCLE Inter

                 IF(ix.eq.0.or.ix.eq.nlcx1.or.iy.eq.0 &
                      .or.iy.eq.nlcy1.or.iz.eq.nlcz1) GOTO 332
                 
                 IF(.NOT.ICOND(ic)) THEN
                    T_a(ic) = T_a(ic) + T_l(ic)*RON(ic)
                    d_a(ic) = d_a(ic) + rden(ic)*RON(ic)
                    R_a(ic) = R_a(ic) + RON(ic)
                 ENDIF

332              Neigh: DO nb = 1,13

                    jx = ix + lcellx(nb)
                    jy = iy + lcelly(nb)
                    jz = iz + lcellz(nb)

                    jc = jx + nlcx2*(jy + nlcy2*jz)

                    IF(NVK(jc).eq.0) CYCLE Neigh

                    IF(jx.ge.0.and.jx.le.nlcx1 &
                         .and.jy.ge.0.and.jy.le.nlcy1) THEN

                       IF(.NOT.ICOND(jc)) THEN
                          RDENS = RON(jc)
                          T_a(jc) = T_a(jc) + T_l(ic)*RON(ic)
                          d_a(jc) = d_a(jc) + rden(ic)*RON(ic)
                          R_a(jc) = R_a(jc) + RON(ic)
                       ELSE
                          RDENS = RODENO
                       ENDIF   
                       
                       T_a(ic) = T_a(ic) + T_l(jc)*RDENS
                       d_a(ic) = d_a(ic) + rden(jc)*RDENS
                       R_a(ic) = R_a(ic) + RDENS

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
                    T_l(ic)  = T_a(ic)/R_a(ic)
                    rden(ic) = d_a(ic)/R_a(ic)
                 ENDIF
                    
              ENDDO
           ENDDO
        ENDDO

      END Subroutine Smooth_TlD
      
