      Subroutine Smooth_Q1V()
        INCLUDE 'common.h'
        INCLUDE 'commonTTM.h'

        REAL*8 Q_a(0:nlc2-1),R_a(0:nlc2-1)
        INTEGER lcellx(13),lcelly(13),lcellz(13)

        DATA lcellx/-1,0,1,-1,0,1,-1,0,1,-1,0,1,-1/
        DATA lcelly/-1,-1,-1,0,0,0,1,1,1,-1,-1,-1,0/
        DATA lcellz/-1,-1,-1,-1,-1,-1,-1,-1,-1,0,0,0,0/

        Q_a(0:nlc2-1) = 0.0d0
        R_a(0:nlc2-1) = 0.0d0

        DO iz = 1,nlcz1
           DO iy = 0,nlcy1
              Internal: DO ix = 0,nlcx1

                 ic = ix + nlcx2*(iy + nlcy2*iz)

                 IF(Q1V(ic).eq.0.0d0) CYCLE Internal

                 IF(ix.eq.0.or.ix.eq.nlcx1.or.iy.eq.0 &
                      .or.iy.eq.nlcy1.or.iz.eq.nlcz1) GOTO 331

                 Q_a(ic) = Q_a(ic) + Q1V(ic)*RON(ic)
                 R_a(ic) = R_a(ic) + RON(ic)
                 
331              Neigh: DO nb = 1,13

                    jx = ix + lcellx(nb)
                    jy = iy + lcelly(nb)
                    jz = iz + lcellz(nb)

                    jc = jx + nlcx2*(jy + nlcy2*jz)

                    IF(Q1V(jc).EQ.0.0d0) CYCLE Neigh

                    IF(jx.ge.0.and.jx.le.nlcx1 &
                         .and.jy.ge.0.and.jy.le.nlcy1) THEN

                       Q_a(ic) = Q_a(ic) + Q1V(jc)*RON(jc)
                       R_a(ic) = R_a(ic) + RON(jc)

                       Q_a(jc) = Q_a(jc) + Q1V(ic)*RON(ic)
                       R_a(jc) = R_a(jc) + RON(ic)
                       
                    ENDIF

                 ENDDO Neigh

              ENDDO Internal
           ENDDO
        ENDDO

        DO iz = 1,nlcz
           DO iy = 1,nlcy
              DO ix = 1,nlcx

                 ic = ix + nlcx2*(iy + nlcy2*iz)

                 IF(Q1V(ic).ne.0.0d0) Q1V(ic) = Q_a(ic)/R_a(ic)
                 
              ENDDO
           ENDDO
        ENDDO

      END Subroutine Smooth_Q1V
      
