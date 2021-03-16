!     This subroutine is aimed on imulation of a substrate presence
!     at some position ZCRIT. Will work with Z components only.
      SUBROUTINE F_SUB()
        INCLUDE 'common.h'
        INCLUDE 'commonTTM.h'

        PARAMETER (cigma =  2.28896d0)    !Ni-Ni at 300 K
        PARAMETER (epsil = 0.5197d0/ENUNIT/2.0d0)   !Ni-SiO2

!       Above are approximate LJ parameters for Ni-Silica interaction
!       We take sigma as an equlibrium distance between substrate and
!       metal. Since substrate is immovable here, we merly take the metal-
!       substrate distance equal to the metal-metal equlibrium distance at
!       the initial temperature. The choice of epsil parameter however is
!       more important and can be astimated based on adhesion energy between
!       two materials. Here we take it, by Chichkov's suggestion, to be 10
!       times weaker that for the case of metal-metal interaction. It might
!       not be a perfect choice, but we assume that at least we get the
!       correct order of magnitude.

        izoffset = mynodez*nlcz
        izg = INT((ZCRIT - zlstart)*nzsign*nlczg/ZL + 1)
        izc = izg - izoffset

!        print *,mypid,zlstart,izc,nlcz
        
!       Want to cut-off all cells that are beyond of the interaction range
        IF(izc.GE.1.AND.izc.LE.nlcz) THEN
!           DO iz = 1,nlcz
           DO iz = izc-1,izc
              DO iy = 1,nlcy
                 Loop: DO ix = 1,nlcx
                    ic = ix + nlcx2*(iy + nlcy2*iz)

                    i = ltop(ic)
                 
311                 CONTINUE

                    IF(i.eq.0) CYCLE Loop 
                 
                    DZ = XD(3,i) - ZCRIT
                    
                    IF(ABS(DZ).lt.2.5d0*cigma) THEN
                       DZ6 = (cigma/DZ)**6.0d0

                       FSUB = 48.0d0*epsil*DZ6*(DZ6-0.5d0)/DZ
                       USUB = 4.0d0*epsil*DZ6*(DZ6-1.0d0)

                       FD(3,i) = FD(3,i) + FSUB
                       POT(i)  = POT(i)  + USUB
                       STEN(i,3,3) = STEN(i,3,3) + FSUB*ABS(DZ)
                    ENDIF

                    i = link(i)
                    
                    GOTO 311
         
                 ENDDO Loop
              ENDDO
           ENDDO
        ENDIF

        RETURN
      END Subroutine F_SUB
