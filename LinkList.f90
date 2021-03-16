      SUBROUTINE LINKLIST()
        INCLUDE 'common.h'

        ltop(0:nlc2-1)=0

        ixoffset = mynodex*nlcx
        iyoffset = mynodey*nlcy
        izoffset = mynodez*nlcz

        DO I = 1,NN1
!          Determine in which linkcell this atom is situated
           ixg = INT((XD(1,I) - xlstart)*nxsign*nlcxg/XL + 1)
           iyg = INT((XD(2,I) - ylstart)*nysign*nlcyg/YL + 1)
           izg = INT((XD(3,I) - zlstart)*nzsign*nlczg/ZL + 1)

           ix = ixg - ixoffset
           iy = iyg - iyoffset
           iz = izg - izoffset

           ip = ix + nlcx2*(iy + nlcy2*iz)

!          Assign atom I to linkcell ip
           J = ltop(ip)
           ltop(ip) = I
           link(I) = J
        ENDDO

        RETURN
      END SUBROUTINE LINKLIST
