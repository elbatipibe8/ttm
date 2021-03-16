!     For the analisys of the results we give the microscopic water variable
!     only for the writing the output data. Following the writing procedure,
!     all these variable must be nullified with Fluid-Void.f90 procedure
      SUBROUTINE Fluid_Set()
        INCLUDE 'common.h'
        INCLUDE 'commonTTM.h'

        RON(0:nlc2 -1) = 0.0d0
        QINK(0:nlc2-1) = 0.0d0
        TPSX(0:nlc2-1) = 0.0d0
        Q1V(0:nlc2 -1) = 0.0d0
        NPR(0:nlc1 -1) = 0
        idlc(0:nlc2-1) = 0
        ippl = 0

        DO iz = 1,nlcz
           DO iy = 1,nlcy
              Loop: DO ix = 1,nlcx
                 ip = ix + nlcx2*(iy + nlcy2*iz)

                 IF(NVK(ip).EQ.1.OR.IBULK(ip)) CYCLE Loop

                 mcount = 0
                 XMASSW = 0.0d0

                 i = ltop(ip)
                 IF(i.eq.0) GOTO 26

27               IF(KTYPE(I).EQ.KWT)THEN
                    mcount = mcount + 1
                    Q1V(ip) = Q1V(ip) + Q1D(3,i)*XMASS(KTYPE(I))/DELTA
                    XMASSW=XMASSW+XMASS(KTYPE(I))
                 ENDIF

                 i = link(i)
                 IF(i.ne.0)GOTO 27

                 IF(mcount.GE.RMINH*DENWAT) THEN
                    ippl = ippl + 1
                    idlc(ippl) = ip
                    RON(ip) = DBLE(mcount)
                    Q1V(ip)= Q1V(ip)/XMASSW
                 ENDIF

26               CONTINUE
              ENDDO Loop
           ENDDO
        ENDDO

        DO ii = 1,ippl
           ip = idlc(ii)    
           i = ltop(ip)

           IF(i.eq.0) GOTO 28

29         TPSX(ip) = TPSX(ip) + &
                (STEN(i,1,1)+STEN(i,2,2)+STEN(i,3,3))*ENUNIT/3.0d0
                    
           QINK(ip)=QINK(ip)+ &
                0.5*XMASS(KTYPE(I))*Q1D(1,i)*Q1D(1,i)/(DELTA*DELTA)+ &
                0.5*XMASS(KTYPE(I))*Q1D(2,i)*Q1D(2,i)/(DELTA*DELTA)+ &
                0.5*XMASS(KTYPE(I))*(Q1D(3,i)/DELTA - &
                Q1V(ip))*(Q1D(3,i)/DELTA-Q1V(ip))


           i = link(i)
           IF(i.ne.0)GOTO 29

           T_l(ip) = QINK(ip)*ENUNIT*2.0d0/BK/DBLE(NDIM)/RON(ip)
           TPSX(ip) = TPSX(ip) + RON(ip)*T_l(ip)*BK
           RON(ip) = RON(ip)*RODENO/DENWAT
                 
28         CONTINUE
        ENDDO
                 
      END SUBROUTINE Fluid_Set
