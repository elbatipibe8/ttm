      Program pot
        INTEGER, PARAMETER :: N = 8192
        REAL*8 V(N),dvdr(N),den(N),dndr(N),F(N),dfdn(N)
        REAL*8 r(N),rho(N)
        INTEGER I

        OPEN(UNIT = 1,FILE = 'Vral.dat')
        OPEN(UNIT = 2,FILE = 'nral.dat')
        OPEN(UNIT = 3,FILE = 'Fnal.dat')
        OPEN(UNIT = 4,FILE = 'allu.data')

        DO I =1,N
           READ(1,*) r(I),v(I),dvdr(I)
           READ(2,*) r(I),den(I),dndr(I)
           READ(3,*) rho(I),F(I),dfdn(I)
           IF(I.NE.0) WRITE(4,444) r(I),V(I),dvdr(I),den(I),dndr(I), &
                rho(I),F(I),dfdn(I)
        ENDDO
        
444     FORMAT(8(E24.16,1x))
      END Program pot
