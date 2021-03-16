      Program table
        IMPLICIT NONE
        REAL*8 T_e(1000),G_e(1000),C_e(1000),R_e(1000)
        INTEGER N,I
        REAL*8 T,G,Ce

        OPEN(UNIT = 1, FILE = 'G_Ce_R_Al.data')
        OPEN(UNIT = 2, FILE = 'G_Ni.dat')
        OPEN(UNIT = 3, FILE = 'Ce_Ni.dat')
        OPEN(UNIT = 4, FILE = 'G_Ce_R_Ni.data')
        
        READ(1,*) N
        READ(1,*) (T_e(i),G,Ce,R_e(i),I=1,N)
        READ(2,*) (T,G_e(i),i = 1,N)
        READ(3,*) (T,C_e(i),i = 1,N)

        WRITE(4,443) N
        WRITE(4,444) (T_e(i),G_e(i)*1.0d+17,C_e(i)*1.0d+5,R_e(i),I=1,N)

443     FORMAT(I4)
444     FORMAT(4(E16.10,1x))

        CLOSE(UNIT = 1)
        CLOSE(UNIT = 2)
        CLOSE(UNIT = 3)
        CLOSE(UNIT = 4)

      END Program table
