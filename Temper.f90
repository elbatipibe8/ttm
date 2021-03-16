!     Calculation of temperature and averaged energies
!     Energies in eV/particle, Temperature in K.
!     You can modify this to calculate temperature 
!     for atoms with KHIST(J)=2 only.
!     Pressure calculation added, April 2002
!     Leonid Zhigilei, 2000
      SUBROUTINE TEMPER(QINT,POTT,TEMPTR)
        INCLUDE 'common.h'
        
        POTT=0.0d0
        QINT=0.0d0                         
        ent_loop: DO I=1,NN1
           POTT=POTT+POT(I)
           QINT=QINT+QIN(I)
        END DO ent_loop
        
        QINT=QINT*ENUNIT
        POTT=POTT*ENUNIT

        TEMPTR=QINT*2.0d0/BK/DBLE(NDIM)/DBLE(NAN)
        
        RETURN
      END SUBROUTINE TEMPER
