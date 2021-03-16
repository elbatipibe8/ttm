       SUBROUTINE NbList_MPI()
         INCLUDE 'common.h'

         INTEGER lcellx(13),lcelly(13),lcellz(13)

!        13 neighbores are enough to construct the full NbList with
!        number of operations proportional NN1/2
         DATA lcellx/-1,0,1,-1,0,1,-1,0,1,-1,0,1,-1/
         DATA lcelly/-1,-1,-1,0,0,0,1,1,1,-1,-1,-1,0/
         DATA lcellz/-1,-1,-1,-1,-1,-1,-1,-1,-1,0,0,0,0/

         NNG(1:NN1) = 0

         DO iz = 1,nlcz1
           DO iy = 0,nlcy1
               DO ix = 0,nlcx1

!                 Define ID of a host link cell
                  ic = ix + nlcx2*(iy + nlcy2*iz)
!                 Find the entry point in the host link cell
                  i = ltop(ic)
                  
!                 If it is zero - the cell is empty
115               IF(i.EQ.0) GOTO 111

                  ktypei = KTYPE(I)

!                 Do not need to consider a particle's neigbores within 
!                 the same link cellif it is belong to the skin
                  IF(ix.eq.0.or.ix.eq.nlcx1 &
                       .and.iy.eq.0.or.iy.eq.nlcy1.or.iz.eq.nlcz1) GOTO 112
                  
                  j = i

!                 Peak next particle                  
110               j = link(j)
                  
                  IF(j.eq.0) GOTO 112

                  ij = ijindex(ktypei,KTYPE(J))
  
                  Rij = RList(ij)
                  R2ij = RList2(ij)

                  DX=XD(1,i)-XD(1,j) 
                  IF (ABS(DX).LT.Rij) THEN
                     DY=XD(2,i)-XD(2,j)  
                     IF (ABS(DY).LT.Rij) THEN
                        DZ=XD(3,i)-XD(3,j)
                        IF (ABS(DZ).LT.Rij) THEN               
                           RR=DX*DX+DY*DY+DZ*DZ
                           IF (RR.LT.R2ij) THEN

                              IF(i.LE.nn1) THEN
                                 NNG(i)=NNG(i)+1
                                 NNNG(NNG(i),i)=j
                              ENDIF
                              IF(j.LE.nn1)THEN
                                 NNG(j)=NNG(j)+1 
                                 NNNG(NNG(j),j)=i
                              ENDIF
                              
                           ENDIF
                        ENDIF
                     ENDIF
                  ENDIF
                  
!                 Will consider next particle if any in the chain
                  GOTO 110

!                 Consider particles from 13 neighbores
!                 Will count neighbor cells within the node's boundaries only
112               DO nb = 1,13
                     
                     jx = ix + lcellx(nb)
                     jy = iy + lcelly(nb)
                     jz = iz + lcellz(nb)
                     
                     IF(jx.ge.0.and.jx.le.nlcx1 &
                          .and.jy.ge.0.and.jy.le.nlcy1) THEN
                        
                        jc = jx + nlcx2*(jy + nlcy2*jz)
                        j = ltop(jc)

114                     IF(j.EQ.0) GOTO 113

                        ij = ijindex(ktypei,KTYPE(J))

!                       Ramp for BS model inside the host cell
                        Rij = RList(ij)
                        R2ij = RList2(ij)

                        DX=XD(1,i)-XD(1,j) 
                        IF (ABS(DX).LT.Rij) THEN
                           DY=XD(2,i)-XD(2,j)  
                           IF (ABS(DY).LT.Rij) THEN
                              DZ=XD(3,i)-XD(3,j)
                              IF (ABS(DZ).LT.Rij) THEN               
                                 RR=DX*DX+DY*DY+DZ*DZ
                                 IF (RR.LT.R2ij) THEN

                                    IF(i.LE.nn1)THEN
                                       NNG(i)=NNG(i)+1
                                       NNNG(NNG(i),i)=j
                                    ENDIF
                                    IF(j.LE.nn1)THEN
                                       NNG(j)=NNG(j)+1 
                                       NNNG(NNG(j),j)=i
                                    ENDIF
      
                                 ENDIF
                              ENDIF
                           ENDIF
                        ENDIF
                        
                        j = link(j)
                        GOTO 114

                     ENDIF
                     
113                  CONTINUE
                     ENDDO
                  
!                 Take next particle in the chain of the host link cell
                  i = link(i)
                  GOTO 115
                  
111               CONTINUE
               ENDDO
            ENDDO
         ENDDO

         NNNMAX = 0
         NNN = 0
         DO I=1,NN1
            IF(NNG(I).GE.NNNMAX) THEN
               NNNMAX = NNG(I)
               NNN = I
            ENDIF
         ENDDO

         IF(NNNMAX.GT.MAXNNB) print *,NNNMAX, " number of neighbores for",NNN, &
              " atom in",mypid," processor is greater than MAXNNB"

         RETURN
       END SUBROUTINE NbList_MPI
