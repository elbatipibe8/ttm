!     Calculation of the local order parameter 
!     for identification of liquid and crystal regions
!     To use it,need to have alat defined in Eftab for example
      SUBROUTINE Liquid()
        INCLUDE 'common.h'
        INCLUDE 'mpif.h'
        DIMENSION QX(6),QY(6),QZ(6)
        REAL*8, DIMENSION (:), ALLOCATABLE :: WAL
        INTEGER status(MPI_STATUS_SIZE)
        REAL*8, DIMENSION (:), ALLOCATABLE :: smsgbuf,rmsgbuf        

        ALLOCATE(WAL(LPMX),smsgbuf(LPMX),rmsgbuf(LPMX))

        CALL Temper(QINT,POTT,TEMPTR)
        CALL MPI_REDUCE(TEMPTR,TEMPTRG,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
             MPI_COMM_WORLD,ierr)
        CALL MPI_BCAST(TEMPTRG,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

        ExpCoef= 1.922962811d-05
        alatT=alat*(1.0d0+ExpCoef*TEMPTRG)  ! approximate lattice parameter

!       Defining a set of six reference vectors (perfect for fcc,
!       but may not be acceptable for all the possible stucters)
        QX=(4*PI/alatT)*(/1.0d0, 0.0d0, 0.0d0, 1.0d0, 0.0d0, 1.0d0/)
        QY=(4*PI/alatT)*(/0.0d0, 1.0d0, 0.0d0, 1.0d0, 1.0d0, 0.0d0/)
        QZ=(4*PI/alatT)*(/0.0d0, 0.0d0, 1.0d0, 0.0d0, 1.0d0, 1.0d0/)

!       Defining a cutoff distance so that the first two neighbor shells
!       in a perfect fcc lattice are included

        RADCUT=1.115d0*alatT
        RCUT2=RADCUT*RADCUT
        WAL(1:NN1) = 0.0d0

!       Finding all the neighbors within the cutoff distance 
!       of each atom using the neighbor list
        outer_loop: DO I=1,NN1
           SUM1=0.0d0
           SUM2=0.0d0
           MXG=0
           inner_loop: DO JJ=1,NNG(I)  ! full list is assumed here
              J=NNNG(JJ,I)
              DZ=XD(3,I)-XD(3,J)
              IF (ABS(DZ).LT.RADCUT) THEN
                 DX=XD(1,I)-XD(1,J)
                 IF (ABS(DX).LT.RADCUT) THEN
                    DY=XD(2,I)-XD(2,J)
                    IF(ABS(DY).LT.RADCUT) THEN
                       RR=DX*DX+DY*DY+DZ*DZ
                       IF(RR.LT.RCUT2) THEN
                          MXG=MXG+1
                          DO K=1,6
                             SUM1=SUM1+COS(QX(K)*DX+QY(K)*DY+QZ(K)*DZ)
                             SUM2=SUM2+SIN(QX(K)*DX+QY(K)*DY+QZ(K)*DZ)
                          ENDDO
                       ENDIF
                    ENDIF
                 ENDIF
              ENDIF
           ENDDO inner_loop
           IF(MXG.NE.0) WAL(I)=(SUM1*SUM1+SUM2*SUM2)/(36.0d0*MXG*MXG)
        ENDDO outer_loop

!       Need to know information about the skin for averaging
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Take the LOP parameters to the skin from all the neighboring nodes

        ivid = nn1

        neigh6: DO nbb = 1,26

           IF(nbb.le.18) THEN
              nb = longarray(nbb)
              nodecheck = longnode
           ELSEIF(nbb.gt.24)THEN
              nb = korarray(nbb-24)
              nodecheck = kornode
           ELSE
              nb = midlarray(nbb-18)
              nodecheck = midlnode
           ENDIF

           l = 1
           smsgbuf(l) = DBLE(nb)
           l = l + 1

           IF(ncom(nb).EQ.1) THEN

              DO icn = 1,nbcell(nb)

                 ic = kbcell(nb,icn)

                 i = ltop(ic)

10               CONTINUE

                 IF(i.GT.0) THEN

                    smsgbuf(l) = WAL(i)

                    l = l + 1

                    i = link(i)

                    GOTO 10

                 ENDIF

              ENDDO

           ENDIF

           msgsend = l - 1

!          Use index "nb" of the loop as a tag fo send/receive pair
           IF(mypid.NE.nbpid(nb)) THEN
              IF(MOD(nodecheck,2).EQ.0) THEN
                 CALL MPI_SEND(smsgbuf,msgsend,MPI_DOUBLE_PRECISION, &
                      nbpid(nb),130+nb,MPI_COMM_WORLD,info)
                 msgrecv = LPMX
                 CALL MPI_RECV(rmsgbuf,msgrecv,MPI_DOUBLE_PRECISION, &
                      nbres(nb),130+nb,MPI_COMM_WORLD,status,info)
                 CALL MPI_GET_COUNT(status,MPI_DOUBLE_PRECISION,icount,info)
                 msgrecv = icount
                 ifrom = status(MPI_SOURCE)
              ELSE
                 msgrecv = LPMX
                 CALL MPI_RECV(rmsgbuf,msgrecv,MPI_DOUBLE_PRECISION, &
                      nbres(nb),130+nb,MPI_COMM_WORLD,status,info)
!                Since there is no exact lenght of the message was passed,
!                need to conut it
                 CALL MPI_GET_COUNT(status,MPI_DOUBLE_PRECISION,icount,info)
                 msgrecv = icount
                 ifrom = status(MPI_SOURCE)
                 CALL MPI_SEND(smsgbuf,msgsend,MPI_DOUBLE_PRECISION, &
                      nbpid(nb),130+nb,MPI_COMM_WORLD,info)
              ENDIF
           ELSE
              msgrecv = msgsend
              rmsgbuf(1:msgrecv) = smsgbuf(1:msgsend)
              ifrom = mypid
           ENDIF

           lim = msgrecv - 1

           l = 1

           nbcome = INT(rmsgbuf(l))
           l = l + 1

           ivid = imstart(ifrom,nbcome)
!          Assighn this information in exactly the same order
!          as the position for this particle was recieved

           DO i = 1,lim
!             Select nexst atom above NN1
              ivid = ivid + 1

              WAL(ivid)  = rmsgbuf(l)

              l = l + 1

           ENDDO

        ENDDO neigh6

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        FD(1,1:NN1) = 0.0d0

!       Averaging over neighboring atoms
        o_loop: DO I=1,NN1
           FD(1,I) = WAL(I)
           MXG=0
           i_loop: DO KK=1,NNG(I)
              J=NNNG(KK,I)
              DZ=XD(3,I)-XD(3,J)
              IF (ABS(DZ).LT.RADCUT) THEN
                 DX=XD(1,I)-XD(1,J)
                 IF (ABS(DX).LT.RADCUT) THEN
                    DY=XD(2,I)-XD(2,J)
                    IF(ABS(DY).LT.RADCUT) THEN
                       RR=DX*DX+DY*DY+DZ*DZ
                       IF (RR.LT.RCUT2) THEN
                          MXG=MXG+1
                          FD(1,I)=FD(1,I) + WAL(J)
                       END IF
                    END IF
                 END IF
              END IF
           END DO i_loop
           FD(1,I)=FD(1,I)/DBLE(MXG+1)
        END DO o_loop

        DEALLOCATE(WAL,smsgbuf,rmsgbuf)

        RETURN
      END SUBROUTINE Liquid
