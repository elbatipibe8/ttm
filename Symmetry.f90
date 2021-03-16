!     Calculation of the Central Symmetry parameter 
!     for identification of crystal structure
!     To use it,need to have alat defined in Eftab for example
      SUBROUTINE Symmetry()
        INCLUDE 'common.h'
        INCLUDE 'mpif.h'
!       Choose 12 neighbores for CSP construction
        INTEGER, PARAMETER :: NCS = 12
        DIMENSION RQ(NSYM)
        INTEGER, DIMENSION (:,:), ALLOCATABLE :: MMMG
        INTEGER, DIMENSION (:), ALLOCATABLE :: MMG
        REAL*8, DIMENSION (:), ALLOCATABLE :: smsgbuf,rmsgbuf
        INTEGER status(MPI_STATUS_SIZE)
        LOGICAL  KKG(NSYM)

        ALLOCATE(smsgbuf(LPMX),rmsgbuf(LPMX))
        ALLOCATE(MMMG(NCS,NN1),MMG(NN1))

        CALL Temper(QINT,POTT,TEMPTR)
        CALL MPI_REDUCE(TEMPTR,TEMPTRG,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
             MPI_COMM_WORLD,ierr)
        CALL MPI_BCAST(TEMPTRG,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

        ExpCoef = 1.922962811d-05
        alatT = alat*(1.0d0+ExpCoef*TEMPTRG)  ! approximate lattice parameter

!       Defining a cutoff distance so that the first two neighbor shells
!       in a perfect fcc lattice are included
        RADCUT=0.855d0*alatT
        RCUT2=RADCUT*RADCUT
        MMMG(1:NCS,1:NN1) = 0
        MMG(1:NN1) = 0
        FD(2,1:NN1) = 0.0d0
        FD(3,1:NN1) = 0.0d0

!       Finding all the neighbors within the cutoff distance 
!       of each atom using the neighbor list
        outer: DO I=1,NN1
           RQ(1:NCS) = RCUT2
           MMM = 0
           inner: DO JJ=1,NNG(I)  ! full list is assumed here
              J=NNNG(JJ,I)
              DZ=XD(3,I)-XD(3,J)
              IF (ABS(DZ).LT.RADCUT) THEN
                 DX=XD(1,I)-XD(1,J)
                 IF (ABS(DX).LT.RADCUT) THEN
                    DY=XD(2,I)-XD(2,J)
                    IF(ABS(DY).LT.RADCUT) THEN
                       RR=DX*DX+DY*DY+DZ*DZ
                       IF(RR.LT.RCUT2) THEN
                          IF(MMM.LT.NCS) MMM = MMM + 1
                          IF(MMM.GT.1) THEN
                             K = 1
                             DO WHILE(K.LE.MMM)
                                IF(RR.LT.RQ(K)) THEN
725                                RRT = RQ(K)
                                   JT = MMMG(K,I)
                                   RQ(K) = RR
                                   MMMG(K,I) = J
                                   RR = RRT
                                   J = JT
                                   K = K + 1
                                   IF(K.LE.MMM) GOTO 725
                                ELSE
                                   K = K + 1
                                ENDIF
                             ENDDO
                          ELSE
                             RQ(1) = RR
                             MMMG(1,I) = J
                          ENDIF
                       ENDIF
                    ENDIF
                 ENDIF
              ENDIF
           ENDDO inner
           IF(MMM.LT.NCS) THEN
              MMG(I) = (MMM/2)*2
           ELSE
              MMG(I) = NCS
           ENDIF
        ENDDO outer

!       Constructing CSP out of MMG neighbores
        i_loop:DO I = 1,NN1
           NLIM = MMG(I)
           IF(NLIM.EQ.0) THEN
              FD(2,I) = 0.0d0
              CYCLE i_loop
           ENDIF
           IF(NLIM.EQ.1) THEN
              FD(2,I) = 1.0d0
              CYCLE i_loop
           ENDIF
           KKG(1:NLIM) = .FALSE.
           j_loop:DO k = 1,NLIM-1
              IF(KKG(k)) CYCLE j_loop
              J = MMMG(k,I)
              RSMIN2 = 4.0d0*RCUT2 !Max possible difference of two RCUT2 vectros
              XCSJ = XD(1,J) - XD(1,I)
              YCSJ = XD(2,J) - XD(2,I)
              ZCSJ = XD(3,J) - XD(3,I)
              FD(3,I) = FD(3,I) + XCSJ*XCSJ + YCSJ*YCSJ + ZCSJ*ZCSJ
              jj_loop:DO m = k + 1,NLIM
                 IF(KKG(m)) CYCLE jj_loop
                 JJ = MMMG(m,I)
                 XCSJJ = XD(1,JJ) - XD(1,I)
                 YCSJJ = XD(2,JJ) - XD(2,I)
                 ZCSJJ = XD(3,JJ) - XD(3,I)
                 RSC2=(XCSJ+XCSJJ)*(XCSJ+XCSJJ) + &
                      (YCSJ+YCSJJ)*(YCSJ+YCSJJ) + &
                      (ZCSJ+ZCSJJ)*(ZCSJ+ZCSJJ)
                 IF(RSC2.LT.RSMIN2) THEN
                    RSMIN2 = RSC2
                    RSMAX2 = XCSJJ*XCSJJ + YCSJJ*YCSJJ + ZCSJJ*ZCSJJ
                    LMIN = m
                 ENDIF
              ENDDO jj_loop
              FD(2,I) = FD(2,I) + RSMIN2
              FD(3,I) = FD(3,I) + RSMAX2
              KKG(k)    = .TRUE.
              KKG(LMIN) = .TRUE.
           ENDDO j_loop
           IF(FD(3,I).NE.0.0d0) THEN
              FD(2,I) = 0.5d0*FD(2,I)/FD(3,I)
           ELSE
              FD(2,I) = 0.0d0
              print *,"Attention to CSP in ",mypid,"proc.", &
                   "Squre distanse is zero. CSP of atom",I," will be set to zero."
           ENDIF
        ENDDO i_loop

        CALL MPI_BARRIER(MPI_COMM_WORLD,ipass1)
        
!       Need to know information about the skin for averaging
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       Take the CSP parameters to the skin from all the neighboring nodes

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
                    smsgbuf(l) = FD(2,i)
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
              FD(2,ivid)  = rmsgbuf(l)
              l = l + 1
           ENDDO
        ENDDO neigh6

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!       Averaging over neighboring atoms
        o_loop: DO I=1,NN1
           FD(3,I) = FD(2,I)
           in_loop: DO KK=1,MMG(I)
              J=MMMG(KK,I)
              FD(3,I)=FD(3,I)+FD(2,J)
           ENDDO in_loop
           FD(3,I)=FD(3,I)/DBLE(MMG(I)+1)
        ENDDO o_loop

        DEALLOCATE(MMMG,MMG,smsgbuf,rmsgbuf)

        CALL MPI_BARRIER(MPI_COMM_WORLD,ipass2)

        RETURN
      END SUBROUTINE Symmetry
