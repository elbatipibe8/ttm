!     This subroutine is responsible for sharing of the EAM density information
!     for particles that are in the skin of the current preocessor
      SUBROUTINE SHARE_EAM()
        INCLUDE 'common.h'
        INCLUDE 'mpif.h'
        INTEGER status(MPI_STATUS_SIZE)
        REAL*8, DIMENSION(:), ALLOCATABLE :: smsgbuf,rmsgbuf

        ALLOCATE(smsgbuf(LPMX),rmsgbuf(LPMX))

        ivid = nn1

        neigh3: DO nbb = 1,26

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

!          Prevent communication based on PBC
           IF(ncom(nb).EQ.1) THEN

              DO icn = 1,nbcell(nb)

                 ic = kbcell(nb,icn)

                 i = ltop(ic)

10               CONTINUE

                 IF(i.GT.0) THEN

                    smsgbuf(l) = RhoDen(i)

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
                      nbpid(nb),52+nb,MPI_COMM_WORLD,info)
                 msgrecv = LPMX
                 CALL MPI_RECV(rmsgbuf,msgrecv,MPI_DOUBLE_PRECISION, &
                      nbres(nb),52+nb,MPI_COMM_WORLD,status,info)
!                Need to count the lenght since its value is not known
                 CALL MPI_GET_COUNT(status,MPI_DOUBLE_PRECISION,icount,info)
                 msgrecv = icount
                 ifrom = status(MPI_SOURCE)
              ELSE
                 msgrecv = LPMX
                 CALL MPI_RECV(rmsgbuf,msgrecv,MPI_DOUBLE_PRECISION, &
                      nbres(nb),52+nb,MPI_COMM_WORLD,status,info)
!                Need to count the lenght since its value is not known
                 CALL MPI_GET_COUNT(status,MPI_DOUBLE_PRECISION,icount,info)
                 msgrecv = icount
                 ifrom = status(MPI_SOURCE)
                 CALL MPI_SEND(smsgbuf,msgsend,MPI_DOUBLE_PRECISION, &
                      nbpid(nb),52+nb,MPI_COMM_WORLD,info)
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
              
!          In the skin we do not construct link cell chaines, therefore
!          Make sure that the atomes recieve their density value in the
!          same order as their position was received
           ivid = imstart(ifrom,nbcome)

           DO i = 1,lim

!             Select next atom in the skin
              ivid = ivid + 1

              RhoDen(ivid)  = rmsgbuf(l)

              l = l + 1

           ENDDO
           
        ENDDO neigh3
        
        DEALLOCATE(smsgbuf,rmsgbuf)

        RETURN
      END SUBROUTINE SHARE_EAM
