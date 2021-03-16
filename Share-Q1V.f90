!     Send/recieve the temperatutre and so on information into the skin
!     layer of the related node
      SUBROUTINE SHARE_Q1V()
        INCLUDE 'common.h'
        INCLUDE 'commonTTM.h'
        INCLUDE 'mpif.h'
        INTEGER status(MPI_STATUS_SIZE)
        REAL*8, DIMENSION(:), ALLOCATABLE :: smsgbuf,rmsgbuf

        ALLOCATE(smsgbuf(MSGTTM3),rmsgbuf(MSGTTM3))

        ixoffset = mynodex*nlcx
        iyoffset = mynodey*nlcy
        izoffset = mynodez*nlcz

        neigh5: DO nbb = 1,26

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
                 
                 smsgbuf(l)     = Xsh(ic) + ximage(nb)
                 smsgbuf(l + 1) = Ysh(ic) + yimage(nb)
                 smsgbuf(l + 2) = Zsh(ic) + zimage(nb)
                 smsgbuf(l + 3) = Q1V(ic)
                 smsgbuf(l + 4) = RON(ic)
                 IF(ICOND(ic)) THEN
                    smsgbuf(l + 5) = 1.0d0
                 ELSE
                    smsgbuf(l + 5) = 0.0d0
                 ENDIF
                 
                 l = l + ITTM3
                 
              ENDDO

           ENDIF

           msgsend = l - 1

!           print *,mypid,msgsend,MSGTTM3
!           CALL MPI_FINALIZE(info)
!           STOP

!          Use index "nb" of the loop as a tag fo send/receive pair
           IF(mypid.NE.nbpid(nb)) THEN
              IF(MOD(nodecheck,2).EQ.0) THEN
                 CALL MPI_SEND(smsgbuf,msgsend,MPI_DOUBLE_PRECISION, &
                      nbpid(nb),104+nb,MPI_COMM_WORLD,info)
                 msgrecv = MSGTTM3
                 CALL MPI_RECV(rmsgbuf,msgrecv,MPI_DOUBLE_PRECISION, &
                      nbres(nb),104+nb,MPI_COMM_WORLD,status,info)
                 CALL MPI_GET_COUNT(status,MPI_DOUBLE_PRECISION,icount,info)
                 msgrecv = icount
                 ifrom = status(MPI_SOURCE)
              ELSE
                 msgrecv = MSGTTM3
                 CALL MPI_RECV(rmsgbuf,msgrecv,MPI_DOUBLE_PRECISION, &
                      nbres(nb),104+nb,MPI_COMM_WORLD,status,info)
                 CALL MPI_GET_COUNT(status,MPI_DOUBLE_PRECISION,icount,info)
                 msgrecv = icount
                 ifrom = status(MPI_SOURCE)
                 CALL MPI_SEND(smsgbuf,msgsend,MPI_DOUBLE_PRECISION, &
                      nbpid(nb),104+nb,MPI_COMM_WORLD,info)
              ENDIF
           ELSE
              msgrecv = msgsend
              rmsgbuf(1:msgrecv) = smsgbuf(1:msgsend)
              ifrom = mypid
           ENDIF

           lim = (msgrecv - 1)/ITTM3

           l = 1

           nbcome = int(rmsgbuf(l))

           l = l + 1

           DO i = 1,lim

!             Define the position of the recieved linkcell
              ixg = INT((rmsgbuf(l  ) - xlstart)*nxsign*nlcxg/XL + 1)
              iyg = INT((rmsgbuf(l+1) - ylstart)*nysign*nlcyg/YL + 1)
              izg = INT((rmsgbuf(l+2) - zlstart)*nzsign*nlczg/ZL + 1)

              ix = ixg - ixoffset
              iy = iyg - iyoffset
              iz = izg - izoffset

              ip = ix + nlcx2*(iy + nlcy2*iz)

              Q1V(ip) = rmsgbuf(l + 3)
              RON(ip) = rmsgbuf(l + 4)
              IF(rmsgbuf(l + 5).EQ.1.0d0) THEN
                 ICOND(ip) = .TRUE.
              ELSE
                 ICOND(ip) = .FALSE.
              ENDIF
              
              l = l + ITTM3

           ENDDO
              
        ENDDO neigh5

        DEALLOCATE(smsgbuf,rmsgbuf)

        RETURN
      END SUBROUTINE SHARE_Q1V
