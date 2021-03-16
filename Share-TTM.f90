!     In this subroutine we performe the information exchange between
!     the processors subject to know temperature and some other
!     parameters in the skin layer
      SUBROUTINE SHARE_TTM()
        INCLUDE 'common.h'
        INCLUDE 'commonTTM.h'
        INCLUDE 'mpif.h'
        INTEGER status(MPI_STATUS_SIZE)
        REAL*8, DIMENSION(:), ALLOCATABLE :: smsgbuf,rmsgbuf

        ALLOCATE(smsgbuf(MSGTTM),rmsgbuf(MSGTTM))

        ixoffset = mynodex*nlcx
        iyoffset = mynodey*nlcy
        izoffset = mynodez*nlcz
        
        neigh4: DO nbb = 1,26

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
                 smsgbuf(l + 3) = T_l(ic)
                 smsgbuf(l + 4) = FSOL(ic)
                 smsgbuf(l + 5) = DBLE(NPR(ic))
                 smsgbuf(l + 6) = DBLE(NVK(ic))
                 smsgbuf(l + 7) = T_e(ic)

                 l = l + ITTM

              ENDDO
              
           ENDIF
           
           msgsend = l - 1

!          Use index "nb" of the loop as a tag fo send/receive pair
           IF(mypid.NE.nbpid(nb)) THEN
              IF(MOD(nodecheck,2).EQ.0) THEN
                 CALL MPI_SEND(smsgbuf,msgsend,MPI_DOUBLE_PRECISION, &
                      nbpid(nb),78+nb,MPI_COMM_WORLD,info)
                 msgrecv = MSGTTM
                 CALL MPI_RECV(rmsgbuf,msgrecv,MPI_DOUBLE_PRECISION, &
                      nbres(nb),78+nb,MPI_COMM_WORLD,status,info)
!                Need to count since we do not know the exact lenght
!                of the recieved message
                 CALL MPI_GET_COUNT(status,MPI_DOUBLE_PRECISION,icount,info)
                 msgrecv = icount
                 ifrom = status(MPI_SOURCE)
              ELSE
                 msgrecv = MSGTTM
                 CALL MPI_RECV(rmsgbuf,msgrecv,MPI_DOUBLE_PRECISION, &
                      nbres(nb),78+nb,MPI_COMM_WORLD,status,info)
                 CALL MPI_GET_COUNT(status,MPI_DOUBLE_PRECISION,icount,info)
                 msgrecv = icount
                 ifrom = status(MPI_SOURCE)
                 CALL MPI_SEND(smsgbuf,msgsend,MPI_DOUBLE_PRECISION, &
                      nbpid(nb),78+nb,MPI_COMM_WORLD,info)
              ENDIF
           ELSE
              msgrecv = msgsend
              rmsgbuf(1:msgrecv) = smsgbuf(1:msgsend)
              ifrom = mypid
           ENDIF

           lim = (msgrecv - 1)/ITTM

           l = 1

           nbcome = INT(rmsgbuf(l))

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

              T_l(ip) = rmsgbuf(l + 3)
              FSOL(ip)= rmsgbuf(l + 4)
              NPR(ip) = INT(rmsgbuf(l + 5))
              NVK(ip) = INT(rmsgbuf(l + 6))              
              T_e(ip) = rmsgbuf(l + 7)

              l = l + ITTM
           ENDDO
              
        ENDDO neigh4
        
        DEALLOCATE(smsgbuf,rmsgbuf)

        RETURN
      END SUBROUTINE SHARE_TTM
