!     In this subroutine, a processor is provided with
!     information about its skin layer. This information
!     is being passed from releative neighbor processors
      SUBROUTINE SHARE()
        INCLUDE 'common.h'
        INCLUDE 'mpif.h'
        INTEGER status(MPI_STATUS_SIZE)
        REAL*8, DIMENSION(:), ALLOCATABLE :: smsgbuf,rmsgbuf

        ALLOCATE(smsgbuf(MSGSHA),rmsgbuf(MSGSHA))

        ixoffset = mynodex*nlcx
        iyoffset = mynodey*nlcy
        izoffset = mynodez*nlcz

        ivid = NN1

!       Loop through all 26 neighbores to get the skin information
        neigh2: DO nbb = 1,26

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

                    smsgbuf(l)     = XD(1,i) + ximage(nb)
                    smsgbuf(l + 1) = XD(2,i) + yimage(nb)
                    smsgbuf(l + 2) = XD(3,i) + zimage(nb)
                    smsgbuf(l + 3) = DBLE(KHIST(i))
                    smsgbuf(l + 4) = DBLE(KTYPE(i))

                    l = l + ISHARE

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
                      nbpid(nb),26+nb,MPI_COMM_WORLD,info)
                 msgrecv = MSGSHA
                 CALL MPI_RECV(rmsgbuf,msgrecv,MPI_DOUBLE_PRECISION, &
                      nbres(nb),26+nb,MPI_COMM_WORLD,status,info)
!                Need to count the lenght since its value is not known
                 CALL MPI_GET_COUNT(status,MPI_DOUBLE_PRECISION,icount,info)
                 msgrecv = icount
                 ifrom = status(MPI_SOURCE)
              ELSE
                 msgrecv = MSGSHA
                 CALL MPI_RECV(rmsgbuf,msgrecv,MPI_DOUBLE_PRECISION, &
                      nbres(nb),26+nb,MPI_COMM_WORLD,status,info)
!                Need to count the lenght since its value is not known
                 CALL MPI_GET_COUNT(status,MPI_DOUBLE_PRECISION,icount,info)
                 msgrecv = icount
                 ifrom = status(MPI_SOURCE)
                 CALL MPI_SEND(smsgbuf,msgsend,MPI_DOUBLE_PRECISION, &
                      nbpid(nb),26+nb,MPI_COMM_WORLD,info)
              ENDIF
           ELSE
              msgrecv = msgsend
              rmsgbuf(1:msgrecv) = smsgbuf(1:msgsend)
              ifrom = mypid
           ENDIF

           lim = (msgrecv - 1)/ISHARE

           l = 1

           nbcome = int(rmsgbuf(l))

           imstart(ifrom,nbcome) = ivid

           l = l + 1

           DO i = 1,lim

!             Select nexst atom above NN1.
!             After this subroutine, ivid - is the current number of particles
!             in the processor including skin area
              ivid = ivid + 1

              XD(1,ivid)  = rmsgbuf(l)
              XD(2,ivid)  = rmsgbuf(l +  1)
              XD(3,ivid)  = rmsgbuf(l +  2)
              KHIST(ivid) = INT(rmsgbuf(l + 3))
              KTYPE(ivid) = INT(rmsgbuf(l + 4))
              l = l + ISHARE

              IF(MOD(ISTEP,NLREN).EQ.0) THEN

!             LINK to the appropriate Linkcell
              ixg = INT((XD(1,ivid) - xlstart)*nxsign*nlcxg/XL + 1)
              iyg = INT((XD(2,ivid) - ylstart)*nysign*nlcyg/YL + 1)
              izg = INT((XD(3,ivid) - zlstart)*nzsign*nlczg/ZL + 1)

              ix = ixg - ixoffset
              iy = iyg - iyoffset
              iz = izg - izoffset
              
              ip = ix + nlcx2*(iy + nlcy2*iz)

!             For SW potential (KEYBS = 5) we need to pass 3-body part of
!             the forces later. Therefore, we requaree the same order of
!             particles in the skin, thought this way is a bit longer
!             In principle, kvid myght be a function of cell's ip, but since
!             we pass the information cell by cell, kvid of the current
!             cell will not overlap with another one until we are done
!             with the current cell
              IF(KEYBS.EQ.5) THEN
                 IF(ltop(ip).eq.0) THEN
                    ltop(ip) = ivid
                    kvid = ivid
                    link(ivid) = 0
                 ELSE
                    link(kvid) = ivid
                    kvid = ivid
                    link(ivid) = 0
                 ENDIF
              ELSE
                 j = ltop(ip)
                 ltop(ip) = ivid
                 link(ivid) = j
              ENDIF

           ENDIF

           ENDDO

        ENDDO neigh2

!       Current nubmer of particles in the processor including skin layer
        NNF = ivid

        IF(NNF.GT.LPMX) THEN 
           print *,"Number of particles NNF = ",NNF," in ",mypid, &
                "processor exceeded LPMX =",LPMX,"due to Share"
!           call Flushout(6)
222        CONTINUE
           GOTO 222
        ENDIF

        DEALLOCATE(smsgbuf,rmsgbuf)

        RETURN
      END SUBROUTINE SHARE
