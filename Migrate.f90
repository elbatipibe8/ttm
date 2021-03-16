!     This subroutine is responsible for passing of the information
!     for particles that go beyond the baoundary of the processor
      SUBROUTINE MIGRATE()
        INCLUDE 'common.h'
        INCLUDE 'mpif.h'
        REAL*8, DIMENSION(:), ALLOCATABLE :: smsgbuf,rmsgbuf
        INTEGER, DIMENSION(:), ALLOCATABLE :: left
        INTEGER status(MPI_STATUS_SIZE)
        INTEGER IEMPTY
        LOGICAL leave
        PARAMETER (IEMPTY = -1)

        ALLOCATE(smsgbuf(MSGMIG),rmsgbuf(MSGMIG))
        ALLOCATE(left(LPMZ))

        ixoffset = mynodex*nlcx
        iyoffset = mynodey*nlcy
        izoffset = mynodez*nlcz

!       If an atom goes byond the boundary of the node it leaves
!       a hole in an array. Set the index hole counter.
        lv = 0
        left(1:LPMZ) = 0

        neigh1: DO nbb = 1,26

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

!          Loop through skin cells
           DO icn = 1,nscell(nb)
              ic = kscell(nb,icn)
              i = ltop(ic)

10            CONTINUE

              IF(i.gt.0) THEN

!                Define wether an atom can leave the system in this node
                 leave = ((LIDX.EQ.0.AND.DABS(XD(1,i)-XCENTR).GT.XLHALF) &
                      .OR.(LIDY.EQ.0.AND.DABS(XD(2,i)-YCENTR).GT.YLHALF) &
                      .OR.(LIDZ.EQ.0.AND.DABS(XD(3,i)-ZCENTR).GT.ZLHALF))

!                Send real atoms (indexes from 1 to NN1) but
!                ignor visitors (indexes above NN1 from skin)
                 IF(i.LE.nn1) THEN
                    IF(.NOT.LEAVE) THEN

                       smsgbuf(l)      = XD(1,i) + ximage(nb)
                       smsgbuf(l +  1) = XD(2,i) + yimage(nb)
                       smsgbuf(l +  2) = XD(3,i) + zimage(nb)
                       smsgbuf(l +  3) = Q1D(1,i)
                       smsgbuf(l +  4) = Q1D(2,i)
                       smsgbuf(l +  5) = Q1D(3,i)
                       smsgbuf(l +  6) = Q2(3*i-2)
                       smsgbuf(l +  7) = Q2(3*i-1)
                       smsgbuf(l +  8) = Q2(3*i  )
                       smsgbuf(l +  9) = Q3(3*i-2)
                       smsgbuf(l + 10) = Q3(3*i-1)
                       smsgbuf(l + 11) = Q3(3*i  )
                       smsgbuf(l + 12) = Q4(3*i-2)
                       smsgbuf(l + 13) = Q4(3*i-1)
                       smsgbuf(l + 14) = Q4(3*i  )
                       smsgbuf(l + 15) = Q5(3*i-2)
                       smsgbuf(l + 16) = Q5(3*i-1)
                       smsgbuf(l + 17) = Q5(3*i)                     
                       smsgbuf(l + 18) = DBLE(KHIST(i))
                       smsgbuf(l + 19) = DBLE(KTYPE(i))

                       l = l + IMIGRATE
                    ENDIF

!                   Set identifier of the atomes leaving the node
                    KHIST(i) = IEMPTY
                    lv = lv + 1
                    left(lv) = i
                    
                 ENDIF

                 i = link(i)
                 GOTO 10
                 
              ENDIF

              ltop(ic) = 0

           ENDDO

           msgsend = l - 1

!          Use blocking MPI to essure for a reliable passing of the information
!          Use index "nb" of the loop as a tag fo send/receive pair
           IF(mypid.NE.nbpid(nb)) THEN
              IF(MOD(nodecheck,2).EQ.0) THEN
                 CALL MPI_SEND(smsgbuf,msgsend,MPI_DOUBLE_PRECISION, &
                      nbpid(nb),nb,MPI_COMM_WORLD,info)
                 msgrecv = MSGMIG
                 CALL MPI_RECV(rmsgbuf,msgrecv,MPI_DOUBLE_PRECISION, &
                      nbres(nb),nb,MPI_COMM_WORLD,status,info)
!                Need to count the lenght since it vvalue is not known
                 CALL MPI_GET_COUNT(status,MPI_DOUBLE_PRECISION,icount,info)
                 msgrecv = icount
                 ibelong = status(MPI_SOURCE)
              ELSE
                 msgrecv = MSGMIG
                 CALL MPI_RECV(rmsgbuf,msgrecv,MPI_DOUBLE_PRECISION, &
                      nbres(nb),nb,MPI_COMM_WORLD,status,info)
!                Need to count the lenght since it value is not known 
                 CALL MPI_GET_COUNT(status,MPI_DOUBLE_PRECISION,icount,info)
                 msgrecv = icount
                 ibelong = status(MPI_SOURCE)
                 CALL MPI_SEND(smsgbuf,msgsend,MPI_DOUBLE_PRECISION, &
                      nbpid(nb),nb,MPI_COMM_WORLD,info)
              ENDIF
           ELSE
              msgrecv = msgsend
              rmsgbuf(1:msgrecv) = smsgbuf(1:msgsend)
              ibelong = mypid
           ENDIF

           l = 1

           DO i = 1,msgrecv/IMIGRATE

              IF(lv.gt.0) THEN
                 idx = left(lv)
                 lv = lv - 1
              ELSE
                 NN1 = NN1 + 1
                 idx = NN1
              ENDIF

              XD(1,idx)  = rmsgbuf(l)
              XD(2,idx)  = rmsgbuf(l +  1)
              XD(3,idx)  = rmsgbuf(l +  2)
              Q1D(1,idx) = rmsgbuf(l +  3)
              Q1D(2,idx) = rmsgbuf(l +  4)
              Q1D(3,idx) = rmsgbuf(l +  5)
              Q2(3*idx-2) = rmsgbuf(l +  6)
              Q2(3*idx-1) = rmsgbuf(l +  7)
              Q2(3*idx  ) = rmsgbuf(l +  8)
              Q3(3*idx-2) = rmsgbuf(l +  9)
              Q3(3*idx-1) = rmsgbuf(l + 10)
              Q3(3*idx  ) = rmsgbuf(l + 11)
              Q4(3*idx-2) = rmsgbuf(l + 12)
              Q4(3*idx-1) = rmsgbuf(l + 13)
              Q4(3*idx  ) = rmsgbuf(l + 14)
              Q5(3*idx-2) = rmsgbuf(l + 15)
              Q5(3*idx-1) = rmsgbuf(l + 16)
              Q5(3*idx  ) = rmsgbuf(l + 17)
              KHIST(idx) = INT(rmsgbuf(l + 18))
              KTYPE(idx) = INT(rmsgbuf(l + 19))

              l = l + IMIGRATE

!             Link the atom to the correct linkcell

              ixg = INT((XD(1,idx) - xlstart)*nxsign*nlcxg/XL + 1)
              iyg = INT((XD(2,idx) - ylstart)*nysign*nlcyg/YL + 1)
              izg = INT((XD(3,idx) - zlstart)*nzsign*nlczg/ZL + 1)

              ix = ixg - ixoffset
              iy = iyg - iyoffset
              iz = izg - izoffset

              ip = ix + nlcx2*(iy + nlcy2*iz)

!             Assign atom i to linkcell ip
              j = ltop(ip)
              ltop(ip) = idx
              link(idx) = j
           ENDDO
        ENDDO neigh1

        IF(NN1.GT.LPMZ) THEN 
           print *,"Number of particles",NN1," in",mypid,"processor", &
                "exceeded LPMZ =",LPMZ,"in Migrate"
!           call Flushout(6)
           STOP
        ENDIF

!       Back-fill any remaining index holes with higher-numbered atoms,
!       decrementing atom count after each move

50      IF(lv.GT.0) THEN
           
           idx = left(lv)
           
60         IF(KHIST(NN1).EQ.IEMPTY) THEN
              NN1 = NN1 - 1
              IF(NN1.EQ.0) THEN
                 lv = 0
                 GOTO 50
              ENDIF
              GOTO 60
           ENDIF

           IF(idx.LT.NN1) THEN

              ixg = INT((XD(1,NN1) - xlstart)*nxsign*nlcxg/XL + 1)
              iyg = INT((XD(2,NN1) - ylstart)*nysign*nlcyg/YL + 1)
              izg = INT((XD(3,NN1) - zlstart)*nzsign*nlczg/ZL + 1)

              ix = ixg - ixoffset
              iy = iyg - iyoffset
              iz = izg - izoffset

              ip = ix + nlcx2*(iy + nlcy2*iz)

!             Find atom and swap its local ID in to the linked list
              j = ltop(ip)

65            IF(j.EQ.NN1) THEN
                 ltop(ip) = idx
                 link(idx) = link(nn1)
                 link(nn1) = 0
              ELSEIF (j.EQ.0) THEN
                 print *,"atom is not found, Migrate",mypid
!                 call Flushout(6)
                 stop
              ELSE

70               jj = j

              j = link(j)

              IF(j.EQ.NN1) THEN
                 link(jj) = idx
                 link(idx) = link(nn1)
                 link(nn1) = 0
              ELSEIF (j.EQ.0) THEN
                 print *,"an error, migrate",mypid
!                 call Flushout(6)
                 stop
              ELSE
                 GOTO 70
              ENDIF
              
           ENDIF

           XD(1,idx)  =  XD(1,nn1)
           XD(2,idx)  =  XD(2,nn1)
           XD(3,idx)  =  XD(3,nn1)
           Q1D(1,idx) = Q1D(1,nn1)
           Q1D(2,idx) = Q1D(2,nn1)
           Q1D(3,idx) = Q1D(3,nn1)
           Q2(3*idx-2) = Q2(3*nn1-2)
           Q2(3*idx-1) = Q2(3*nn1-1)
           Q2(3*idx  ) = Q2(3*nn1  )
           Q3(3*idx-2) = Q3(3*nn1-2)
           Q3(3*idx-1) = Q3(3*nn1-1)
           Q3(3*idx  ) = Q3(3*nn1  )
           Q4(3*idx-2) = Q4(3*nn1-2)
           Q4(3*idx-1) = Q4(3*nn1-1)
           Q4(3*idx  ) = Q4(3*nn1  )
           Q5(3*idx-2) = Q5(3*nn1-2)
           Q5(3*idx-1) = Q5(3*nn1-1)
           Q5(3*idx  ) = Q5(3*nn1  )           
           KHIST(idx) = KHIST(nn1)
           KTYPE(idx) = KTYPE(nn1)

           NN1 = NN1 - 1
           
        ENDIF
        
        lv = lv - 1

        GOTO 50
        
     ENDIF

     DEALLOCATE(smsgbuf,rmsgbuf)
     DEALLOCATE(left)

     RETURN
   END SUBROUTINE MIGRATE
