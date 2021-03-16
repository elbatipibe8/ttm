!     Perform initial settings for bulk targets. Processors are given a
!     specific ID whether the belong to MD or purely TTM part. The same is
!     done for linkcells that are different in numbers for TTM and Md part.
!     In particular, TTM operations are approximatelyt 10 times faster than
!     those of MD for each MD step. Therefore, the size of TTM part is mdev=10
!     times bigger in the direction of NRB.
      SUBROUTINE SetBulk()
        INCLUDE 'common.h'
        INCLUDE 'commonTTM.h'
        INCLUDE 'mpif.h'
        INTEGER status1(MPI_STATUS_SIZE),status2(MPI_STATUS_SIZE)

        ixoffset = mynodex*nlcx
        iyoffset = mynodey*nlcy
        izoffset = mynodez*nlcz

        IF(TIME.EQ.0.0d0) THEN
           zdist = 0
           
           DO i = 1,NN1
              IF(KHIST(i).GT.10) THEN
                 zget = ZCENTR-XD(3,i)
                 IF(zget.gt.zdist) zdist = zget
              ENDIF
           ENDDO
           
           IF(mypid.ne.0)CALL MPI_SEND(zdist,1,MPI_DOUBLE_PRECISION, &
                0,0,MPI_COMM_WORLD,ierror)
           IF(mypid.eq.0)  THEN
              DO nrec = 1,nnodes-1
                 CALL MPI_RECV(zradii,1,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,0, &
                      MPI_COMM_WORLD,status1,ierror)
                 IF(zradii.ge.zdist) zdist = zradii
              ENDDO
           ENDIF
           CALL MPI_BCAST(zdist,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        ENDIF

!       Redefine the system boundaries so that the space belonging to the bulk
!       would comprise much greater space to provide
!       thermal conduction mechanism as it is in bulk. Depending on the
!       material the bulk space can very in size significantly. The minimal
!       size of the bulk margines is defined in TTMInit as "TTMBULK" variable.
!       Broadcast the coordinate of the beggining of pure bulk section
!       Redefine the link cell positions for pure bulk section in geometrical
!       progression. This definition results in less memory consumtion

        critza = 0.0d0
        DO iz = 1,nlcz
           DO iy = 1,nlcy
              DO ix = 1,nlcx
                 ip = ix + nlcx2*(iy + nlcy2*iz)

                 izt=INT((ZCENTR-ZED(ip)-zdist)*nzsign*nlczg/ZL+1)-marbulk

                 IF(izt.ge.1) THEN
                    izts = ABS(ZED(ip) - ZCENTR)/(ZED(ip) - ZCENTR)
                    difz(ip) = adifz*bulkin**izt
                    ZED(ip) = ZED(ip)+ 1.0d+10*izts*(-adifz*(izt+0.5d0) + &
                         adifz*(bulkin**izt-1.0d0)/(bulkin-1.0d0) + &
                         0.5d0*adifz*bulkin**izt)
                    IF(ABS(ZED(ip)-ZCENTR).GE.critza) critza=ABS(ZED(ip)-ZCENTR)
                 ENDIF

              ENDDO
           ENDDO
        ENDDO
        
        IF(mypid.ne.0)CALL MPI_SEND(critza,1,MPI_DOUBLE_PRECISION, &
             0,0,MPI_COMM_WORLD,ierror)
        IF(mypid.eq.0)  THEN
           DO nrec = 1,nnodes-1
              CALL MPI_RECV(critz,1,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,0, &
                   MPI_COMM_WORLD,status2,ierror)
              IF(critz.ge.critza.and.critz.ne.0.0d0) critza = critz
           ENDDO
        ENDIF
        CALL MPI_BCAST(critza,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

        IF(critza.LT.TTMBULK) THEN
           IF(mypid.eq.0) THEN
              print *,"The size of ttm part in bulk",critza,"section"
              print *,"is not enough to build up space for bulk"
              print *,"simulation",TTMBULK,"[A] Resize the system or increase"
              print *,"'bulkin' parameter in TTMInit.f90. The job is halted"
!              call Flushout(6)
           ENDIF
           CALL MPI_FINALIZE(info)
           STOP
        ENDIF

        IF(critza.GT.4.0d0*TTMBULK) THEN
           IF(mypid.eq.0) THEN
              print *,"The size of ttm part in bulk section ",critza
              print *,"is too big as compare to its min value",TTMBULK,"[A]"
              print *,"Resize the system or decrease 'bulkin' parameter in"
              print *,"TTMInit.f90 The job is halted now."
!              call Flushout(6)
           ENDIF
           CALL MPI_FINALIZE(info)
           STOP
        ENDIF

        IF(mypid.eq.0) THEN
           print *,"The size of ttm part is",critza,"[A] in Z"
           print *," "
!           call Flushout(6)
        ENDIF

        RETURN
      END SUBROUTINE SetBulk
