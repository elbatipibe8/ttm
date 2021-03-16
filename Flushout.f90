!     flushout flushes the buffer of iunit to force printing/writing
!     to a file.  It is useful for debugging 
      subroutine flushout (iunit)
        character(LEN=255):: filename
        LOGICAL:: ex, nmd
        INTEGER, INTENT(IN):: iunit
        INTEGER:: ios
!       get the name of the external file
        inquire (unit=iunit, iostat=ios, err=905, &
             exist=ex, named=nmd, name=filename)
        if (ex .and. nmd) then
!          close unit to flush buffer
           close (unit=iunit, iostat=ios, err=905)
!          reconnect to external file
           open (unit=iunit, iostat=ios, err=905, &
                &         status='unknown', position='APPEND', file=filename)

!          Nov 2003, this positioning has been removed from sputs flush
!          move pointer to the bottom of the file
!105       continue
!          read (unit=iunit, iostat=ios, err=905, &
!          fmt='(a)', end=905) str
!          goto 105
        endif
905     continue
        if (ios.gt.0) write(*,*) ' flush: unit=',iunit,', iostat=',ios
        return
      end  subroutine flushout
