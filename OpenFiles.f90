!     THIS SUBROUTINE OPENS MOST OF THE INPUT/OUTPUT FILES 
!     Leonid Zhigilei, 2000
      SUBROUTINE OpenFiles(openf)
        INCLUDE 'common.h'
        character*14 str1, str2
        character*4 str
        character*8 openf
        
        OPEN (UNIT = 10,FILE=openf)
        
        opener_loop: DO 
           read(10,1001,end=1004) str,iunit,str1,str2
1001       FORMAT(A4,1x,I2,1x,A14,1x,A7) 
           IF (str(1:1).eq.'#') CYCLE opener_loop
           IF (str(1:1).eq.'*') THEN
              ldr=leng(str1)
              IF(ldr.eq.0) Then
                 WRITE(*,*) 'Wrong name of directory',str1
                 STOP
              ENDIF
              drname(1:ldr)=str1(1:ldr)
              CYCLE opener_loop
           ENDIF

           lf1=leng(str1)
           IF (lf1.eq.0) Then
              WRITE(*,*) 'Wrong format in the file ',openf
              STOP
           ENDIF
           lf2=leng(str2)
           OPEN (iunit,file=str1(1:lf1),status=str2(1:lf2),err=255)
           CYCLE opener_loop
255        WRITE(*,*) 'Error while reading file ',str1
           STOP
1004       EXIT opener_loop
        END DO opener_loop

        RETURN
      END SUBROUTINE OpenFiles
