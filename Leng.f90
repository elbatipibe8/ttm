      function leng( str )
        character str*(*)

        do leng = len(str), 1, -1
           if ( str(leng:leng).ne.' ' ) return
        end do

        return
      end function leng
