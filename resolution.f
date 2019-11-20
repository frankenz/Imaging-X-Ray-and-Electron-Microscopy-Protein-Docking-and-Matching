C=============================================================
      SUBROUTINE RESOLUTION (np,surf,xresl)
c       FORMER subroutine reslser(np,surf,xresl)
C=============================================================
C
c     sub for the changing the amplitude of various resolutions
c     this sub is partly based on sub resl23 

      INCLUDE 'param.h'
c      INCLUDE 'blk2_3.h'
    
      complex*8 surf(np,np,np)
c
c      print *,'xresl in reslser1=',xresl
c      print *,'np in reslser1=',np
      xresl2=xresl*xresl
      np2=np/2
      do i=1,np
         do j=1,np
            do k=1,np
               
c     origin shift correction
               ii=i-np2
               jj=j-np2
               kk=k-np2
               if (ii.lt.0) then 
                  ii=i
               else 
                  ii=i-np
               end if
               if (jj.lt.0) then 
                  jj=j
               else 
                  jj=j-np
               end if
               if (kk.lt.0) then 
                  kk=k
               else 
                  kk=k-np
               end if
c               
               dis=(float((ii)**2+(jj)**2+(kk)**2))
c               
               if (dis.gt.xresl2) then
                  surf (i,j,k)=(0.,0.)
               end if
c               
            end do
         end do
      end do
c       
c
      RETURN
      END











