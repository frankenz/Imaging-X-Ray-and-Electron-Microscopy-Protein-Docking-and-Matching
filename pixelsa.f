C------------------------------------------------------------------
       SUBROUTINE PIXELSA (NP,SURF,XR,XI)
C------------------------------------------------------------------
       include 'param.h'
       COMPLEX SURF(np,np,np),TEST
C
       TEST=CMPLX(XR,XI)
       NPXL=0
       DO I1=1,NP
         DO I2=1,NP
           DO I3=1,NP
            tr=real(SURF(I1,I2,I3))
            ti=aimag(surf(i1,i2,i3))
            IF (abs(tr-xr).lt.0.004999) then
            if (abs(ti-xi).lt.0.004999) NPXL=NPXL+1
             end if
           END DO
         END DO
       END DO
       PRINT 10,XR,XI,NPXL
  10   FORMAT(2X,'THE NUMBER OF PIXELS (',F6.2,',',F6.2,') IS ',I8)
       
       RETURN
       END

