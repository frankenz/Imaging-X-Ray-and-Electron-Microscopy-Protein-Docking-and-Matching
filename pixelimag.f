C------------------------------------------------------------------
       SUBROUTINE PIXELimag(NP,nmol,SURF)
C------------------------------------------------------------------
       include 'param.h'
       COMPLEX SURF(np,np,np)
C
       NPXL=0
       DO I1=1,NP
         DO I2=1,NP
           DO I3=1,NP
             IF (aimag(SURF(I1,I2,I3)).NE.0.0) then
             NPXL=NPXL+1
c            print 999,aimag(SURF(I1,I2,I3)),i1,i2,i3
c999   format ('Imaginary value =',f10.5,' at',3i4)
c            if (nmol.eq.1) then
c             if (i1.eq.57.and.i2.eq.50.and.i3.eq.46) then
c              ti=aimag(SURF(I1,I2,I3))
c              surf(i1,i2,i3)=cmplx(0.0,ti)
c            print 999,aimag(SURF(I1,I2,I3)),i1,i2,i3
c999   format ('Imaginary value =',f10.5,' at',3i4)
c             else
c              tr=real(surf(I1,I2,I3))
c              SURF(I1,I2,I3)=cmplx(0.0,0.0)
c             end if
c            end if
c            if (nmol.eq.2) then
c             if (i1.eq.51.and.i2.eq.31.and.i3.eq.55) then
c              ti=aimag(SURF(I1,I2,I3))
c              surf(i1,i2,i3)=cmplx(0.0,ti)
c            print 999,aimag(SURF(I1,I2,I3)),i1,i2,i3
c             else
c              tr=real(surf(I1,I2,I3))
c              SURF(I1,I2,I3)=cmplx(0.0,0.0)
c             end if
c            end if
c            else
c            SURF(I1,I2,I3)=cmplx(0.0,0.0)
             END IF
           END DO
         END DO
       END DO
       PRINT 10,NPXL
  10   FORMAT(2X,'THE NUMBER OF NON-0 IMAGINARY PIXELS IS ',I8)
       RETURN
       END

