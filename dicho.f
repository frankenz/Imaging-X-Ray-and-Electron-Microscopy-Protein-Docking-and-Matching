**********************************************************************
      SUBROUTINE DICHO (np,dt_range,rangeu,SURF)
C**********************************************************************
      INCLUDE 'param.h'
c     INCLUDE 'blk2_3.h'

c      Common /wsp/ wspace(maxwsp)
c      COMMON /wfft/ iwl(1038),iwm(1038),iwn(1038),iwork(256)
      COMMON /wfft/ iwl(1742),iwm(1742),iwn(1742),iwork(432)

C
      COMPLEX*8 SURF(NP,NP,NP)
      INTEGER iord, iopt,ierr
      INTEGER NP3
      REAL dt_range
C
C     initialize parameters for the mfft routines.
      iord = 1
      iopt =0
C
      NP3 = NP**3
C
C     INVERSE FOURIER TRANSFORM
 20   call newfft(surf,np,np,np,np,iwl,iwm,iwn,iopt,1,
     1              iord,iwork,ierr)
C
C Check that the fft functions worked well, if not stop.
      IF (ierr.eq.1) THEN
         print *,'ERROR, MFFT DATA DIMENTIONS ARE NOT CORRECT.'
         stop
      END IF
      IF (ierr.eq.3) THEN
         PRINT*, 'ERROR, FFT TABLES ARE NOT CORRECTLY INITIALIZED'
         stop
      END IF
C
c      print *,'in bfft rangel,rengeu=',rangel,rangeu
c deviding in np^3 and dichotomizing into 1 and 0
c according to the limits set by rangel and rangeu

      do i=1,np
         do j=1,np
            do k=1,np
               surf(i,j,k)=surf(i,j,k)/(NP3)
               if (REAL(surf (i,j,k)) .LT. dt_range) then
                  surf (i,j,k)=(0.,0.) 
               else
                  surf (i,j,k)=(1.,0.)
               end if
               
            end do
         end do
      end do
C
C
      RETURN 
      END
      
