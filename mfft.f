C=========================================================================
      SUBROUTINE MFFT (ij,Z,L1,L2,L3,IER)
C=========================================================================
      INCLUDE 'param.h'
C     L1=L2=L3=NP
      COMPLEX*8 Z(L1,L2,L3)
      complex*8 wspace
C     common /wsp/ wspace(maxwsp)
      common /wfft/ iwl(1742),iwm(1742),iwn(1742),iwork(432)
c      common /wfft/ iwl(1038),iwm(1038),iwn(1038),iwork(256)
      INTEGER iord, iopt
      iord = 1
      iopt = 0

C
      np=l1
C Initialization of FFT tables
 10   IF (ij.eq.0) THEN 
         call newfft(z,np,np,np,np,iwl,iwm,iwn,iopt,0,iord,
     1        iwork,ierr)
C
C check that the fft functions worked well, if not stop.
         IF (ierr.ne.0) THEN
            print*,'****************************'
            IF (ierr.eq.1) THEN
               print *,'ERROR, MFFT DATA DIMENTIONS ARE NOT CORRECT.'
            ELSE IF (ierr.eq.2) THEN
               PRINT*,'ERROR, NP HAS PRIME FACTORS DIFFERENT FROM 2,
     1              3, 5'
            ELSE
               PRINT*, 'ERROR, FFT TABLES ARE NOT CORRECTLY INITIALIZE'
            END IF
            stop
         END IF
      END IF
C  Do a direct FFT
C     call cfft3d (1,L1,L2,L3,z,l1,l2,wspace)
 20   CALL newfft(z,np,np,np,np,iwl,iwm,iwn,iopt,-1,iord,
     1            iwork,ierr)
c
c
C check that the fft functions worked well, if not stop.
      IF (ierr.eq.1) THEN
         print *,'ERROR, MFFT DATA DIMENTIONS ARE NOT CORRECT.'
         stop
      END IF
      IF (ierr.eq.3) THEN
         PRINT*, 'ERROR, FFT TABLES ARE NOT CORRECTLY INITIALIZED'
         stop
      END IF
C
C
      RETURN
      END
