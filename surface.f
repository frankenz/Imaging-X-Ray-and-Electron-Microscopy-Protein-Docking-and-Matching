C**********************************************************************
      SUBROUTINE SURFACE (SURFA,SURFB,XYZPOT,ISURF,XA,IHDA,NMOL)
C**********************************************************************
      INCLUDE 'param.h'
      INCLUDE 'blk2_3.h'
C
      COMPLEX*8 SURFA(NP,NP,NP),SURFB(NP,NP,NP)
      INTEGER   ISURF(NP,NP,NP), ifft
      DIMENSION XA(MAXATM,5)
      DIMENSION XYZPOT(MAXPOT,4)
      DIMENSION IHDA(MAXATM)
      CHARACTER*3 CRES

c
C PROJECT AND EDGE (ONLY FOR MOLECULE A)
C
      NP2=INT(NP/2)
c the first projection derermines the inner surface layer
      CALL PROJECT (SURFA,XA,NMOL,0.)
c       print *,'after 1st project of A'
c       call pixels (np,surfa,1.0,0.0)
c       call pixels (np,surfa,rho,0.0)
c       call pixels (np,surfa,0.0,0.0)
c
C  preparing a single edge with rho inside
      CALL EDGE (NP,RHO,SURFA)
c      print *,'after edge of A'
c       call pixels (np,surfa,1.0,0.0)
c       call pixels (np,surfa,rho,0.0)
c       call pixels (np,surfa,0.0,0.0)

c       CALL PRTMOL (SURFA,RHO,STEPSZ,NP,NP2,1)
c
C CREATE A 'MASK' TO IMPROVE THE SURFACE LAYER OF MOLECULE A
C this done first on the grid of B. Then SURFB use in the rest
C by add a water aurface i.e. 1.4 Amstrong to the atomic radiuos
c
      PRINT*
      PRINT*,'MOLECULE A OUTER SURFACE LAYER EXTENDED 
     1BY',XXM,'TIMES OF GRID STEP SIZE' 

      CALL PROJECT (SURFB,XA,NMOL,XXM)
c      print *,'after 2nd project of A'
c      print *,' Surface of A extended by ',xxm
c      call pixels (np,surfb,1.0,0.0)
c      call pixels (np,surfb,rho,0.0)
c      call pixels (np,surfb,0.0,0.0)
C
C fill empty spots (bag grid fill) with value
      CALL FILL (NP,SURFB,ISURF,RHO)
c      print *,'after fill'
c      call pixels (np,surfb,1.0,0.0)
c      call pixels (np,surfb,rho,0.0)
c      call pixels (np,surfb,0.0,0.0)
C
C The ifft flag is needed in order to know if there is a need to
C do the initialization to the fft or not. The initialization is done
C only once, before that first time we call the fft routines. So if
C the RESL option is on, we have to do the fft initialization before
C the fft calls and we have to flag that we already did fft, so we
C will not do the initialization again at the end of the routine.
      ifft=0
C

C APPLY THE MASK
C
C copy the intreior of molecule A from grid B (creation place)
C to its origin source. this for use on less memory
      CALL APMASK (NP,SURFA,SURFB)
c      print *,'after apmask of A'
c       call pixels (np,surfa,1.0,0.0)
c       call pixels (np,surfa,rho,0.0)
c       call pixels (np,surfa,0.0,0.0)

c
      IF (IOMTA.EQ.1.OR.KCUT.GT.0)THEN
	 CALL PROJECTOUT (SURFA,XA,1,10.)
      ENDIF

C   read delphi potential and prepare potential surface for molecule A
c      print *,'after projectout'
c      call pixels (np,surfa,1.0,0.0)
c      call pixels (np,surfa,rho,0.0)
c      call pixels (np,surfa,0.0,0.0)
C
c        CALL PRTMOL (SURFA,RHO,STEPSZ,NP,NP2,1)
C Projection of hydrophobicity for molecule A
      IF (IHDRF .EQ. 1) THEN
C     print *,'in surface hd=',hd        
         CALL PROJECTHD(SURFA,XA,IHDA,NMOL,10.)
      ENDIF
C
C Projection of electrostatics for molecule A
      IF (IPOT.EQ.1) THEN
         CALL POTENTIAL (SURFA,SURFB,ISURF,XYZPOT)
c     CALL PRTPOT (SURFA,STEPSZ,NP,NP2,1)
      END IF
C
c     call pixelreal (np,surfa,1.0)
c     call pixelreal (np,surfa,rho)
c     call pixelreal (np,surfa,0.0)
c     call pixelimag (np,1,surfa)
c
C   FFT FOR MOLECULE A
C
c         print *,'ifft=',ifft  
C If RESL was not operated then ifft = 0, which signals the MFFT
C routine to do the initialization of the arrays before doing the
C fft. If RESL was on, then ifft hold the value of 1, which flags
C MFFT routine NOT to do initialization of the fft arrays.
 300  CALL MFFT (ifft,SURFA,NP,NP,NP,IER)
C     
      RETURN
      END
