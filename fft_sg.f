C=========================================================================
      SUBROUTINE fft_sg (ij,Z,L1,L2,L3,IER)
C=========================================================================
      INCLUDE 'param.h'
c     L1=L2=L3=NP
      COMPLEX*8 Z(L1,L2,L3)
      complex*8 wspace
      common /wsp/ wspace(maxwsp)
C
      if (ij.eq.0) call cfft3di(L1,L2,L3,wspace,maxwsp)
C
C
      call cfft3d (1,L1,L2,L3,z,l1,l2,wspace)
C
      RETURN
      END
