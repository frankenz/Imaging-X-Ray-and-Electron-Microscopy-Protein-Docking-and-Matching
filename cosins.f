C************************************************************
      SUBROUTINE COSINS (R,RAD)
C************************************************************
      INCLUDE 'cosin.h'
C
      DIMENSION R(3)
C
      THEANG=R(1)*RAD
      PHIANG=R(2)*RAD
      PSIANG=R(3)*RAD
      COPHI=COS(PHIANG)
      SIPHI=SIN(PHIANG)
      COTHE=COS(THEANG)
      SITHE=SIN(THEANG)
      COPSI=COS(PSIANG)
      SIPSI=SIN(PSIANG)
      RETURN
      END
