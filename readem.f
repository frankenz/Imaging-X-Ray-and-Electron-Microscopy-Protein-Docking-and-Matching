C-----------------------------------------------------------------------
      SUBROUTINE READEM (INAT,NMOL,IVOX,X,IHD) 
C-----------------------------------------------------------------------

C     ROUTINE TO READ AND DEFINE ATOMS BY THEIR GRID INTERVAL
      INCLUDE 'param.h'      
      CHARACTER*1 Skp
      CHARACTER*1 CI
      REAL DENSCUT,DENS,GIF,GIM,GIS,RA
      INTEGER S,M,F,NS,NM,NF,IVOX,OUT,NACI

      DIMENSION X(MAXATM,5),IHD(MAXATM)
      IVOX=0
      NACI=1
      CI='A'

      DO I=1,MAXATM
         X(I,5)=0.
      END DO

      PRINT*,'MOLECULE',NMOL,' DATA:'

	IF(NMOL.EQ.1) THEN
	OUT=40
	END IF
c
	IF(NMOL.EQ.2) THEN
	OUT=45
	END IF

      IC1=0
      NAT=0
      NRES=1

      READ (INAT,*) NF
      READ (INAT,*) NM
      READ (INAT,*) NS
      READ (INAT,*) DENSCUT
      READ (INAT,*) GIF
      READ (INAT,*) GIM
      READ (INAT,*) GIS
c

      IF ((GIF.NE.GIM).OR.(GIF.NE.GIS).OR.(GIM.NE.GIS)) THEN
      PRINT 33
      stop
      END IF

      RA= sqrt(GIF*GIF+GIM*GIM+GIS*GIS)/2.
c
      DO S=1,NS
         DO M=1,NM
            DO F=1,NF
               READ (INAT,*) DENS
               IF(DENS.GE.DENSCUT) THEN
                  IVOX=IVOX+1
                  X(IVOX,1)=(F-1)*GIF
                  X(IVOX,2)=(M-1)*GIM
                  X(IVOX,3)=(S-1)*GIS
                  X(IVOX,4)=RA
                  X(IVOX,5)=DENS

                  NAT=NAT+1
                  IC1=IC1+1
                  IF (IC1.EQ.51) THEN
                     NRES=NRES+1
                     IC1=0
                  END IF
                  WRITE (OUT,111) NAT,NRES,X(IVOX,1),X(IVOX,2),
     1                 X(IVOX,3),DENS

               END IF
            END DO
         END DO
      END DO
    
 111  FORMAT('HETATM',I5,2X,'H1',2X,'HT1',1X,A1,I4,4X,3F8.3,2X,'1.00',
     1       F6.1)

      
 33   FORMAT ('*** NOTE VOXEL SIZE VALUES ARE NOT EQUAL ***')
      PRINT*,'GRID POINTS FAST ',NF
      PRINT*,'GRID POINTS MEDIUM ',NM
      PRINT*,'GRID POINTS SLOW ',NS
      PRINT*,'DENSITY CUTOFF ',DENSCUT      
      PRINT*,'VOXEL SIZE FAST ',GIF
      PRINT*,'VOXEL SIZE MEDIUM ',GIM
      PRINT*,'VOXEL SIZE SLOW ',GIS
      PRINT*,'VIRTUAL ATOM RADII ',RA

 
      PRINT 999,IVOX,NMOL
      PRINT*

C      GENERAL PRINTING MESSAGES
 999  FORMAT (I6,' VOXELS READ FOR MOLECULE ',I1)

      RETURN
      END




      
