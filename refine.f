c**********************************************************
      SUBROUTINE REFINE (SURFA,SURFB,ISURF,
     *                   JP,kxmax,kymax,kzmax,RT,XB,scmax)
c**********************************************************
      INCLUDE 'param.h'
      INCLUDE 'blk2_3.h'
      INCLUDE 'cosin.h'
C
      COMPLEX*8 SURFA(NP,NP,NP),SURFB(NP,NP,NP)
      INTEGER   ISURF(NP,NP,NP)
C
      DIMENSION RT(3),SCORES(7,7,7),XXB(MAXATM,5),XB(MAXATM,5)
      dimension itrn1(7,7,7),itrn2(7,7,7),itrn3(7,7,7)
      dimension ntrn1(3,3,3),ntrn2(3,3,3),ntrn3(3,3,3),scores1(3,3,3)
C
      DATA ANG1/-0.104719758/,ang2/-0.104719758/,ang3/-0.104719758/
      DATA dang/0.0349065848/
C

      nsteps=7
c
      DO J=1,NB
        DO I=1,5
          XXB(J,I)=XB(J,I)
        END DO
      END DO

      kxin=kxmax
      kyin=kymax
      kzin=kzmax
C
      do i=1,nsteps
         do j=1,nsteps
            do k=1,nsteps
               scores(i,j,k)=0.
            end do
         end do
      end do 
c
c     print *,'angles=',ang1/rad,ang2/rad,ang3/rad,dang/rad
c
      do i1=1,nsteps
         the=ang1+(i1-1)*dang
         cothe=cos(the)
         sithe=sin(the)
         do i2=1,nsteps
            phi=ang2+(i2-1)*dang
            cophi=cos(phi)
            siphi=sin(phi)
            do i3=1,nsteps
               psi=ang3+(i3-1)*dang
               copsi=cos(psi)
               sipsi=sin(psi)
               if (iirot.eq.2) call rotxyz (xxb,xb,nb)
c               if (iirot.lt.2) call rotate (xxb,xb,nb,itest)
c     print *,'the,phi,psi=',the/rad,phi/rad,psi/rad
               if (ipot.gt.1) then
                  print *,'refinement with elec not on yet'
               end if
               call project (surfb,xb,2,1.)
               IF (IOMTB.EQ.1.OR.KCUT.GT.0) THEN
                  CALL PROJECTOUT (SURFB,XB,2,1.)
               END IF
               call productnew (surfa,surfb,stepsz,0,np)
               CALL PEAKS(12,THE,PHI,PSI,HPEAK,SURFB,ISURF)
               backspace (12)
               read (12,999) Idum,Tdum,Pdum,Pdum,itrn1(i1,i2,i3),
     1              itrn2(i1,i2,i3),itrn3(i1,i2,i3),scores(i1,i2,i3)
c           print *,'trans=',itrn1(i1,i2,i3),itrn2(i1,i2,i3),
c    1               itrn3(i1,i2,i3),scores(i1,i2,i3)
            end do
         end do
      end do
 999  format (I10,1X,3F10.2,2X,3I4,1X,F15.2)

      scmax=-9999
      do i1=1,nsteps
         do i2=1,nsteps
            do i3=1,nsteps
               if (scores(i1,i2,i3).gt.scmax) then
                  imax1=i1
                  imax2=i2
                  imax3=i3
                  scmax=scores(i1,i2,i3)
                  kxmax=itrn1(i1,i2,i3)
                  kymax=itrn2(i1,i2,i3)
                  kzmax=itrn3(i1,i2,i3)
               end if
            end do
         end do
      end do
      rt(1)=(ang1+dang*(imax1-1))/rad
      rt(2)=(ang2+dang*(imax2-1))/rad
      rt(3)=(ang3+dang*(imax3-1))/rad
      kxmax=kxmax-np2-1
      kymax=kymax-np2-1
      kzmax=kzmax-np2-1

      if (iabs(kxmax-kxin).gt.9.or.iabs(kymax-kyin).gt.9.or.
     1     iabs(kzmax-kzin).gt.9) then
         print 998,jp,kxmax,kymax,kzmax
 998     format ('Stopped refinement for solution ',i3,'. The position o
     1f the highest score (',3i4,') is too far from the starting positio
     2n')
         return
      end if
c     
c check if maximum score was obtained at the edges of the search region
c
c     print *,'imax values=',imax1,imax2,imax3
      if (imax1.eq.1.or.imax2.eq.1.or.imax3.eq.1.
     1  or.imax1.eq.nsteps.or.imax2.eq.nsteps.or.imax3.eq.nsteps) then
         continue
      else
         print *,'Finished refinement of solution',jp
         return
      end if
c
 1000 continue
      do i1=1,3
         do i2=1,3
            do i3=1,3
               scores1(i1,i2,i3)=0.
            end do
         end do
      end do
c
c define limits of data to be copied to next-step matrix of scores
c
      if (imax1.eq.1) then
         ist1=1
         id1=1
      else
         ist1=imax1-1
         id1=-(imax1-2)        
      end if
      if (imax2.eq.1) then
         ist2=1
         id2=1
      else
         ist2=imax2-1
         id2=-(imax2-2)
      end if
      if (imax3.eq.1) then
         ist3=1
         id3=1
      else
         ist3=imax3-1
         id3=-(imax3-2)
      end if
      if (imax1.eq.nsteps) then
         ils1=nsteps
      else
         ils1=imax1+1
      end if
      if (imax2.eq.nsteps) then
         ils2=nsteps
      else
         ils2=imax2+1
      end if
      if (imax3.eq.nsteps) then
         ils3=nsteps
      else
         ils3=imax3+1
      end if
c     print *,'ist1,ils1,id1=',ist1,ils1,id1
c     print *,'ist2,ils2,id2=',ist2,ils2,id2
c     print *,'ist3,ils3,id3=',ist3,ils3,id3
c
c copying the relevant scores to the new 3x3x3 matrix
c   
      do i1=ist1,ils1
         do i2=ist2,ils2
            do i3=ist3,ils3
c           print *,'copying',i1,i2,i3,scores(i1,i2,i3),' to',
c    1      i1+id1,i2+id2,i3+id3
               scores1(i1+id1,i2+id2,i3+id3)=scores(i1,i2,i3)
               ntrn1(i1+id1,i2+id2,i3+id3)=itrn1(i1,i2,i3)
               ntrn2(i1+id1,i2+id2,i3+id3)=itrn2(i1,i2,i3)
               ntrn3(i1+id1,i2+id2,i3+id3)=itrn3(i1,i2,i3)
            end do
         end do
      end do
c     print *,'scores1 before completion'
c     print *,((scores1(1,i2,i3),i2=1,3),i3=1,3)
c     print *,((scores1(2,i2,i3),i2=1,3),i3=1,3)
c     print *,((scores1(3,i2,i3),i2=1,3),i3=1,3)
c
c completing the missing data in the 3x3x3 matrix
c
      thest=rt(1)*rad-dang
      phist=rt(2)*rad-dang
      psist=rt(3)*rad-dang
c     print *,thest/rad,phist/rad,psist/rad
      do i1=1,3
         the=thest+dang*(i1-1)
         cothe=cos(the)
         sithe=sin(the)
         do i2=1,3
            phi=phist+dang*(i2-1)
            cophi=cos(phi)
            siphi=sin(phi)
            do i3=1,3
               psi=psist+dang*(i3-1)
               copsi=cos(psi)
               sipsi=sin(psi)
               if(scores1(i1,i2,i3).eq.0) then
c             print *,'completing for i1,i2,i3=',i1,i2,i3
c             print *,the/rad,phi/rad,psi/rad
                  if (iirot.eq.2) call rotxyz (xxb,xb,nb)
c                  if (iirot.lt.2) call rotate (xxb,xb,nb,itest)
                  call project (surfb,xb,2,1.)
                  call productnew (surfa,surfb,stepsz,0,np)
                  CALL PEAKS(12,THE,PHI,PSI,HPEAK,SURFB,ISURF)
                  backspace (12)
                  read (12,999) Idum,Tdum,Pdum,Pdum,ntrn1(i1,i2,i3),
     1                 ntrn2(i1,i2,i3),ntrn3(i1,i2,i3),scores1(i1,i2,i3)
               end if
            end do
         end do
      end do
c     print *,'scores1 after completion'
c     print *,((scores1(1,i2,i3),i2=1,3),i3=1,3)
c     print *,((scores1(2,i2,i3),i2=1,3),i3=1,3)
c     print *,((scores1(3,i2,i3),i2=1,3),i3=1,3)
c
c
c find maximum in the 3x3x3 matrix of orientations
c
      scmax=-9999
      nsteps=3
      do i1=1,nsteps
         do i2=1,nsteps
            do i3=1,nsteps
               if (scores1(i1,i2,i3).gt.scmax) then
                  imax1=i1
                  imax2=i2
                  imax3=i3
                  scmax=scores1(i1,i2,i3)
                  kxmax=ntrn1(i1,i2,i3)-np2-1
                  kymax=ntrn2(i1,i2,i3)-np2-1
                  kzmax=ntrn3(i1,i2,i3)-np2-1
               end if
            end do
         end do
      end do
      rt(1)=(thest+dang*(imax1-1))/rad
      rt(2)=(phist+dang*(imax2-1))/rad
      rt(3)=(psist+dang*(imax3-1))/rad
c     print *,'scmax after 2nd round=',scmax
c     print *,'kxmax...',kxmax,kymax,kzmax
c     print *,rt
c     print *,'current imax values',imax1,imax2,imax3
      if (iabs(kxmax-kxin).gt.9.or.iabs(kymax-kyin).gt.9.or.
     1     iabs(kzmax-kzin).gt.9) then
         print 998,jp,kxmax,kymax,kzmax
         return
      end if
c     
      if (imax1.eq.1.or.imax2.eq.1.or.imax3.eq.1.
     1   or.imax1.eq.nsteps.or.imax2.eq.nsteps.or.imax3.eq.nsteps) then
         continue
      else
         print *,'Finished refinement of solution',jp
         return
      end if
c
      do i1=1,3
         do i2=1,3
            do i3=1,3
               scores(i1,i2,i3)=scores1(i1,i2,i3)
               itrn1(i1,i2,i3)=ntrn1(i1,i2,i3)
               itrn2(i1,i2,i3)=ntrn2(i1,i2,i3)
               itrn3(i1,i2,i3)=ntrn3(i1,i2,i3)
            end do
         end do
      end do
      go to 1000            
C
      RETURN
      END
