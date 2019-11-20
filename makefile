SRC = apmask.f      exchng.f            productnew.f \
blkdta000.f         project.f           pixels.f \
renum_pdb.f         unpeck.f		readem.f\
projectout.f\
fill.f              rotxyz.f            projpot.f \
pixelreal.f         pixelimag.f \
fold.f              reg.f               prtmol.f \
prtpot.f \
rotate.f            rotatepot.f         rotxyzpot.f \
cosins.f \
srtp.f              pixelsa.f \
srt1.f \
surface.f \
edge.f \
edgeout.f \
edgein.f            molsize.f           unfold.f \
origin.f \
peaklst.f \
peaks.f             potential.f \
prtmolpot.f \
wrtc.f \
exchgi.f            readData.f \
mfft.f \
newfft.f \
mfftp.f             mfftz0.f            mfftov.f \
mfftom.f            mfftiv.f            mfftis.f \
mfftdv.f            mfftdm.f            mfftds.f \
mfftim.f            mfftp1.f            mfftp2.f \
mfftp4.f            mffta5.f            mfftb5.f \
mfftc5.f            mffta9.f            mfftb9.f \
mfftc9.f            mffta4.f            mfftb4.f \
mfftc4.f            mffta6.f            mfftb6.f \
mfftc6.f            mffta8.f            mfftb8.f \
mfftc8.f            mffta7.f            mfftb7.f \
mfftc7.f            mfftp3.f \
charToReal.f        findIdx.f           findNP.f\
main.f              residueID.f         projecthd.f \
resolution.f        dicho.f             refine.f

          
 
 
OBJS =  apmask.o     exchng.o            productnew.o \
blkdta000.o         project.o           pixels.o \
renum_pdb.o         unpeck.o 		readem.o\
projectout.o\
fill.o              rotxyz.o            projpot.o \
pixelreal.o         pixelimag.o \
fold.o              reg.o               prtmol.o \
prtpot.o \
rotate.o            rotatepot.o         rotxyzpot.o \
cosins.o \
srtp.o              pixelsa.o \
srt1.o \
surface.o \
edge.o \
edgeout.o \
edgein.o            molsize.o           unfold.o \
origin.o \
peaklst.o \
peaks.o             potential.o \
prtmolpot.o \
wrtc.o \
exchgi.o            readData.o \
mfft.o \
newfft.o \
mfftp.o             mfftz0.o            mfftov.o \
mfftom.o            mfftiv.o            mfftis.o \
mfftdv.o            mfftdm.o            mfftds.o \
mfftim.o            mfftp1.o            mfftp2.o \
mfftp4.o            mffta5.o            mfftb5.o \
mfftc5.o            mffta9.o            mfftb9.o \
mfftc9.o            mffta4.o            mfftb4.o \
mfftc4.o            mffta6.o            mfftb6.o \
mfftc6.o            mffta8.o            mfftb8.o \
mfftc8.o            mffta7.o            mfftb7.o \
mfftc7.o            mfftp3.o \
charToReal.o        findIdx.o           findNP.o\
main.o              residueID.o         projecthd.o \
resolution.o        dicho.o             refine.o

FC = ifort
FFLAGS = -O2
LDFLAGS = 

all: FitEM2EM.exe

FitEM2EM.exe: $(OBJS)
	$(FC) $(FFLAGS) -o FitEM2EM.exe $(OBJS)
