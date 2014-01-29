OBJECTS = Nucleus.o Mass.o Chart.o Yrast.o TlArray.o  LevelDensity.o Angle.o Nuclide.o LightP.o Evap.o AngleDist.o Random.o TlBarDist.o Scission.o  Weight.o   SigCharged.o SigBarDist.o Fus.o 


ALLOBJECTS := $(patsubst %.cpp,%.o,$(wildcard *.cpp))
FOBJECTS:=$(patsubst %.f,%.o,$(wildcard *.f))
CFLAGS= -c -Wall -W -O3 -ggdb
COMPILER= c++ 
FC=gfortran

all:  testDecay testDecayROOT testTheWidth testFusion testGemini testWidth

testDecay:  testDecay.o gemini.a
	$(COMPILER) -o testDecay testDecay.o gemini.a

testDecayROOT:  testDecayROOT.o gemini.a
	$(COMPILER) -o testDecayROOT testDecayROOT.o gemini.a $(shell root-config --cflags --libs) 

testTheWidth:  testTheWidth.o gemini.a
	$(COMPILER) -o testTheWidth testTheWidth.o gemini.a

testFusion:  testFusion.o gemini.a
	$(COMPILER) -o testFusion testFusion.o gemini.a  

testGemini:: testGemini.o gemini.o gemini.a
	$(FC) -o testGemini testGemini.o gemini.o gemini.a -lstdc++

testWidth:: testWidth.o gemini.o gemini.a
	$(FC) -o testWidth testWidth.o gemini.o gemini.a -lstdc++

gemini.a: $(OBJECTS)
	ar rcs gemini.a $(OBJECTS)
	ranlib gemini.a

$(ALLOBJECTS): %.o : %.cpp
	$(COMPILER) $(CFLAGS) $< -o $@

$(FOBJECTS): %.o : %.f
	$(FC) $(CFLAGS) $< -o $@


clean:
	rm -f *.o gemini.a

