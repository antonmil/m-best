objs=HungarianBPMex.o BPSolver.o DualSolver.o PiSBPSolver.o PSBPSolver.o HungarianBP.o



all:HungarianBPMex.mexa64 CluComputeObjMex.mexa64

HungarianBPMex.mexa64:${objs}
	mex -v -O -output HungarianBPMex ${objs}
CluComputeObjMex.mexa64:CluComputeObjMex.cpp
	mex -v -O CluComputeObjMex.cpp

%.o:%.cpp
	mex -v -O -c $^

clean:
	rm -rf *.mexa64
	rm -rf *.o
