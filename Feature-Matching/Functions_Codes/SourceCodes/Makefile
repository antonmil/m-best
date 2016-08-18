objs=HungarianBPMex.o BPSolver.o DualSolver.o PiSBPSolver.o PSBPSolver.o HungarianBP.o
#MATLAB=/usr/local/MATLAB/R2016a/
MEX=mex


all:HungarianBPMex CluComputeObjMex AdditionOfPhiMex

HungarianBP:${objs}
	${MEX} -v -O -output HungarianBP ${objs}
CluComputeObjMex:CluComputeObjMex.o
	${MEX} -v -O -output CluComputeObj CluComputeObjMex.o
AdditionOfPhiMex:AdditionOfPhiMex.o
	${MEX} -v -O -output AdditionOfPhiMex AdditionOfPhiMex.o

%.o:%.cpp
	${MEX} -c $^ -O

clean:
	rm -rf *.mex*
	rm -rf *.o
