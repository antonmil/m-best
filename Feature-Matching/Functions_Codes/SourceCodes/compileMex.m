%% 
mex -c BPSolver.cpp
mex -c DualSolver.cpp
mex -c HungarianBP.cpp
mex -c PiSBPSolver.cpp
mex -c PSBPSolver.cpp
mex -c HungarianBPMex.cpp
mex -c CluComputeObjMex.cpp

mex HungarianBP.o HungarianBPMex.o BPSolver.o DualSolver.o PiSBPSolver.o PSBPSolver.o -output HungarianBP

mex  CluComputeObjMex.o -output CluComputeObj
