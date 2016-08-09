/** HungarianBPMex.cpp --- 
 *
 * Copyright (C) 2015 Zhen Zhang
 *
 * Author: Zhen Zhang <zzhang@desktop.zzhang.org>
 *
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 * 
 * 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 * 
 * 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */



#include <cstring>
#include "PSBPSolver.h"
#include "PiSBPSolver.h"
#include "HungarianBP.h"
#include "mex.h"


void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
     if(nrhs != 9)
     {
	  mexErrMsgIdAndTxt("HungarianBP:InvalidNumverOfInputs",
			    "Need four inputs NofNodes, StatesofNodes, C, thetac,  outerIters, innerIters, LB");
     }
     int NofNodes = mxGetScalar(prhs[0]);

     double *tNofStates = mxGetPr(prhs[1]);
     int *NofStates = new int[NofNodes];
     int TotalNodedim = 0;
     for(int i = 0; i < NofNodes; i++)
     {
	  NofStates[i] = tNofStates[i];
	  TotalNodedim += NofStates[i];
     }
     
     mwSize nofClusters = mxGetNumberOfElements(prhs[2]);
     if(mxGetNumberOfElements(prhs[2]) != mxGetNumberOfElements(prhs[3]))
	  mexErrMsgIdAndTxt("MRFInference:InvalidInput",
			    "The dimension of c and thetac must be the same");
     
     std::vector< std::vector<int> >Clusters(nofClusters);
     std::vector< boost::shared_array<double> > Potentials(nofClusters);
     std::vector< std::vector<boost::shared_array<double> > > Phi;
     //mexPrintf("NofClusters: %d\n", nofClusters);
     for(int index = 0; index < nofClusters; index++)
     {
	  const mxArray *clu_ptr = mxGetCell(prhs[2], index);
	  int clustersize = mxGetM(clu_ptr) * mxGetN(clu_ptr);
	  double *pcluele = mxGetPr(clu_ptr);
	  std::vector<int> Cluster(clustersize);
	  for(int i = 0; i < clustersize; i++)
	  {
	       Cluster[i] = pcluele[i];
	  }
	  const mxArray *pot_ptr = mxGetCell(prhs[3], index);
	  int potsize = mxGetM(pot_ptr) * mxGetN(pot_ptr);
	  boost::shared_array<double> Potential(new double[potsize]);
	  memcpy(Potential.get(), mxGetPr(pot_ptr), sizeof(double) * potsize);
	  Clusters[index] = Cluster;
	  Potentials[index] = Potential;
     }

     int outerIters = mxGetScalar(prhs[6]);
     int innerIters = mxGetScalar(prhs[7]);
     double LB = mxGetScalar(prhs[8]);

     
     zzhang::CHungarianBP solver;
     solver.ConstructProblemWithConstraints( NofNodes,
					     NofStates,
					     Clusters,
					     Potentials,
					     Phi, LB);

     memcpy(solver.GetU(), mxGetPr(prhs[4]), sizeof(double) * NofNodes);
     memcpy(solver.GetV(), mxGetPr(prhs[5]), sizeof(double) * NofNodes);

     solver.BPClear(outerIters, innerIters);
     
     std::vector<int> Decode;
     solver.GetDecode(Decode);
     plhs[0] = mxCreateDoubleMatrix(NofNodes, 1, mxREAL);
     double* decodeout = mxGetPr(plhs[0]);
     memset(decodeout, 0, sizeof(double) * NofNodes);
     for(int i = 0; i < Decode.size(); i++) decodeout[i] = Decode[i];
     std::vector<int> TotalDim;
     solver.GetClustersAndPotentials(Clusters, Potentials, TotalDim);
     mwSize ndim = 1;
     mwSize dims = Clusters.size();
     plhs[1] = mxCreateCellArray(1, &dims);
     plhs[2] = mxCreateCellArray(1, &dims);
     plhs[3] = mxCreateDoubleMatrix(1,1,mxREAL);
     plhs[4] = mxCreateDoubleMatrix(1,NofNodes,mxREAL);
     plhs[5] = mxCreateDoubleMatrix(1,NofNodes,mxREAL);
     double *dualptr = mxGetPr(plhs[3]);
     solver.GetDual(*dualptr);
     for(int i = 0; i < dims; i++)
     {
	  mxArray *clu = mxCreateDoubleMatrix(1, Clusters[i].size(), mxREAL);
	  mxArray *pot = mxCreateDoubleMatrix(1, TotalDim[i], mxREAL);
	  double *pclu = mxGetPr(clu);
	  for(int j = 0; j < Clusters[i].size(); j++)
	  {
	       pclu[j]=Clusters[i][j];
	  }
	  memcpy(mxGetPr(pot), Potentials[i].get(), sizeof(double) * TotalDim[i]);
	  mxSetCell(plhs[1], i, clu);
	  mxSetCell(plhs[2], i, pot);
     }
     memcpy(mxGetPr(plhs[4]), solver.GetU(), sizeof(double) * NofNodes);
     memcpy(mxGetPr(plhs[5]), solver.GetV(), sizeof(double) * NofNodes);
     delete []NofStates;
}
