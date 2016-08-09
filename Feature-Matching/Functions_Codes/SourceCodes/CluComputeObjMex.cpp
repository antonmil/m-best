/** MRFMex.cpp --- 
 *
 * Copyright (C) 2015 Zhen Zhang
 *
 * Author: Zhen Zhang <zzhang@desktop.zzhang.org>
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see http://www.gnu.org/licenses/.
 */



#include <cstring>
#include "mex.h"

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
     if(nrhs != 4)
     {
	  mexErrMsgIdAndTxt("CluComputeObjMex:InvalidNumverOfInputs",
			    "Need four inputs Decodes, NofNodes, StatesofNodes, C, thetac");
     }
     double *decodes = mxGetPr(prhs[0]);
     int NofNodes = mxGetM(prhs[1]) * mxGetN(prhs[1]);;
     double v = 0;
     double *tNofStates = mxGetPr(prhs[1]);
     int *NofStates = new int[NofNodes];
     int TotalNodedim = 0;
     for(int i = 0; i < NofNodes; i++)
     {
	  NofStates[i] = tNofStates[i];
	  TotalNodedim += NofStates[i];
     }
     
     mwSize nofClusters = mxGetNumberOfElements(prhs[3]);
     if(mxGetNumberOfElements(prhs[2]) < mxGetNumberOfElements(prhs[3]))
	  mexErrMsgIdAndTxt("MRFInference:InvalidInput",
			    "The dimension of c and thetac must be the same");
     for(int index = 0; index < nofClusters; index++)
     {
	  const mxArray *clu_ptr = mxGetCell(prhs[2], index);
	  int clustersize = mxGetM(clu_ptr) * mxGetN(clu_ptr);
	  double *pcluele = mxGetPr(clu_ptr);

	  const mxArray *pot_ptr = mxGetCell(prhs[3], index);
	  int potsize = mxGetM(pot_ptr) * mxGetN(pot_ptr);
	  if(potsize == 0)
	       continue;
	  int idx = 0;
	  for(int i = 0; i < clustersize; i++)
	  {
	       int ii = static_cast<int>(pcluele[i]);
	       idx *= NofStates[ii];
	       idx += static_cast<int>(decodes[ii]);
	  }
	  v += mxGetPr(pot_ptr)[idx];
     }
     
  
     plhs[0] = mxCreateDoubleMatrix(1, 1,  mxREAL);
     * mxGetPr(plhs[0]) = v;
     delete []NofStates;
}
