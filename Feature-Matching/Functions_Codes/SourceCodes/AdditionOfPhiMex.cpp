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
     if(nrhs != 3)
     {
	  mexErrMsgIdAndTxt("AdditionOfPhi:InvalidNumverOfInputs",
			    "Need four inputs Potentials, Phi, gamma");
     }
     double gamma = mxGetScalar(prhs[2]);
   
     
     mwSize nofClusters = mxGetNumberOfElements(prhs[0]);
     if(mxGetNumberOfElements(prhs[0]) < mxGetNumberOfElements(prhs[1]))
	  mexErrMsgIdAndTxt("AdditionOfPhi:InvalidInput",
			    "The dimension of c and thetac must be the same");
     plhs[0] = mxCreateCellArray(1, &nofClusters);
     mwSize nofPhis = mxGetNumberOfElements(prhs[1]);
     for(int index = 0; index < nofClusters; index++)
     {
	  const mxArray *theta_ptr = mxGetCell(prhs[0], index);
	  int thetasize = mxGetM(theta_ptr) * mxGetN(theta_ptr);
	  double *ptheta = mxGetPr(theta_ptr);
	  mxArray *pot = mxCreateDoubleMatrix(1, thetasize, mxREAL);
	  mxSetCell(plhs[0], index, pot);

	  memcpy(mxGetPr(pot), ptheta, sizeof(double) * thetasize);
	  double *nptheta = mxGetPr(pot);
	  if(index < nofPhis)
	  {
	  
	       const mxArray *phi_ptr = mxGetCell(prhs[1], index);
	       int phisize = mxGetM(phi_ptr) * mxGetN(phi_ptr);
	       double *pphi = mxGetPr(phi_ptr);
	       if(phisize == 0)
		    continue;
	       for(int i = 0; i < thetasize; i++)
	       {
		    nptheta[i] -= gamma * pphi[i];
	       }
	  }
     }
}
