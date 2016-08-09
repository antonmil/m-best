/** BPSolver.cpp --- 
 *
 * Copyright (C) 2015 Zhen Zhang
 *
 * Author: Zhen Zhang <zhen@zzhang.org>
 *
 * All rights reserved.
 * Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 * 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer. 
 * 
 * 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution. 
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <map>
#include <algorithm>
#include <fstream>
#include <cstring>
#include <iostream>
#include <cmath>
#include <list>
#include <random>
#include <queue>
#include <utility>
#include <iterator>
#include "BPSolver.h"

//#include <quadmath.h>
#ifdef MATLAB_MEX_FILE
#include "mex.h"
#endif
//#define OUTPUT_POTENTIAL
bool zzhang::CBPDualSolver::SmoothBP(int maxIter, double epsilon)
{
     m_epsilon = epsilon;
     SmoothBPInit(epsilon);
     for(int iter = 0; iter < maxIter; iter++)
     {
	  double lastdual = m_CurrentUB;
	  SmoothBPOneIter();
	  SmoothComputeObj();
	  printf("iter = %d, dual = %12.7f, non-smooth dual = %12.7f, primal = %12.7f, diff = %12.7f\n", iter, m_CurrentUB,  m_SmoothDual2, m_BestLB, lastdual - m_CurrentUB);
	  if(fabs(m_CurrentUB-lastdual) < 1e-6)
	       break;
     }
     return true;
}
void zzhang::CBPDualSolver::SmoothBPOneIter()
{
     double inv_epsilon = 1.0 / m_epsilon;
     double inv_si = 0.0;
     for(int ci = 0; ci < m_EClusterVec.size(); ci++)
     {
	  if(m_SubEClusters[ci].size() == 0) continue;
	  double *bc = m_Beliefs[ci].get();
	  double maxval = -1e20;
	  
	  for(int si = 0; si < m_SubEClusters[ci].size(); si++)
	  {
	       int sidx = m_SubEClusters[ci][si];
	       double *bs =m_Beliefs[sidx].get();
	       double *ps = new double[m_TotalDim[sidx]];
	       int *cvttable = m_IdxConvertingTable[ci][si].get();
	       for(int xc = 0; xc < m_TotalDim[ci]; xc++)
	       {
	    	    bc[xc] += bs[cvttable[xc]];
	       }
	       for(int xs = 0; xs < m_TotalDim[sidx]; xs++)
	       {
		    bs[xs] = -1e20;
		    ps[xs] = 0;
	       }
	       for(int xc = 0; xc < m_TotalDim[ci]; xc++)
	       {
		    if(bc[xc] > bs[cvttable[xc]])
			 bs[cvttable[xc]] = bc[xc];
	       }
	       for(int xc = 0; xc < m_TotalDim[ci]; xc++)
	       {
		    int xs = cvttable[xc];
		    ps[xs] += exp(inv_epsilon * (bc[xc] - bs[xs]));
	       }
	       for(int xs = 0; xs < m_TotalDim[sidx]; xs++)
	       {
		    bs[xs] += m_epsilon * log(ps[xs]);
		    bs[xs] *= 0.5;
	       }
	       for(int xc = 0; xc < m_TotalDim[ci]; xc++)
	       {
	    	    bc[xc] -= bs[cvttable[xc]];
	       }
	       delete []ps;
	  }
     }
}
bool zzhang::CBPDualSolver::SmoothBPInit(double epsilon)
{
     m_epsilon = epsilon;
     double max = -1e20;
     double min = 1e20;
     m_offset = 0.0;
     if(epsilon <= 0) return false;
     else{
	  double inv_epsilon = 1.0 / epsilon;
	  m_PseudoDistribution = std::vector< boost::shared_array<long double> > (m_Beliefs.size());
	  for(int i = 0; i < m_Beliefs.size(); i++)
	  {
	       //First step, normalise

	       m_PseudoDistribution[i] = boost::shared_array<long double>(new long double[m_TotalDim[i]]);
	       long double *pxc = m_PseudoDistribution[i].get();
	       double *bxc = m_Beliefs[i].get();
	       double maxv = -1e20;
	       for(int xc = 0; xc < m_TotalDim[i]; xc++)
	       {
		    if(bxc[xc] >= maxv) maxv = bxc[xc];
	       }
	       int alldim = m_TotalDim[i];
	       m_offset += maxv;
	       for(int xc = 0; xc < alldim; xc++)
	       {
		    long double t1 =  bxc[xc] - maxv;
		    
		    pxc[xc] = exp(inv_epsilon * t1) + 1e-40;
	       }

	  }
     }
	 return true;
}
void zzhang::CBPDualSolver::SmoothComputeObj()
{
     double inv_epsilon = 1.0 / m_epsilon;
     double dual = 0.0;
     double primal = 0.0;
     m_SmoothDual2 = 0.0;
     boost::shared_array<int> CurrentDecode(new int[m_NodeSize]);
     std::for_each(CurrentDecode.get(), CurrentDecode.get() + m_NodeSize,
		   [](int& p) {p = -1;});
     
     for(int ci = 0; ci < m_Beliefs.size(); ci++)
     {
	  double *bc = m_Beliefs[ci].get();
	  double maxval = bc[0];
	  m_LocalMaximum[ci] = 0;
	  double sum = 0.0;
	  
	  for(int xc = 1; xc < m_TotalDim[ci]; xc++)
	  {
	       if(bc[xc] > maxval)
	       {
		    maxval = bc[xc];
		    m_LocalMaximum[ci] = xc;
	       }
	  }
	  for(int xc = 0; xc < m_TotalDim[ci]; xc++)
	  {
	       sum += exp(inv_epsilon * (bc[xc] - maxval));
	  }
	  dual += m_epsilon * log(sum) + maxval;
	  m_SmoothDual2 += maxval;
	  if(m_EClusterVec[ci].size() == 1){
	       
	       CurrentDecode.get()[m_EClusterVec[ci][0]] = m_LocalMaximum[ci];
	  }
     }
     for (int *p = CurrentDecode.get(); p < CurrentDecode.get() + m_NodeSize;
	  p++) {
	  if (*p == -1) {
	       int ni = p - CurrentDecode.get();
	       int i = m_NodeToECluIdx[ni][0];
	       int t = m_LocalMaximum[i];
	       for (int jdim = m_EClusterVec[i].size() - 1; jdim >= 0; jdim--) {
		    int jidx = t % m_NodeStates.get()[m_EClusterVec[i][jdim]];
		    if (m_EClusterVec[i][jdim] == ni) {
			 CurrentDecode.get()[m_EClusterVec[i][jdim]] = jidx;
			 break;
		    }
		    t /= m_NodeStates.get()[m_EClusterVec[i][jdim]];

	       }
	  }
     }
   
     double primal2 = 0.0;
     for(int ci = 0; ci < m_EClusterVec.size(); ci++)
     {
	  int idx = 0;
	  for (int j = 0; j < m_EClusterVec[ci].size(); j++) {
	       if (j == 0)
		    idx += CurrentDecode.get()[m_EClusterVec[ci][j]];
	       else
		    idx = idx * m_NodeStates.get()[m_EClusterVec[ci][j]]
			 + CurrentDecode.get()[m_EClusterVec[ci][j]];
	  }
	  double t = m_PseudoDistribution[ci][idx];
	  primal2 += m_Beliefs[ci][idx];
     }
     m_CurrentUB = dual;
     if(primal2 > m_BestLB)
     {
	  m_BestLB = primal2;
	  m_DecodedSolution = CurrentDecode;
     }
}

void zzhang::CBPDualSolver::GenPartionGraph()
{
     
     std::vector<bool> FactorQueue(m_EClusterVec.size(), false);
     int nsize = 0;
     int cpart = 0;
     while(nsize < m_EClusterVec.size()){
	  std::vector<int> CurrentPart;
	  std::vector<bool> IsUsed(m_NodeSize, false);
	  for(int i = 0; i < m_EClusterVec.size(); i++)
	  {
	       bool canadd = true;
	       if(FactorQueue[i] == true) continue;
	       if(m_SubEClusters[i].size() == 0)
	       {
		    FactorQueue[i] = true;
		    nsize ++;
		    continue;
	       }
	       for(int j = 0; j < m_EClusterVec[i].size(); j++)
	       {
		    if(IsUsed[m_EClusterVec[i][j]])
		    {
			 canadd = false;
			 break;
		    }
	       }
	       if(canadd)
	       {
		    
		    for(int j = 0; j < m_EClusterVec[i].size(); j++)
		    {
			 IsUsed[m_EClusterVec[i][j]] = true;
		    }
		    FactorQueue[i] = true;
		    CurrentPart.push_back(i);
	       }
	  }
	  nsize += CurrentPart.size();
	  m_PartionGraph.push_back(CurrentPart);
	  cpart ++;
     }
     
}
void zzhang::CBPDualSolver::PBPOneIter()
{
     for(int pi = 0; pi < m_PartionGraph.size(); pi++)
     {
#pragma omp for 
	  for(int j = 0; j < m_PartionGraph[pi].size(); j++)
	  {
	       int i = m_PartionGraph[pi][j];
	       double tcinverce = 1.0 / m_SubEClusters[i].size();
	       int dim = m_TotalDim[i];
	       double *bigdat = m_Beliefs[i].get();
	       for (int si = 0; si < m_SubEClusters[i].size(); si++) {
		    for (int vi = 0; vi < m_TotalDim[i]; vi++) {
			 int sidx = m_SubEClusters[i][si];
			 double *smalldata = m_Beliefs[sidx].get();
			 bigdat[vi] += *(smalldata
					 + m_IdxConvertingTable[i][si].get()[vi]);
		    }
	       }
	       bigdat = m_Beliefs[i].get();
	       double LocalMax = -Huge;
	       m_LocalMaximum[i] = 0;
	       for (int vi = 0; vi < m_TotalDim[i]; vi++) {
		    if (*(bigdat + vi) > LocalMax) {
			 LocalMax = *(bigdat + vi);
			 m_LocalMaximum[i] = vi;
		    }
	       }
	       for (int si = 0; si < m_SubEClusters[i].size(); si++) {
		    int sidx = m_SubEClusters[i][si];
		    m_LocalMaximum[sidx] =
			 m_IdxConvertingTable[i][si].get()[m_LocalMaximum[i]];
		    for (int supi = 0; supi < m_SupEClustersIdx[sidx].size();
			 supi++) {
			 int supidx = m_SupEClustersIdx[sidx][supi];
			 int suppos = m_SupEClustersPos[sidx][supi];
			 if (m_LocalMaximum[supidx] == -1)
			      continue;
		    }
		    double *smalldata = m_Beliefs[sidx].get();
		    std::for_each(smalldata, smalldata + m_TotalDim[sidx],
				  [](double& p) {p=-2*Huge;});
		    for (int vi = 0; vi < m_TotalDim[i]; vi++) {
			 if (*(smalldata + m_IdxConvertingTable[i][si].get()[vi])
			     < *(bigdat + vi)) {
			      *(smalldata + m_IdxConvertingTable[i][si].get()[vi]) =
				   *(bigdat + vi);
			 }
		    }
		    std::for_each(smalldata, smalldata + m_TotalDim[sidx],
				  [&](double& p) {p*=tcinverce;});
	       }
	       for (int si = 0; si < m_SubEClusters[i].size(); si++) {
		    int sidx = m_SubEClusters[i][si];
		    double *bc = m_Beliefs[i].get();
		    double *smalldata = m_Beliefs[sidx].get();
		    for (int vi = 0; vi < m_TotalDim[i]; vi++) {
			 bc[vi] -= *(smalldata
				     + m_IdxConvertingTable[i][si].get()[vi]);
		    }
	       }
	  }
     }
}
bool zzhang::CBPDualSolver::Init() {
	srand(456);
	clock_t init_start = clock();
	std::vector<std::vector<int> > TNodeToEclusterIdx(m_NodeSize);
	std::vector<int> HasSuperSet(m_EClusterVec.size());
	std::map<std::pair<int, int>, bool> IsVisited;
	for (int i = 0; i < m_EClusterVec.size(); i++) {
		HasSuperSet[i] = 0;
		for (int j = 0; j < m_EClusterVec[i].size(); j++) {
			TNodeToEclusterIdx[m_EClusterVec[i][j]].push_back(i);
		}
	}
	
	for (int i = 0; i < m_NodeSize; i++) {
		for (int j = 0; j < TNodeToEclusterIdx[i].size(); j++) {
			int cj = TNodeToEclusterIdx[i][j];
			for (int k = j + 1; k < TNodeToEclusterIdx[i].size(); k++) {
				int ck = TNodeToEclusterIdx[i][k];
				if (cj == ck)
					continue;
				std::pair<int, int> currentjk(cj, ck);
				if (IsVisited.find(currentjk) != IsVisited.end())
					continue;
				IsVisited[currentjk] = true;
				IsVisited[std::pair<int, int>(ck, cj)] = true;
				if (m_EClusterSet[cj].size() == m_EClusterSet[ck].size())
				continue;
				if (m_EClusterSet[cj].size() < m_EClusterSet[ck].size()) {
					int temp = cj;
					cj = ck;
					ck = temp;
				}

				if (std::includes(m_EClusterSet[cj].begin(),
						m_EClusterSet[cj].end(), m_EClusterSet[ck].begin(),
						m_EClusterSet[ck].end())) {
					HasSuperSet[ck] = 1;
				}

			}
		}
	}
	IsVisited.clear();
	int maxsize = HasSuperSet.size();
	for (int i = 0; i < maxsize; i++) {
		if (HasSuperSet[i]) {
			m_SubEClusters[i].clear();
			continue;
		}
		for (int ni = 0; ni < m_EClusterVec[i].size(); ni++) {
			int nidx = m_EClusterVec[i][ni];
			for (int k = 0; k < TNodeToEclusterIdx[nidx].size(); k++) {
				int cj = i, ck = TNodeToEclusterIdx[nidx][k];
				if (cj == ck)
					continue;
				if (HasSuperSet[ck])
					continue;
				std::pair<int, int> currentjk(cj, ck);
				if (IsVisited.find(currentjk) != IsVisited.end())
					continue;
				IsVisited[currentjk] = true;
				IsVisited[std::pair<int, int>(ck, cj)] = true;
				std::set<int> intersection;
				std::set_intersection(m_EClusterSet[cj].begin(),
						m_EClusterSet[cj].end(), m_EClusterSet[ck].begin(),
						m_EClusterSet[ck].end(),
						std::inserter(intersection, intersection.begin()));
				if (m_EClusterIdx.find(intersection) != m_EClusterIdx.end())
					continue;
				std::vector<int> intvec(intersection.begin(),
						intersection.end());
				AddECluster(intvec, -1);
				//std::cout << "Added New intersection\n";
				for (int j = 0; j < intvec.size(); j++) {
					TNodeToEclusterIdx[intvec[j]].push_back(
							m_EClusterVec.size() - 1);

				}
				HasSuperSet.push_back(1);
			}
		}
	}

	IsVisited.clear();
	for (int i = 0; i < m_EClusterSet.size(); i++) {
		if (HasSuperSet[i])
			continue;
		std::map<int, bool> isvisited;
		for (int ni = 0; ni < m_EClusterSet[i].size(); ni++) {
			int nidx = m_EClusterVec[i][ni];
			for (int k = 0; k < TNodeToEclusterIdx[nidx].size(); k++) {
				int cj = i, ck = TNodeToEclusterIdx[nidx][k];
				if (cj == ck)
					continue;
				if (isvisited.find(ck) != isvisited.end())
					continue;
				isvisited[ck] = true;
				if (m_EClusterSet[cj].size() == m_EClusterSet[ck].size())
					continue;

				if (std::includes(m_EClusterSet[cj].begin(),
						m_EClusterSet[cj].end(), m_EClusterSet[ck].begin(),
						m_EClusterSet[ck].end())) {
					m_SubEClusters[cj].push_back(ck);
				}
			}

		}
	}
	clock_t time_explore_graph_end = clock();
//	printf("Explore Graph Structure Spends: %f seconds\n",
//	       (time_explore_graph_end - init_start) * 1.0 / CLOCKS_PER_SEC);

	CacheIndex();
	
	m_Potentials.clear();
	m_LocalMaximum = std::vector<int>(m_EClusterVec.size());
	m_DecodedSolution = boost::shared_array<int>(new int[m_NodeSize]);
	m_NodeToECluIdx = TNodeToEclusterIdx;
	clock_t cache_idx_cvt = clock();
	printf("Cache index converting table spends: %f seconds\n",
	       (cache_idx_cvt - time_explore_graph_end) * 1.0 / CLOCKS_PER_SEC);

	return true;
}
bool zzhang::CBPDualSolver::CacheIndex()
{
     for (int i = 0; i < m_EClusterSet.size(); i++) {
	  m_IdxConvertingTable.push_back(std::vector<boost::shared_array<int> >( m_SubEClusters[i].size()));
	  if(m_EClusterSet[i].size() == 2)
	  {
	       int rii = m_EClusterVec[i][0];
	       int rjj = m_EClusterVec[i][1];
	       for (int si = 0; si < m_SubEClusters[i].size(); si++) {
		    boost::shared_array<int> cvttable(new int[m_TotalDim[i]]);
		    int sidx = m_SubEClusters[i][si];
		    int ii = m_EClusterVec[sidx][0];
		    int *NofStates = m_NodeStates.get();
		    int *pcvttable = cvttable.get();
		    if(ii == rii)
		    {
			 for(int xi = 0; xi < NofStates[rii]; xi++)
			 {
			      for(int xj = 0; xj < NofStates[rjj]; xj++)
			      {
				   int xij = xi * NofStates[rjj] + xj;
				   pcvttable[xij] = xi;
			      }
			 }
		    }
		    else{
			 for(int xi = 0; xi < NofStates[rii]; xi++)
			 {
			      for(int xj = 0; xj < NofStates[rjj]; xj++)
			      {
				   int xij = xi * NofStates[rjj] + xj;
				   pcvttable[xij] = xj;
			      }
			 }
		    }
		    m_IdxConvertingTable[i][si] = cvttable;
		    m_SupEClustersIdx[sidx].push_back(i);
		    m_SupEClustersPos[sidx].push_back(si);
	       }
	  }
	  else
	  {
	       for (int si = 0; si < m_SubEClusters[i].size(); si++) {
		    boost::shared_array<int> cvttable;
		    int sidx = m_SubEClusters[i][si];
		    this->IdxConverting(m_EClusterVec[i], m_EClusterVec[sidx],
					cvttable);
		    m_IdxConvertingTable[i][si] = cvttable;
		    m_SupEClustersIdx[sidx].push_back(i);
		    m_SupEClustersPos[sidx].push_back(si);
	       }
	  }
	  int pidx = m_EClusterPotIdx[i];
	  if (pidx == -1) {
	       boost::shared_array<double> belief(new double[m_TotalDim[i]]);
	       std::for_each(belief.get(), belief.get() + m_TotalDim[i],
			     [](double& p) {p=0.0;});
	       m_Beliefs.push_back(belief);
	  } else {
	       m_Beliefs.push_back(m_Potentials[pidx]);
	  }

     }
	 return true;
}
void zzhang::CBPDualSolver::BPClear(int maxOuterIter, int maxInnerIter)
{
     for(int oiter = 0; oiter < maxOuterIter; oiter++)
     {
			 
	  for(int iiter = 0; iiter < maxInnerIter; iiter++)
	  {
	       for(int i = 0; i < m_EClusterVec.size(); i++)
	       {
		    m_LocalMaximum[i] = -1;
	       }
	       double lastdual = m_CurrentUB;
	       BPOneIter();
	       ComputeObj();
	       clock_t current_clk = clock();
	       printf(
		    "Iter=%d Dual=%16.10f NFPrimal=%16.10f Primal=%16.10f Diff=%16.10f IntGap=%16.10f Time=%16.8f\n",
		    iiter, m_CurrentUB, m_primal2, m_BestLB, lastdual - m_CurrentUB,
		    m_BestLB - m_CurrentUB,
		    1.0 * (current_clk - m_starttime) / CLOCKS_PER_SEC);
	       if(m_CurrentUB < m_LBBab)
		    break;
	       if (std::fabs(lastdual - m_CurrentUB) < 1e-6) {
		    break;
	       }
	       fflush(stdout);
	  }
	  if(m_CurrentUB < m_LBBab)
	       break;
	  if(std::fabs(m_CurrentUB - m_primal2) < 1e-5)
	       break;
	  if(maxOuterIter == 1) break;
	  if(oiter == maxOuterIter) break;
	  clock_t tighten_start = clock();
	  TightenCluster(15);
	  clock_t tighten_end = clock();
	  printf("Time Spend For tighten %12.5f\n", (1.0 * tighten_end - 1.0 * tighten_start) / CLOCKS_PER_SEC);
	  fflush(stdout);
     }
}
void zzhang::CBPDualSolver::DualCoordinateDescend(int MaxIter)
{
     m_starttime=clock();
     for(int jj = 0; jj < MaxIter; jj++)
     {
	  double lastdual = m_CurrentUB;
	  BPOneIter();
          for(int kki = 0; kki < 1; kki ++)
          {
	  int kstart = rand();
	  for(int kk = 0; kk < m_Phik.size(); kk++)
	  {
	       int k = (kstart + kk) % m_Phik.size();
	       DCDUpdateGamma(k);
	  }}
	  ComputeObj();
	  clock_t current_clk = clock();
#ifdef MATLAB_MEX_FILE
	  #if 0
	  mexPrintf(
	       "DCDIter=%d UB=%16.10f LB=%16.10f Diff=%16.10f IntGap=%16.10f Time=%16.8f\n",
	       jj, m_CurrentUB, m_BestLB, lastdual - m_CurrentUB,
	       m_BestLB - m_CurrentUB,
	       1.0 * (current_clk - m_starttime) / CLOCKS_PER_SEC);
	  #endif
#else
	  printf(
	       "Iter=%d UB=%16.10f LB=%16.10f Diff=%16.10f IntGap=%16.10f Time=%16.8f\n",
	       jj, m_CurrentUB, m_BestLB, lastdual - m_CurrentUB,
	       m_BestLB - m_CurrentUB,
	       1.0 * (current_clk - m_starttime) / CLOCKS_PER_SEC);
#endif
	  if (std::fabs(lastdual - m_CurrentUB) < 1e-6) {
	       break;
	  }
     }

}
void zzhang::CBPDualSolver::SGBP(int MaxIter)
{
     double m_bestdual = 1e20;
     m_starttime=clock();
     for( int sgIter = 0; sgIter < 200; sgIter++)
     {
	  double lastdual_sg = m_CurrentUB;

	  for(int i = 0; i < MaxIter; i++)
	  {
	       double lastdual = m_CurrentUB;
	       BPOneIter();
	       ComputeObj();
	       if(m_bestdual > m_CurrentUB) m_bestdual = m_CurrentUB;
	       clock_t current_clk = clock();
#ifdef MATLAB_MEX_FILE
	       #if 0
	        mexPrintf(
		     "SGIter=%d BestDual = %16.10f, Current Dual=%16.10f, Best Decoded Primal=%16.10f Diff=%16.10f IntGap=%16.10f Time=%16.8f\n",
		     i, m_bestdual, m_CurrentUB, m_BestLB, lastdual - m_CurrentUB,
		     m_BestLB - m_CurrentUB,
		     1.0 * (current_clk - m_starttime) / CLOCKS_PER_SEC);
		#endif
#else
	       printf(
		    "Iter=%d BestDual = %16.10f, Current Dual=%16.10f, Best Decoded Primal=%16.10f Diff=%16.10f IntGap=%16.10f Time=%16.8f\n",
		    i, m_bestdual, m_CurrentUB, m_BestLB, lastdual - m_CurrentUB,
		    m_BestLB - m_CurrentUB,
		    1.0 * (current_clk - m_starttime) / CLOCKS_PER_SEC);
#endif
	       if (std::fabs(lastdual - m_CurrentUB) < 1e-6) {
		    break;
	       }
		 if (std::fabs(m_bestdual - m_BestLB) < 1e-6) {
		break;
	  }
	
	  }
	  if (std::fabs(lastdual_sg - m_CurrentUB) < 1e-6) {
	       break;
	  }
          if (std::fabs(m_bestdual - m_BestLB) < 1e-6) {
		break;
	  }
	  std::vector<double> SG(m_Phik.size(), 0.0);
	  double l2norm = 0.0;
	  for(int k = 0; k < m_Phik.size(); k++)
	  {
	       double sg = 0;
	       for(int cii = 0; cii < ActiveConstraints[k].size(); cii++)
	       {
		    int ci = ActiveConstraints[k][cii];
		    double *phic = m_Phik[k][ci].get();
		    sg -= phic[m_LocalMaximum[ci]];
	       }
	       SG[k] = sg;
	       printf("Current SG[%d]: %f\n", k, sg);
	       l2norm += sg * sg;
	  }
	  double step = 0.0;
	  if(fabs(l2norm) >= 1e-6)
               step = 10.0 / (sgIter + 1)/sqrt(l2norm);
          else 
               step = 1.0;   
	  for(int k = 0; k < m_Phik.size(); k++)
	  {
	       double oldgamma = m_Gamma[k];
	       m_Gamma[k] -=  step * SG[k] ;
	       if(m_Gamma[k] < 0) m_Gamma[k] = 0;
	       double DeltaGamma = (m_Gamma[k] - oldgamma);
	       for(int cii = 0; cii < ActiveConstraints[k].size(); cii++)
	       {
		    int ci = ActiveConstraints[k][cii];
		    double *bc = m_Beliefs[ci].get();
		    double *phic = m_Phik[k][ci].get();
		    for(int vi = 0; vi < m_TotalDim[ci]; vi++)
		    {
			 bc[vi] -= DeltaGamma * phic[vi];
		    }
	       }
	  }
     }
}
void zzhang::CBPDualSolver::DCDUpdateGamma(int k)
{
     double cgamma = 0.1;
     double L = 0.0;
     double R = cgamma;
     if(m_Gamma[k] !=0)
     {
	  L = 0;
	  R = m_Gamma[k] + 0.1;
	  cgamma = R;
     }
     double lsg = -m_Gamma[k]; double rsg = 100;
     clock_t s1 = clock();
     while(1)
     {
	  double sg = 0;
	  double DeltaGamma = (cgamma - m_Gamma[k]);
	  m_Gamma[k] = cgamma;
	  for(int cii = 0; cii < ActiveConstraints[k].size(); cii++)
	  {
	       int ci = ActiveConstraints[k][cii];
	       double *bc = m_Beliefs[ci].get();
	       double *phic = m_Phik[k][ci].get();
	       double maxv = -1e20;
	       int maxidx = -1;
	       for(int vi = 0; vi < m_TotalDim[ci]; vi++)
	       {
		    bc[vi] -= DeltaGamma * phic[vi];
		    if(bc[vi]  > maxv)
		    {
			 maxv = bc[vi];
			 maxidx = vi;
		    }
	       }
	       sg -= phic[maxidx];
	  }
	  rsg = sg;
	  if(sg >= 0)
	       break;
	  else
	  {
	       lsg = sg;
	       L = cgamma;
	       cgamma += cgamma;
	       R = cgamma;
	  }
     }
#ifdef MATLAB_MEX_FILE
//     mexPrintf("Find Right bound spend: %12.7f", (clock() - s1) * 1.0 / CLOCKS_PER_SEC);
     clock_t s2 = clock();
#endif
     while(abs(L - R) > 1e-6)
     {
	  double sg = 0.0;
	  cgamma = (L * (fabs(rsg) + 100.0) + R * (fabs(lsg) + 100.0)) / (fabs(lsg) + fabs(rsg) + 200);
	  //cgamma = (L + R) / 2;
	  double DeltaGamma = (cgamma - m_Gamma[k]);
	  m_Gamma[k] = cgamma;
	  for(int cii = 0; cii < ActiveConstraints[k].size(); cii++)
	  {
	       int ci = ActiveConstraints[k][cii];
	       double *bc = m_Beliefs[ci].get();
	       double *phic = m_Phik[k][ci].get();
	       double maxv = -1e20;
	       int maxidx = -1;
	       for(int vi = 0; vi < m_TotalDim[ci]; vi++)
	       {
		    bc[vi] -= DeltaGamma * phic[vi];
		    if(bc[vi]  > maxv)
		    {
			 maxv = bc[vi];
			 maxidx = vi;
		    }
	       }
	       sg -= phic[maxidx];
	       m_LocalMaximum[ci] = maxidx;
	  }
	  if(fabs(sg) < 1e-6)
	       sg = 0.0;
	  if(sg == 0)
	       break;
	  if(sg > 0)
	  {
	       rsg = sg;
	       R = cgamma;
	  }
	  else
	  {
	       lsg = sg;
	       L = cgamma;
	  }
     }
#ifdef MATLAB_MEX_FILE
//     mexPrintf("Find gamma bound spend: %12.7f", (clock() - s2) * 1.0 / CLOCKS_PER_SEC);
#endif
}
void zzhang::CBPDualSolver::ConstrainedBP(int maxOuterIter, int maxInnerIter)
{
     int NofConstraints = m_Phik.size();
     double gamma0 = 1;
     m_starttime = clock();
     double up_scale = 4;
     for(int cditer = 0; cditer < 10; cditer++)
     {
	  int rstart = rand();
	  for(int rk = 0; rk < NofConstraints; rk++)
	  {
	       int k = (rk + rstart) % NofConstraints;
	       double L = -m_Gamma[k];
	       double R = gamma0;
	       for(int j = 0; j < m_Phik[k].size(); j++)
	       {
		    if(m_Phik[k][j] == NULL)
			 continue;
		    double *phic = m_Phik[k][j].get();
		    double *bc = m_Beliefs[j].get();
		    for(int xi = 0; xi < m_TotalDim[j]; xi++)
		    {
			 bc[xi] -= R * phic[xi];
		    }
	       }
	       m_Gamma[k] = R;
	       m_primal2 = -1e20;
	       double lsg = -1.0;
	       double rsg = 1.0;
	       while(1)
	       {
		    m_primal2 = -1e20;
		    BPClear(maxOuterIter, maxInnerIter);
		    double sg = -IsFeasible[k];
		    rsg = sg;
		    printf("Current SG: %12.5f\n", sg);
		    if(sg < 0)
		    {
			 for(int j = 0; j < m_Phik[k].size(); j++)
			 {
			      if(m_Phik[k][j] == NULL)
				   continue;
			      double *phic = m_Phik[k][j].get();
			      double *bc = m_Beliefs[j].get();
			      for(int xi = 0; xi < m_TotalDim[j]; xi++)
			      {
				   bc[xi] -= (up_scale - 1) * R * phic[xi];
			      }
			 }
			 lsg = sg;
			 L = R;
			 R = R * up_scale;
			 m_Gamma[k] = R;

		    }
		    else{
			 break;
		    }
	       }
	       double LastGamma = R;
	       while(fabs(L-R) > 1e-10)
	       {
		    double gamma = (L * fabs(rsg) + R * fabs(lsg)) / (fabs(rsg) + fabs(lsg));
		    double deltaGamma = (gamma - LastGamma);
		    m_Gamma[k] = gamma;
		    LastGamma = gamma;
		    for(int j = 0; j < m_Phik[k].size(); j++)
		    {
			 if(m_Phik[k][j] == NULL)
			      continue;
			 double *phic = m_Phik[k][j].get();
			 double *bc = m_Beliefs[j].get();
			 for(int xi = 0; xi < m_TotalDim[j]; xi++)
			 {
			      bc[xi] -= deltaGamma * phic[xi];
			 }
		    }
		    m_primal2 = -1e20;
		    BPClear(maxOuterIter, maxInnerIter);
		    if(fabs(m_CurrentUB - m_BestLB) < 1e-6)
			 break;
		    double sg = -IsFeasible[k];
		    
		    printf("Current Gamma[%d]: %12.5f\n", k, gamma);
		    if((fabs(m_CurrentUB - m_primal2) < 1e-6) && (fabs(sg) < 1e-6))
			 break;
		    if(sg < 0)
		    {
			 lsg = sg;
			 L = gamma;
		    }
		    else if(sg > 0)
		    {
			 rsg = sg;
			 R = gamma;
		    }
		    else{
			 break;
		    }
	       }
	  }
	  if(NofConstraints == 1)
	       break;
     }
}

void zzhang::CBPDualSolver::BPOneIter()
{
     int start = rand();
     for (int ii = m_EClusterVec.size() - 1; ii >= 0; ii--) {
	  int i = (ii + start) % m_EClusterVec.size();
	  if (m_SubEClusters[i].size() == 0)
	       continue;

	  double tcinverce = 1.0 / m_SubEClusters[i].size();
	  int dim = m_TotalDim[i];
	  double *bigdat = m_Beliefs[i].get();
	  for (int si = 0; si < m_SubEClusters[i].size(); si++) {
	       for (int vi = 0; vi < m_TotalDim[i]; vi++) {
		    int sidx = m_SubEClusters[i][si];
		    double *smalldata = m_Beliefs[sidx].get();
		    bigdat[vi] += *(smalldata
				    + m_IdxConvertingTable[i][si].get()[vi]);
	       }
	  }
	  bigdat = m_Beliefs[i].get();
	  double LocalMax = -Huge;
	  m_LocalMaximum[i] = 0;
	  for (int vi = 0; vi < m_TotalDim[i]; vi++) {
	       if (*(bigdat + vi) > LocalMax) {
		    LocalMax = *(bigdat + vi);
		    m_LocalMaximum[i] = vi;
	       }
	  }
	  for (int si = 0; si < m_SubEClusters[i].size(); si++) {
	       int sidx = m_SubEClusters[i][si];
	       m_LocalMaximum[sidx] =
		    m_IdxConvertingTable[i][si].get()[m_LocalMaximum[i]];
	       for (int supi = 0; supi < m_SupEClustersIdx[sidx].size();
		    supi++) {
		    int supidx = m_SupEClustersIdx[sidx][supi];
		    int suppos = m_SupEClustersPos[sidx][supi];
		    if (m_LocalMaximum[supidx] == -1)
			 continue;
	       }
	       double *smalldata = m_Beliefs[sidx].get();
	       std::for_each(smalldata, smalldata + m_TotalDim[sidx],
			     [](double& p) {p=-2*Huge;});
	       for (int vi = 0; vi < m_TotalDim[i]; vi++) {
		    if (*(smalldata + m_IdxConvertingTable[i][si].get()[vi])
			< *(bigdat + vi)) {
			 *(smalldata + m_IdxConvertingTable[i][si].get()[vi]) =
			      *(bigdat + vi);
		    }
	       }
	       std::for_each(smalldata, smalldata + m_TotalDim[sidx],
			     [&](double& p) {p*=tcinverce;});
	  }
	  for (int si = 0; si < m_SubEClusters[i].size(); si++) {
	       int sidx = m_SubEClusters[i][si];
	       double *bc = m_Beliefs[i].get();
	       double *smalldata = m_Beliefs[sidx].get();
	       for (int vi = 0; vi < m_TotalDim[i]; vi++) {
		    bc[vi] -= *(smalldata
				+ m_IdxConvertingTable[i][si].get()[vi]);
	       }
	  }
     }

     // #if 0
     //NodeBased BP
     start = rand();
     for(int ii = 0; ii < m_EClusterVec.size(); ii++)
     {
	  int i = (ii + start) % m_EClusterVec.size();
	  if(m_SupEClustersIdx[i].size() == 0)
	       continue;
	  double *bi_tmp = new double[m_TotalDim[i]];
	  
	  double scale = 1.0 / (m_SupEClustersIdx[i].size() + 1);
	  double *bs = m_Beliefs[i].get();
	  for(int supi = 0; supi < m_SupEClustersIdx[i].size(); supi++)
	  {
	       int sup_idx = m_SupEClustersIdx[i][supi];
	       int sub_idx = m_SupEClustersPos[i][supi];
	       for(int xi = 0; xi < m_TotalDim[i]; xi++)
		    bi_tmp[xi] = -1e20;
	       double *bc = m_Beliefs[sup_idx].get();
	       int *cvttable = m_IdxConvertingTable[sup_idx][sub_idx].get();
	       for(int xi = 0; xi < m_TotalDim[sup_idx]; xi++)
	       {
		    if(bc[xi] >= bi_tmp[cvttable[xi]])
			 bi_tmp[cvttable[xi]] = bc[xi];
	       }
	       for(int xi = 0; xi < m_TotalDim[sup_idx]; xi++)
	       {
		    bc[xi] -= bi_tmp[cvttable[xi]];
	       }
	       for(int xi = 0; xi < m_TotalDim[i]; xi++)
		    bs[xi] += bi_tmp[xi];
	  }
	  double cmax = -1e20;
	  for(int xi = 0; xi < m_TotalDim[i]; xi++)
	  {
	       if(bs[xi] >= cmax)
	       {
		    cmax = bs[xi];
		    m_LocalMaximum[i] = xi;
	       }
	       bs[xi] *= scale;
	  }
	  for(int supi = 0; supi < m_SupEClustersIdx[i].size(); supi++)
	  {
	       cmax = -1e20;
	       int sup_idx = m_SupEClustersIdx[i][supi];
	       int sub_idx = m_SupEClustersPos[i][supi];
	       double *bc = m_Beliefs[sup_idx].get();
	       int *cvttable = m_IdxConvertingTable[sup_idx][sub_idx].get();
	       for(int xi = 0; xi < m_TotalDim[sup_idx]; xi++)
	       {
		    bc[xi] += bs[cvttable[xi]];
		    if(bc[xi] >= cmax)
		    {
			 cmax = bc[xi];
			 m_LocalMaximum[sup_idx] = xi;
		    }
	       }
	  }
	  delete []bi_tmp;
     }

     //#endif
}
bool zzhang::CBPDualSolver::ParallelBP(int maxIter)
{
     	m_starttime = clock();
	int ri = 0;
	std::for_each(m_LocalMaximum.begin(), m_LocalMaximum.end(),
			[](int &p) {p=-1;});
	int tighten = 0;
	int veryhigh = 0;
	while (std::fabs(m_CurrentUB - m_BestLB) > 1e-6) {
	     m_PartionGraph.clear();
	     GenPartionGraph();

		static int first = 0;
		int it = first;
		for (it = 0; it < maxIter; it++) {
			double lastdual = m_CurrentUB;
			PBPOneIter();
			ComputeObj();	       
			clock_t current_clk = clock();
			printf(
			     "Iter=%d Dual=%16.10f Primal=%16.10f Diff=%16.10f IntGap=%16.10f Time=%16.8f\n",
			     it, m_CurrentUB, m_BestLB, lastdual - m_CurrentUB,
			     m_BestLB - m_CurrentUB,
			     1.0 * (current_clk - m_starttime) / CLOCKS_PER_SEC);
			fflush(stdout);
			if (std::fabs(lastdual - m_CurrentUB) < 1e-8) {
			     break;
			}
		}


		int ncluster = 20;
		double bound;
		static bool isfirst = true;
		if(std::fabs(m_CurrentUB - m_BestLB) < 1e-6)
			break;
		if (isfirst) {
			isfirst = false;
			if (abs(m_BestLB - m_CurrentUB) < 1)
				maxIter = 600;
			else
				maxIter = 20;
		}
		if (abs(m_BestLB - m_CurrentUB) < 1)
			maxIter = 600;
		clock_t tighten_start = clock();
		TightenCluster(ncluster);
		veryhigh = 0;
		clock_t tighten_end = clock();
		printf("Time Spend For tighten %12.5f\n", (1.0 * tighten_end - 1.0 * tighten_start) / CLOCKS_PER_SEC);
	        
	}
	return true;
	
}
bool zzhang::CBPDualSolver::BeliefPropagation(int maxIter) {
	m_starttime = clock();
	int ri = 0;
	std::for_each(m_LocalMaximum.begin(), m_LocalMaximum.end(),
			[](int &p) {p=-1;});
	int tighten = 0;
	int veryhigh = 0;
	while (std::fabs(m_CurrentUB - m_BestLB) > 1e-6) {
		static int first = 0;
		int it = first;
		for (it = 0; it < maxIter; it++) {
			double lastdual = m_CurrentUB;
			BPOneIter();
			ComputeObj();	       
#ifdef OUTPUT_POTENTIAL
			if(it == 100)
			{
				FILE *fp = fopen("Output.UAI.LG","w");
				fprintf(fp, "MARKOV\n");
				fprintf(fp, "%d", m_NodeSize);
				for(int ni = 0; ni < m_NodeSize; ni++) fprintf(fp, " %d", m_NodeStates.get()[ni]);
				fprintf(fp, "\n%d", m_EClusterVec.size());
				for(int ci = 0; ci < m_EClusterVec.size(); ci++)
				{
					fprintf(fp, "\n%d", m_EClusterVec[ci].size());
					for(int ni = 0; ni < m_EClusterVec[ci].size(); ni++)
					fprintf(fp, " %d", m_EClusterVec[ci][ni]);
				}
				for(int pi = 0; pi < m_Beliefs.size(); pi++)
				{
					fprintf(fp, "\n%d", m_TotalDim[pi]);
					for(int vi = 0; vi < m_TotalDim[pi]; vi++)
					{
						fprintf(fp, " %f", m_Beliefs[pi].get()[vi]);
					}
				}
				fclose(fp);
			}
#endif
			clock_t current_clk = clock();
			printf(
			     "Iter=%d Dual=%16.10f Primal=%16.10f Diff=%16.10f IntGap=%16.10f Time=%16.8f\n",
			     it, m_CurrentUB, m_BestLB, lastdual - m_CurrentUB,
			     m_BestLB - m_CurrentUB,
			     1.0 * (current_clk - m_starttime) / CLOCKS_PER_SEC);
			fflush(stdout);
			if (std::fabs(lastdual - m_CurrentUB) < 1e-8) {
			     break;
			}
		}


		int ncluster = 20;
		double bound;
		static bool isfirst = true;
		if(std::fabs(m_CurrentUB - m_BestLB) < 1e-6)
			break;

		SmoothBP(1000, 0.000001);
		
		if (isfirst) {
			isfirst = false;
			if (abs(m_BestLB - m_CurrentUB) < 1)
				maxIter = 600;
			else
				maxIter = 20;
		}
		if (abs(m_BestLB - m_CurrentUB) < 1)
			maxIter = 600;
		clock_t tighten_start = clock();
		TightenCluster(ncluster);
		veryhigh = 0;
		clock_t tighten_end = clock();
		printf("Time Spend For tighten %12.5f\n", (1.0 * tighten_end - 1.0 * tighten_start) / CLOCKS_PER_SEC);
	        
	}
	return true;
}
int zzhang::CBPDualSolver::TightenCluster(int nclus_to_add) {
	std::vector<std::vector<int> > PossibleClusters;
	std::vector<std::set<int> > PossibleClustersSet;
	std::vector<std::vector<int> > pSubEclusters;
	std::vector<std::vector<boost::shared_array<int> > > pIdxConvertingTable;
	std::map<std::set<int>, bool> IsVisited;
	std::vector<std::pair<int, double> > pClusterWeight;
	std::vector<int> pClusterDim;
	//Find all possible extended clusters
	std::set< std::pair<int,int> > PossibleUnion;

	for(int i = 0; i < m_NodeSize; i++)
	{
	     for(int j = 0; j < m_NodeToECluIdx[i].size(); j++)
	     {
		  int jidx = m_NodeToECluIdx[i][j];
		  if(m_EClusterVec[jidx].size() == 1) continue;
		  for(int k = j + 1; k < m_NodeToECluIdx[i].size(); k++)
		  {
		       int kidx = m_NodeToECluIdx[i][k];
		       if(m_EClusterVec[kidx].size() == 1) continue;
		       PossibleUnion.insert(std::pair<int,int>(jidx,kidx));
		  }
	     }
	}

	for(std::set< std::pair<int,int> >::iterator it = PossibleUnion.begin();
	    it != PossibleUnion.end();
	    it++)
	{
	     int jidx = it->first;
	     int kidx = it->second;
	     int ci = it->first;
	     int cj = it->second;
	     int locali = m_LocalMaximum[jidx];
	     int localj = m_LocalMaximum[kidx];
	     bool willadd = false;
	     for(int si = 0; si < m_SubEClusters[jidx].size(); si++)
	     {
		  int sidxi = m_SubEClusters[jidx][si];
		  for(int sk = 0; sk < m_SubEClusters[kidx].size(); sk++)
		  {
		       int sidxk = m_SubEClusters[kidx][sk];
		       if(sidxi == sidxk )
			    if(m_IdxConvertingTable[jidx][si].get()[locali]
			       != m_IdxConvertingTable[kidx][sk].get()[localj])
			    {
				 willadd = true;
				 std::set<int> unionset;
				 std::set_union(m_EClusterSet[ci].begin(),
						m_EClusterSet[ci].end(), m_EClusterSet[cj].begin(),
						m_EClusterSet[cj].end(),
						std::inserter(unionset, unionset.begin()));
				 if (unionset.size() > 8)
				      continue;
				 if (m_EClusterIdx.find(unionset) != m_EClusterIdx.end())
				      continue;
				 if (IsVisited.find(unionset) != IsVisited.end())
				      continue;
				 IsVisited[unionset] = true;
				 std::vector<int> unionvec(unionset.begin(), unionset.end());
				 PossibleClusters.push_back(unionvec);
				 PossibleClustersSet.push_back(unionset);
				 int dim = 1;
				 for (int ni = 0; ni < unionvec.size(); ni++) {
				      dim *= m_NodeStates.get()[unionvec[ni]];
				 }
				 pClusterDim.push_back(dim);
				 break;
			    }
		       if(willadd) break;
		  }
	     }
	}

	#if 0
	for (int i = 0; i < m_EClusterSet.size(); i++) {
		for (int supi = 0; supi < m_SupEClustersIdx[i].size(); supi++) {

			for (int supj = supi + 1; supj < m_SupEClustersIdx[i].size();
					supj++) {
				int ci = m_SupEClustersIdx[i][supi];
				int cj = m_SupEClustersIdx[i][supj];
				int locali = m_LocalMaximum[ci];
				int localj = m_LocalMaximum[cj];
				int si = m_SupEClustersPos[i][supi];
				int sj = m_SupEClustersPos[i][supj];
				
				
				if (m_IdxConvertingTable[ci][si].get()[locali]
						!= m_IdxConvertingTable[cj][sj].get()[localj]) {
					std::set<int> unionset;
					std::set_union(m_EClusterSet[ci].begin(),
							m_EClusterSet[ci].end(), m_EClusterSet[cj].begin(),
							m_EClusterSet[cj].end(),
							std::inserter(unionset, unionset.begin()));
					if (unionset.size() > 8)
						continue;
					if (m_EClusterIdx.find(unionset) != m_EClusterIdx.end())
						continue;
					if (IsVisited.find(unionset) != IsVisited.end())
						continue;
					IsVisited[unionset] = true;
					std::vector<int> unionvec(unionset.begin(), unionset.end());
					PossibleClusters.push_back(unionvec);
					PossibleClustersSet.push_back(unionset);
					int dim = 1;
					for (int ni = 0; ni < unionvec.size(); ni++) {
						dim *= m_NodeStates.get()[unionvec[ni]];
					}
					pClusterDim.push_back(dim);
				}
			}
		}
	}


#endif
	pIdxConvertingTable = std::vector< std::vector<boost::shared_array<int> > >(PossibleClustersSet.size());
	pSubEclusters = std::vector< std::vector<int> > (PossibleClustersSet.size());

#pragma omp for 
	//Find sub-clusters for possible extended clusters
	for (int i = 0; i < PossibleClustersSet.size(); i++) {
	     //	pIdxConvertingTable.push_back(std::vector<boost::shared_array<int> >());
	     //	pSubEclusters.push_back(std::vector<int>());
		std::map<int, bool> IsAdded;
		std::set<int> PossibleSubClusters;
		std::vector<int> CurrentFactor(m_NodeSize, 0);
		for (int ni = 0; ni < PossibleClusters[i].size(); ni++) {
			int nidx = PossibleClusters[i][ni];
			CurrentFactor[nidx] = 1;
			for (int k = 0; k < m_NodeToECluIdx[nidx].size(); k++) {
				int ck = m_NodeToECluIdx[nidx][k];
				PossibleSubClusters.insert(ck);
/*				if (IsAdded.find(ck) != IsAdded.end())
					continue;
				IsAdded[ck] = true;
				if (std::includes(PossibleClustersSet[i].begin(),
						PossibleClustersSet[i].end(), m_EClusterSet[ck].begin(),
						m_EClusterSet[ck].end())) {
					pSubEclusters[i].push_back(ck);
					boost::shared_array<int> cvttable;
					IdxConverting(PossibleClusters[i], m_EClusterVec[ck],
							cvttable);
					pIdxConvertingTable[i].push_back(cvttable);
					}*/
			}
		}
		for(std::set<int>::iterator it = PossibleSubClusters.begin();
		    it != PossibleSubClusters.end();
		    it++)
		{
		     int sidx = *it;
		     bool IsSub = true;
		     for(int ni = 0; ni < m_EClusterVec[sidx].size(); ni++)
		     {
			  if(CurrentFactor[m_EClusterVec[sidx][ni]] == 0)
			  {
			       IsSub = false;
			       continue;
			  }
		     }
		     if(IsSub)
		     {
			  pSubEclusters[i].push_back(sidx);
			  boost::shared_array<int> cvttable;
			  IdxConverting(PossibleClusters[i], m_EClusterVec[sidx],
					cvttable);
			  pIdxConvertingTable[i].push_back(cvttable);
		     }
		}
			 
		
	}
	pClusterWeight = std::vector< std::pair<int, double> >(PossibleClustersSet.size());

	#pragma omp for
	//Compute weight for all possible exteneded clusters
	for (int i = 0; i < PossibleClustersSet.size(); i++) {
	        boost::shared_array<double> all(new double[pClusterDim[i]]);
		memset(all.get(), 0, sizeof(double) * pClusterDim[i]);
		double *bigdat = all.get();
		double local_ub = 0.0;
		double local_exact = -Huge;
		for (int si = 0; si < pSubEclusters[i].size(); si++) {
			int sidx = pSubEclusters[i][si];
			double *smalldat = m_Beliefs[sidx].get();
			local_ub += smalldat[m_LocalMaximum[sidx]];
		}
		for (int vi = 0; vi < pClusterDim[i]; vi++) {
			for (int si = 0; si < pSubEclusters[i].size(); si++) {
				int sidx = pSubEclusters[i][si];
				double *smalldat = m_Beliefs[sidx].get();
				*bigdat += *(smalldat + pIdxConvertingTable[i][si].get()[vi]);
			}
			if (local_exact < *bigdat)
				local_exact = *bigdat;
			bigdat++;
		}
		pClusterWeight[i] = 
		     std::pair<int, double>(i, local_ub - local_exact);
	}
	int rclusters = std::min(nclus_to_add, (int) PossibleClustersSet.size());
	std::partial_sort(pClusterWeight.begin(),
			pClusterWeight.begin() + rclusters, pClusterWeight.end(),
			[](std::pair<int, double>& p, std::pair<int, double>& q) {
				return p.second > q.second;});
	if (rclusters > 0) {
		printf("Max decrease : %f\n", pClusterWeight[0].second);
	}

	for (int i = 0; i < rclusters; i++) {
		int cidx = pClusterWeight[i].first;
		AddECluster(PossibleClusters[cidx], -1);
		printf("Cluster Added: ");
		for(int i = 0; i < PossibleClusters[cidx].size(); i++)
		{
		     printf(" %d", PossibleClusters[cidx][i]);
		     m_NodeToECluIdx[PossibleClusters[cidx][i]].push_back(m_EClusterVec.size() - 1);
		}
		printf("\n");
		m_IdxConvertingTable.push_back(std::vector<boost::shared_array<int> >());
		int rcidx = m_EClusterVec.size() - 1;
		m_SubEClusters[rcidx] = pSubEclusters[cidx];
		m_IdxConvertingTable[rcidx] = pIdxConvertingTable[cidx];
	        boost::shared_array<double> belief(new double[m_TotalDim[rcidx]],
				std::default_delete<double[]>());
		std::for_each(belief.get(), belief.get() + m_TotalDim[rcidx],
				[](double& p) {p=0.0;});
		m_Beliefs.push_back(belief);
		for (int si = 0; si < m_SubEClusters[rcidx].size(); si++) {
			int sidx = m_SubEClusters[rcidx][si];
			m_SupEClustersIdx[sidx].push_back(rcidx);
			m_SupEClustersPos[sidx].push_back(si);
		}
		m_LocalMaximum.push_back(-1);
	}

	return rclusters;
}
