/** HungarianBP.cpp --- 
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
#include <cmath>
#include <iterator>
#include "HungarianBP.h"
#ifdef MATLAB_MEX_FILE
#include "mex.h"
#endif
#ifndef EPS
#define EPS 1e-2
#endif

bool zzhang::CHungarianBP::Init()
{
     int N = m_NodeStates.get()[0];
     for(int i = 0; i < m_NodeSize; i++)
     {
	  std::vector<int> clu(1);
	  clu[0] = i ;
	  AddECluster(clu, -1);
     }
#ifndef WIN32
     unsigned char IsEdges[N * N];
	 memset(IsEdges, 0, sizeof(char) * N * N);
#else
	 boost::shared_array<unsigned char> IsEdges(new unsigned char[N * N]);
	 memset(IsEdges.get(), 0, sizeof(char) * N * N);
#endif
     
     bool res = zzhang::CPiSBPDualSolver::Init();
     if(!res) return res;
     m_bi = std::vector< boost::shared_array<double> > (m_NodeSize);
     for(int i = 0; i < m_EClusterVec.size(); i++)
     {
	  if(m_EClusterVec[i].size() == 1)
	  {
	       int nidx = m_EClusterVec[i][0];
	       m_bi[nidx] = m_Beliefs[i];
	  }
	  if(m_EClusterVec[i].size() == 2)
	  {
	       int ni = m_EClusterVec[i][0];
	       int nj = m_EClusterVec[i][1];
	       IsEdges[ni * N + nj] = 1;
	       IsEdges[nj * N + ni] = 1;
	  }
     }
     m_ui = boost::shared_array<double>(new double[N]);
     m_vi = boost::shared_array<double>(new double[N]);
     m_CurrentFdecode = boost::shared_array<int>(new int[m_NodeSize]);
     m_DecodedSolution = boost::shared_array<int>(new int[m_NodeSize]);
     memset(m_ui.get(), 0 , sizeof(double) * N);
     memset(m_vi.get(), 0 , sizeof(double) * N);
     for(int i = 0; i < m_NodeSize; i++)
     {
	  m_CurrentFdecode.get()[i] = i;
     }
     m_sum_uv = 0.0;
     return true;
}

void zzhang::CHungarianBP::ComputeObj()
{
     double dual = 0.0;
     int NofConstraints = m_Phik.size();
     IsFeasible =      std::vector<double>(NofConstraints, 0.0);
     for (int i = 0; i < m_Beliefs.size(); i++) {
	  if (m_EClusterVec[i].size() == 1)
	  {
	       int idx = m_EClusterVec[i][0];
	       
	       m_LocalMaximum[i] = m_CurrentFdecode[idx];
	  }
	  if (m_LocalMaximum[i] == -1) {
	       double v = -2 * Huge;
	       double* pstart = m_Beliefs[i].get();
	       for (double *p = pstart; p < pstart + m_TotalDim[i]; p++) {
		    if (*p > v) {
			 v = *p;
			 m_LocalMaximum[i] = static_cast<int>(p - pstart);
		    }
	       }
	  }
	  dual += m_Beliefs[i].get()[m_LocalMaximum[i]];
     }
     dual += m_sum_uv;
     if(m_CurrentUB > dual) m_CurrentUB = dual;
     double primal = 0;
     
 
     for (int i = 0; i < m_EClusterVec.size(); i++) {
	  if(m_EClusterVec[i].size() == 1)
	  {
	       continue;
	  }
	  int idx = 0;
	  for (int j = 0; j < m_EClusterVec[i].size(); j++) {
	       if (j == 0)
		    idx += m_CurrentFdecode.get()[m_EClusterVec[i][j]];
	       else
		    idx = idx * m_NodeStates.get()[m_EClusterVec[i][j]]
			 + m_CurrentFdecode.get()[m_EClusterVec[i][j]];
	  }
	  primal += m_Beliefs[i].get()[idx];

     }
     primal += m_sum_uv;
     if (primal - m_BestLB > 1e-10) {
	  m_BestLB = primal;
	  memcpy(m_DecodedSolution.get(), m_CurrentFdecode.get(), sizeof(int) * m_NodeSize);
     }
}

void zzhang::CHungarianBP::UpdateUV()
{
     int dim = m_NodeStates.get()[0];
     cost **assigncost,  lapcost;
     cost *nu, *nv;
     nu = new cost[dim];
     nv = new cost[dim];
     row i, *colsol;
     col j, *rowsol;
     assigncost = new cost*[dim];
     for (i = 0; i < dim; i++)
	  assigncost[i] = new cost[dim];
     
     rowsol = new col[dim];
     colsol = new row[dim];
     double * u = m_ui.get();
     double * v = m_vi.get();
     for (i = 0; i < m_NodeSize; i++)
     {
	  double * bi = m_bi[i].get();
	  for(int j = 0; j < dim; j++)
	  {
	       bi[j] += u[i] + v[j];
	       assigncost[i][j] = -bi[j] * 10000;
	  }
     }
     for (i = m_NodeSize; i < dim; i++)
     {
	  for(int j = 0; j < dim; j++)
	  {
	       assigncost[i][j] = 0;
	  }
     }
     memset(nu,0,sizeof(cost) * dim);
     memset(nv,0,sizeof(cost) * dim);
     lapcost = lap(dim, assigncost, rowsol, colsol, nu, nv);
     m_sum_uv = 0;
     for(int i = 0; i < dim; i++)
     {
	  u[i] = -double(nu[i]) / 10000; v[i] = - double(nv[i]) / 10000;
	  m_sum_uv += u[i] + v[i];
     }
     int *Cdecode = m_CurrentFdecode.get();
     for(int i = 0; i < m_NodeSize; i++)
     {
	  double * bi = m_bi[i].get();
	  Cdecode[i] = rowsol[i];
	  for(int j = 0; j < dim; j++)
	  {
	       bi[j] -= u[i] + v[j];
	  }
     }
     for(i = 0; i < dim; i++)
     {
	  
	  delete [] assigncost[i];
     }
     
     delete [] rowsol;
     delete [] colsol;
     delete [] assigncost;
     delete [] nu;
     delete [] nv;
}

void zzhang::CHungarianBP::TightenTriplet()
{
     int nclus_to_add = 5;
     std::vector<std::vector<int> > PossibleClusters;
     std::vector<std::set<int> > PossibleClustersSet;
     std::vector<std::vector<int> > pSubEclusters;
     std::vector<std::vector<boost::shared_array<int> > > pIdxConvertingTable;
     std::map<std::set<int>, bool> IsVisited;
     std::vector<std::pair<int, double> > pClusterWeight;
     std::vector<int> pClusterDim;
     //Find all possible extended clusters
     for (int i = 0; i < m_EClusterSet.size(); i++) {
	  if(m_EClusterSet[i].size() !=1) continue;
	  for (int supi = 0; supi < m_SupEClustersIdx[i].size(); supi++) {
	       int ci = m_SupEClustersIdx[i][supi];
	       
	       for (int supj = supi + 1; supj < m_SupEClustersIdx[i].size();
		    supj++) {
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

     //Find sub-clusters for possible extended clusters
     for (int i = 0; i < PossibleClustersSet.size(); i++) {
	  pIdxConvertingTable.push_back(std::vector<boost::shared_array<int> >());
	  pSubEclusters.push_back(std::vector<int>());
	  std::map<int, bool> IsAdded;
	  for (int ni = 0; ni < PossibleClusters[i].size(); ni++) {
	       int nidx = PossibleClusters[i][ni];

	       for (int k = 0; k < m_NodeToECluIdx[nidx].size(); k++) {
		    int ck = m_NodeToECluIdx[nidx][k];
		    if (IsAdded.find(ck) != IsAdded.end())
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
		    }
	       }
	  }
     }

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
	  pClusterWeight.push_back(
	       std::pair<int, double>(i, local_ub - local_exact));
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
	       printf(" %d", PossibleClusters[cidx][i]);
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

}
void zzhang::CHungarianBP::BPClear(int maxOuterIter, int maxInnerIter)
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
	       UpdateUV();
	       ComputeObj();
	       clock_t current_clk = clock();
#ifdef MATLAB_MEX_FILE
	       #if 0
	       mexPrintf(
		    "Iter=%d Dual=%16.10f Primal=%16.10f Diff=%16.10f IntGap=%16.10f Time=%16.8f\n",
		    iiter, m_CurrentUB, m_BestLB, lastdual - m_CurrentUB,
		    m_BestLB - m_CurrentUB,
		    1.0 * (current_clk - m_starttime) / CLOCKS_PER_SEC);
	       #endif
#else
	       printf(
		    "Iter=%d Dual=%16.10f Primal=%16.10f Diff=%16.10f IntGap=%16.10f Time=%16.8f\n",
		    iiter, m_CurrentUB, m_BestLB, lastdual - m_CurrentUB,
		    m_BestLB - m_CurrentUB,
		    1.0 * (current_clk - m_starttime) / CLOCKS_PER_SEC);

#endif	       
	       if(m_CurrentUB < m_LBBab)
		    break;
//	       if (std::fabs(lastdual - m_CurrentUB) < 1e-6) {
	       //		    break;
//	       }
	  }
	  if(m_CurrentUB < m_LBBab)
	       break;
	  if(std::fabs(m_CurrentUB - m_BestLB) < 1e-5)
	       break;
	  if(maxOuterIter == 1) break;
	  if(oiter == maxOuterIter) break;
	  clock_t tighten_start = clock();
	  TightenCluster(15);
	  clock_t tighten_end = clock();
//	  printf("Time Spend For tighten %12.5f\n", (1.0 * tighten_end - 1.0 * tighten_start) / CLOCKS_PER_SEC);
     }
}

zzhang::CHungarianBP::~CHungarianBP()
{
}

int lap(int dim, 
        cost **assigncost,
        col *rowsol, 
        row *colsol, 
        cost *u, 
        cost *v)

// input:
// dim        - problem size
// assigncost - cost matrix

// output:
// rowsol     - column assigned to row in solution
// colsol     - row assigned to column in solution
// u          - dual variables, row reduction numbers
// v          - dual variables, column reduction numbers

{
     bool unassignedfound;
     row  i, imin, numfree = 0, prvnumfree, f, i0, k, freerow, *pred, *free;
     col  j, j1, j2, endofpath, last, low, up, *collist, *matches;
     cost min, h, umin, usubmin, v2, *d;

     free = new row[dim];       // list of unassigned rows.
     collist = new col[dim];    // list of columns to be scanned in various ways.
     matches = new col[dim];    // counts how many times a row could be assigned.
     d = new cost[dim];         // 'cost-distance' in augmenting path calculation.
     pred = new row[dim];       // row-predecessor of column in augmenting/alternating path.

     // init how many times a row will be assigned in the column reduction.
     for (i = 0; i < dim; i++)  
	  matches[i] = 0;

     // COLUMN REDUCTION 
     for (j = dim-1; j >= 0; j--)    // reverse order gives better results.
     {
	  // find minimum cost over rows.
	  min = assigncost[0][j]; 
	  imin = 0;
	  for (i = 1; i < dim; i++)  
	       if (assigncost[i][j] < min) 
	       { 
		    min = assigncost[i][j]; 
		    imin = i;
	       }
	  v[j] = min; 

	  if (++matches[imin] == 1) 
	  { 
	       // init assignment if minimum row assigned for first time.
	       rowsol[imin] = j; 
	       colsol[j] = imin; 
	  }
	  else
	       colsol[j] = -1;        // row already assigned, column not assigned.
     }

     // REDUCTION TRANSFER
     for (i = 0; i < dim; i++) 
	  if (matches[i] == 0)     // fill list of unassigned 'free' rows.
	       free[numfree++] = i;
	  else
	       if (matches[i] == 1)   // transfer reduction from rows that are assigned once.
	       {
		    j1 = rowsol[i]; 
		    min = BIG;
		    for (j = 0; j < dim; j++)  
			 if (j != j1)
			      if (assigncost[i][j] - v[j] < min) 
				   min = assigncost[i][j] - v[j];
		    v[j1] = v[j1] - min;
	       }

     // AUGMENTING ROW REDUCTION 
     int loopcnt = 0;           // do-loop to be done twice.
     do
     {
	  loopcnt++;

	  // scan all free rows.
	  // in some cases, a free row may be replaced with another one to be scanned next.
	  k = 0; 
	  prvnumfree = numfree; 
	  numfree = 0;             // start list of rows still free after augmenting row reduction.
	  while (k < prvnumfree)
	  {
	       i = free[k]; 
	       k++;

	       // find minimum and second minimum reduced cost over columns.
	       umin = assigncost[i][0] - v[0]; 
	       j1 = 0; 
	       usubmin = BIG;
	       for (j = 1; j < dim; j++) 
	       {
		    h = assigncost[i][j] - v[j];
		    if (h < usubmin)
			 if (h >= umin) 
			 { 
			      usubmin = h; 
			      j2 = j;
			 }
			 else 
			 { 
			      usubmin = umin; 
			      umin = h; 
			      j2 = j1; 
			      j1 = j;
			 }
	       }

	       i0 = colsol[j1];
	       if (usubmin - umin > EPS ) 
		    // change the reduction of the minimum column to increase the minimum
		    // reduced cost in the row to the subminimum.
		    v[j1] = v[j1] - (usubmin - umin);
	       else                   // minimum and subminimum equal.
		    if (i0 >= 0)         // minimum column j1 is assigned.
		    { 
			 // swap columns j1 and j2, as j2 may be unassigned.
			 j1 = j2; 
			 i0 = colsol[j2];
		    }

	       // (re-)assign i to j1, possibly de-assigning an i0.
	       rowsol[i] = j1; 
	       colsol[j1] = i;

	       if (i0 >= 0)           // minimum column j1 assigned earlier.
		    if (usubmin - umin > EPS ) 
			 // put in current k, and go back to that k.
			 // continue augmenting path i - j1 with i0.
			 free[--k] = i0; 
		    else 
			 // no further augmenting reduction possible.
			 // store i0 in list of free rows for next phase.
			 free[numfree++] = i0; 
	  }
     }
     while (loopcnt < 2);       // repeat once.

     // AUGMENT SOLUTION for each free row.
     for (f = 0; f < numfree; f++) 
     {
	  freerow = free[f];       // start row of augmenting path.

	  // Dijkstra shortest path algorithm.
	  // runs until unassigned column added to shortest path tree.
	  for (j = 0; j < dim; j++)  
	  { 
	       d[j] = assigncost[freerow][j] - v[j]; 
	       pred[j] = freerow;
	       collist[j] = j;        // init column list.
	  }

	  low = 0; // columns in 0..low-1 are ready, now none.
	  up = 0;  // columns in low..up-1 are to be scanned for current minimum, now none.
	  // columns in up..dim-1 are to be considered later to find new minimum, 
	  // at this stage the list simply contains all columns 
	  unassignedfound = false;
	  do
	  {
	       if (up == low)         // no more columns to be scanned for current minimum.
	       {
		    last = low - 1; 

		    // scan columns for up..dim-1 to find all indices for which new minimum occurs.
		    // store these indices between low..up-1 (increasing up). 
		    min = d[collist[up++]]; 
		    for (k = up; k < dim; k++) 
		    {
			 j = collist[k]; 
			 h = d[j];
			 if (h <= min)
			 {
			      if (h < min)     // new minimum.
			      { 
				   up = low;      // restart list at index low.
				   min = h;
			      }
			      // new index with same minimum, put on undex up, and extend list.
			      collist[k] = collist[up]; 
			      collist[up++] = j; 
			 }
		    }

		    // check if any of the minimum columns happens to be unassigned.
		    // if so, we have an augmenting path right away.
		    for (k = low; k < up; k++) 
			 if (colsol[collist[k]] < 0) 
			 {
			      endofpath = collist[k];
			      unassignedfound = true;
			      break;
			 }
	       }

	       if (!unassignedfound) 
	       {
		    // update 'distances' between freerow and all unscanned columns, via next scanned column.
		    j1 = collist[low]; 
		    low++; 
		    i = colsol[j1]; 
		    h = assigncost[i][j1] - v[j1] - min;

		    for (k = up; k < dim; k++) 
		    {
			 j = collist[k]; 
			 v2 = assigncost[i][j] - v[j] - h;
			 if (v2 < d[j])
			 {
			      pred[j] = i;
			      if (v2 == min)   // new column found at same minimum value
				   if (colsol[j] < 0) 
				   {
					// if unassigned, shortest augmenting path is complete.
					endofpath = j;
					unassignedfound = true;
					break;
				   }
			      // else add to list to be scanned right away.
				   else 
				   { 
					collist[k] = collist[up]; 
					collist[up++] = j; 
				   }
			      d[j] = v2;
			 }
		    }
	       } 
	  }
	  while (!unassignedfound);

	  // update column prices.
	  for (k = 0; k <= last; k++)  
	  { 
	       j1 = collist[k]; 
	       v[j1] = v[j1] + d[j1] - min;
	  }

	  // reset row and column assignments along the alternating path.
	  do
	  {
	       i = pred[endofpath]; 
	       colsol[endofpath] = i; 
	       j1 = endofpath; 
	       endofpath = rowsol[i]; 
	       rowsol[i] = j1;
	  }
	  while (i != freerow);
     }

     // calculate optimal cost.
     cost lapcost = 0;
     for (i = 0; i < dim; i++)  
     {
	  j = rowsol[i];
	  u[i] = assigncost[i][j] - v[j];
	  lapcost = lapcost + assigncost[i][j]; 
     }

     // free reserved memory.
     delete[] pred;
     delete[] free;
     delete[] collist;
     delete[] matches;
     delete[] d;

     return lapcost;
}
/**
 *@param dim Dimensions of the problem.
 *@param AssignMent 
 */
static double SecondBestLap(int dim,
			    const boost::shared_array<int>& CurrentAssignMent,
			    boost::shared_array<int> &AssignMent,
			    cost **assigncost)
{
     col * rowsol = new col[dim];
     row * colsol = new row[dim];
     cost * u = new cost[dim];
     cost * v = new cost[dim];
     double Huge = 1e8;
     double min_cost = 1e20;
     std::vector<bool> IsVisited(dim, false);
     AssignMent = boost::shared_array<int> (new int[dim]);
     for(int i = 0; i < dim; i++)
     {
	  int CurrentLabel = CurrentAssignMent[i];
	  double TmpCost = assigncost[i][CurrentLabel];
	  assigncost[i][CurrentLabel] = Huge;
	  double lapcost = lap(dim, assigncost, rowsol, colsol, u, v);
	  assigncost[i][CurrentLabel] = TmpCost;
	  if(lapcost < min_cost)
	  {
	       min_cost = lapcost;
	       for(int j = 0; j < dim; j++)
	       {
		    AssignMent[j] = rowsol[j];
	       }
	  }
     }
     delete [] rowsol;
     delete [] colsol;
     delete [] u;
     delete [] v;
     return min_cost;
}
bool MBestLap(int dim, int N, std::vector< boost::shared_array<int> > & AssignMents, cost ** assigncost)
{
     col * rowsol = new col[dim];
     row * colsol = new row[dim];
     cost * u = new cost[dim];
     cost * v = new cost[dim];
     std::vector< boost::shared_array<int> > Xs(N);
     std::vector< boost::shared_array<int> > Ys(N);
     std::vector< double > Yvs(N);
     std::vector< double > Xvs(N);
     std::vector< std::pair<int, int> > Additional_Constraints;

     for(int i = 0; i < N; i++)
     {
	  if(i == 1)
	  {
	       boost::shared_array<int> BestDecoding(new int[dim]);
	       double lapcost = lap(dim, assigncost, BestDecoding.get(), colsol, u, v);
	       Ys[0] = BestDecoding;
	       Yvs[0] = lapcost;
	  }
	  else{
	       int lidx = -1;
	       
	  }
     }
	 return true;
}
