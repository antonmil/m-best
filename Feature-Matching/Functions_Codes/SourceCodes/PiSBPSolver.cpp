/** PiSBPSolver.cpp--
 *
 * Copyright (C) 2016 Zhen Zhang
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


#include "PiSBPSolver.h"
#include <cstdlib>
#include <random>
#include <iostream>
#include <algorithm>
#include <iterator>
#ifdef MATLAB_MEX_FILE
#include "mex.h"
#endif
bool zzhang::CPiSBPDualSolver::Init() {

     srand(456);
     clock_t init_start = clock();

     std::vector<std::vector<int> > TNodeToEclusterIdx(m_NodeSize);
     std::map<std::pair<int, int>, bool> IsVisited;
     for (int i = 0; i < m_EClusterVec.size(); i++) {
	  for (int j = 0; j < m_EClusterVec[i].size(); j++) {
	       TNodeToEclusterIdx[m_EClusterVec[i][j]].push_back(i);
	  }
     }
     for(int i = 0; i < m_EClusterVec.size(); i++)
     {
	  if(m_EClusterVec[i].size() != 2)
	       continue;
	  int ii = m_EClusterVec[i][0]; int jj = m_EClusterVec[i][1];
	  std::set<int> clu1; clu1.insert(ii);
	  std::set<int> clu2; clu2.insert(jj);
	  if(m_EClusterIdx.find(clu1) == m_EClusterIdx.end())
	  {
	       std::vector<int> clu1vec(1); clu1vec[0] = ii;
	       AddECluster(clu1vec, -1);
	  }
	  if(m_EClusterIdx.find(clu2) == m_EClusterIdx.end())
	  {
	       std::vector<int> clu2vec(1); clu2vec[0] = ii;
	       AddECluster(clu2vec, -1);
	  }
     }

     int OriginalSize;
     do {
	  OriginalSize = m_EClusterVec.size();
	  for(int i = 0; i < m_NodeSize; i++)
	  {
	       for (int j = 0; j < TNodeToEclusterIdx[i].size(); j++) {
		    int cj = TNodeToEclusterIdx[i][j];
		    if(m_EClusterVec[cj].size() <= 2)
			 continue;
		    for(int k =j + 1; k < TNodeToEclusterIdx[i].size(); k++)
		    {
			 int ck = TNodeToEclusterIdx[i][k];
			 if(m_EClusterVec[cj].size() <= 2)
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
			 for (int j = 0; j < intvec.size(); j++) {
			      TNodeToEclusterIdx[intvec[j]].push_back(
				   m_EClusterVec.size() - 1);

			 }
		    }
	       }
	  }

     } while (OriginalSize != m_EClusterVec.size());

     IsVisited.clear();
     for(int i = 0; i < m_EClusterVec.size(); i++)
     {
	  if(m_EClusterVec[i].size() == 2)
	  {
	       int ii = m_EClusterVec[i][0]; int jj = m_EClusterVec[i][1];
	       std::set<int> clu1; clu1.insert(ii);
	       std::set<int> clu2; clu2.insert(jj);
	       int idx1 = m_EClusterIdx[clu1];
	       int idx2 = m_EClusterIdx[clu2];
	       m_SubEClusters[i].push_back(idx1);
	       m_SubEClusters[i].push_back(idx2);
	  }
     }
     
     for (int i = 0; i < m_EClusterSet.size(); i++) {
	  if(m_EClusterVec[i].size() <= 2)
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
     
     std::vector<std::vector<int> > RCSubEClusters;
     for (int i = 0; i < m_EClusterSet.size(); i++) {
	  std::map<int, bool> NeedNotToBeadded;

	  RCSubEClusters.push_back(std::vector<int>());
	  for (int si = 0; si < m_SubEClusters[i].size(); si++) {
	       int sidx = m_SubEClusters[i][si];
	       for (int ssi = 0; ssi < m_SubEClusters[sidx].size(); ssi++) {
		    int ssidx = m_SubEClusters[sidx][ssi];
		    NeedNotToBeadded[ssidx] = true;
	       }
	  }
	  for (int si = 0; si < m_SubEClusters[i].size(); si++) {
	       int sidx = m_SubEClusters[i][si];
	       if (NeedNotToBeadded.find(sidx) == NeedNotToBeadded.end()) {
		    RCSubEClusters[i].push_back(sidx);
	       }
	  }
     }

     m_SubEClusters = RCSubEClusters;
     clock_t time_explore_graph_end = clock();
     // printf("Explore Graph Structure Spends: %f seconds\n",
//	    (time_explore_graph_end - init_start) * 1.0 / CLOCKS_PER_SEC);

     CacheIndex();
     
     m_Potentials.clear();
     m_LocalMaximum = std::vector<int>(m_EClusterVec.size());
     m_DecodedSolution = boost::shared_array<int>(new int[m_NodeSize]);
     m_NodeToECluIdx = TNodeToEclusterIdx;
     clock_t cache_idx_cvt = clock();
     //   printf("Cache index converting table spends: %f seconds\n",
//	    (cache_idx_cvt - time_explore_graph_end) * 1.0 / CLOCKS_PER_SEC);

     
#ifdef OUTPUT_ORDER
     for (int i = 0; i < m_Beliefs.size(); i++) {
	  if (m_EClusterVec[i].size() < 4)
	       continue;
	  //std::cout << "Here" << std::endl;
	  int dim = m_TotalDim[i];
	  std::vector<std::pair<double, int> > PotentialWithID(dim);
	  for (int j = 0; j < dim; j++) {
	       PotentialWithID[j] = std::pair<double, int>(m_Beliefs[i].get()[j],
							   j);
	  }
	  std::sort(PotentialWithID.begin(), PotentialWithID.end(),
		    GreaterPair<double>);
	  m_Order[i] = std::vector<int>(dim);
	  m_Error[i] = std::vector<double>(dim);
	  int z = 11;
	  for (int j = 0; j < dim; j++) {
	       m_Order[i][j] = PotentialWithID[j].second;
	       m_Error[i][j] = PotentialWithID[j].first - PotentialWithID[z].first;
	       if (m_Error[i][j] > 0)
		    m_Error[i][j] = 0;
	  }
     }
#endif
     
     return true;

}


bool zzhang::CPiSBPDualSolverDropHOP::BeliefPropagation(int maxIter)
{
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
	       int start = 0; //rand();
	       tighten++;

	       for (int ii = m_EClusterVec.size() - 1; ii >= 0; ii--) {
		    int i = ii;
		    if (m_SubEClusters[i].size() == 0)
			 continue;
		    if (tighten > 50 && m_EClusterVec[i].size() > 2) {
			 continue;
		    }
		    if (veryhigh > 10 && m_EClusterVec[i].size() > 4) {
			 //continue;
		    }
		    double tcinverce = 1.0 / m_SubEClusters[i].size();
		    int dim = m_TotalDim[i];
		    double *bigdat = m_Beliefs[i].get();
		    for (int si = 0; si < m_SubEClusters[i].size(); si++) {
			 for (int vi = 0; vi < dim; vi++) {
			      int sidx = m_SubEClusters[i][si];
			      double *smalldata = m_Beliefs[sidx].get();
			      bigdat[vi] += *(smalldata
					      + m_IdxConvertingTable[i][si].get()[vi]);
			 }

		    }
		    bigdat = m_Beliefs[i].get();
		    double LocalMax = -Huge;
		    for (int vi = 0; vi < dim; vi++) {
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
			 for (int vi = 0; vi < dim; vi++) {
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
			 for (int vi = 0; vi < dim; vi++) {
			      bc[vi] -= *(smalldata
					  + m_IdxConvertingTable[i][si].get()[vi]);
			 }
		    }
	       }
	       ComputeObj();
	       veryhigh++;
	       if (tighten == 200)
		    tighten = 0;

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
		    "Iter=%d UB=%16.10f LB=%16.10f Diff=%16.10f IntGap=%16.10f Time=%16.8f\n",
		    it, m_CurrentUB, m_BestLB, lastdual - m_CurrentUB,
		    m_BestLB - m_CurrentUB,
		    1.0 * (current_clk - m_starttime) / CLOCKS_PER_SEC);
	       fflush(stdout);
	       if (std::fabs(lastdual - m_CurrentUB) < 1e-8) {
		    if (tighten > 50) {
			 tighten = 0;
			 continue;
		    }
		    tighten = 0;
		    veryhigh = 0;
		    first = it + 1;
		    break;
	       }
	  }

	  int ncluster = 20;
	  double bound;
	  static bool isfirst = true;
	  if (isfirst) {
	       isfirst = false;
	       if (abs(m_BestLB - m_CurrentUB) < 1)
		    maxIter = 600;
	       else
		    maxIter = 100;
	  }
	  if (abs(m_BestLB - m_CurrentUB) < 1)
	       maxIter = 600;
	  TightenCluster(ncluster);
	  veryhigh = 0;
	  // TightenCycle(ncluster, bound, 1);
	  // if(bound < CLUSTER_THR)
	  // {
	  //      TightenCycle(ncluster, bound, 2);
	  // }
     }
     return true;
}
