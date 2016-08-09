/** DualSolver.cpp --- 
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

#include <iostream>
#include <fstream>
#include <cstring>
#include <algorithm>
#include "DualSolver.h"
bool zzhang::CDualSolver::ReadUAIBin( const std::string& FileName)
{
     clock_t start = clock();
     FILE *fp = fopen(FileName.c_str(), "rb");
     if (fp == NULL) {
	  perror("Fatal error");
	  return false;
     }
     fread(&m_NodeSize, sizeof(int), 1, fp);
     if (m_NodeSize <= 0) {
	  fprintf(stderr, "Fatal error: the amount of nodes must be positive\n");
	  return false;
     }
     m_NodeStates = boost::shared_array<int>(new int[m_NodeSize]);
     fread(m_NodeStates.get(), sizeof(int), m_NodeSize, fp);
     int ClusterCnt;
     fread( &ClusterCnt, sizeof(int), 1, fp);
     int *clu_buffer = new int[1024];
     std::vector<std::vector<int> > ClustersToBeAdded(ClusterCnt);
     std::vector<int> CluDims(ClusterCnt);
     for (int i = 0; i < ClusterCnt; i++) {
	  int ClusterSize;
	  fread( &ClusterSize, sizeof(int), 1, fp);
	  fread(clu_buffer, sizeof(int), ClusterSize, fp);
	  ClustersToBeAdded[i] = std::vector<int>(clu_buffer, clu_buffer + ClusterSize);
	  CluDims[i] = 1;
	  for(int j = 0; j < ClusterSize; j++)
	       CluDims[i] *= m_NodeStates.get()[clu_buffer[j]];
     }
     m_Potentials = std::vector< boost::shared_array<double> >(ClusterCnt);
     for(int i = 0; i < ClusterCnt; i++)
     {
	  m_Potentials[i] = boost::shared_array<double>(new double[CluDims[i]]);
	  fread(m_Potentials[i].get(), sizeof(double), CluDims[i], fp);
     }
     fclose(fp);
     delete clu_buffer;
     for (int i = 0; i < ClusterCnt; i++) {
	  this->AddECluster(ClustersToBeAdded[i], i);
     }
     clock_t end = clock();
     printf("Time spend: %f \n", (end - start) * 1.0 / CLOCKS_PER_SEC);
     //fflush(stdout);
     return Init();
	
}
bool zzhang::CDualSolver::ReadUAIFile(const std::string& FileName) {
     clock_t start = clock();
     std::ifstream fin(FileName);
     if (!fin.good()) {
	  return false;
     }
     m_FName = FileName;
     const static std::string LOG_SUFFIX = ".LG";
     bool isLog = true;
     for (int i = 0; i < LOG_SUFFIX.size(); i++) {
	  if (FileName[FileName.size() - 3 + i] != LOG_SUFFIX[i]) {
	       isLog = false;
	       break;
	  }
     }
     const static std::string BIN_SUFFIX = ".BIN";
     bool isBin = true;
     for (int i = 0; i < BIN_SUFFIX.size(); i++) {
	  if (FileName[FileName.size() - 4 + i] != BIN_SUFFIX[i]) {
	       isBin = false;
	       break;
	  }
     }
     if(isBin)
	  return ReadUAIBin(FileName);
     FILE *fp = fopen(FileName.c_str(), "rb");
     if (fp == NULL) {
	  perror("Fatal error");
	  return false;
     }
     static char header_buffer[1024] = { '\0' };
     fscanf(fp, "%s", header_buffer);
     if (strcmp(header_buffer, "MARKOV") != 0) {
	  fprintf(stderr, "Fatal error: file format error\n");
	  return false;
     }
     fscanf(fp, "%d", &m_NodeSize);
     if (m_NodeSize <= 0) {
	  fprintf(stderr, "Fatal error: the amount of nodes must be positive\n");
	  return false;
     }
     m_NodeStates = boost::shared_array<int>(new int[m_NodeSize]);
     for (int *p = m_NodeStates.get(); p < m_NodeStates.get() + m_NodeSize;
	  p++) {
	  fscanf(fp, "%d", p);
     }
     int ClusterCnt;
     fscanf(fp, "%d", &ClusterCnt);
     std::vector<std::vector<int> > ClustersToBeAdded(ClusterCnt);
     for (int i = 0; i < ClusterCnt; i++) {
	  int ClusterSize;
	  fscanf(fp, "%d", &ClusterSize);
	  for (int j = 0; j < ClusterSize; j++) {
	       int cnode;
	       fscanf(fp, "%d", &cnode);
	       ClustersToBeAdded[i].push_back(cnode);
	  }
     }
     for (int i = 0; i < ClusterCnt; i++) {
	  int dim = 0;
	  fscanf(fp, "%d", &dim);
	  if (dim == -1) {
	       std::cerr << "Not Implemented" << std::endl;
	       return false;
	  }
	  m_Potentials.push_back(
	       boost::shared_array<double>(new double[dim]));
	  for (int j = 0; j < dim; j++) {
	       fscanf(fp, "%lf", (m_Potentials[i].get() + j));
	       if (!isLog) {
		    if (m_Potentials[i].get()[j] > 1e-6)
			 m_Potentials[i].get()[j] = log(m_Potentials[i].get()[j]);
		    else {
			 m_Potentials[i].get()[j] = -Huge;
		    }
	       }
	  }
     }
     for (int i = 0; i < ClusterCnt; i++) {
	  this->AddECluster(ClustersToBeAdded[i], i);
     }
     clock_t end = clock();
     printf("Time spend Readfile: %f \n", (end - start) * 1.0 / CLOCKS_PER_SEC);
     fflush(stdout);
     start = clock();
     bool res = Init();
     end = clock();
     printf("Time spend Init: %f \n", (end - start) * 1.0 / CLOCKS_PER_SEC);
     return res;
}

bool zzhang::CDualSolver::ConstructProblemWithConstraints(
     int NofNodes,
     int *NodeStates,
     const std::vector< std::vector<int> >& Clusters,
     const std::vector< boost::shared_array<double> >& Potentials,
     const std::vector< std::vector< boost::shared_array< double > > >& Phiki,
     double LBBab
     )
{
     m_LBBab = LBBab;
     m_NodeSize = NofNodes;
     m_NodeStates = boost::shared_array<int>(new int[m_NodeSize],
					 std::default_delete<int[]>());
     memcpy(m_NodeStates.get(), NodeStates, sizeof(int) * m_NodeSize);
     m_Potentials = Potentials;
     m_Phik = Phiki;
     int NofConstraints = m_Phik.size();
     //printf("Number of Constraints: %d\n", NofConstraints);
     m_Gamma = std::vector<double>(NofConstraints, 0.0);
     for(int i = 0; i < Clusters.size(); i++)
     {
	  this->AddECluster(Clusters[i], i);
     }
     m_starttime = clock();
     bool res = Init();
     ActiveConstraints = std::vector< std::vector<int> >( m_Phik.size() );
     for(int k = 0; k < m_Phik.size(); k++)
     {
	  for(int ci = 0; ci < m_Phik[k].size(); ci++)
	  {
	       if(m_Phik[k][ci].get() != NULL)
		    ActiveConstraints[k].push_back(ci);
	  }
     }
     return false;
}

bool zzhang::CDualSolver::IdxConverting(int cluidx, int hidx,
					std::vector<int> & lidx) const {
     if (hidx < 0 || hidx > m_TotalDim[cluidx]) {
	  return false;
     }
     int t = hidx;
     lidx = std::vector<int>(m_EClusterVec[cluidx].size());
     for (int j = m_EClusterVec[cluidx].size() - 1; j >= 0; j--) {
	  lidx[j] = t % m_NodeStates.get()[m_EClusterVec[cluidx][j]];
	  t /= m_NodeStates.get()[m_EClusterVec[cluidx][j]];
     }
     return true;
}
/**
 * Generating index converting table between two clusters, where ecvec2 should be subset of ecvec1.
 * @name  - IdxConverting
 * @return If succeed return true, otherwise false.
 */
bool zzhang::CDualSolver::IdxConverting(const std::vector<int>& ecvec1,
					const std::vector<int>& ecvec2, boost::shared_array<int>& CvtTable) const {
#ifndef WIN32
     int States[ecvec1.size()];
#else
						boost::shared_array<int> States(new int[ecvec1.size()]);
#endif
     std::vector<int> node_pos(m_NodeSize, -1);
#ifndef WIN32
     int N2InN1[ecvec2.size()];
#else
	 boost::shared_array<int> N2InN1(new int[ecvec2.size()]);
#endif
     for (int i = 0; i < ecvec2.size(); i++) {
	  node_pos[ecvec2[i]] = i;
     }
     for (int i = 0; i < ecvec1.size(); i++)
     {
	  int nidx = ecvec1[i];
	  if(node_pos[nidx] != -1)
	  {
	       N2InN1[node_pos[nidx]] = i;
	  }
     }
     int totaldim = 1;
     
     for (int i = 0; i < ecvec1.size(); i++) {
	  States[i] = m_NodeStates.get ()[ecvec1[i]];
	  totaldim *= m_NodeStates.get()[ecvec1[i]];
     }
     CvtTable = boost::shared_array<int>(new int[totaldim],
				     std::default_delete<int[]>());
     int *nstates = m_NodeStates.get();
     int *cvttable = CvtTable.get();
#ifndef WIN32
     int AssignMent[ecvec1.size()];
#else
	 boost::shared_array<int> AssignMent(new int[ecvec1.size()]);
#endif
     for(int i = 0; i < ecvec1.size(); i++)
     {
	  AssignMent[i] = 0;
     }
     int mn = ecvec1.size() - 1;
     for (int i = 0; i < totaldim; i++) {
	  int c = 0;
	  int ni = ecvec1.size() - 1;
	  int t = i;
	  int cvtidx = 0;

	  for (int j = 0; j < ecvec2.size(); j++) {
	       cvtidx = cvtidx * nstates[ecvec2[j]] + AssignMent[N2InN1[j]];
	  }
	  cvttable[i] = cvtidx;
	  switch(ni)
	  {
	  case 0:
	       AssignMent[0] ++;
	       break;
	  case 1:
	       AssignMent[1] ++;
	       if(AssignMent[1] == States[1])
	       {
		    AssignMent[1] = 0;
		    AssignMent[0] ++;
	       }
	       break;
	  case 2:
	       AssignMent[2] ++;
	       if(AssignMent[2] == States[2])
	       {
		    AssignMent[2] = 0;
		    AssignMent[1] ++;
		    if(AssignMent[1] == States[1])
		    {
			 AssignMent[1] = 0;
			 AssignMent[0] ++;
		    }
	       }
	       break;
	  case 3:
	       AssignMent[3] ++;
	       if(AssignMent[3] == States[3])
	       {
		    AssignMent[3] = 0;
		    AssignMent[2] ++;
		    if(AssignMent[2] == States[2])
		    {
			 AssignMent[2] = 0;
			 AssignMent[1] ++;
			 if(AssignMent[1] == States[1])
			 {
			      AssignMent[1] = 0;
			      AssignMent[0] ++;
			 }
		    }
	       }
	       break;
	  default:
	       for(int j = ni; j >=0; j--)
	       {
		    AssignMent[j] ++;
		    if(AssignMent[j] < States[j])
		    {
			 break;
		    }
		    AssignMent[j] = 0;
	       }
	       break;
	  }
     }
     return true;
}

bool zzhang::CDualSolver::AddECluster(std::vector<int> cluster, int pfuncid) {
     std::set<int> clu_set;
     for (int i = 0; i < cluster.size(); i++) {
	  clu_set.insert(cluster[i]);
     }
     std::map<std::set<int>, int>::iterator it = m_EClusterIdx.find(clu_set);
     if (it == m_EClusterIdx.end()) {
	  int totalDim = 1;
	  for (int i = 0; i < cluster.size(); i++) {
	       totalDim *= m_NodeStates.get()[cluster[i]];
	  }
	  this->m_EClusterSet.push_back(clu_set);
	  this->m_EClusterVec.push_back(cluster);
	  while (m_EClusterPotIdx.size() < m_EClusterSet.size()) {
	       m_EClusterPotIdx.push_back(-1);
	  }
	  this->m_TotalDim.push_back(totalDim);
	  this->m_EClusterPotIdx[m_EClusterSet.size() - 1] = pfuncid;
	  this->m_EClusterIdx[clu_set] = m_EClusterSet.size() - 1;
	  this->m_SubEClusters.push_back(std::vector<int>());
	  this->m_SupEClustersIdx.push_back(std::vector<int>());
	  this->m_SupEClustersPos.push_back(std::vector<int>());
     } else {
	  if (pfuncid == -1)
	       return true;
	  if (this->m_EClusterPotIdx[it->second] == -1) {
	       this->m_EClusterPotIdx[it->second] = pfuncid;
	       return true;
	  }
	  boost::shared_array<int> cvtTable;
	  this->IdxConverting(this->m_EClusterVec[it->second], cluster, cvtTable);
	  for (int i = 0; i < m_TotalDim[it->second]; i++) {
	       this->m_Potentials[this->m_EClusterPotIdx[it->second]].get()[i] +=
		    this->m_Potentials[pfuncid].get()[cvtTable.get()[i]];
	  }
     }
     return true;
}

void zzhang::CDualSolver::ComputeObj() {
     static bool isfirst = true;
     double dual = 0.0;
     boost::shared_array<int> CurrentDecode(new int[m_NodeSize]);
     std::for_each(CurrentDecode.get(), CurrentDecode.get() + m_NodeSize,
		   [](int& p) {p = -1;});
     int NofConstraints = m_Phik.size();
     IsFeasible =      std::vector<double>(NofConstraints, 0.0);
     for (int i = 0; i < m_Beliefs.size(); i++) {
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
	  if (m_EClusterSet[i].size() == 1) {
	       CurrentDecode.get()[m_EClusterVec[i][0]] = m_LocalMaximum[i];
	  }
	  dual += m_Beliefs[i].get()[m_LocalMaximum[i]];
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
     m_CurrentUB = dual;
     double primal = 0;
     double primal2 = 0;
 
     for (int i = 0; i < m_EClusterVec.size(); i++) {
	  int idx = 0;
	  for (int j = 0; j < m_EClusterVec[i].size(); j++) {
	       if (j == 0)
		    idx += CurrentDecode.get()[m_EClusterVec[i][j]];
	       else
		    idx = idx * m_NodeStates.get()[m_EClusterVec[i][j]]
			 + CurrentDecode.get()[m_EClusterVec[i][j]];
	  }
	  primal += m_Beliefs[i].get()[idx];
	  primal2 += m_Beliefs[i].get()[idx];
	  for(int k = 0; k < NofConstraints; k++)
	  {
	       if(m_Phik[k].size() > i)
	       {
		    if(m_Phik[k][i].get() != NULL)
		    {
			 IsFeasible[k] += m_Phik[k][i].get()[idx];
			 primal += m_Gamma[k] * m_Phik[k][i].get()[idx];
			 
		    }

	       }
	  }
     }
     for(int k = 0; k < NofConstraints; k++)
     {
	  if(IsFeasible[k] >= 1e-6)
	  {
	       primal = -1e20;
	       break;
	  }
     }
     if (primal - m_BestLB > 1e-10) {
	  m_BestLB = primal;
	  m_DecodedSolution = CurrentDecode;
     }
     if( primal2 > m_primal2)
	  m_primal2 = primal2;
     m_DecodedSolution = CurrentDecode;

     return;
}
