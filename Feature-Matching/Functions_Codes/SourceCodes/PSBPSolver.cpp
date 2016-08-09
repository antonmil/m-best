/*
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

#include "PSBPSolver.h"
#include <cstdlib>
#include <random>
#include <algorithm>
bool zzhang::CPSBPDualSolver::Init() {
	srand(456);
	clock_t init_start = clock();

	std::vector<std::vector<int> > TNodeToEclusterIdx(m_NodeSize);
	std::map<std::pair<int, int>, bool> IsVisited;
	for (int i = 0; i < m_EClusterVec.size(); i++) {
		for (int j = 0; j < m_EClusterVec[i].size(); j++) {
			TNodeToEclusterIdx[m_EClusterVec[i][j]].push_back(i);
		}
	}

	for (int i = 0; i < m_EClusterVec.size(); i++) {
		std::set<int> ncluster;
		for (int ni = 0; ni < m_EClusterVec[i].size(); ni++) {
			ncluster.clear();
			for (int j = 0; j < m_EClusterVec[i].size(); j++) {
				if (j != ni)
					ncluster.insert(m_EClusterVec[i][j]);
			}
			std::map<std::set<int>, int>::iterator it;
			it = m_EClusterIdx.find(ncluster);
			int sidx;
			if (it == m_EClusterIdx.end()) {
				std::vector<int> ncluvec(ncluster.begin(), ncluster.end());
				AddECluster(ncluvec, -1);
				sidx = m_EClusterVec.size() - 1;
				for (int ncluit = 0; ncluit < ncluvec.size(); ncluit++) {
					TNodeToEclusterIdx[ncluvec[ncluit]].push_back(sidx);
				}
			} else {
				sidx = it->second;
			}
			m_SubEClusters[i].push_back(sidx);
		}
	};
	clock_t time_explore_graph_end = clock();
	printf("Explore Graph Structure Spends: %f seconds\n",
	       (time_explore_graph_end - init_start) * 1.0 / CLOCKS_PER_SEC);
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
