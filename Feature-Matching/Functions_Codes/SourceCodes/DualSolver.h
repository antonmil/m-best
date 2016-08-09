/* DualSolver.h --- 
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

#ifndef DUALSOLVER_H
#define DUALSOLVER_H 1
#include <set>
#include <map>
#include <vector>
#include <memory>
#include <string>
#include <list>
#include <ctime>
#include <boost/shared_array.hpp>
#include "common.h"

namespace zzhang {
/**
 *Base class for dual map solver
 */
     const double Huge = 10e50;
     typedef std::map<std::pair<int, int>, int> mapType;
     typedef std::vector<std::vector<std::pair<int, double> > > adj_type;
     const int RandomSize = 25;

     class DLLAPI CDualSolver {
     protected:
	  std::vector< std::vector<int> > ActiveConstraints;
	  /**
	   * Specifying that how many nodes a factor graph has.
	   */
	  int m_NodeSize;
	  /**
	   * Specifying that how many states each node takes.
	   */
	  boost::shared_array<int> m_NodeStates;
	  /**
	   *Extended Clusters stored in vector formulation.
	   */
	  std::vector<std::vector<int> > m_EClusterVec;
	  /**
	   *Extended Clusters stored in set formulation.
	   */
	  std::vector<std::set<int> > m_EClusterSet;
	  /**
	   *Sub-Extended-Clusters of a clusters
	   */
	  std::vector<std::vector<int> > m_SubEClusters;
	  /**
	   *Sup-Extended-Clusters of a clusters, Idx
	   */
	  std::vector<std::vector<int> > m_SupEClustersIdx;
	  /**
	   *Sup-Extended-Clusters of a clusters, Pos
	   */
	  std::vector<std::vector<int> > m_SupEClustersPos;
	  /**
	   *Beliefs
	   */
	  std::vector<boost::shared_array<double> > m_Beliefs;
	  std::vector< std::vector<int> > m_UsefulIdxs;
	  std::vector< std::vector<boost::shared_array<double> > > m_Phik;
	  std::vector< double > m_Gamma;
	  /**
	   *Potentials
	   */
	  std::vector<boost::shared_array<double> > m_Potentials;
	  /**
	   *Dimensions of each belief
	   */
	  std::vector<int> m_TotalDim;
	  /**
	   *Converting table for ids in belief (Used for partial max or partial sum)
	   */
	  std::vector<std::vector<boost::shared_array<int> > > m_IdxConvertingTable;
	  /**
	   *Local maximum for each extended cluster
	   */
	  std::vector<int> m_LocalMaximum;
	  /**
	   *The decoded solution
	   */
	  boost::shared_array<int> m_DecodedSolution;
	  /**
	   *Used to store which clusters a nodes belongs to.
	   */
	  std::vector<std::vector<int> > m_NodeToECluIdx;
	  /**
	   *Used to store which potentials an extended cluster corresponds to.
	   */
	  std::vector<int> m_EClusterPotIdx;
	  /**
	   * Used to tag if two edges are connected
	   */
	  std::map<int, int> m_IsConnected;
	  /**
	   *Mapping an extended cluster to its position in the extended cluster table.
	   */
	  std::map<std::set<int>, int> m_EClusterIdx;
	  /**
	   *Current upper bound for MAP Inference.
	   */
	  double m_CurrentUB;
	  /**
	   *Current lower bound for MAP Inference.
	   */
	  double m_BestLB;
	  /**
	   * Output file name.
	   */
	  std::string m_FName;
	  /**
	   *Lower bound for branch and bound
	   */
	  double m_LBBab;
	  /**
	   *the time we start to compute
	   */
	  clock_t m_starttime;
	  double m_primal2;
	  std::vector<double> IsFeasible;

     public:
	  /**
	   * @name  - CDualSolver
	   * Empty constructor.
	   */
	  CDualSolver() :
	       m_CurrentUB(Huge), m_BestLB(-Huge), m_LBBab(-Huge) {
	  }
	  virtual ~CDualSolver() {
	  }
	  /**
	   * Reading a factor graph from a UAI file.
	   * @name  - ReadUAIFile
	   * @return If succeed return true, otherwise false.
	   * 
	   */
	  bool ReadUAIFile(const std::string& FileName);
	  bool ConstructProblemWithConstraints(
	       int NofNodes,
	       int *NodeStates,
	       const std::vector< std::vector<int> >& Clusters,
	       const std::vector< boost::shared_array<double> >& Potentials,
	       const std::vector< std::vector< boost::shared_array< double > > >& Phiki,
	       double LBBab = -1e20
	       );
	  bool GetDecode(std::vector<int>& Decode)
	  {
	       Decode = std::vector<int>(m_DecodedSolution.get(),
					 m_DecodedSolution.get() + m_NodeSize);
	       return true;
	  }
	  void GetClustersAndPotentials(std::vector< std::vector<int> >& Clusters,
					std::vector< boost::shared_array<double> >& Potentials,
					std::vector<int>& TotalDim
	       ){
	       Clusters = m_EClusterVec;
	       Potentials = m_Beliefs;
	       TotalDim = m_TotalDim;
	  }
/**
 * @name  - ReadUAIBin
 * Read uai model stored in bianry files
 * @return 
 */
	  bool ReadUAIBin(const std::string& FileName);
     protected:
	  virtual bool Init() {
	       m_primal2 = -1e20;
	       return true;
	  }
	  ;
	  /**
	   * Generating index converting table between two clusters, where ecvec2 should be subset of ecvec1.
	   * @name  - IdxConverting
	   * @return If succeed return true, otherwise false.
	   */
	  bool IdxConverting(const std::vector<int>& ecvec1,
			     const std::vector<int>& ecvec2,
			     boost::shared_array<int>& CvtTable) const;
	  /**
	   *Convert an high dimension encoded idx of a cluster to a low dimension decoded one
	   * @name IdxConverting
	   * @return If succeed return true, otherwise false.
	   */
	  bool IdxConverting(int cluidx, int hidx, std::vector<int> & lidx) const;
	  /**
	   * Adding a cluster to the extended cluster table, where pfuncid is the idx of potentials, and \f$-1\f$ no potential corresponds to the cluster.
	   * @name  - AddECluster
	   * @return  If succeed return true, otherwise false.
	   * 
	   */
	  bool AddECluster(std::vector<int> cluster, int pfuncid);
	  virtual std::string GetStoreSuffix()= 0;
	  /**
	   * Compute the dual objective and decoded primal.
	   * @name  - ComputeObj
	   * 
	   */
	  virtual void ComputeObj();
	  
     };
}

#endif // DUALSOLVER_H

// Local Variables:
// mode: c++
// End:
