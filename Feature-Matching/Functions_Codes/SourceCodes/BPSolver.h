/** BPSolver.h --- 
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

#ifndef BPSolver_h
#define BPSolver_h 1
#include <set>
#include <map>
#include <vector>
#include <memory>
#include <string>
#include <list>
#include <ctime>
#include "DualSolver.h"

//#define OUTPUT_ORDER

#define Inf 9999999999.9
namespace zzhang {
#ifdef OUTPUT_ORDER
template<typename T>
bool GreaterPair(const std::pair<T, int> &a, const std::pair<T, int> &b) {
	return a.first > b.first;
}
#endif
/**
 *Approximately solve MAP inference problem using belief propagation dual solver.
 */
class DLLAPI CBPDualSolver: public CDualSolver {
#ifdef OUTPUT_ORDER
protected:
	std::map<int, std::vector<int> > m_Order;
	std::map<int, std::vector<double> > m_Error;
#endif
private:
	double m_SmoothDual2;
	/**
	 * Pseudo distribution used in smoothing BP.
	 */
	std::vector< boost::shared_array<long double> > m_PseudoDistribution;
	double m_offset;
	/**
	 * The parameter espilon used in smoothing BP.
	 */
	double m_epsilon;
	void SmoothComputeObj();
	void SmoothBPOneIter();
	bool SmoothBPInit(double epsilon=1);
public:
	bool SmoothBP(int maxIter, double epsilon = 0.1);
public:
	/**
	 * Using belief propagation to solve the MAP inference problem, where maxIter is the number of max iterations.
	 * @name  - BeliefPropagation
	 * @return If succeed return true, otherwise false.
	 *
	 */
	bool BeliefPropagation(int maxIter);
	bool ParallelBP(int maxIter);
	void BPOneIter();
	void BPClear(int maxOuterIter, int maxInnerIter);
	void ConstrainedBP(int maxOuterIter, int maxInnerIter);
	void SGBP(int maxIter);
	void DualCoordinateDescend(int MaxIter);
	void DCDUpdateGamma(int i);
	std::vector<double> &GetGamma(){
	     return m_Gamma;
	}
	
protected:
	int TightenCluster(int nclus_to_add);
	void PBPOneIter();
	void GenPartionGraph();
	std::vector< std::vector<int> > m_PartionGraph;
	/**
	 * Initialisation of the solver.
	 * @name  - Init
	 * @return If succeed return true, otherwise false.
	 * 
	 */
	virtual bool Init();
	bool CacheIndex();
	virtual std::string GetStoreSuffix() {
		return std::string(".MI");
	}

};
}
#endif // BPSolver_h
// Local Variables:
// mode: c++
// End:
