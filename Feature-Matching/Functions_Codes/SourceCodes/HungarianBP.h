/** HungarianBP.h --- 
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


#ifndef HBP_H
#define HBP_H 1

#include "PiSBPSolver.h"


//The below part is taken from LPJAV
/************************************************************************
*
*  lap.h
   version 1.0 - 21 june 1996
   author  Roy Jonker, MagicLogic Optimization Inc.
   
   header file for LAP
*
**************************************************************************/

/*************** CONSTANTS  *******************/

  #define BIG 100000000

/*************** TYPES      *******************/

  typedef int row;
  typedef int col;
  typedef int cost;

/*************** FUNCTIONS  *******************/

extern int lap(int dim, 
	       cost **assigncost,
	       col *rowsol, 
	       row *colsol, 
	       cost *u, 
	       cost *v);

bool MBestLap(int dim,
	      int N,
	      std::vector< boost::shared_array<int> > & AssignMents,
	      cost ** assigncost);
namespace zzhang{
     /**
      * Hungarian BP
      */
     class  DLLAPI CHungarianBP : public CPiSBPDualSolver{
     public:
	  /**
	   * Hungarian-BP Iterations
	   */
	  void BPClear(int maxOuterIter, int maxInnerIter);
     protected:
	  /**
	   * Node reparametrisations
	   */
	  std::vector< boost::shared_array<double> > m_bi;
	  /**
	   * Dual variables ui
	   */
	  boost::shared_array<double> m_ui;
	  /**
	   * Dual variables vi
	   */
	  boost::shared_array<double> m_vi;
	  /**
	   * Current Feasible Decode
	   */
	  boost::shared_array<int> m_CurrentFdecode;
	  /**
	   * Current sum of u and v
	   */
	  double m_sum_uv;
	  /**
	   * Triplets;
	   */
	  std::vector< std::vector<int> > Triplet;
	  /**
	   * Initialization
	   */
	  virtual bool Init();
	  /**
	   * The suffix for storing results;
	   */
	  virtual std::string GetStoreSuffix() {
	       return std::string(".HuBP");
	  }
	  /**
	   * Compute current objective value;
	   */
	  virtual void ComputeObj();
	  /**
	   * Updating variable u and v.
	   */
	  void UpdateUV();
	  /**
	   * Tighten the relaxation by adding triplet.
	   */
	  void TightenTriplet();
     public:
	  /**
	   * Get pointers point to U.
	   */
	  double *GetU(){
	       return m_ui.get();
	  }
	  /**
	   * Get pointers point to V.
	   */
	  double *GetV(){
	       return m_vi.get();
	  }
	  /**
	   * Set u and v.
	   */
	  bool setUV(double *u, double *v)
	  {
	       memcpy(m_ui.get(), u, sizeof(double) * m_NodeSize);
	       memcpy(m_vi.get(), v, sizeof(double) * m_NodeSize);
	  }
     public:
	  virtual ~CHungarianBP();
	  /**
	   * Get dual objective
	   */
	  void GetDual(double &dual){dual = m_CurrentUB;};
	 
     };

}



#endif // HBP_H
