#ifndef FALAISE_CIMRMAN_SOLUTION_H
#define FALAISE_CIMRMAN_SOLUTION_H

// Standard headers
#include <iostream>
#include <memory>

#include "tkrec/PreclusterSolution.h"

#include <datatools/logger.h>

namespace tkrec {

  class Solution
  {
  private:
		
    std::vector<PreclusterSolutionHdl> precluster_solutions;
    
  public:
    
    datatools::logger::priority verbosity = datatools::logger::PRIO_FATAL;
    
  public:
    
    Solution() = default;	
    Solution(const std::vector<PreclusterSolutionHdl> & precluster_solutions);	
    virtual ~Solution() = default;    

    std::vector<PreclusterSolutionHdl> & get_precluster_solutions();
    std::vector<ConstPreclusterSolutionHdl> get_precluster_solutions() const;
    
    std::vector<ConstTrajectoryHdl> get_trajectories() const;

    void print(std::ostream & out_ = std::cout) const;
  };

  typedef std::shared_ptr<Solution> SolutionHdl;
  typedef std::shared_ptr<const Solution> ConstSolutionHdl;

} //  end of namespace tkrec

#endif // FALAISE_CIMRMAN_SOLUTION_H
