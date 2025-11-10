// Cimrman headers
#include "tkrec/Solution.h"

// Bayeux:
#include <datatools/exception.h>

// ClassImp(tkrec::Solution);

namespace tkrec {

  Solution::Solution(const std::vector<PreclusterSolutionHdl> & _precluster_solutions)
  {
    precluster_solutions = _precluster_solutions;
  }
  
  std::vector<PreclusterSolutionHdl> & Solution::get_precluster_solutions()
  {
    return precluster_solutions;
  }  
  
  std::vector<ConstPreclusterSolutionHdl> Solution::get_precluster_solutions() const
  {
    std::vector<ConstPreclusterSolutionHdl> pr_solutions;
  	for(const auto & pr_solution : precluster_solutions)
  	{
  	  pr_solutions.push_back(pr_solution);
  	}
    return pr_solutions;
  }

  std::vector<ConstTrajectoryHdl> Solution::get_trajectories() const
  {
    std::vector<ConstTrajectoryHdl> all_trajectories;
    for(ConstPreclusterSolutionHdl precluster_solution : precluster_solutions)
    {
      std::vector<ConstTrajectoryHdl> trajectories = precluster_solution->get_trajectories();
      all_trajectories.insert(all_trajectories.end(), trajectories.begin(), trajectories.end());
    }
    return all_trajectories;
  }


  void Solution::print(std::ostream & out_) const
  {
    out_ <<"Solution: " << std::endl;
    out_ << "	" << precluster_solutions.size() << " precluster solutions: " << std::endl;
    for(auto i = 0u; i < precluster_solutions.size(); i++)
    {
      out_ << "	" << i << ". ";
      precluster_solutions[i]->print();
    }
  }

} //  end of namespace tkrec
