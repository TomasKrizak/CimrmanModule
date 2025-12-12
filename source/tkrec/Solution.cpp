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

  void Solution::set_num_of_trajectories(unsigned int no_trajectories)
  {
    num_of_trajectories = no_trajectories;
  }
  
  void Solution::set_num_of_segments(unsigned int no_segments)
  {
    num_of_segments = no_segments;
  }
  
  void Solution::set_num_of_unclustered_hit(unsigned int no_unclustered_hit)
  {
    num_of_unclustered_hit = no_unclustered_hit;
  }
  
  void Solution::set_total_summed_chi2(double chi2_sum)
  {
    total_summed_chi2 = chi2_sum;
  }
  
  unsigned int Solution::get_num_of_trajectories() const
  {
    return num_of_trajectories;
  }
  
  unsigned int Solution::get_num_of_segments() const
  {
    return num_of_segments;
  }
  
  unsigned int Solution::get_num_of_unclustered_hit() const
  {
    return num_of_unclustered_hit;
  }
  
  double Solution::get_total_summed_chi2() const
  {
    return total_summed_chi2;
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
