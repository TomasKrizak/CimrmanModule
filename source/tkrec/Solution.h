#ifndef FALAISE_CIMRMAN_SOLUTION_H
#define FALAISE_CIMRMAN_SOLUTION_H

// Standard headers
#include <iostream>
#include <memory>
#include <limits>
#include <vector>

// Bayeux headers
#include <datatools/logger.h>

namespace tkrec {

  // Forward declaration 
  class PreclusterSolution;
  using PreclusterSolutionHdl = std::shared_ptr<PreclusterSolution>;
  using ConstPreclusterSolutionHdl = std::shared_ptr<const PreclusterSolution>;
  
  class Trajectory;
  using TrajectoryHdl = std::shared_ptr<Trajectory>;
  using ConstTrajectoryHdl = std::shared_ptr<const Trajectory>;



  class Solution
  {
  private:
		
    std::vector<PreclusterSolutionHdl> precluster_solutions;
    unsigned int num_of_trajectories = 0;
    unsigned int num_of_segments = 0;
    unsigned int num_of_unclustered_hit = 0;
    double total_summed_chi2 = std::numeric_limits<double>::max();
    
  public:
    
    datatools::logger::priority verbosity = datatools::logger::PRIO_FATAL;
    
  public:
    
    Solution() = default;	
    Solution(const std::vector<PreclusterSolutionHdl> & precluster_solutions);	
    virtual ~Solution() = default;  
   
    std::vector<PreclusterSolutionHdl> & get_precluster_solutions();
    std::vector<ConstPreclusterSolutionHdl> get_precluster_solutions() const;
    
    std::vector<ConstTrajectoryHdl> get_trajectories() const;

    void set_num_of_trajectories(unsigned int no_trajectories);
    void set_num_of_segments(unsigned int no_segments);
    void set_num_of_unclustered_hit(unsigned int no_unclustered_hit);
    void set_total_summed_chi2(double chi2_sum);
    
    unsigned int get_num_of_trajectories() const;
    unsigned int get_num_of_segments() const;
    unsigned int get_num_of_unclustered_hit() const;
    double get_total_summed_chi2() const;

    void print(std::ostream & out_ = std::cout) const;
  };

  typedef std::shared_ptr<Solution> SolutionHdl;
  typedef std::shared_ptr<const Solution> ConstSolutionHdl;

} //  end of namespace tkrec

#endif // FALAISE_CIMRMAN_SOLUTION_H
