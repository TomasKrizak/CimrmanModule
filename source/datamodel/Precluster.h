#ifndef FALAISE_CIMRMAN_PRECLUSTER_H
#define FALAISE_CIMRMAN_PRECLUSTER_H

// Standard headers
#include <iostream>
#include <memory>
#include <vector>

// Bayeux headers
#include "datatools/logger.h"

namespace cimrman::datamodel {

  // forward declarations
  class TrackerHit;
  using TrackerHitHdl = std::shared_ptr<TrackerHit>;
  using ConstTrackerHitHdl = std::shared_ptr<const TrackerHit>;
  
  class Cluster;
  using ClusterHdl = std::shared_ptr<Cluster>;
  using ConstClusterHdl = std::shared_ptr<const Cluster>;
  
  class PreclusterSolution;
  using PreclusterSolutionHdl = std::shared_ptr<PreclusterSolution>;
  using ConstPreclusterSolutionHdl = std::shared_ptr<const PreclusterSolution>;


  class Precluster
  {
  private:
		
    bool prompt = true;
    int side = -1;
    
    std::vector<TrackerHitHdl> unclustered_tracker_hits;
    std::vector<ClusterHdl> linear_clusters;
    std::vector<PreclusterSolutionHdl> precluster_solutions;
    
  public:
    
    datatools::logger::priority verbosity = datatools::logger::PRIO_FATAL;
    
  public:
    
    Precluster() = default;	
    Precluster(const std::vector<TrackerHitHdl> & tracker_hits, bool _prompt, int _side);	
    virtual ~Precluster() = default;    
    bool is_prompt() const;
    bool is_delayed() const;
    int get_side() const;
    void set_side(int _side);

    std::vector<TrackerHitHdl> & get_unclustered_tracker_hits();
    std::vector<ConstTrackerHitHdl> get_unclustered_const_tracker_hits() const;
    
    std::vector<ClusterHdl> & get_clusters();
    std::vector<ConstClusterHdl> get_clusters() const;
    
    std::vector<PreclusterSolutionHdl> & get_precluster_solutions();
    std::vector<ConstPreclusterSolutionHdl> get_precluster_solutions() const;

    void print(std::ostream & out_ = std::cout) const;
  };

  typedef std::shared_ptr<Precluster> PreclusterHdl;
  typedef std::shared_ptr<const Precluster> ConstPreclusterHdl;

} //  end of namespace cimrman::datamodel

#endif // FALAISE_CIMRMAN_PRECLUSTER_H
