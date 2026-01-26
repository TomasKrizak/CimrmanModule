#ifndef FALAISE_CIMRMAN_EVENT_H
#define FALAISE_CIMRMAN_EVENT_H

// Standard headers
#include <vector>
#include <iostream>

// Bayeux headers
#include <bayeux/datatools/logger.h>
#include <bayeux/datatools/properties.h>

namespace tkrec {

  // Forward declarations
  class OMHit;
  using OMHitHdl = std::shared_ptr<OMHit>;
  using ConstOMHitHdl = std::shared_ptr<const OMHit>;
  
  class TrackerHit;
  using TrackerHitHdl = std::shared_ptr<TrackerHit>;
  using ConstTrackerHitHdl = std::shared_ptr<const TrackerHit>;
  
  class Cluster;
  using ClusterHdl = std::shared_ptr<Cluster>;
  using ConstClusterHdl = std::shared_ptr<const Cluster>;
  
  class Precluster;
  using PreclusterHdl = std::shared_ptr<Precluster>;
  using ConstPreclusterHdl = std::shared_ptr<const Precluster>;
  
  class Solution;
  using SolutionHdl = std::shared_ptr<Solution>;
  using ConstSolutionHdl = std::shared_ptr<const Solution>;
  
  class Visu;
  

  // Model of a tracker reconstructed event
  class Event
  {
  private:
	
    int run_number   = -1; ///< Run number
    int event_number = -1; ///< Event number

    std::vector<OMHitHdl> OM_hits; ///< List of calo hits (only for visualization)	
    std::vector<TrackerHitHdl> tracker_hits; ///< List of tracker hits	
    std::vector<PreclusterHdl> preclusters; ///< List of preclusters		
    std::vector<SolutionHdl> solutions; ///< List of reconstruction solutions
    
    // things excluded by Cimrman, stored to pass them along to other modules
    std::vector<TrackerHitHdl> invalid_tracker_hits; // prompt tracker hits without drift radius - unusable
    std::vector<ClusterHdl> unfitted_clusters; ///< List of unfitted clusters (not completely implemented - not needed)	
		
  public:
	
    // basic functionality section:
    Event() = default;
    Event(int _run_number, int _event_number);
    virtual ~Event() = default;

    void set_event_ids(int _run_number, int _event_number);
    /// Reset the event internal data
    void reset();

    std::vector<OMHitHdl> & get_OM_hits(); 		
    std::vector<ConstOMHitHdl> get_OM_hits() const; 		
    
    std::vector<TrackerHitHdl> & get_tracker_hits(); 		
    std::vector<ConstTrackerHitHdl> get_tracker_hits() const; 		

    std::vector<PreclusterHdl> & get_preclusters();		
    std::vector<ConstPreclusterHdl> get_preclusters() const;	
    
    std::vector<SolutionHdl> & get_solutions();		
    std::vector<ConstSolutionHdl> get_solutions() const;		

    std::vector<TrackerHitHdl> & get_invalid_tracker_hits();
    std::vector<ConstTrackerHitHdl> get_invalid_tracker_hits() const;

    std::vector<ClusterHdl> & get_unfitted_clusters();		
    std::vector<ConstClusterHdl> get_unfitted_clusters() const;	
    
    int get_run_number() const;
    int get_event_number() const;

    void add_OM_hit(const OMHitHdl & OMhit);
    void add_tracker_hit(const TrackerHitHdl & trhit);

    void add_precluster(const std::vector<TrackerHitHdl> & tracker_hits, bool is_prompt, int side);

    void print(std::ostream & out_ = std::clog) const;

    friend class Visu; ///< Private access from visualization engine

  };

} //  end of namespace tkrec

#endif // FALAISE_CIMRMAN_EVENT_H
