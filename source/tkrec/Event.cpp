// Cimrman headers
#include "tkrec/Event.h"
#include "tkrec/OMHit.h"
#include "tkrec/TrackerHit.h"
#include "tkrec/Cluster.h"
#include "tkrec/Precluster.h"
#include "tkrec/Solution.h"

// Standard headers
#include <cmath>
#include <iomanip>

// Bayeux:
#include "bayeux/datatools/clhep_units.h"

// ClassImp(tkrec::Event);

namespace tkrec {
    
  Event::Event(int _run_number ,int _event_number)
  {
    set_event_ids(_event_number, _run_number);
    return;
  }

  void Event::reset()
  {
    event_number = -1;
    run_number = -1;

    OM_hits.clear();
    tracker_hits.clear();
    invalid_tracker_hits.clear();
    preclusters.clear();
    solutions.clear();
  }

  void Event::set_event_ids(int _run_number, int _event_number)
  {
    event_number = _event_number;
	  run_number = _run_number;
  }
 
  std::vector<OMHitHdl> & Event::get_OM_hits()
  {
    return OM_hits;
  }

  std::vector<ConstOMHitHdl> Event::get_OM_hits() const
  {
    std::vector<ConstOMHitHdl> hits;
    for (const auto & hit : OM_hits)
    {
      hits.push_back(hit);
    }
    return hits;
  }

  std::vector<TrackerHitHdl> & Event::get_tracker_hits()
  {
    return tracker_hits;
  }

  std::vector<ConstTrackerHitHdl> Event::get_tracker_hits() const
  {
    std::vector<ConstTrackerHitHdl> hits;
    for (const auto& h : tracker_hits)
    {
      hits.push_back(h);
    }
    return hits;
  }

    
  std::vector<TrackerHitHdl> & Event::get_invalid_tracker_hits()
  {
    return invalid_tracker_hits;
  }

  std::vector<ConstTrackerHitHdl> Event::get_invalid_tracker_hits() const
  {
    std::vector<ConstTrackerHitHdl> hits;
    for (const auto& h : invalid_tracker_hits)
    {
      hits.push_back(h);
    }
    return hits;
  }


  std::vector<ClusterHdl> & Event::get_unfitted_clusters()
  {
    return unfitted_clusters;
  }

  std::vector<ConstClusterHdl> Event::get_unfitted_clusters() const
  {
    std::vector<ConstClusterHdl> clusters;
    for (const auto & cl : unfitted_clusters)
    {
      clusters.push_back(cl);
    }
    return clusters;
  }


  std::vector<PreclusterHdl> & Event::get_preclusters()
  {
    return preclusters;
  }
  
  std::vector<ConstPreclusterHdl> Event::get_preclusters() const
  {
    std::vector<ConstPreclusterHdl> precl;
    for (const auto & p : preclusters)
    {
      precl.push_back(p);
    }
    return precl;
  }


  std::vector<SolutionHdl> & Event::get_solutions()
  {
    return solutions;
  }

  std::vector<ConstSolutionHdl> Event::get_solutions() const
  {
    std::vector<ConstSolutionHdl> sol;
    for (const auto & s : solutions)
    {
      sol.push_back(s);
    }
    return sol;
  }

  int Event::get_run_number() const
  {
    return run_number;
  }

  int Event::get_event_number() const
  {	
    return event_number;
  }

  void Event::add_OM_hit(const OMHitHdl & OMhit)
  {
    OM_hits.push_back(OMhit);
  }	
  
  void Event::add_tracker_hit(const TrackerHitHdl & trhit)
  {
    if( !trhit->has_valid_R() && trhit->is_prompt())
    {
      invalid_tracker_hits.push_back(trhit);
    }
    else
    {
      tracker_hits.push_back(trhit);    
    }
  }	

  void Event::print(std::ostream & out_) const
  {
    out_ << std::endl;
    out_ << "RUN " << run_number << " | EVENT " << event_number << std::endl << std::endl;
    out_ << "Collection of OM hits: " << std::endl;

    for(const auto & OM_hit : OM_hits) 
    {
	  OM_hit->print(out_);
    }
	
    out_ << std::endl;
    out_ << "Collection of tracker hits: " << std::endl;
	
    for(const auto & tr_hit : tracker_hits) 
    {
	  tr_hit->print(out_);
    }

    out_ << std::endl;
    out_ << "Collection of preclusters: " << std::endl;
	
    for(const auto & precluster : preclusters)
    {	
      precluster->print(out_);
    }
    out_ << std::endl;
    out_ << "Collection of solutions: " << std::endl;
	
    for(const auto & solution : solutions)
    {	
      solution->print(out_);
    }
    out_ << std::endl;
  }

  void Event::add_precluster(const std::vector<TrackerHitHdl> & tracker_hits, bool is_prompt, int side)
  {
    preclusters.emplace_back(std::make_shared<Precluster>(tracker_hits, is_prompt, side));
  }

} // end of namespace tkrec  
