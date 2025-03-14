// TK headers
#include "tkrec/Event.h"

// Standard headers
#include <cmath>
#include <iomanip>

// - Bayeux:
#include <bayeux/datatools/clhep_units.h>

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
    tracker_hits.push_back(trhit);
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

  void Event::add_precluster(const std::vector<ConstTrackerHitHdl> & tracker_hits, bool is_prompt, int side)
  {
    preclusters.emplace_back(std::make_shared<Precluster>(tracker_hits, is_prompt, side));
  }


/*

  std::vector<TKtrhitHdl> TKEvent::filter_side(const std::vector<TKtrhitHdl>& _hits, int side)
  {
    vector<TKtrhitHdl> hits;	
    for(auto& hit : _hits)
      {
	if( side == hit->get_SRL('s'))
	  {
	    hits.push_back( hit );
	  }
      }
    return hits;
  }

  std::vector<TKtrhitHdl> TKEvent::filter_usable(const std::vector<TKtrhitHdl>& _hits)
  {
    vector<TKtrhitHdl> hits;	
    for(auto& hit : _hits)
      {
	// not using broken or too big (incorrectly associated) tracker hits
	if( hit->get_r() != -1.0 && hit->get_r() < 35.0 && hit->get_r() > 2.0 )
	  {
	    hits.push_back( hit );
	  }
      }
    return hits;
  }

  std::vector<TKtrhitHdl> TKEvent::filter_unassociated(const std::vector<TKtrhitHdl>& _hits)
  {
    vector<TKtrhitHdl> hits;	
    for(auto& hit : _hits)
      {
	if( not hit->has_associated_track() )
	  {
	    hits.push_back( hit );
	  }
      }
    return hits;
  }

  std::vector<TKtrhitHdl> TKEvent::filter_distant(const std::vector<TKtrhitHdl>& _hits)
  {
    double distance = 3.0; // in cells: 1 == 44mm
    vector<TKtrhitHdl> hits;
    for(auto i = 0u; i < _hits.size(); i++)
      {
	bool close = false;			
	int RL[2] = {_hits[i]->get_SRL('R'),_hits[i]->get_SRL('L')};
	for(auto j = 0u; j < _hits.size(); j++)
	  {
	    if(i == j) continue;
	    if(pow(RL[0] - _hits[i]->get_SRL('R'), 2) + pow(RL[1] - _hits[i]->get_SRL('L'), 2) <= distance*distance)
	      {
		close = true;
		// continue or break?	
	      }		
	  }
	if( close )
	  { 
	    hits.push_back(_hits[i]);
	  } 
      }
    return hits;
  }

  std::vector<TKtrhitHdl> TKEvent::filter_unclustered(const std::vector<TKtrhitHdl>& _hits, const TKEvent & event_)
  {
    vector<TKtrhitHdl> hits;
    for(auto& hit : _hits)
      {
	int SRL[3] = {hit->get_SRL('S'), hit->get_SRL('R'), hit->get_SRL('L')};
	bool clustered = false;			
	for(auto j = 0u; j < event_.clusters.size(); j++)
	  {	
	    for(auto k = 0u; k < event_.clusters[j]->get_tr_hits().size(); k++)
	      {	
		if(event_.clusters[j]->get_tr_hits()[k]->get_SRL('S') == SRL[0] &&
		   event_.clusters[j]->get_tr_hits()[k]->get_SRL('R') == SRL[1] &&
		   event_.clusters[j]->get_tr_hits()[k]->get_SRL('L') == SRL[2])
		  {
		    clustered = true;
		  }
	      }
	  }
	if( clustered  == false )
	  {
	    hits.push_back( hit );
	  }
      }
    return hits;
  }

  std::vector<TKtrhitHdl> TKEvent::filter_close_hits(const std::vector<TKtrhitHdl>& _hits,
						     double phi,
						     double r,
						     double distance_limit)
  {
    vector<TKtrhitHdl> hits;	
    for(auto& hit : _hits)
      {
	double R = hit->get_r();
	double x = hit->get_xy('x');
	double y = hit->get_xy('y');
	double distance = std::fabs(r - x*sin(phi) + y*cos(phi)) - R;
	if( std::fabs(distance) <= distance_limit )
	  {
	    hits.push_back( hit );
	  }
      }
    return hits;
  }

  void TKEvent::set_r(std::string drift_model, std::string association_mode)
  {	
    for(auto tr_hit = 0u; tr_hit < tr_hits.size(); tr_hit++)
      {
	double r = std::numeric_limits<double>::quiet_NaN();
	if(tr_hits[tr_hit]->get_tsp('0') != -1)
	  {
	    double min_time = std::numeric_limits<double>::infinity();
	    // associates tracker hits to OM with minimal time difference
	    if( association_mode == "time" ) 
	      {
		//double calo_hit = (calo_tdc * 6.25) - 400.0 + (400.0 * peak_cell / 1024.0);
		for(auto om_hit = 0u; om_hit < OM_hits.size(); om_hit++)
		  {
		    int64_t TDC_diff = 2*tr_hits[tr_hit]->get_tsp('0') - OM_hits[om_hit]->get_OM_TDC() + 44;
		    if(       TDC_diff  < 800 && 
			      TDC_diff  > 0   &&
			      double(TDC_diff) < min_time)
		      {
			min_time = 6.25 * TDC_diff;
			tr_hits[tr_hit]->set_associated_OMhit(OM_hits[om_hit]);
		      }
		  }
	      }

	    // associates tracker hits to OM with minimal distance
	    else if( association_mode == "distance" )
	      {
		double min_distance = std::numeric_limits<double>::infinity();
		for(auto om_hit = 0u; om_hit < OM_hits.size(); om_hit++)
		  {
		    double delta_x = tr_hits[tr_hit]->get_xy('x') - OM_hits[om_hit]->get_xyz('x');
		    double delta_y = tr_hits[tr_hit]->get_xy('y') - OM_hits[om_hit]->get_xyz('y');
		    double delta_z = tr_hits[tr_hit]->get_h()     - OM_hits[om_hit]->get_xyz('z');
		    double distance = sqrt( delta_x*delta_x + delta_y*delta_y + delta_z*delta_z );
		    if( distance < min_distance ) 
		      {
			int64_t TDC_diff = 2*tr_hits[tr_hit]->get_tsp('0') - OM_hits[om_hit]->get_OM_TDC() + 44;
			if( TDC_diff <= 0 ) 
			  {
			    continue;
			  }
			min_time = 6.25 * TDC_diff;
			tr_hits[tr_hit]->set_associated_OMhit(OM_hits[om_hit]);
			min_distance = distance;
		      }
					
		  }
	      }
	    else clog << "invalid association model: choose \"distance\" or \"time\"." << endl;

	    if( drift_model == "Manchester" )
	      {
		const double A1 = 0.570947153108633;
		const double B1 = 0.580148313540993;
		const double C1 = 1.6567483468611;
		const double A2 = 1.86938462695651;
		const double B2 = 0.949912427483918;

		const double t_usec = min_time / 1000.0;
		const double ut = 10. * t_usec;
		
		r = A1 * ut / (std::pow(ut, B1) + C1);
		if (r > A2 * ut / (std::pow(ut, B2))) 
		  {
		    r = A2 * ut / (std::pow(ut, B2));
		  }	

		r *= 10.0;	
	      }
	    // modification of Betsy's model
	    else if( drift_model == "Betsy" )
	      {
		const double a1 = 0.828;
		const double b1 = -0.907;
		const double a2 = 0.402;
		const double b2 = -1.955;
		
		min_time = min_time / 1000.0;
		if(min_time < 3.0845 && min_time > 0.0)
		  {
		    r = pow(min_time/a1, 1.0/(1.0-b1)) * 10.0;
		  }
		else
		  {
		    r = pow(min_time/a2, 1.0/(1.0-b2)) * 10.0;
		  }
	      }
	    else cerr << "invalid drift time model: choose \"Betsy\" or \"Manchester\"." << endl;
	  }
	else
	  {
	    r = -1.0;
	  }
	tr_hits[tr_hit]->set_r(r);
	tr_hits[tr_hit]->set_sigma_R();
      }
  }
*/

} // end of namespace tkrec  
