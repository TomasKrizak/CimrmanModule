// Interface from Falaise
#include "tkrec/Algos.h"

// Standard headers
#include <iomanip>
#include <algorithm>
#include <vector>
#include <array>

// Boost:
#include <boost/multi_array.hpp>

// Bayeux:
#include <bayeux/datatools/exception.h>
#include <bayeux/datatools/clhep_units.h>

// Root:
#include <TH2F.h>
#include <TCanvas.h>


namespace tkrec {
 
  EventRecMode event_recmode_from_label(const std::string & label_)
  {
    if (label_ == "electron_kinked") return EventRecMode::electron_kinked;
    if (label_ == "electron_straight") return EventRecMode::electron_straight;
    return EventRecMode::undefined;
  }
 
  void TKEventRecConfig::parse(const datatools::properties & config_)
  {
    if (config_.has_key("verbosity")) {
      std::string verbosityLabel = config_.fetch_string("verbosity");
      auto parsedVerb = datatools::logger::get_priority(verbosityLabel);
      DT_THROW_IF(parsedVerb == datatools::logger::PRIO_UNDEFINED,
                  std::logic_error,
                  "Undefined verbosity label " << std::quoted(verbosityLabel));
      this->verbosity = parsedVerb;
    }

    if (config_.has_key("mode")) {
      std::string evRecModeLabel = config_.fetch_string("mode");
      auto evRecMode = event_recmode_from_label(evRecModeLabel);
      DT_THROW_IF(evRecMode == EventRecMode::undefined,
                  std::logic_error,
                  "Undefined event reconstruction mode " << std::quoted(evRecModeLabel));
      this->mode = evRecMode;
    }
    
    if (config_.has_flag("visualization")) {
      this->visualization = true;
    }
  
    if (config_.has_flag("save_sinograms")) {
      this->save_sinograms = true;
    }
  
    if (config_.has_flag("force_default_sigma_r")) {
      this->force_default_sigma_r = true;
    }
  
    if (config_.has_key("default_sigma_r")) {
      auto value = config_.fetch_real_with_explicit_dimension("default_sigma_r",
							        "length");
      DT_THROW_IF(value <= 1.0e-3 * CLHEP::mm or value >= 50. * CLHEP::mm,
                  std::logic_error,
                  "Invalid default sigma for drift radius"); 
      this->default_sigma_r = value / CLHEP::mm;
    }
  
    if (config_.has_key("clustering_max_distance")) {
      this->clustering_max_distance =
        config_.fetch_real_with_explicit_dimension("clustering_max_distance", 
                                                   "length");
      DT_THROW_IF(this->clustering_max_distance < 44.0 * CLHEP::mm,
		  std::logic_error,
		  "Invalid clustering_max_distance value");
    }

    if (config_.has_key("clustering_hit_association_distance")) {
      this->clustering_hit_association_distance =
        config_.fetch_real_with_explicit_dimension("clustering_hit_association_distance",
                                                   "length");
    }
    
    if (config_.has_key("clustering_no_iterations")) {
      this->clustering_no_iterations =
        config_.fetch_integer_scalar("clustering_no_iterations"); //??
    }
    
    if (config_.has_key("clustering_resolution_phi")) {
      this->clustering_resolution_phi =
        config_.fetch_integer_scalar("clustering_resolution_phi"); //??
    }
    
    if (config_.has_key("clustering_resolution_r")) {
      this->clustering_resolution_r =
        config_.fetch_integer_scalar("clustering_resolution_r"); // fetch_integer_scalar??
    }

    if (config_.has_key("clustering_max_initial_precision_r")) {
      this->clustering_max_initial_precision_r =
        config_.fetch_real_with_explicit_dimension("clustering_max_initial_precision_r",
                                                   "length");
    }
    
    if (config_.has_key("clustering_zoom_factor")) {
      this->clustering_zoom_factor =
        config_.fetch_dimensionless_real("clustering_zoom_factor");
    }
    
    if (config_.has_key("clustering_uncertainty")) {
      this->clustering_uncertainty =
        config_.fetch_real_with_explicit_dimension("clustering_uncertainty",
                                                   "length");
    }
    
    if (config_.has_key("chi_square_threshold")) {
      this->chi_square_threshold =
        config_.fetch_dimensionless_real("chi_square_threshold");
      DT_THROW_IF(this->chi_square_threshold < 1.0,
		  std::logic_error,
		  "Invalid chi_square_threshold value");
    }
    
    if (config_.has_key("polylines_max_vertical_distance")) {
      this->polylines_max_vertical_distance =
        config_.fetch_real_with_explicit_dimension("polylines_max_vertical_distance",
                                                   "length");
      DT_THROW_IF(this->polylines_max_vertical_distance < 0.0 * CLHEP::mm,
		  std::logic_error,
		  "Invalid polylines_max_vertical_distance value");
    }
    
    if (config_.has_key("polylines_max_extention_distance")) {
      this->polylines_max_extention_distance =
        config_.fetch_real_with_explicit_dimension("polylines_max_extention_distance",
                                                   "length");
      DT_THROW_IF(this->polylines_max_extention_distance < 0.0,
		  std::logic_error,
		  "Invalid polylines_max_extention_distance value");
    }
    
    if (config_.has_key("polylines_min_tracker_hits_distance")) {
      this->polylines_min_tracker_hits_distance =
        config_.fetch_real_with_explicit_dimension("polylines_min_tracker_hits_distance",
                                                   "length");
      DT_THROW_IF(this->polylines_min_tracker_hits_distance < 0.0,
		  std::logic_error,
		  "Invalid polylines_min_tracker_hits_distance value");
    }
    
    if (config_.has_key("polylines_max_kink_angle")) {
      this->polylines_max_kink_angle =
        config_.fetch_dimensionless_real("polylines_max_kink_angle");
      DT_THROW_IF(this->polylines_max_kink_angle < 0.0 || this->polylines_max_kink_angle > 180.0,
		  std::logic_error,
		  "Invalid polylines_max_kink_angle value");
    }
    
    if (config_.has_key("polylines_max_trajectories_middlepoint_distance")) {
      this->polylines_max_trajectories_middlepoint_distance =
        config_.fetch_real_with_explicit_dimension("polylines_max_trajectories_middlepoint_distance",
                                                   "length");
      DT_THROW_IF(this->polylines_max_trajectories_middlepoint_distance < 0.0,
		  std::logic_error,
		  "Invalid polylines_max_trajectories_middlepoint_distance value");
    }
    
    if (config_.has_key("polylines_max_trajectory_endpoints_distance")) {
      this->polylines_max_trajectory_endpoints_distance =
        config_.fetch_real_with_explicit_dimension("polylines_max_trajectory_endpoints_distance",
                                                   "length");
      DT_THROW_IF(this->polylines_max_trajectory_endpoints_distance < 0.0,
		  std::logic_error,
		  "Invalid polylines_max_trajectory_endpoints_distance value");
    }
    
    if (config_.has_key("polylines_max_trajectory_connection_angle")) {
      this->polylines_max_trajectory_connection_angle =
        config_.fetch_dimensionless_real("polylines_max_trajectory_connection_angle");
      DT_THROW_IF(this->polylines_max_trajectory_connection_angle < 0.0 || this->polylines_max_trajectory_connection_angle > 180.0,
		  std::logic_error,
		  "Invalid polylines_max_trajectory_connection_angle value");
    }
    
    if (config_.has_key("polylines_min_distance_from_foil")) {
      this->polylines_min_distance_from_foil =
        config_.fetch_real_with_explicit_dimension("polylines_min_distance_from_foil",
                                                   "length");
      DT_THROW_IF(this->polylines_min_distance_from_foil < 0.0,
		  std::logic_error,
		  "Invalid polylines_min_distance_from_foil value");   
    }
    
    return;
  }
 
  Algos::Algos(const Geometry & geom_)
    : _geom_(geom_) 
  {
    return;
  }
 
  Algos::~Algos()
  {
    if (is_initialized()) {
      reset();
    }
    return;
  }

  void Algos::set_event(Event & event_)
  {
    _event_ = &event_;
    if (_visu_) {
      _visu_->set_event(*_event_);
    }
    return;
  }

  bool Algos::has_event() const
  {
    return _event_ != nullptr;
  }

  bool Algos::is_initialized() const
  {
    return _config_.mode != EventRecMode::undefined;
  }
  
  void Algos::initialize(const TKEventRecConfig & evrecconf_)
  {
    _config_ = evrecconf_;
    DT_THROW_IF(_config_.mode == EventRecMode::undefined,
                std::logic_error,
                "Undefined reconstruction mode");

    if (_config_.visualization) {
      _visu_ = std::make_unique<Visu>(_geom_);
    }
    
    return;
  }
  
  void Algos::reset()
  {
    _event_ = nullptr;
    _config_ = TKEventRecConfig();
    return;
  }
  
  void Algos::process(Event & event_)
  {
    set_event(event_);
    
    // common preclustering to identify identify prompt and delayed preclusters and distant groups of tracker hits
    precluster();
     
    // electron reconstruction (creates prompt precluster solutions)
    if (_config_.mode == EventRecMode::electron_kinked) {
      DT_LOG_DEBUG(_config_.verbosity, "Processing 'electron polyline reconstruction' mode...");
      _process_electron_kinked_();
    } else if (_config_.mode == EventRecMode::electron_straight) {
      DT_LOG_DEBUG(_config_.verbosity, "Processing 'electron line reconstruction' mode...");
      _process_electron_straight_();
    }
	
	  // alpha reconstruction (creates delayed precluster solutions)
    if (_config_.reconstruct_alphas) {
      DT_LOG_DEBUG(_config_.verbosity, "Processing 'alpha reconstruction' mode...");
      _process_alpha_();
    } 
  
    // combines prompt and delayed precluster solutions into common solutions
    create_solutions();

    if (_visu_) {
      _visu_->make_top_projection();
      _visu_->build_event();
    }
    return;
  }

  void Algos::_process_electron_kinked_()
  {
    Legendre_transform_cluster_finder();
    make_MLM_fits();
    combine_into_precluster_solutions();
    create_polyline_trajectories();
    refine_trajectories();
    return;
  }
  
  void Algos::_process_electron_straight_()
  {
    Legendre_transform_cluster_finder();
    make_MLM_fits();
    combine_into_precluster_solutions();
    create_line_trajectories();
    //refine_trajectories();
    return;
  }

  void Algos::_process_alpha_()
  {
    // ??
    return;
  }
  
  // step 1: preclustering
  void Algos::precluster()
  {
    // basic preclustering into delayed/prompt hits and into halves of tracker
    std::vector<ConstTrackerHitHdl> prompt_hits_side0;
    std::vector<ConstTrackerHitHdl> prompt_hits_side1;
    std::vector<ConstTrackerHitHdl> delayed_hits_side0;
    std::vector<ConstTrackerHitHdl> delayed_hits_side1;
    
    for(auto & hit : _event_->get_tracker_hits())
    {
      if(hit->is_prompt())
      {
        if(!hit->has_valid_R())
        {
         continue;        
        }
        
        switch(hit->get_SRL()[0])
        {
        case 0:
          prompt_hits_side0.push_back(hit);
          break;
        case 1:
          prompt_hits_side1.push_back(hit);
          break;
        }
      }
      else 
      {
        switch(hit->get_SRL()[0])
        {
        case 0:
          delayed_hits_side0.push_back(hit);
          break;
        case 1:
          delayed_hits_side1.push_back(hit);
          break;
        }
      }
    }
    
    // additional spatial clustering 
    std::vector<std::vector<ConstTrackerHitHdl>> temp_preclusters = separate_hits(prompt_hits_side0);
    for(auto & prec : temp_preclusters)
    {
      _event_->add_precluster(prec, true, 0);
    }
    temp_preclusters = separate_hits(prompt_hits_side1);
    for(auto & prec : temp_preclusters)
    {
       _event_->add_precluster(prec, true, 1);
    }
    temp_preclusters = separate_hits(delayed_hits_side0);    
    for(auto & prec : temp_preclusters)
    {
      _event_->add_precluster(prec, false, 0);
    }
    temp_preclusters = separate_hits(delayed_hits_side1);
    for(auto & prec : temp_preclusters)
    {
       _event_->add_precluster(prec, false, 1);
    }
    
    return;
  }
  
  
  // step 2: clustering
  void Algos::Legendre_transform_cluster_finder()
  {
    for(auto & precluster : _event_->get_preclusters())
    {
      // clusterizes separatelly every prompt cluster
      if( precluster->is_prompt() )
      {        
        std::vector<ConstTrackerHitHdl> & unclustered_hits = precluster->get_unclustered_tracker_hits(); 
        std::vector<ClusterHdl> & clusters = precluster->get_clusters();

        // recursively applies Legendre transform and spatial clustering to find linear clusters
        // creates clusters and removes its tracker hits from "unclustered_hits"
        clusterize(unclustered_hits, clusters);
      }
    }
    return;
  }
  
  /*
  // step 2: clustering
  void Algos::Legendre_transform_cluster_finder()
  {
    const double association_distance = _config_.clustering_hit_association_distance;
    for(auto & precluster : _event_->get_preclusters())
    {
      // clusterizes separatelly every prompt cluster
      if( precluster->is_prompt() )
      {
        std::vector<ConstTrackerHitHdl> & hits = precluster->get_unclustered_tracker_hits(); 
        while( hits.size() > 2u )
        {
          // Legendre transform on hits to find phi and r candidates
          double phi_estimate, r_estimate;
          find_cluster_Legendre( hits, phi_estimate, r_estimate );
          
          // identifying which tracker hits belong to phi, r line candidate
          std::vector<ConstTrackerHitHdl> cluster_hits;
          separate_close_hits_to_line(hits, cluster_hits, phi_estimate, r_estimate, association_distance);
          
          // additional spatial separation to filter coincidental associations to the cluster          
          std::vector<std::vector<ConstTrackerHitHdl>> sub_clusters = separate_hits(cluster_hits);
          
          // finding the biggest group
          int largest = 0;
          auto size = 0u;
          for(auto i = 0u; i < sub_clusters.size(); i++)
          {
            if(sub_clusters[i].size() > size)
            {
              size = sub_clusters[i].size();
              largest = i; 
            }
          }

          // creating a cluster from the biggest groups if it has at least 3 hits
          // putting the rest back into unclustered hits  
          bool cluster_found = false;         
          for(auto j = 0u; j < sub_clusters.size(); j++)
          {
            if(j == largest && sub_clusters[j].size() > 2u)
            {
              // create and add new cluster to the precluster
              ClusterHdl cluster = std::make_shared<Cluster>( sub_clusters[j], phi_estimate, r_estimate, false );
              precluster->get_clusters().push_back( cluster );
              cluster_found = true;
            }
            else
            {
              // putting the separated hits back into unclustered hits
              hits.insert(hits.end(), sub_clusters[j].begin(), sub_clusters[j].end());
            }
          }
          
          // no meaningful candidate found (repeating the proccess would not give a different result)
          if( not cluster_found ) 
          {
            break;
          }
        }
      }
    }
    return;
  }*/
  
  
  void Algos::clusterize(std::vector<ConstTrackerHitHdl> & tracker_hits, std::vector<ClusterHdl> & clusters)
  {
    if( tracker_hits.size() < 3u ) return;
    
    const double association_distance = _config_.clustering_hit_association_distance;
    
    // Legendre transform on hits to find phi and r candidates
    double phi_estimate, r_estimate;
    find_cluster_Legendre( tracker_hits, phi_estimate, r_estimate );
    
    // identifying which tracker hits belong to phi, r line candidate
    std::vector<ConstTrackerHitHdl> cluster_hits;
    separate_close_hits_to_line(tracker_hits, cluster_hits, phi_estimate, r_estimate, association_distance);
    
    // additional spatial separation to filter coincidental associations to the cluster          
    std::vector<std::vector<ConstTrackerHitHdl>> sub_clusters = separate_hits(cluster_hits);
    
    // finding the biggest group
    int largest = 0;
    auto size = 0u;
    for(auto i = 0u; i < sub_clusters.size(); i++)
    {
      if(sub_clusters[i].size() > size)
      {
        size = sub_clusters[i].size();
        largest = i; 
      }
    }
    
    // creating a cluster from the biggest groups if it has at least 3 hits
    // putting the rest back into unclustered hits  
    bool cluster_found = false; 
    for(auto j = 0u; j < sub_clusters.size(); j++)
    {
      if(j == largest && sub_clusters[j].size() > 2u)
      {
        // create and add new cluster to the precluster
        ClusterHdl cluster = std::make_shared<Cluster>( sub_clusters[j], phi_estimate, r_estimate, false );
        clusters.push_back( cluster );
        cluster_found = true;
      }
      else
      {
        // putting the separated hits back into unclustered hits
        tracker_hits.insert(tracker_hits.end(), sub_clusters[j].begin(), sub_clusters[j].end());
      }
    }
    
    // meaningful candidate found (repeating the proccess on the rest of hits)
    if( cluster_found ) 
    {
      std::vector<std::vector<ConstTrackerHitHdl>> sub_groups = separate_hits(tracker_hits); 
      tracker_hits.clear();
      for(auto & sub_group : sub_groups)
      {
        clusterize(sub_group, clusters);
        tracker_hits.insert(tracker_hits.end(), sub_group.begin(), sub_group.end());
      }
    }
  }
  
  
  void Algos::separate_close_hits_to_line(std::vector<ConstTrackerHitHdl> & hits,
                                          std::vector<ConstTrackerHitHdl> & hits_separated,
                                          const double phi, const double r, const double distance_threshold)
  {
    auto it = hits.begin();
    while(it != hits.end())
    {
      auto & hit = *it;
      double R = hit->get_R();
      double x = hit->get_x();
      double y = hit->get_y();
      double distance = std::abs(r - x * std::sin(phi) + y * std::cos(phi)) - R;
      
      if( std::abs(distance) <= distance_threshold )
      {
        hits_separated.push_back( std::move(hit) );
        it = hits.erase(it);
      }
      else
      {
        it++;
      }
    }
    return;
  }
  
  void Algos::find_cluster_Legendre(const std::vector<ConstTrackerHitHdl> & hits, double & phi_estimate, double & r_estimate) const
  {
    bool save_sinograms           = _config_.save_sinograms;
    const double zoom_factor      = _config_.clustering_zoom_factor;
    const uint32_t iterations     = _config_.clustering_no_iterations;
    const uint32_t resolution_phi = _config_.clustering_resolution_phi;
    const uint32_t resolution_r   = _config_.clustering_resolution_r;
    const double max_precision_r  = _config_.clustering_max_initial_precision_r;
    const double sigma            = _config_.clustering_uncertainty;

    // TODO: think through the precision logic
    //int resolution_r = std::min(int(delta_R / max_precision_r), resolution_r);
    
    // enclosing tracker hits in a smallest possible rectangle (min_x, max_y) x (min_y, max_y)
    std::pair<std::vector<ConstTrackerHitHdl>::const_iterator,
              std::vector<ConstTrackerHitHdl>::const_iterator> minmax_X, minmax_Y;
    minmax_X = std::minmax_element(hits.begin(), hits.end(), [](const ConstTrackerHitHdl & hit1, const ConstTrackerHitHdl & hit2)
    { 
      return hit1->get_x() <= hit2->get_x();
    });
    minmax_Y = std::minmax_element(hits.begin(), hits.end(), [](const ConstTrackerHitHdl & hit1, const ConstTrackerHitHdl & hit2)
    { 
      return hit1->get_y() <= hit2->get_y(); 
    });
    
    double min_x = hits[minmax_X.first - hits.begin()]->get_x() - _geom_.tc_radius;
    double max_x = hits[minmax_X.second - hits.begin()]->get_x() + _geom_.tc_radius;
    double min_y = hits[minmax_Y.first - hits.begin()]->get_y() - _geom_.tc_radius;
    double max_y = hits[minmax_Y.second - hits.begin()]->get_y() + _geom_.tc_radius;
    
    // center of the box enclosing the tracker hits
    double center_X = (max_x + min_x) / 2.0;
    double center_Y = (max_y + min_y) / 2.0;
    
    // size of region (delta_phi x delta_R) to be investigated
    double delta_phi = M_PI;
    double delta_R = std::sqrt(std::pow(max_x - min_x, 2.0) + std::pow(max_y - min_y, 2.0));

    // peak_phi, peak_R store information about peak candidate
    double peak_phi = M_PI / 2.0;
    double peak_R = 0.0;

    for(int iter = 0; iter < iterations; iter++)
    {
      double r_min = peak_R - (delta_R / 2.0);
      double r_max = peak_R + (delta_R / 2.0);
      double phi_min = peak_phi - (delta_phi / 2.0);
      double phi_max = peak_phi + (delta_phi / 2.0);
      double offset = delta_phi / (2.0 * resolution_phi);

      TH2F sinograms("sinograms", "sinograms; phi; r",
         resolution_phi,
         phi_min + offset,
         phi_max + offset,
         resolution_r,
         r_min,
         r_max);

      for(auto & hit : hits)
      {
      // double sigma = hit->get_sigma_R(); // does not work well - better to have global higher sigma
      for(int k = 0; k <= resolution_phi; k++)
      {
        double phi = phi_min + ( k * delta_phi / double(resolution_phi) );
        // r - legendre transform of a center of a circle (Hough transform)
        double r = (hit->get_x() - center_X) * std::sin(phi) - (hit->get_y() - center_Y) * std::cos(phi);
        for(int half = 0; half < 2; half++)
        {	
          // mu - legendre transform of half circle (+R/-R)
          double mu = (r + (2.0 * half - 1.0) * hit->get_R());	

          // gauss is calculated only for -3 to 3 sigma region to cut time							
          double r1 = mu - 3.0 * sigma;
          double r2 = mu + 3.0 * sigma;

          // bin numbers coresponding to r1 and r2 values
          int bin1 = (double(resolution_r) * (r1 - r_min) / (r_max - r_min)) + 1;
          int bin2 = (double(resolution_r) * (r2 - r_min) / (r_max - r_min)) + 1;

          // real values of r coresponding to each bin 
          double r_j1 = r_min + (r_max - r_min) * double(bin1) / double(resolution_r);
          double r_j2;
          for(int binj = bin1; binj < bin2 + 1; binj++)
          {
            r_j2 = r_j1 + (r_max - r_min) / double(resolution_r);
            
            // average probability density in a bin given by gauss distribution with mean in mu 
            double weight = ( std::erf( (r_j2 - mu)/(std::sqrt(2.0)*sigma) ) - std::erf( (r_j1 - mu)/(std::sqrt(2.0)*sigma) ) ) 
                          / (2.0 * delta_R / double(resolution_r));
            
            // result is 2D histogram of several sinusoid functions f(phi) in convolution with gauss in r
            sinograms.Fill( phi, (r_j2 + r_j1) / 2.0, weight );
            r_j1 = r_j2;			
            }								
          }			
        }	
      }										

      // Get bin number of maximum value
      int maxBin = sinograms.GetMaximumBin();

      // Get X and Y values corresponding to the maximum bin
      int bin_phi, bin_R, bin_Z;
      sinograms.GetBinXYZ(maxBin, bin_phi, bin_R, bin_Z);
      peak_phi = sinograms.GetXaxis()->GetBinCenter(bin_phi);
      peak_R = sinograms.GetYaxis()->GetBinCenter(bin_R);

      delta_phi = delta_phi / zoom_factor;
      delta_R = delta_R / zoom_factor;

      if( save_sinograms ) 
      {
        TCanvas c2("sinograms", "sinograms", 1000, 800);
        c2.cd();
        sinograms.SetStats(0);
        sinograms.SetContour(100);
        sinograms.Draw("COLZ");
        c2.SaveAs(Form("Events_visu/clustering-run-%d_event-%d_side-%d_iter-%d.png",
                        _event_->get_run_number(),
                        _event_->get_event_number(),
                        hits.front()->get_SRL()[0],
                        iter));
        c2.Close();
      }
    }
    
    // result should be between -pi and pi
    if( peak_phi > M_PI/2.0 )
    {
      peak_phi -= M_PI;
      peak_R *= -1.0;
    }
    
    phi_estimate = peak_phi;
    r_estimate  = peak_R + center_X * std::sin(peak_phi) - center_Y * std::cos(peak_phi);
    return;
  }
  
  // step 3: MLM line fitting + ambiguity checking and solving
  void Algos::make_MLM_fits()
  {
    for(auto & precluster : _event_->get_preclusters())
    {
      // fits every prompt cluster
      if( precluster->is_prompt() )
      {
        for(auto & cluster : precluster->get_clusters())
        {
          // LinearFit object is created (with phi_estimate, r_estimate)
          LinearFitHdl primary_fit = std::make_shared<LinearFit>(cluster);
       		cluster->get_linear_fits().push_back(primary_fit);
       
          // check and detects the type of ambiguity of a cluster
          detect_ambiguity_type(cluster);
          
          // creates mirror fit obejct in case of ambiguous clusters with its own likelihood support structure
          create_mirror_fit(cluster);
          // for both fits (in case of ambiguity) finds ML estimate for (phi, r, theta, h)       
          for(auto & fit : cluster->get_linear_fits()) 
          {
            make_ML_estimate(fit);
          }
        }
      }
    }
    return;
  }
  
  void Algos::make_ML_estimate(LinearFitHdl fit)
  {
    // Newton's method for finding a root of a function applied on the derivative of likelihood to obtain optimal phi.
    // Likelihood is 4D function but 3 variables are solved analytically as a function of phi. 
    Likelihood & lik = fit->likelihood;
    double phi_0 = fit->phi;
    double h = 0.000001;
    for(int i = 0; i < 3; i++)
    {
      double derivative = lik.log_likelihood_derivative(phi_0);
      phi_0 += -h * derivative / (lik.log_likelihood_derivative(phi_0 + h) - derivative);
    }
    double r_0 = (lik.Rr + lik.Rx * std::sin(phi_0) - lik.Ry * std::cos(phi_0)) / lik.R;
    
    fit->phi = phi_0;
    fit->r = r_0;
    
    fit->a = std::tan(phi_0);
    fit->b = -r_0 / std::cos(phi_0);
    
    // calculating chi squared
    double min_likelihood = -lik.log_likelihood_value(phi_0);
    double no_of_measurements = lik.no_R + lik.no_Z;
    double chi_squared = min_likelihood / no_of_measurements;
    fit->chi_squared = chi_squared;
   
    // if at least 2 tracker hits have usable Z position ML fit is calculated
    if( lik.no_Z > 1 )
    {
      double denominator = lik.Cov_Zyy * std::pow(std::sin(phi_0), 2.0) 
                          + 2.0*(lik.Cov_Zxy) * std::sin(phi_0)*std::cos(phi_0) 
                          + lik.Cov_Zxx * std::pow(std::cos(phi_0), 2.0);
      
      if( denominator != 0.0)
      {
        double tan_theta = (lik.Cov_Zzy * std::sin(phi_0) + lik.Cov_Zzx * std::cos(phi_0)) / denominator;
        double h_temp = (lik.Zz / lik.Z) - (lik.Zx * std::cos(phi_0) + lik.Zy * std::sin(phi_0)) * tan_theta / lik.Z;
        
        fit->h = h_temp;
        fit->theta = std::atan(tan_theta);

        fit->c = tan_theta / std::cos(phi_0);
        fit->d = h_temp - r_0 * std::tan(phi_0) * tan_theta;
      }
    }
    // in case of only 1 tracker hit Z position the fit goes horizontally at the height of the one tracker hit
    else if( lik.no_Z == 1 ) 
    {
      fit->h = lik.Zz / lik.Z;
      
      fit->d = lik.Zz / lik.Z;
    } 
    
    return;
  }
  
  void Algos::detect_ambiguity_type(ClusterHdl cluster)
  {
    std::vector<ConstTrackerHitHdl> hits = cluster->get_tracker_hits();
    
    // detecting ambiguity type
    double x0 = hits.front()->get_x();
    double y0 = hits.front()->get_y();
    bool ambiguous;
    
    // type 1 == mirror image along line x = x0 
    ambiguous = true;
    for(auto i = 1u; i < hits.size(); i++)
    {	
      if( x0 != hits[i]->get_x() )
      {
        ambiguous = false;
        break;
      }
    }
    if( ambiguous == true )
    {
      cluster->set_ambiguity_type( 1 );
      return;
    }
    
    // type 2 == mirror image along line y = y0 
    ambiguous = true;
    for(auto i = 1u; i < hits.size(); i++)
    {	
      if( y0 != hits[i]->get_y() )
      {
        ambiguous = false;
        break;
      }
    }
    if( ambiguous == true )
    {
    cluster->set_ambiguity_type( 2 );
    	return;
    }
    
    // type 3 == mirror image along line y = x + (y0-x0) 
    ambiguous = true;
    for(auto i = 1u; i < hits.size(); i++)
    {	
      if( y0 - hits[i]->get_y() != x0 - hits[i]->get_x() )
      {
        ambiguous = false;
        break;
      }
    }
    if( ambiguous == true )
    {
      cluster->set_ambiguity_type( 3 );
      return;
    }
    
    // type 4 == mirror image along line y = -x + (y0-x0) 
    ambiguous = true;
    for(auto i = 1u; i < hits.size(); i++)
    {	
      if( y0 - hits[i]->get_y() != hits[i]->get_x() - x0 )
      {
        ambiguous = false;
       break;
      }
    }
    if( ambiguous == true )
    {
      cluster->set_ambiguity_type( 4 );
    	return;
    }

    cluster->set_ambiguity_type( 0 );
    return;
  }
  
  void Algos::create_mirror_fit(ClusterHdl cluster)
  {  
    if( cluster->get_ambiguity_type() == 0 )
    {
      return;
    }
    if( cluster->get_linear_fits().size() != 1 ) 
    {
      return;
    }
    
    ConstLinearFitHdl fit = cluster->get_linear_fits().front();
    const std::vector<ConstTrackerHitHdl> & hits = cluster->get_tracker_hits();
    double x0 = hits.front()->get_x();
    double y0 = hits.front()->get_y();
    
    // new copied LinearFit object 
    LinearFitHdl mirror_fit = std::make_shared<LinearFit>(fit); 
    cluster->get_linear_fits().push_back(mirror_fit);
    
    double a0 = fit->a;
    double b0 = fit->b;
    double mirror_a, mirror_b;
    
    switch(cluster->get_ambiguity_type())
    {
    case 1:
      mirror_a = -a0;
      mirror_b = 2.0 * a0 * x0 + b0;
      break;
    case 2:
      mirror_a = -a0;
      mirror_b = 2.0 * y0 - b0;
      break;
    case 3:
      mirror_a = 1.0 / a0;
      mirror_b = y0 - x0 - (b0 / a0) + (y0 - x0) / a0;
      break;
    case 4:
      mirror_a = 1.0 / a0;
      mirror_b = y0 + x0 + (b0 / a0) - (y0 + x0) / a0;
      break;
    }  
    
    double mirror_phi = std::atan(mirror_a);
    double mirror_r = -mirror_b / std::sqrt(mirror_a * mirror_a + 1.0);
	  
    mirror_fit->a = mirror_a;
    mirror_fit->b = mirror_b;
    mirror_fit->phi = mirror_phi;
    mirror_fit->r = mirror_r;
    
    cluster->set_phi_estimate( mirror_phi );
    cluster->set_r_estimate( mirror_r );
    
    // with the exception of type 1, mirrot fits correspond to the opposite halves of tracker hits
    // this changes the sign of the sums in the likelihood that contain drift radii
    if( cluster->get_ambiguity_type() != 1 )
    {  
      mirror_fit->likelihood.Rr *= -1.0;
      mirror_fit->likelihood.Rrx *= -1.0;
      mirror_fit->likelihood.Rry *= -1.0;
      mirror_fit->likelihood.Cov_Rrx *= -1.0;
      mirror_fit->likelihood.Cov_Rry *= -1.0;
    }
    
    return;
  }
  
  // step 4: Linear fits are associated to tracker hits and combined into a precluster solutions
  void Algos::combine_into_precluster_solutions()
  {
    // proccess is done separately for each precluster
    for(auto & precluster : _event_->get_preclusters())
    {
      // counting the number (N) of ambiguous clusters 
      int no_ambiguities = 0;
      for(auto & cluster : precluster->get_clusters())
      {
        if( cluster->get_ambiguity_type() )
        {
          no_ambiguities++;
        }
      }
      // TODO : limit the number of ambiguities? What is too much?
      // 2^N precluster solutions should exist (each combination of ambiguities)
      int no_precluster_solutions = 1 << no_ambiguities;
      std::vector<PreclusterSolutionHdl> & precluster_solutions = precluster->get_precluster_solutions();
      precluster_solutions.reserve( no_precluster_solutions );
      // creating empty solutions
      for(int i = 0; i < no_precluster_solutions; i++)
      { 
        PreclusterSolutionHdl precluster_solution = std::make_shared<PreclusterSolution>();
        precluster_solutions.push_back(precluster_solution);
        precluster_solution->get_unclustered_tracker_hits() = precluster->get_unclustered_tracker_hits();
      }
      
      // filling the precluster solutions with different track combinations
      int N1 = 1;
      for(auto & cluster : precluster->get_clusters())
      {
        // for each unambiguous cluster its single track is added to all solutions
        if(cluster->get_linear_fits().size() == 1u)
        {
          for(auto & precluster_solution : precluster_solutions)
          {  
            TrackHdl track = std::make_shared<Track>(cluster->get_linear_fits().front(), cluster->get_tracker_hits());
            track->sort_associations();
            precluster_solution->add_track(track);
          }
        }
        
        // for each ambiguous cluster its two tracks are added to different halves of solutions
        // example scheme for 3 ambiguities
        //  [0,0,0,0,1,1,1,1]
        //  [0,0,1,1,0,0,1,1]
        //  [0,1,0,1,0,1,0,1]
        else if(cluster->get_linear_fits().size() == 2u)
        {
          int k = 0;
          for(int i = 0; i < N1; i++)
          {
            int N2 = no_precluster_solutions / (N1 * 2);
            for(int j = 0; j < N2; j++)
            {  
              TrackHdl track = std::make_shared<Track>(cluster->get_linear_fits()[0], cluster->get_tracker_hits());
              track->sort_associations();
              precluster_solutions[k]->add_track(track);
              k++;
            }
            for(int j = 0; j < N2; j++)
            {  
              TrackHdl track = std::make_shared<Track>(cluster->get_linear_fits()[1], cluster->get_tracker_hits());
              track->sort_associations();
              precluster_solutions[k]->add_track(track);
              k++;
            }
          }
          N1 *= 2;
        }
      }
    }
    return;
  }
  
  // step 5 (straight track) alternative: Track -> Trajectories
  void Algos::create_line_trajectories()
  {
    // proccess is done separately for each precluster
    for(auto & precluster : _event_->get_preclusters())
    {
      // also separately for each precluster solution
      for(auto & precluster_solution : precluster->get_precluster_solutions())
      {
        std::vector<TrackHdl> & untrajectorized_tracks = precluster_solution->get_unprocessed_tracks();
        std::vector<TrajectoryHdl> & trajectories = precluster_solution->get_trajectories();
        for(auto & trajectory_candidate : untrajectorized_tracks)
        {
          // new Trajectory object for each found trajectory candidate
          TrajectoryHdl trajectory = std::make_shared<Trajectory>(trajectory_candidate);
          trajectories.push_back(trajectory);
          create_line_trajectory_points(trajectory);  
        }       
        untrajectorized_tracks.clear();
      }
    }
    return;
  }
  
  // step 5: Kink finding and connecting into polyline trajectories
  void Algos::create_polyline_trajectories()
  {
    // proccess is done separately for each precluster
    for(auto & precluster : _event_->get_preclusters())
    {
      // also separately for each precluster solution
      for(auto & precluster_solution : precluster->get_precluster_solutions())
      {
      
        std::vector<TrackHdl> unprocessed_tracks = precluster_solution->get_unprocessed_tracks(); 
        std::vector<std::vector<TrackHdl>> trajectory_candidates = find_polyline_candidates(unprocessed_tracks,
        																																				  precluster->get_side());
        // remove_fake_segments( trajectory_candidates ); // TODO needed??
        
        for(auto & trajectory_candidate : trajectory_candidates)
        {
          // new Trajectory object for each found trajectory candidate
          TrajectoryHdl trajectory = std::make_shared<Trajectory>(trajectory_candidate);
          precluster_solution->get_trajectories().push_back(trajectory);
          
          // creating trajectory points // TODO this could be one function (probably works)
          if(trajectory->has_kink())
          {
            //create_polyline_trajectory_points(trajectory);
            build_polyline_trajectory(precluster_solution, trajectory);
          }
          // line trajectories are more simple
          else
          {
            create_line_trajectory_points(trajectory);  
          }
          
          // remove trajectorized tracks from untrajectorized tracks (removing duplicates)       
          auto it = unprocessed_tracks.begin();
          while(it != unprocessed_tracks.end())
          {
            bool erased = false;   
            for(auto & track2 : trajectory_candidate)
            {
              if(*it == track2)
              {
                it = unprocessed_tracks.erase(it);
                erased = true;
                break;
              }
            }
            if(not erased)
            {
              it++;
            }
          }
        }       
      }
    }
    return;
  }
  
  std::vector<std::vector<TrackHdl>> Algos::find_polyline_candidates(std::vector<TrackHdl> & tracks, const int side) const
  {
    const double vertical_threshold = _config_.polylines_max_vertical_distance;
    const double min_distance = _config_.polylines_min_tracker_hits_distance;
    const double min_distance_from_foil = _config_.polylines_min_distance_from_foil;
  
    std::vector<std::vector<TrackHdl>> trajectory_candidates;
    int no_tracks = (int)tracks.size();
    if(no_tracks == 0)
    {
      return trajectory_candidates;
    }
    if(no_tracks == 1)
    {
       trajectory_candidates.push_back( std::vector<TrackHdl>(1, tracks.front()) );
       return trajectory_candidates;
    }
		
    // array stores the kink candidates - stores "true" for pairs (i, j)
    // of linear tracks that are candidates for connection into kinked trajectory  
    typedef boost::multi_array<bool, 2> array_type;
    typedef array_type::index index;
    array_type connections(boost::extents[no_tracks][no_tracks]);
    // Assign values to the elements:
    for(index ia = 0; ia != no_tracks; ++ia) 
      for(index ja = 0; ja != no_tracks; ++ja)
      {
    connections[ia][ja] = false;
      }

    // connection candidate (each pair (i, j) of tracks) are filtered based on several criteria 
    const snemo::geometry::gg_locator & ggLocator = _geom_.geo_loc->geigerLocator();
    // x coordinate of the outside border of tracker (side 1) 
    double tracker_x_max = ggLocator.getXCoordOfLayer(1, 8) + _geom_.tc_radius;
    // x coordinate of the inside border (source foil gap) of tracker (side 1)
    double tracker_x_min = std::max(ggLocator.getXCoordOfLayer(1, 0) - _geom_.tc_radius, 
    				     min_distance_from_foil);
    // y coordinate of the border (source foil gap) of tracker (side 1)
    double tracker_y_max = ggLocator.getYCoordOfRow(1, 112) + _geom_.tc_radius;
    
    
    for(int i = 0; i < no_tracks; i++)
    {
      TrackHdl track1 = tracks[i];
      for(int j = i+1; j < no_tracks; j++)
      {
        TrackHdl track2 = tracks[j];

        // for each pair of tracks calculate the (x,y) coordinates of the tracks intersection in 2D (kink point candidates)
        Point kink_point = get_intersection(track1, track2);
        
      // 1. kink position cuts
        // (checking if the (x, y) intersection is outside of the given side of tracker)
        if(kink_point.y < -tracker_y_max || tracker_y_max < kink_point.y) continue;
        if(kink_point.x < -tracker_x_max || tracker_x_max < kink_point.x) continue;
        if(side == 0 && kink_point.x > -tracker_x_min) continue;
        if(side == 1 && kink_point.x < tracker_x_min) continue;		
        
      // 2. vertical distance cut         
        // (calculating the vertical distance in the intersection point)		
        double z1 = track1->get_c() * kink_point.x + track1->get_d();
        double z2 = track2->get_c() * kink_point.x + track2->get_d();

        // track_connection_vertical_threshold
        if(std::abs(z2 - z1) > vertical_threshold) continue;			
        
        // z1 = 0 or z2 = 0 means one or both tracks have failed vertical reconstruction (missing timestamps,...) 
        if(z1 == 0 || z2 == 0) continue;
        
      // 3. associated tracker hit cuts
        // (connection is fake if no associated tracker hits are near the candidate kink point)
        
        PointHdl kink_point_hdl = std::make_shared<Point>(kink_point);
        // detecting close hits on track1
        bool match1_found = false;
        for(auto & association : track1->get_associations())
        {
          if( distance_2D( kink_point_hdl, association.point) < min_distance )
          {
            match1_found = true;
            break;
          }
        }     
        if( not match1_found ) continue;
        
        // detecting close hit on track2
        bool match2_found = false;
        for(auto & association : track2->get_associations())
        {
          if( distance_2D( kink_point_hdl, association.point) < min_distance )
          {
            match2_found = true;
            break;
          }
        }  
        if( not match2_found ) continue;
        
        connections[i][j] = true;
        connections[j][i] = true;			
      }
    }
    
    // connection_counter stores the number of kink candidates for each linear track
    std::vector<int> connection_counter(no_tracks, 0);
    for(int i = 0; i < no_tracks; i++)
      for(int j = 0; j < no_tracks; j++)
      {
        if(connections[j][i])
        {
          connection_counter[i]++;
        }
      }
    
    // connecting the tracks = "trajectorizing"
    std::vector<bool> trajectorized(no_tracks, false);
    for(int i = 0; i < no_tracks; i++)
    {
      // tracks with no kink candidates can be finalized as a trajectory
      if(connection_counter[i] == 0)
      {
        trajectory_candidates.push_back( std::vector<TrackHdl>(1, tracks[i]) );
        trajectorized[i] = true;
      }
      
      // we look for untrajectorized track with exactly 1 connection = beggining or end of the polyline track
      else if(connection_counter[i] == 1 && not trajectorized[i])
      {
        std::vector<TrackHdl> composite_track;
        composite_track.push_back(tracks[i]);
        trajectorized[i] = true;

        int next_index = i;
        // for loop ensure that the algorithm stops after no_tracks-1 steps in case that "connection_counter[next_index] == 1" fails
        for(int count = 0; count < no_tracks-1; count++)
        {
          for(int j = 0; j < no_tracks; j++)
          {
            if(connections[next_index][j] && not trajectorized[j])
            {
              next_index = j;
              composite_track.push_back(tracks[next_index]);
                    trajectorized[next_index] = true;
                  break;
            }
          }
          // checking for the end od polyline trajectory
          if(connection_counter[next_index] == 1)
          {
            break;
          }
        }
        trajectory_candidates.push_back(composite_track);
      }
    } 
    return trajectory_candidates;
  }
 
 
  // for a simple trajectory without kinks only the endpoint are created    
  void Algos::create_line_trajectory_points(TrajectoryHdl & trajectory)
  {
    DT_THROW_IF(trajectory->get_segments().size() != 1, std::logic_error, "No line trajectory for creating trajectory points");

    TrackHdl & track = trajectory->get_segments().front();
    std::vector<PointHdl> & trajectory_points = trajectory->get_trajectory_points();
        
    ConstPointHdl point1 = track->get_associations().front().point;
    ConstPointHdl point2 = track->get_associations().back().point;
    
    trajectory_points.push_back(std::make_shared<Point>(point1));
    trajectory_points.push_back(std::make_shared<Point>(point2));
     
    return;
  }
  
  
  // calculates the angle between the three point  
  double calculate_angle(const Point & point1, const Point & point2, const Point & point3)
  {
    Point vec1 = {point1.x - point2.x,
                  point1.y - point2.y,
                  point1.z - point2.z };
    Point vec2 = {point3.x - point2.x,
                  point3.y - point2.y,
                  point3.z - point2.z };
    double angle = vec1.x * vec2.x + vec1.y * vec2.y + vec1.z * vec2.z;
    double norm1 = std::sqrt(std::pow(vec1.x, 2) + pow(vec1.y, 2) + pow(vec1.z, 2));
    double norm2 = std::sqrt(std::pow(vec2.x, 2) + pow(vec2.y, 2) + pow(vec2.z, 2));
    angle /= (norm1 * norm2);
    angle = std::max(-1.0, std::min(1.0, angle));
    angle = std::acos(angle);
    angle = angle * 180.0 / M_PI;
    return angle;
  } 
  
  // takes a trajectory with multiple track segments and calculates the kink points and end points 
  void Algos::build_polyline_trajectory(PreclusterSolutionHdl & precluster_solution, TrajectoryHdl & trajectory)
  {
    DT_THROW_IF(trajectory->get_segments().size() < 2, std::logic_error, "No polyline trajectory for creating trajectory points");
    
    const double max_angle = _config_.polylines_max_kink_angle; 
    
    std::vector<PointHdl> & tr_points = trajectory->get_trajectory_points();
    std::vector<TrackHdl> & segments = trajectory->get_segments();

    // finding the kink points
    std::vector<PointHdl> kink_points;
    for(auto i = 0u; i < segments.size() - 1; i++)
    {
      Point intersection = get_intersection(segments[i], segments[i+1]);
      kink_points.push_back( std::make_shared<Point>(intersection) );
    }

    // first segment
    PointHdl kink_point = kink_points.front();
    PointHdl current_point;
    PointHdl alternative_point;
    PointHdl next_alternative_point;
    PointHdl next_kink_point;
    
    TrackHdl segment = segments.front();
    double t_kink = segment->calculate_t(*kink_point);
    int hit_count = segment->hit_split_counter(t_kink); 
    // TODO what if the counts are equal? (case hit_count == 0)
    std::vector<Association> & associations = segment->get_associations();
    if(hit_count > 0)
    {
      current_point = associations.back().point;
      alternative_point = associations.front().point;
    }
    else
    {
      current_point = associations.front().point;
      alternative_point = associations.back().point;
    }
    tr_points.push_back( std::make_shared<Point>(current_point) );


    // middle segments
    for(auto i = 1u; i < segments.size() - 1; i++)
    {
      segment = segments[i];
      next_kink_point = kink_points[i];
      double t_kink_next = segment->calculate_t(*next_kink_point);
      bool fake_next_connection = false;
       
      t_kink = segment->calculate_t(*kink_point);
      hit_count = segment->hit_split_counter(t_kink); 
      if(hit_count > 0)
      {
        next_alternative_point = segment->get_associations().back().point;
        if(t_kink_next < t_kink)
        {
          fake_next_connection = true;
        }
      }
      else
      {
        next_alternative_point = segment->get_associations().front().point;
        if(t_kink_next > t_kink)
        {
          fake_next_connection = true;
        }
      }
      
      double angle;
      if( not fake_next_connection )
      {
        angle = 180.0 - calculate_angle( *current_point, *kink_point, *next_kink_point );     
      }
      else 
      {
        angle = 180.0 - calculate_angle( *current_point, *kink_point, *next_alternative_point ); // TODO what if the hit does not have height??
      }
      
      if( angle < max_angle )
      {
        tr_points.push_back( std::make_shared<Point>(kink_point) );
        
        if( not fake_next_connection )
        {
          current_point = kink_point;
          kink_point = next_kink_point;  
          alternative_point = next_alternative_point;   
        }
        else
        {
          tr_points.push_back( std::make_shared<Point>(next_alternative_point) );
          
          std::vector<TrackHdl> new_trajectory_candidate(segments.begin() + i + 1, segments.end());
          TrajectoryHdl new_trajectory = std::make_shared<Trajectory>(new_trajectory_candidate);
          precluster_solution->get_trajectories().push_back( new_trajectory );
          
          segments.resize(i + 1);
          
          // recursion on the rest of the segments
          if(new_trajectory->get_segments().size() > 1)
          {
            new_trajectory->mark_as_kinked();
            build_polyline_trajectory( precluster_solution, new_trajectory );
          }
          else
          {
            create_line_trajectory_points( new_trajectory );
          }

          return;
        }                
      }
      else
      {
        tr_points.push_back( std::make_shared<Point>(alternative_point) );
        
        std::vector<TrackHdl> new_trajectory_candidate(segments.begin() + i, segments.end());
        TrajectoryHdl new_trajectory = std::make_shared<Trajectory>(new_trajectory_candidate);
        precluster_solution->get_trajectories().push_back( new_trajectory );
        
        segments.resize(i);
        if( i == 1 )
        {
          trajectory->mark_as_straight();
        }
        
        // recursion on the rest of the segments
        if(new_trajectory->get_segments().size() > 1)
        {
          new_trajectory->mark_as_kinked();
          build_polyline_trajectory( precluster_solution, new_trajectory );
        }
        else
        {
          create_line_trajectory_points( new_trajectory );
        }

        return;
      }
    }
    
    // last segment
    segment = segments.back();
    kink_point = kink_points.back();
    t_kink = segment->calculate_t(*kink_point);
    hit_count = segment->hit_split_counter(t_kink);
    
    PointHdl last_point;
    if(hit_count > 0)
    {
      last_point = segment->get_associations().back().point;
    }
    else
    {
      last_point = segment->get_associations().front().point;
    }
    
    double angle = 180.0 - calculate_angle( *current_point, *kink_point, *last_point );   
    if( angle < max_angle )
    {
      tr_points.push_back( std::make_shared<Point>(kink_point) );
      tr_points.push_back( std::make_shared<Point>(last_point) );
    }
    else
    {
      tr_points.push_back( std::make_shared<Point>(alternative_point) );
        
      TrajectoryHdl new_trajectory = std::make_shared<Trajectory>( segments.back());
      precluster_solution->get_trajectories().push_back( new_trajectory );
      
      segments.pop_back();
      if( segments.size() == 1 )
      {
        trajectory->mark_as_straight();
      }
      create_line_trajectory_points( new_trajectory );
    }
    return;
  }
  
  // step 6: Trajectory refinement
  //  - clustering refinement, trajectory connecting, fit quality calculation
  // (Potentially maximum likelihood refinement of the polyline fits)
  void Algos::refine_trajectories()
  {
  	// proccess is done separately for each precluster
    for(auto & precluster : _event_->get_preclusters())
    {
      // also separately for each precluster solution
      for(auto & precluster_solution : precluster->get_precluster_solutions())
      {
        for(auto & trajectory : precluster_solution->get_trajectories())
        {
          // polyline trajectories are handled differently
          if(trajectory->has_kink())
          {
            remove_wrong_hits_associations(precluster_solution, trajectory);
          }
        }
        
        // trying to add unclustered tracker hits to existing trajectories
        refine_clustering(precluster_solution);
        connect_close_trajectories(precluster_solution);
        
        // calculating all fit quality metrics
        for(auto & trajectory : precluster_solution->get_trajectories())
        {
          evaluate_trajectory(trajectory);
          /*
          std::cout << "	MSE: " << trajectory->get_MSE() << std::endl;
          std::cout << "	MSE_R: " << trajectory->get_MSE_R() << std::endl;
          std::cout << "	MSE_Z: " << trajectory->get_MSE_Z() << std::endl;
          std::cout << "	chi: " << trajectory->get_chi_squared() << std::endl;
          std::cout << "	chi_R: " << trajectory->get_chi_squared_R() << std::endl;
          std::cout << "	chi_Z: " << trajectory->get_chi_squared_Z() << "\n" << std::endl;*/
        }
      }
    }
    return;
  }
  
  // calculates mean square errors and chi squares of trajectories
  void Algos::evaluate_trajectory(TrajectoryHdl & trajectory)
  {   
    double chi_squared = 0.0;
    double chi_squared_R = 0.0;
    double chi_squared_Z = 0.0;
    
    double MSE = 0.0;
    double MSE_R = 0.0;
    double MSE_Z = 0.0;
  
    // number of measurements of R and Z used in the fit
    int no_R = 0;
    int no_Z = 0;
    
    // degrees of freedom
    int DOF_R = 0;
    int DOF_Z = 1;
    
    for(auto & segment : trajectory->get_segments())
    {
      segment->evaluate();
      
      const auto & associations = segment->get_associations();
      int no_R_seg = std::count_if(associations.begin(),
                           associations.end(),
                           [](const Association & association){
                             return association.tracker_hit->has_valid_R();
                           });
      no_R += no_R_seg;   
      int no_Z_seg = std::count_if(associations.begin(),
                           associations.end(),
                           [](const Association & association){
                             return association.tracker_hit->has_valid_Z();
                           });  
      no_Z += no_Z_seg;
      
      MSE_R += segment->get_square_error_R();
      MSE_Z += segment->get_square_error_Z();  
                    
      chi_squared_R += segment->get_chi_squared_R();
      chi_squared_Z += segment->get_chi_squared_Z();
      
      DOF_R += 2;
      DOF_Z += 1;
    }
    MSE = (MSE_R + MSE_Z);
    chi_squared = (chi_squared_R + chi_squared_Z);
    
    // mean square error is the sum of square errors of all segments divided by the number of measurements
    MSE_R = MSE_R / double(no_R);
    if(no_Z == 0)
    {
    	MSE_Z = datatools::invalid_real();
    }
    else
    {
      MSE_Z = MSE_Z / double(no_Z);
    }
    MSE = MSE / double(no_R + no_Z);
    
    // chi_squared is the sum of chi_squares of all segments divided by number of measurement - DOF
    if(no_R - DOF_R > 0)
    {
      chi_squared_R = chi_squared_R / double(no_R - DOF_R);    
    }
    else
    {
      chi_squared_R = datatools::invalid_real();
    }
    
    if(no_Z - DOF_Z > 0)
    {
      chi_squared_Z = chi_squared_Z / double(no_Z - DOF_Z);
    }
    else
    {
      chi_squared_Z = datatools::invalid_real();
    }
    
    if(no_R + no_Z - DOF_R - DOF_Z > 0)
    {
      chi_squared = chi_squared / double(no_R + no_Z - DOF_R - DOF_Z);    
    }
    else
    {
      chi_squared_Z = datatools::invalid_real();
    }

    trajectory->set_MSE( MSE );
    trajectory->set_MSE_R( MSE_R );
    trajectory->set_MSE_Z( MSE_Z );
    
    trajectory->set_chi_squared( chi_squared );
    trajectory->set_chi_squared_R( chi_squared_R );
    trajectory->set_chi_squared_Z( chi_squared_Z );
    return;
  }
  
  
  Point get_middle_point(const Point & point1, const Point & point2)
  {
    Point middle_point;
    middle_point.x = (point1.x + point2.x) / 2.0;
    middle_point.y = (point1.y + point2.y) / 2.0;
    middle_point.z = (point1.z + point2.z) / 2.0;
    return middle_point;
  } 
  
  
  // connecting very close trajectories with small angles
  void Algos::connect_close_trajectories(PreclusterSolutionHdl & precluster_solution)
  {
    const double distance_threshold = _config_.polylines_max_trajectories_middlepoint_distance;
    const double trajectory_connection_distance = _config_.polylines_max_trajectory_endpoints_distance;
    const double max_angle = _config_.polylines_max_trajectory_connection_angle;
    auto i = 0u;
    auto & trajectories = precluster_solution->get_trajectories(); 
    while(i < trajectories.size())
    {       
      bool connection_found = false;
      auto & trajectory_points1 = trajectories[i]->get_trajectory_points();
      for(auto j = i + 1; j < trajectories.size(); j++)
      {
        auto & trajectory_points2 = trajectories[j]->get_trajectory_points();
        
        // case A: front - front
        PointHdl point1 = trajectory_points1.front();
        PointHdl point2 = trajectory_points2.front();
        double distance = distance_3D(point1, point2); 
        if( distance < trajectory_connection_distance )
        {
          Point middle_point = get_middle_point(*point1, *point2);
          double angle = 180.0 - calculate_angle(*trajectory_points1[1],
                                                 middle_point,
                                                 *trajectory_points2[1]);
          if(angle < max_angle)
          {
            TrackHdl segment1 = trajectories[i]->get_segments().front();
            double distance1 = segment1->horizontal_distance_to_line(middle_point); 
            TrackHdl segment2 = trajectories[j]->get_segments().front();
            double distance2 = segment2->horizontal_distance_to_line(middle_point);             
            if( distance1 < distance_threshold && distance2 < distance_threshold )
            {
              std::reverse(trajectory_points1.begin(), trajectory_points1.end());
              trajectory_points1.pop_back();
              trajectory_points2.front() = std::make_shared<Point>(middle_point);
              
              auto & segments1 = trajectories[i]->get_segments();
              auto & segments2 = trajectories[j]->get_segments();
              std::reverse(segments1.begin(), segments1.end());

              trajectory_points1.insert(trajectory_points1.end(),
                                        trajectory_points2.begin(),
                                        trajectory_points2.end());           
              segments1.insert(segments1.end(), 
                               segments2.begin(), 
                               segments2.end());
              
              trajectories.erase(trajectories.begin() + j);
              trajectories[i]->mark_as_kinked();
              connection_found = true;
              break;
            }
          }
        }
        
        // case B: front - back
        point1 = trajectory_points1.front();
        point2 = trajectory_points2.back();
        distance = distance_3D(point1, point2); 
        if( distance < trajectory_connection_distance )
        {
          Point middle_point = get_middle_point(*point1, *point2);
          double angle = 180.0 - calculate_angle(*trajectory_points1[1] ,
                                                 middle_point,
                                                 *trajectory_points2[trajectory_points2.size()-1]);
          if(angle < max_angle)
          {
            TrackHdl segment1 = trajectories[i]->get_segments().front();
            double distance1 = segment1->horizontal_distance_to_line(middle_point); 
            
            TrackHdl segment2 = trajectories[j]->get_segments().back();
            double distance2 = segment2->horizontal_distance_to_line(middle_point);        
            if( distance1 < distance_threshold && distance2 < distance_threshold )
            {
              trajectory_points1.front() = std::make_shared<Point>(middle_point);
              trajectory_points2.pop_back();
              std::reverse(trajectory_points1.begin(), trajectory_points1.end());
              std::reverse(trajectory_points2.begin(), trajectory_points2.end());
              trajectory_points1.insert(trajectory_points1.end(),
                                        trajectory_points2.begin(),
                                        trajectory_points2.end());           
              
              auto & segments1 = trajectories[i]->get_segments();
              auto & segments2 = trajectories[j]->get_segments();
              std::reverse(segments1.begin(), segments1.end());
              std::reverse(segments2.begin(), segments2.end());
              segments1.insert(segments1.end(), 
                               segments2.begin(), 
                               segments2.end());
              
              trajectories.erase(trajectories.begin() + j);
              trajectories[i]->mark_as_kinked();
              connection_found = true;
              break;
            }
          }
        }
        
        // case C: back - front
        point1 = trajectory_points1.back();
        point2 = trajectory_points2.front();
        distance = distance_3D(point1, point2); 
        if( distance < trajectory_connection_distance )
        {
         Point middle_point = get_middle_point(*point1, *point2);
          double angle = 180.0 - calculate_angle(*trajectory_points1[trajectory_points1.size()-1] ,
                                                 middle_point,
                                                 *trajectory_points2[1]);
          if(angle < max_angle)
          {
            TrackHdl segment1 = trajectories[i]->get_segments().back();
            double distance1 = segment1->horizontal_distance_to_line(middle_point); 
            
            TrackHdl segment2 = trajectories[j]->get_segments().front();
            double distance2 = segment2->horizontal_distance_to_line(middle_point);   
            if( distance1 < distance_threshold && distance2 < distance_threshold )
            {
              trajectory_points1.back() = std::make_shared<Point>(middle_point);
              
              auto & segments2 = trajectories[j]->get_segments();
              auto & segments1 = trajectories[i]->get_segments();

              trajectory_points1.insert(trajectory_points1.end(),
                                        trajectory_points2.begin() + 1,
                                        trajectory_points2.end());       
                                            
              segments1.insert(segments1.end(), 
                               segments2.begin(), 
                               segments2.end());
              
              trajectories.erase(trajectories.begin() + j);
              trajectories[i]->mark_as_kinked();
              connection_found = true;
              break;
            }
          }
        }
        
        // case D: back - back
        point1 = trajectory_points1.back();
        point2 = trajectory_points2.back();
        distance = distance_3D(point1, point2); 
        if( distance < trajectory_connection_distance )
        {
          Point middle_point = get_middle_point(*point1, *point2);
          double angle = 180.0 - calculate_angle(*trajectory_points1[trajectory_points1.size()-1] ,
                                                  middle_point,
                                                 *trajectory_points2[trajectory_points2.size()-1]);
          if(angle < max_angle)
          {
            TrackHdl segment1 = trajectories[i]->get_segments().back();
            double distance1 = segment1->horizontal_distance_to_line(middle_point); 
            
            TrackHdl segment2 = trajectories[j]->get_segments().back();
            double distance2 = segment2->horizontal_distance_to_line(middle_point);
            if( distance1 < distance_threshold && distance2 < distance_threshold )
            {
              trajectory_points1.back() = std::make_shared<Point>(middle_point);
              trajectory_points2.pop_back();
              std::reverse(trajectory_points2.begin(), trajectory_points2.end());
              
              auto & segments2 = trajectories[j]->get_segments();
              auto & segments1 = trajectories[i]->get_segments();
              std::reverse(segments2.begin(), segments2.end());

              trajectory_points1.insert(trajectory_points1.end(),
                                        trajectory_points2.begin(),
                                        trajectory_points2.end());           
              segments1.insert(segments1.end(), 
                               segments2.begin(), 
                               segments2.end());
              
              trajectories.erase(trajectories.begin() + j);
              trajectories[i]->mark_as_kinked();
              connection_found = true;
              break;
            }
          }
        }
      }
      if(not connection_found)
      {
        i++;
      }
      else
      {
        trajectories[i]->update_segments();
      }

    }
    return;
  }
  
  
  // checks whether all tracker hits are associated to correct segment inside a polyline trajectory
  // if not - associates it to the correct segment or removes it from the trajectory
  void Algos::remove_wrong_hits_associations(PreclusterSolutionHdl & precluster_solution, TrajectoryHdl & trajectory)
  {
    const double tolerance = 0.01;
    std::vector<TrackHdl> & segments = trajectory->get_segments();
    std::vector<PointHdl> & traj_points = trajectory->get_trajectory_points();
    std::vector<ConstTrackerHitHdl> & unclustered_tracker_hits = precluster_solution->get_unclustered_tracker_hits();
    
    for(auto i = 0u; i < segments.size(); i++)
    {
      ConstPointHdl start = traj_points[i];
      ConstPointHdl end   = traj_points[i+1];
      double t_start = segments[i]->calculate_t(start);
      double t_end = segments[i]->calculate_t(end);
      if(t_start > t_end)
      {      
        std::swap(t_start, t_end);
      }
     
      std::vector<Association> & associations = segments[i]->get_associations(); 
      auto j = 0u;
      while(j < associations.size())
      {
        ConstTrackerHitHdl hit = associations[j].tracker_hit;
        
        double t_j = associations[j].parameter_t;
        // checking if the point is in bounds of its segments (correct association)
        if( t_start <= t_j + tolerance && t_j - tolerance <= t_end ) 
        {
          j++;
          continue;
        }
        // erasing the hit association point from the trajectory
        associations.erase(associations.begin() + j);
        
        // marking the hit as unclustered
        unclustered_tracker_hits.push_back(hit);
      }
    }
    return;
  }   
    
  // TODO improve the association criteria for hits near kinks!!!
  // goes through unclustered tracker hits of a precluster solution and checks
  // if it can be associated to some segment of some trajectory
  void Algos::refine_clustering(PreclusterSolutionHdl & precluster_solution)
  {
    const double distance_threshold = _config_.clustering_hit_association_distance;
    const double max_distance = _config_.polylines_max_extention_distance;
    
    auto & unclustered_tracker_hits = precluster_solution->get_unclustered_tracker_hits();
    auto & trajectories = precluster_solution->get_trajectories();
    
    if(unclustered_tracker_hits.size() == 0 || trajectories.size() == 0) return; 
    
    // computing parametric bounds of segments of trajectories
    std::vector<std::vector<double>> t_starts;
    std::vector<std::vector<double>> t_ends;
    t_starts.reserve(trajectories.size());
    t_ends.reserve(trajectories.size());
    for(auto & trajectory : trajectories)
    {
      auto & traj_points = trajectory->get_trajectory_points(); 
      auto & segments = trajectory->get_segments();
      std::vector<double> trajectory_starts;
      std::vector<double> trajectory_ends;
      trajectory_starts.reserve(segments.size());
      trajectory_ends.reserve(segments.size());

      for(auto i = 0u; i < segments.size(); i++)
      {
        ConstPointHdl start = traj_points[i];
        ConstPointHdl end = traj_points[i+1];
        
        double t_start = segments[i]->calculate_t(start);
        double t_end = segments[i]->calculate_t(end);
        trajectory_starts.push_back(t_start);
        trajectory_ends.push_back(t_end);
      }
      t_starts.push_back(trajectory_starts);
      t_ends.push_back(trajectory_ends);
    }
    
    // looking for good associations
    auto i = 0u;
    while(i < unclustered_tracker_hits.size())
    {
      ConstTrackerHitHdl hit = unclustered_tracker_hits[i];
      double hit_x = hit->get_x();
      double hit_y = hit->get_y();
      double hit_R = hit->get_R();
      bool clustered = false;
      for(auto j = 0u; j < trajectories.size(); j++)
      {
        if( clustered ) break;
        
        auto & segments = trajectories[j]->get_segments();
        auto & traj_points = trajectories[j]->get_trajectory_points(); 
        for(auto k = 0u; k < segments.size(); k++)
        {
          double cos_phi = std::cos(segments[k]->get_phi());
          double sin_phi = std::sin(segments[k]->get_phi());
          
          // checking if the hit is close to the line (segment without bounds)
          //TODO possible to add vertical distance as well
          double distance = std::abs(segments[k]->get_r() - hit_x * sin_phi + hit_y * cos_phi) - hit_R;
          if( std::abs(distance) > distance_threshold )
          { 
                continue;
          }
                  
         // checking if the hit is within bounds of the segment
          double t = hit_x * cos_phi + hit_y * sin_phi;
          bool association_found = false;
          bool trajectory_extention_needed = false;
          auto position = 0u;
          std::vector<Association> & associations_k = segments[k]->get_associations();
          PointHdl end_point;   
       
          // tracker hit is in the middle of a segment
          if(std::min(t_starts[j][k], t_ends[j][k]) <= t &&
             std::max(t_starts[j][k], t_ends[j][k]) >= t)
          {
            association_found = true;
            
            // finding where the hit belong on the line (tracker_hits and tracker_hit_points are sorted vectors)
            while(position < associations_k.size())
            {
              if(t < associations_k[position].parameter_t)
              {
                break;
              }
              else
              {
                position++;           
              }
            }        
          }
          
          // first and last segment has to be handled separately - new hit can prolong the segment
          // first segment case
          else if( k == 0 )
          {
            if( (t_ends[j][k] > t_starts[j][k] && t < t_starts[j][k]) ||
                (t_ends[j][k] < t_starts[j][k] && t >= t_starts[j][k]) )
            {
              double distance = std::abs(t - t_starts[j][k]);
              if(distance < max_distance)
              {
                association_found = true;
                trajectory_extention_needed = true;
                position = 0;
                end_point = traj_points.front();
              }
            }
          }
          
          // last segment case
          else if( k == segments.size() - 1 )
          {
            if( (t_ends[j][k] > t_starts[j][k] && t >= t_ends[j][k]) ||
                (t_ends[j][k] < t_starts[j][k] && t <= t_ends[j][k]) )
            {
              double distance = std::abs(t - t_ends[j][k]);
              if(distance < max_distance)
              {
                association_found = true;
                trajectory_extention_needed = true;
                position = associations_k.size();
                end_point = traj_points.back();
              }
            }
          } 
          
          if( association_found )
          {
            Association new_association(hit);
            new_association.parameter_t = t;
            
            // creating new association point and storing at the right place
            double x = t * cos_phi + segments[k]->get_r() * sin_phi;
            double y = t * sin_phi - segments[k]->get_r() * cos_phi;
            double z = t * std::tan(segments[k]->get_theta()) + segments[k]->get_h();
            new_association.point = std::make_shared<Point>(x, y, z); 
            associations_k.insert(associations_k.begin() + position, new_association);
            
            if(trajectory_extention_needed)
            {
              end_point->x = x;
              end_point->y = y;
              end_point->z = z;
            }

            unclustered_tracker_hits.erase(unclustered_tracker_hits.begin() + i);
            clustered = true;
            break;
          }
        }
      }
      if(not clustered)
      {
        i++;
      }
    }
    
    return;
  }
  

  // step 7: Combining precluster solutions into all solutions
  void Algos::create_solutions()
  {
    // caluclating needed number of solutions
    int no_solutions = 1;
    for(auto & precluster : _event_->get_preclusters())
    {
      int size_of_precluster = precluster->get_precluster_solutions().size(); 
      if( size_of_precluster != 0 )
      {
        no_solutions *= size_of_precluster;
      }
    }
    if(no_solutions == 0) return;
    // TODO: limit the number of solutions?
    
    std::vector<SolutionHdl> & solutions = _event_->get_solutions();
    solutions.reserve( no_solutions );
    solutions.push_back( std::make_shared<Solution>() );
    
    // constructing each solution as a unique combination of precluster solutions
    // each solution points to exactly one precluster solution of each cluster
    for(auto & precluster : _event_->get_preclusters())
    {
      std::vector<PreclusterSolutionHdl> & precluster_solutions = precluster->get_precluster_solutions(); 
      
      int no_copied_solution = solutions.size();
      for(auto i = 1u; i < precluster_solutions.size(); i++)
      {
        for(int j = 0; j < no_copied_solution; j++)
        {
          solutions.push_back(std::make_shared<Solution>( solutions[j]->get_precluster_solutions() ));
        }
      }
      
      for(auto i = 0u; i < precluster_solutions.size(); i++)
      {
        for(int j = 0; j < no_copied_solution; j++)
        {
          auto & prec_sol = solutions[i * no_copied_solution + j]->get_precluster_solutions();
          prec_sol.push_back( precluster_solutions[i] );
        }
      }
    }
    
    // sorts the solutions based on overall reconstruction qualities
    // better solutions are first (preliminary)
    sort_solutions(solutions);

    // TODO: possible to add solution discard criteria 
    
    return;
  }
  
  bool compare_solutions(const SolutionHdl solution1, const SolutionHdl solution2)
  {
    int no_trajectories1 = 0;
    int no_segments1 = 0;
    int no_unclustered_hit1 = 0;
    double chi2_sum1 = 0.0;

    for(auto & precluster_solution : solution1->get_precluster_solutions())
    {
      auto & trajectories = precluster_solution->get_trajectories();
      no_trajectories1 += trajectories.size();
      for(auto & trajectory : trajectories)
      {
        no_segments1 += trajectory->get_segments().size();
        chi2_sum1 += trajectory->get_chi_squared();
      }
      no_unclustered_hit1 += precluster_solution->get_unclustered_tracker_hits().size();
    }
    
    int no_trajectories2 = 0;
    int no_segments2 = 0;
    int no_unclustered_hit2 = 0;
    double chi2_sum2 = 0.0;
    
    for(auto & precluster_solution : solution2->get_precluster_solutions())
    {
      auto & trajectories = precluster_solution->get_trajectories();
      no_trajectories2 += trajectories.size();
      for(auto & trajectory : trajectories)
      {
        no_segments2 += trajectory->get_segments().size();
        chi2_sum2 += trajectory->get_chi_squared();
      }
      no_unclustered_hit2 += precluster_solution->get_unclustered_tracker_hits().size();
    }
    
    //TODO how shoudl this work??
    if( no_segments1 == no_segments2 )
    {
      if( no_trajectories1 == no_trajectories2 )
      {
        return (chi2_sum2 > chi2_sum1); 
      }
      return (no_trajectories2 > no_trajectories1);
    }
    else
    {
      return (no_segments1 < no_segments2);
    }
  }
  
  
  void Algos::sort_solutions(std::vector<SolutionHdl> & solutions)
  {
    std::sort(solutions.begin(), solutions.end(), compare_solutions); 
  }
  
  // separates a vector of tracker hits into spatialy distant subgroups (vectors)
  std::vector<std::vector<ConstTrackerHitHdl>> Algos::separate_hits(const std::vector<ConstTrackerHitHdl> & hits)
  {
    std::vector<std::vector<ConstTrackerHitHdl>> clusters_temp;
    clusters_temp.reserve( hits.size() );
    for(const auto & hit : hits)
    {
      clusters_temp.push_back( std::vector<ConstTrackerHitHdl>(1, hit) );
    }
    
    for(auto & cluster1 : clusters_temp)
    {
      for(auto & cluster2 : clusters_temp)
      {
        if(cluster1 != cluster2 && clusters_close(cluster1, cluster2))
        {
          std::move(cluster2.begin(), cluster2.end(), std::back_inserter(cluster1));
          cluster2.clear(); 
        }
      }
    }
    
    std::vector<std::vector<ConstTrackerHitHdl>> clusters_out;
    for(auto & cluster : clusters_temp)
    {
      if(!cluster.empty())
      {
        clusters_out.push_back(cluster);
      }
    }
    return clusters_out;	
  }
  
  // checks if the minimum distance between triggered anode wires of two clusters is lower than "preclustering_distance_treshold"
  bool Algos::clusters_close(const std::vector<ConstTrackerHitHdl> & cluster1,
                             const std::vector<ConstTrackerHitHdl> & cluster2) const
  {
  	const double distance_treshold = _config_.clustering_max_distance;
    for(const auto & hit1 : cluster1)
    {
      for(const auto & hit2 : cluster2)
      {
        // treshold distance of 3 tracker cells
        if( std::pow(hit1->get_x() - hit2->get_x(), 2.0) + std::pow(hit1->get_y() - hit2->get_y(), 2.0) < std::pow(distance_treshold, 2.0) )
        {
          return true;
        }
      }
    }	
    return false; 
  }



/*
  TKclusterHdl Algos::find_cluster(const std::vector<TKtrhitHdl> & hits)
  {
    TKclusterHdl foundClusterHdl;
    //number of different values of phi among which the cluster is being searched for
    const int bins_phi = 100;
    const auto sideId = hits.front()->get_SRL('s');
    for(auto i = 0u; i < hits.size(); i++)
      {
        if( sideId != hits[i]->get_SRL('s'))
          {
            DT_LOG_DEBUG(_config_.verbosity, "cluster not found - hits from both sides included");
            return foundClusterHdl;
          }
      }

    if( hits.size() < 3 ) 
      {
        DT_LOG_DEBUG(_config_.verbosity, "cluster not found - too few usable hits");
        return foundClusterHdl;
      }
        
    int max_count[bins_phi] = {0};
    double argmax_R[bins_phi];
    int global_max = 0;
    for(int i = 0; i < bins_phi; i++)
      {
        std::vector<double> boundaries;
        std::vector<int> hit_count;
        double theta = i * M_PI / double(bins_phi);             
        double boundary_down, boundary_up;
                
        boundaries.push_back(-10000);
        hit_count.push_back(0);
        boundaries.push_back(10000);
                
        for(auto j = 0u; j < hits.size(); j++)
          {
            const auto & hitPtr = hits.at(j);
            double x = hitPtr->get_xy('x');
            double y = hitPtr->get_xy('y');
            if(i < bins_phi/2)
              {
                boundary_down = (x-22.0) * sin(theta) - (y+22.0) * cos(theta);
                boundary_up   = (x+22.0) * sin(theta) - (y-22.0) * cos(theta);
              }
            else
              {         
                boundary_down = (x-22.0) * sin(theta) - (y-22.0) * cos(theta);
                boundary_up   = (x+22.0) * sin(theta) - (y+22.0) * cos(theta);  
              }
                        
            int k = 0;
            while( boundary_down > boundaries.at(k) )
              {
                k++;
              } 
            if( boundary_down < boundaries.at(k) )
              {
                boundaries.insert(boundaries.begin()+k, boundary_down);
                hit_count.insert(hit_count.begin()+k, hit_count.at(k-1));
              }
            while( boundary_up > boundaries.at(k) )                     
              {
                hit_count.at(k)++;
                k++;
              }
            if( boundary_up < boundaries.at(k) )
              {
                boundaries.insert(boundaries.begin()+k, boundary_up);
                hit_count.insert(hit_count.begin()+k, hit_count.at(k-1)-1);
              }

          }
                
        for(auto j = 0u; j < hit_count.size(); j++)
          {
            if(hit_count.at(j) > max_count[i])
              {
                max_count[i] = hit_count.at(j);
                argmax_R[i] = (boundaries.at(j) + boundaries.at(j+1))/2.0;
              }
          }

        if(max_count[i] > global_max)
          global_max = max_count[i];
      }

    int phi_bin_min = 0;
    int phi_bin_max = 0;        
    double phi_min;
    double phi_max;
        
    if(max_count[0] == global_max)
      {
        while(max_count[phi_bin_max] == global_max)
          {
            phi_bin_max++;
          }
        phi_bin_min = bins_phi-1;
        while(max_count[phi_bin_min] == global_max)
          {
            phi_bin_min--;
          }

        phi_max = phi_bin_max * M_PI / double(bins_phi);
        phi_min = phi_bin_min * M_PI / double(bins_phi) - M_PI;
      }
    else
      {
        while(max_count[phi_bin_min] < global_max)
          {
            phi_bin_min++;
          }             
        phi_bin_max = phi_bin_min;
        while(max_count[phi_bin_max] == global_max)
          {
            phi_bin_max++;
          }

        phi_min = (phi_bin_min-1.0) * M_PI / double(bins_phi);
        phi_max = (phi_bin_max) * M_PI / double(bins_phi);
      }

    double R_0 = argmax_R[phi_bin_max-1];
    double theta_0 = double(phi_bin_max-1) * M_PI / double(bins_phi);
        
    std::vector<TKtrhitHdl> cluster_hits;
    for(auto j = 0u; j < hits.size(); j++)
      {
        auto & hitPtr = hits.at(j);
        double boundary_down, boundary_up;
        double x = hitPtr->get_xy('x');
        double y = hitPtr->get_xy('y');
        if(phi_bin_max-1 < bins_phi/2)
          {
            boundary_down = (x-22.0) * sin(theta_0) - (y+22.0) * cos(theta_0);
            boundary_up   = (x+22.0) * sin(theta_0) - (y-22.0) * cos(theta_0);
          }
        else
          {             
            boundary_down = (x-22.0) * sin(theta_0) - (y-22.0) * cos(theta_0);
            boundary_up   = (x+22.0) * sin(theta_0) - (y+22.0) * cos(theta_0);  
          }
        if(boundary_down <= R_0 && R_0 <= boundary_up)
          {
            cluster_hits.push_back(hitPtr);
          }
      }
    cluster_hits = TKEvent::filter_distant(cluster_hits);
    if(cluster_hits.size() < 3)
      {
        return foundClusterHdl;
      }
    foundClusterHdl = std::make_shared<TKcluster>(cluster_hits, phi_min, phi_max);
    return foundClusterHdl;
  }*/
  
  
  
   /*
   // distance check for trakcer hits - if the kink candidate does 
  // not have neighbouring cells it is probably fake candidate
  void Algos::split_fake_polyline_candidates(std::vector<std::vector<TrackHdl>> & trajectory_candidates )
  {
    const double min_distance = _config_.min_distance_threshold;
    
    for(auto i = 0u; i < trajectory_candidates.size(); i++)
    {
      // finding the kink points
      for(auto j = 0u; j < trajectory_candidates[i].size()-1; j++)
      {
        TrackHdl segment1 = trajectory_candidates[i][j];
        TrackHdl segment2 = trajectory_candidates[i][j+1];
        
        Point point = get_intersection(segment1, segment2);
        PointHdl kink_point = std::make_shared<Point>(point);
      
        bool match1_found = false;
        for(auto & association : segment1->get_associations())
        {
          if( distance_2D( kink_point, association.point) < min_distance )
          {
            match1_found = true;
            break;
          }
        }      
        bool match2_found = false;
        for(auto & association : segment2->get_associations())
        {
          if( distance_2D( kink_point, association.point) < min_distance )
          {
            match2_found = true;
            break;
          }
        }  
        
        if( match1_found && match2_found )
        {
          continue;
        }
        else
        {
          std::vector<TrackHdl> new_trajectory_candidate(std::make_move_iterator(trajectory_candidates[i].begin() + j + 1),
                                                         std::make_move_iterator(trajectory_candidates[i].end()));
          trajectory_candidates[i].resize(j + 1);
          trajectory_candidates.push_back( new_trajectory_candidate );
        }
      } 
    }
    return;
  }
 */
 
 /*
  // checks if polyline trajectory candidates with at least 2 kinks even have 
  // some tracker hits in between the kinks (if not, they are considered as fake) 
  void Algos::remove_fake_segments(std::vector<std::vector<TrackHdl>> & trajectory_candidates )
  {
    for(auto i = 0u; i < trajectory_candidates.size(); i++)
    {
      if( trajectory_candidates[i].size() < 3u) continue;
      
      TrackHdl segment1 = trajectory_candidates[i][0];
      TrackHdl segment2 = trajectory_candidates[i][1];
      Point kink_point1 = get_intersection( segment1, segment2 );
      Point kink_point2;
      for(auto j = 1u; j < trajectory_candidates[i].size()-1; j++)
      {
        segment1 = trajectory_candidates[i][j];
        segment2 = trajectory_candidates[i][j+1];
        kink_point2 = get_intersection( segment1, segment2 );    
  
        double t_start = segment1->calculate_t(kink_point1);
        double t_end = segment1->calculate_t(kink_point2);
        if(t_start > t_end)
        {
          std::swap(t_start, t_end);
        }
        
        bool hit_found = false;
        for(const auto & association : segment1->get_associations())
        {
          if( t_start < association.parameter_t 
             && t_end > association.parameter_t )
          {
            hit_found = true;
            break;
          }
        }  
        std::cout << "middle segment" << std::endl;
        if( not hit_found ) 
        {
          std::cout << "fake found!!!!!!!!!!!!!!!" << std::endl;
        }
          
        std::swap(kink_point1, kink_point2);
      }
    
    }
    return;
  }
  */
  
  /*
   // takes a trajectory with multiple track segments and calculates the kink points and end points 
  void Algos::create_polyline_trajectory_points(TrajectoryHdl & trajectory)
  {
    DT_THROW_IF(trajectory->get_segments().size() < 2, std::logic_error, "No polyline trajectory for creating trajectory points");
    
    std::vector<PointHdl> & tr_points = trajectory->get_trajectory_points();
    std::vector<TrackHdl> & segments = trajectory->get_segments();

    tr_points.reserve(segments.size() + 1);
    tr_points.push_back(std::make_shared<Point>());
    
    // finding the kink points
    for(auto i = 0u; i < segments.size()-1; i++)
    {
      Point point = get_intersection(segments[i], segments[i+1]);
      PointHdl kink_point = std::make_shared<Point>(point);
      tr_points.push_back( kink_point );      
    }
    
    // finding the start-point of the polyline
    std::vector<Association> & associations1 = segments.front()->get_associations();
    PointHdl & kink_point_front = tr_points[1];
    
    double cos_phi = std::cos(segments.front()->get_phi());
    double sin_phi = std::sin(segments.front()->get_phi());  
    double t_kink = kink_point_front->x * cos_phi + kink_point_front->y * sin_phi;
    int hit_count_positive = 0;
    int hit_count_negative = 0;
    for(auto j = 0u; j < associations1.size(); j++)
    {
      if(associations1[j].parameter_t > t_kink)
      {
        hit_count_positive++;
      }
      else
      {
        hit_count_negative++;
      }
    }
    
    if(hit_count_positive > hit_count_negative)
    {
      tr_points.front() = std::make_shared<Point>( associations1.back().point );
    }
    // TODO what if the counts are equal?
    else
    {
      tr_points.front() = std::make_shared<Point>( associations1.front().point );
    }
    
    
    // TODO add kink limiter
    // finding the end-point of the polyline
    std::vector<Association> & associations2 = segments.back()->get_associations();
    PointHdl & kink_point_back = tr_points[segments.size()-1];

    cos_phi = std::cos(segments.back()->get_phi());
    sin_phi = std::sin(segments.back()->get_phi());  
    t_kink = kink_point_back->x * cos_phi + kink_point_back->y * sin_phi;
    hit_count_positive = 0;
    hit_count_negative = 0;
    for(auto j = 0u; j < associations2.size(); j++)
    {
      if(associations2[j].parameter_t > t_kink)
      {
        hit_count_positive++;
      }
      else
      {
        hit_count_negative++;
      }
    }
    
    if(hit_count_positive > hit_count_negative)
    {
      tr_points.push_back( std::make_shared<Point>(associations2.back().point) );
    }
    else
    {
      tr_points.push_back( std::make_shared<Point>(associations2.front().point) );
    }
    
    return;
  }*/
  
  
  
} //  end of namespace tkrec

