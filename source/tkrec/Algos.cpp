// Cimrman headers
#include "tkrec/Algos.h"

// XXX Temporary root headers
#include <TCanvas.h>
#include <TH1F.h>

// Standard headers
#include <algorithm>
#include <stack>
#include <cmath>

namespace tkrec {
 
  ElectronRecMode electron_recmode_from_label(const std::string & label_)
  {
    if (label_ == "polyline") return ElectronRecMode::polyline;
    if (label_ == "line") return ElectronRecMode::line;
    return ElectronRecMode::undefined;
  }
 
  void CimrmanAlgoConfig::parse(const datatools::properties & config_)
  {
    //=======================================================
    //################## general algo config ################
    //=======================================================
  
    if (config_.has_key("verbosity")) {
      std::string verbosityLabel = config_.fetch_string("verbosity");
      auto parsedVerb = datatools::logger::get_priority(verbosityLabel);
      DT_THROW_IF(parsedVerb == datatools::logger::PRIO_UNDEFINED,
                  std::logic_error,
                  "Undefined verbosity label " << std::quoted(verbosityLabel));
      this->verbosity = parsedVerb;
    }

    DT_LOG_DEBUG(this->verbosity, "Loading configuration...");

    if (config_.has_key("electron_mode")) {
      std::string elRecModeLabel = config_.fetch_string("electron_mode");
      auto elRecMode = electron_recmode_from_label(elRecModeLabel);
      DT_THROW_IF(elRecMode == ElectronRecMode::undefined,
                  std::logic_error,
                  "Undefined event reconstruction mode " << std::quoted(elRecModeLabel));
      this->electron_mode = elRecMode;
    }
    
    if (config_.has_flag("visualization_2D")) {
      this->visualization_2D = true;
      DT_LOG_DEBUG(this->verbosity, "Saving 2D visualizations");
    }
    
    if (config_.has_flag("visualization_3D")) {
      this->visualization_3D = true;
      DT_LOG_DEBUG(this->verbosity, "Saving 3D visualizations");
    }
    
    if (config_.has_flag("use_provided_preclustering")) {
      this->use_provided_preclustering = true;
      DT_LOG_DEBUG(this->verbosity, "Using provided TCD bank");
    }
    
    if (config_.has_flag("reconstruct_alphas")) {
      this->reconstruct_alphas = true;
      DT_LOG_DEBUG(this->verbosity, "Alpha tracking enabled");
    }
    
    if (config_.has_flag("keep_only_best_solutions")) {
      this->keep_only_best_solutions = true;
      DT_LOG_DEBUG(this->verbosity, "Keeping only best solutions");
    }
  
    if (config_.has_flag("force_default_sigma_r")) {
      this->force_default_sigma_r = true;
      DT_LOG_DEBUG(this->verbosity, "Default sigma used");
    }
  
    if (config_.has_key("default_sigma_r")) {
      auto value = config_.fetch_real_with_explicit_dimension("default_sigma_r",
							        "length");
      DT_THROW_IF(value <= 1.0e-3 * CLHEP::mm or value >= 50.0 * CLHEP::mm,
                  std::logic_error,
                  "Invalid default sigma for drift radius"); 
      this->default_sigma_r = value / CLHEP::mm;
    }
  
    /*
    if (config_.has_key("chi_square_threshold")) {
      this->chi_square_threshold =
        config_.fetch_dimensionless_real("chi_square_threshold");
      DT_THROW_IF(this->chi_square_threshold < 0.0,
		              std::logic_error,
		              "Invalid chi_square_threshold value");
    }*/
  
    //=======================================================
    //############# Electron clustering config ##############
    //=======================================================
  
    if (config_.has_key("clustering.max_distance")) {
      this->clustering.max_distance =
        config_.fetch_real_with_explicit_dimension("clustering.max_distance", 
                                                   "length");
      DT_THROW_IF(this->clustering.max_distance < 44.0 * CLHEP::mm,
		              std::logic_error,
		              "Invalid clustering.max_distance value");
      DT_LOG_DEBUG(this->verbosity, "clustering.max_distance: " << this->clustering.max_distance << "mm");
    }

    if (config_.has_key("clustering.hit_association_distance")) {
      this->clustering.hit_association_distance =
        config_.fetch_real_with_explicit_dimension("clustering.hit_association_distance",
                                                   "length");
      DT_LOG_DEBUG(this->verbosity, "clustering.hit_association_distance: " 
                    << this->clustering.hit_association_distance << "mm");
    }
    
    if (config_.has_flag("clustering.sinograms.save_sinograms")) {
      this->clustering.save_sinograms = true;
      DT_LOG_DEBUG(this->verbosity, "Clustering - saving sinograms");
    }
    
    if (config_.has_key("clustering.sinograms.iterations")) {
      this->clustering.iterations =
        config_.fetch_integer_scalar("clustering.sinograms.iterations");
      DT_LOG_DEBUG(this->verbosity, "clustering.sinograms.iterations: " 
                                     << this->clustering.iterations);
    }
    
    if (config_.has_key("clustering.sinograms.zoom_factor")) {
      this->clustering.zoom_factor =
        config_.fetch_dimensionless_real("clustering.sinograms.zoom_factor");
      DT_LOG_DEBUG(this->verbosity, "clustering.sinograms.zoom_factor: " 
                                    << this->clustering.zoom_factor);
    }
    
    if (config_.has_key("clustering.sinograms.resolution_phi")) {
      this->clustering.resolution_phi =
        config_.fetch_integer_scalar("clustering.sinograms.resolution_phi");
      DT_LOG_DEBUG(this->verbosity, "clustering.sinograms.resolution_phi: " 
                                    << this->clustering.resolution_phi << " bins"); 
    }
    
    if (config_.has_key("clustering.sinograms.resolution_r")) {
      this->clustering.resolution_r =
        config_.fetch_integer_scalar("clustering.sinograms.resolution_r"); 
      DT_LOG_DEBUG(this->verbosity, "clustering.sinograms.resolution_r: " 
                                    << this->clustering.resolution_r << " bins"); 
    }

    if (config_.has_key("clustering.sinograms.uncertainty")) {
      this->clustering.uncertainty =
        config_.fetch_real_with_explicit_dimension("clustering.sinograms.uncertainty",
                                                   "length");
      DT_LOG_DEBUG(this->verbosity, "clustering.sinograms.uncertainty: " 
                                    << this->clustering.uncertainty << "mm"); 
    }
    
    if (config_.has_key("clustering.sinograms.max_initial_precision_r")) {
      this->clustering.max_initial_precision_r =
        config_.fetch_real_with_explicit_dimension("clustering.sinograms.max_initial_precision_r",
                                                   "length");
      DT_LOG_DEBUG(this->verbosity, "clustering.sinograms.max_initial_precision_r: " 
                                    << this->clustering.max_initial_precision_r << "mm"); 

    }
    
    ///=======================================================
    ///############## alpha clustering config ################
    ///=======================================================
    
    if (config_.has_key("alphas.clustering_resolution_phi")) {
      this->alphas.clustering_resolution_phi =
        config_.fetch_integer_scalar("alphas.clustering_resolution_phi");
      DT_LOG_DEBUG(this->verbosity, "alphas.clustering_resolution_phi: " 
                                    << this->alphas.clustering_resolution_phi << " bins"); 
    }    
    
    if (config_.has_key("alphas.max_possible_drift_time")) {
      this->alphas.max_possible_drift_time =
        config_.fetch_real_with_explicit_dimension("alphas.max_possible_drift_time",
                                                   "time");
		  
		  DT_THROW_IF(this->alphas.max_possible_drift_time < 100.0 * CLHEP::ns,
		              std::logic_error,
		              "Invalid alphas.max_possible_drift_time value - too strict"); 
      DT_LOG_DEBUG(this->verbosity, "alphas.max_possible_drift_time: " 
                                    << this->alphas.max_possible_drift_time << "ns");   
    }
    
    if (config_.has_key("alphas.min_possible_drift_time")) {
      this->alphas.min_possible_drift_time =
        config_.fetch_real_with_explicit_dimension("alphas.min_possible_drift_time",
                                                   "time");
      DT_THROW_IF(this->alphas.min_possible_drift_time < 0.0 * CLHEP::ns,
		              std::logic_error,
		              "Invalid alphas.min_possible_drift_time value - negative time");   
      DT_LOG_DEBUG(this->verbosity, "alphas.min_possible_drift_time: " 
                                    << this->alphas.min_possible_drift_time << "ns"); 
    }
    
    if (config_.has_key("alphas.time_step")) {
      this->alphas.time_step =
        config_.fetch_real_with_explicit_dimension("alphas.time_step",
                                                   "time");
      DT_THROW_IF(this->alphas.time_step < 10.0 * CLHEP::ns,
		              std::logic_error,
		              "Invalid alphas.time_step value - too small steps (too many steps)");   
      DT_LOG_DEBUG(this->verbosity, "alphas.time_step: " 
                                    << this->alphas.time_step << "ns"); 
    }
    
    if (config_.has_flag("alphas.sinograms.save_sinograms")) {
      this->alphas.save_sinograms = true;
      DT_LOG_DEBUG(this->verbosity, "Alpha clustering - saving sinograms"); 
    }
    
    if (config_.has_key("alphas.sinograms.iterations")) {
      this->alphas.iterations =
        config_.fetch_integer_scalar("alphas.sinograms.iterations");
      DT_LOG_DEBUG(this->verbosity, "alphas.sinograms.iterations: " 
                                    << this->alphas.iterations);
    }
    
    if (config_.has_key("alphas.sinograms.zoom_factor")) {
      this->alphas.zoom_factor =
        config_.fetch_dimensionless_real("alphas.sinograms.zoom_factor");
      DT_LOG_DEBUG(this->verbosity, "alphas.sinograms.zoom_factor: " 
                                    << this->alphas.zoom_factor);
    }
    
    if (config_.has_key("alphas.sinograms.phi_step")) {
      this->alphas.phi_step =
        config_.fetch_real_with_explicit_dimension("alphas.sinograms.phi_step",
                                                   "angle");
        
      DT_THROW_IF(this->alphas.phi_step < 0.05 * CLHEP::deg,
                 std::logic_error,
		            "Invalid alphas.sinograms.phi_step value - too small");
      DT_LOG_DEBUG(this->verbosity, "alphas.sinograms.phi_step: " 
                                    << this->alphas.phi_step / CLHEP::deg << " degrees");
    }
    
    if (config_.has_key("alphas.sinograms.delta_r")) {
      this->alphas.delta_r =
        config_.fetch_real_with_explicit_dimension("alphas.sinograms.delta_r",
                                                   "length");
      DT_THROW_IF(this->alphas.delta_r < 0.0 * CLHEP::mm,
		              std::logic_error,
		              "Invalid alphas.sinograms.delta_r value");
      DT_LOG_DEBUG(this->verbosity, "alphas.sinograms.delta_r: " 
                                    << this->alphas.delta_r << "mm");   
    }
    
    if (config_.has_key("alphas.sinograms.resolution_r")) {
      this->alphas.resolution_r =
        config_.fetch_integer_scalar("alphas.sinograms.resolution_r"); 
      DT_LOG_DEBUG(this->verbosity, "alphas.sinograms.resolution_r: " 
                                    << this->alphas.resolution_r << " bins");  
    }
    
    if (config_.has_key("alphas.sinograms.uncertainty")) {
      this->alphas.uncertainty =
        config_.fetch_real_with_explicit_dimension("alphas.sinograms.uncertainty",
                                                   "length");
      DT_THROW_IF(this->alphas.uncertainty < 0.0 * CLHEP::mm,
		              std::logic_error,
		              "Invalid alphas.sinograms.uncertainty value");   
      DT_LOG_DEBUG(this->verbosity, "alphas.sinograms.uncertainty: " 
                                    << this->alphas.uncertainty << "mm");  
    }
    
    ///=======================================================
    ///########## polyline reconstruction config #############
    ///=======================================================
    
    if (config_.has_key("polylines.max_vertical_distance")) {
      this->polylines.max_vertical_distance =
        config_.fetch_real_with_explicit_dimension("polylines.max_vertical_distance",
                                                   "length");
      DT_THROW_IF(this->polylines.max_vertical_distance < 0.0 * CLHEP::mm,
	                std::logic_error,
	                "Invalid polylines.max_vertical_distance value");
      DT_LOG_DEBUG(this->verbosity, "polylines.max_vertical_distance: " 
                                    << this->polylines.max_vertical_distance << "mm"); 
    }
    
    if (config_.has_key("polylines.max_extention_distance")) {
      this->polylines.max_extention_distance =
        config_.fetch_real_with_explicit_dimension("polylines.max_extention_distance",
                                                   "length");
      DT_THROW_IF(this->polylines.max_extention_distance < 0.0 * CLHEP::mm,
		              std::logic_error,
		              "Invalid polylines.max_extention_distance value");
      DT_LOG_DEBUG(this->verbosity, "polylines.max_extention_distance: " 
                                    << this->polylines.max_extention_distance << "mm"); 
    }
    
    if (config_.has_key("polylines.max_tracker_hits_distance")) {
      this->polylines.max_tracker_hits_distance =
        config_.fetch_real_with_explicit_dimension("polylines.max_tracker_hits_distance",
                                                   "length");
      DT_THROW_IF(this->polylines.max_tracker_hits_distance < 0.0 * CLHEP::mm,
		              std::logic_error,
		              "Invalid polylines.max_tracker_hits_distance value");
      DT_LOG_DEBUG(this->verbosity, "polylines.max_tracker_hits_distance: " 
                                    << this->polylines.max_tracker_hits_distance << "mm"); 
    }
    
    if (config_.has_key("polylines.max_kink_angle")) {
      this->polylines.max_kink_angle =
        config_.fetch_real_with_explicit_dimension("polylines.max_kink_angle",
                                                   "angle");
      DT_THROW_IF(this->polylines.max_kink_angle < 0.0 * CLHEP::deg 
                || this->polylines.max_kink_angle > 180.0 * CLHEP::deg,
		            std::logic_error,
		            "Invalid polylines.max_kink_angle value");
      DT_LOG_DEBUG(this->verbosity, "polylines.max_kink_angle: " 
                                    << this->polylines.max_kink_angle / CLHEP::deg << " degrees"); 
    }
    
    if (config_.has_key("polylines.max_trajectories_middlepoint_distance")) {
      this->polylines.max_trajectories_middlepoint_distance =
        config_.fetch_real_with_explicit_dimension("polylines.max_trajectories_middlepoint_distance",
                                                   "length");
      DT_THROW_IF(this->polylines.max_trajectories_middlepoint_distance < 0.0 * CLHEP::mm,
		            std::logic_error,
		            "Invalid polylines.max_trajectories_middlepoint_distance value");
      DT_LOG_DEBUG(this->verbosity, "polylines.max_trajectories_middlepoint_distance: " 
                                    << this->polylines.max_trajectories_middlepoint_distance << "mm"); 
    }
    
    if (config_.has_key("polylines.max_trajectory_endpoints_distance")) {
      this->polylines.max_trajectory_endpoints_distance =
        config_.fetch_real_with_explicit_dimension("polylines.max_trajectory_endpoints_distance",
                                                   "length");
      DT_THROW_IF(this->polylines.max_trajectory_endpoints_distance < 0.0 * CLHEP::mm,
		              std::logic_error,
		              "Invalid polylines.max_trajectory_endpoints_distance value");
      DT_LOG_DEBUG(this->verbosity, "polylines.max_trajectory_endpoints_distance: " 
                                    << this->polylines.max_trajectory_endpoints_distance << "mm"); 
    }
    
    if (config_.has_key("polylines.max_trajectory_connection_angle")) {
      this->polylines.max_trajectory_connection_angle =
        config_.fetch_real_with_explicit_dimension("polylines.max_trajectory_connection_angle",
                                                   "angle");
        
      DT_THROW_IF(this->polylines.max_trajectory_connection_angle < 0.0 * CLHEP::deg 
                || this->polylines.max_trajectory_connection_angle > 180.0 * CLHEP::deg,
		            std::logic_error,
		            "Invalid polylines.max_trajectory_connection_angle value");
      DT_LOG_DEBUG(this->verbosity, "polylines.max_trajectory_connection_angle: " 
                                    << this->polylines.max_trajectory_connection_angle / CLHEP::deg << " degrees");
    }
    
    if (config_.has_key("polylines.min_distance_from_foil")) {
      this->polylines.min_distance_from_foil =
        config_.fetch_real_with_explicit_dimension("polylines.min_distance_from_foil",
                                                   "length");
      DT_THROW_IF(this->polylines.min_distance_from_foil < 0.0 * CLHEP::mm,
		              std::logic_error,
		              "Invalid polylines.min_distance_from_foil value"); 
      DT_LOG_DEBUG(this->verbosity, "polylines.min_distance_from_foil: " 
                                    << this->polylines.min_distance_from_foil << "mm");  
    }
    
    if (config_.has_key("polylines.min_distance_from_main_walls")) {
      this->polylines.min_distance_from_main_walls =
        config_.fetch_real_with_explicit_dimension("polylines.min_distance_from_main_walls",
                                                   "length");
      DT_THROW_IF(this->polylines.min_distance_from_main_walls < 0.0 * CLHEP::mm,
		              std::logic_error,
		              "Invalid polylines.min_distance_from_main_walls value"); 
      DT_LOG_DEBUG(this->verbosity, "polylines.min_distance_from_main_walls: " 
                                    << this->polylines.min_distance_from_main_walls << "mm");   
    }
    
    if (config_.has_key("polylines.min_distance_from_X_walls")) {
      this->polylines.min_distance_from_X_walls =
        config_.fetch_real_with_explicit_dimension("polylines.min_distance_from_X_walls",
                                                   "length");
      DT_THROW_IF(this->polylines.min_distance_from_X_walls < 0.0 * CLHEP::mm,
		              std::logic_error,
		              "Invalid polylines.min_distance_from_X_walls value");  
      DT_LOG_DEBUG(this->verbosity, "polylines.min_distance_from_X_walls: " 
                                    << this->polylines.min_distance_from_X_walls << "mm");  
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
    return _config_.electron_mode != ElectronRecMode::undefined;
  }
  
  void Algos::initialize(const CimrmanAlgoConfig & elrecconf_)
  {
    _config_ = elrecconf_;
    DT_THROW_IF(_config_.electron_mode == ElectronRecMode::undefined,
                std::logic_error,
                "Undefined reconstruction mode");

    if (_config_.visualization_2D || _config_.visualization_3D){
      _visu_ = std::make_unique<Visu>(_geom_);
    }
    
    // initializing prompt sinogram manager
    if(_config_.electron_mode != ElectronRecMode::undefined)
    {
      prompt_sinogram_manager.prompt = true;

      // fetching the setting from electron config
      Sinogram::Settings & set = prompt_sinogram_manager.get_settings();
      
      set.save_sinograms = _config_.clustering.save_sinograms;
      set.resolution_phi = _config_.clustering.resolution_phi;
      set.resolution_r   = _config_.clustering.resolution_r;
      set.iterations     = _config_.clustering.iterations;
      set.zoom_factor    = _config_.clustering.zoom_factor;
      set.sigma          = _config_.clustering.uncertainty;
    
      // reserving space for the array
      prompt_sinogram_manager.get_dual_space().reserve( set.resolution_phi * set.resolution_r );
      
      DT_LOG_DEBUG(_config_.verbosity, "Sinogram manager for prompt hits initialized");
    }

    // initializing delayed sinogram manager
    if(_config_.reconstruct_alphas)
    {
      delayed_sinogram_manager.prompt = false;
      
      // fetching the setting from alpha config
      Sinogram::Settings & set = delayed_sinogram_manager.get_settings();
      
      set.save_sinograms = _config_.alphas.save_sinograms;
      set.resolution_phi = M_PI / _config_.alphas.phi_step;
      set.resolution_r   = _config_.alphas.resolution_r;
      set.iterations     = _config_.alphas.iterations;
      set.zoom_factor    = _config_.alphas.zoom_factor;
      set.sigma          = _config_.alphas.uncertainty;
      
      // reserving space for the array
      delayed_sinogram_manager.get_dual_space().reserve( set.resolution_phi * set.resolution_r );
     
      DT_LOG_DEBUG(_config_.verbosity, "Sinogram manager for delayed hits initialized");
    }
    
    return;
  }

  void Algos::reset()
  {
    _event_ = nullptr;
    _config_ = CimrmanAlgoConfig();
    return;
  }
  
  void Algos::process(Event & event_)
  {
    set_event(event_);
    
    // step 1
    // preclustering to identify identify prompt and delayed preclusters and 
    // distant groups of tracker hits in case that provided TCD was not used
    if( _event_->get_preclusters().empty() )
    {
      DT_LOG_DEBUG(_config_.verbosity, "Step 1: Preclustering");
      precluster();
    }
    else
    {
      DT_LOG_DEBUG(_config_.verbosity, "Skipping Step 1: Preclustering - TCD bank provided");
    }
    
    // step 2
    // prompt tracker hit clustering 
    DT_LOG_DEBUG(_config_.verbosity, "Step 2: Prompt tracker hit clustering");
    Legendre_transform_clustering();
    
    // delayed tracker hit clustering
    if (_config_.reconstruct_alphas)
    {
      DT_LOG_DEBUG(_config_.verbosity, "Step 2: Delayed tracker hit clustering");
      alpha_clustering(); 
    }

    // step 3
    DT_LOG_DEBUG(_config_.verbosity, "Step 3: MLM line fitting + ambiguity checking and solving");
    make_MLM_fits();
    
    // step 4
    DT_LOG_DEBUG(_config_.verbosity, "Step 4: Linear fits are associated to tracker hits and combined into a precluster solutions");
    combine_into_precluster_solutions();
    
    // step 5    
    if (_config_.electron_mode == ElectronRecMode::polyline)
    {
      DT_LOG_DEBUG(_config_.verbosity, "Step 5: Kink finding and connecting into polyline trajectories (electron trajectories)");
      create_polyline_trajectories();
    }
    else if (_config_.electron_mode == ElectronRecMode::line)
    {
      DT_LOG_DEBUG(_config_.verbosity, "Step 5: Creating line trajectories (electron trajectories)");
      create_line_trajectories( electron );      
    }
    
    if (_config_.reconstruct_alphas)
    {
       DT_LOG_DEBUG(_config_.verbosity, "Step 5: Creating line trajectories (alpha trajectories)");
      create_line_trajectories( alpha );
    }
  
    // step 6
    if (_config_.electron_mode == ElectronRecMode::polyline)
    {
      DT_LOG_DEBUG(_config_.verbosity, "Step 6: Trajectory kinks and clustering refinements");
      refinements();
    }
    else
    {
      DT_LOG_DEBUG(_config_.verbosity, "Skipping step 6: not available for line trajectories");
    }

    // step 7
    DT_LOG_DEBUG(_config_.verbosity, "Step 7: Evaluating trajectories");
    evaluate_trajectories();

    // step 8
    DT_LOG_DEBUG(_config_.verbosity, "Step 8: Combining precluster solutions into all solutions");
    create_solutions();

    if (_visu_)
    {
    
      auto & preclusters = _event_->get_preclusters();
      int count_delayed = std::count_if(preclusters.begin(), preclusters.end(), [](const auto & precl){return precl->is_delayed();} );
      
      auto & solutions = _event_->get_solutions();
      bool long_traj = false;
      for(const auto & sol : solutions)
      {
        const auto & traj = sol->get_trajectories();
        auto large_traj = std::find_if(traj.begin(), traj.end(), [](const auto & tr){return tr->get_segments().size() > 6;});
        if( large_traj != traj.end() ) 
        {
          long_traj = true;
          break;
        }
      }
      
      
      if(_config_.visualization_2D /*&& long_traj && count_delayed > 0*/)
      {
        DT_LOG_DEBUG(_config_.verbosity, "Creating and saving 2D visualizations");
        _visu_->make_top_projection();
      }
      
      if(_config_.visualization_3D /*&& long_traj && count_delayed > 0*/)
      {
        DT_LOG_DEBUG(_config_.verbosity, "Creating and saving 3D visualizations");
        _visu_->build_event();
      }
    }
    return;
  }

  
//____________________________________________________________________________
//////////////////////////////////////////////////////////////////////////////
//// step 1: preclustering ///////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
  
  // master function of step 1
  void Algos::precluster()
  {
    // basic preclustering into delayed/prompt hits and into halves of tracker
    std::vector<TrackerHitHdl> prompt_hits_side0;
    std::vector<TrackerHitHdl> prompt_hits_side1;
    std::vector<TrackerHitHdl> delayed_hits_side0;
    std::vector<TrackerHitHdl> delayed_hits_side1;
    
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
    std::vector<std::vector<TrackerHitHdl>> temp_preclusters = separate_hits(prompt_hits_side0);
    DT_LOG_DEBUG(_config_.verbosity, "Prompt preclusters on side 0: " << temp_preclusters.size());
    for(auto & prec : temp_preclusters)
    {
      _event_->add_precluster(prec, true, 0);
    }
    temp_preclusters = separate_hits(prompt_hits_side1);
    DT_LOG_DEBUG(_config_.verbosity, "Prompt preclusters on side 1: " << temp_preclusters.size());
    for(auto & prec : temp_preclusters)
    {
       _event_->add_precluster(prec, true, 1);
    }
    temp_preclusters = separate_hits(delayed_hits_side0);    
    DT_LOG_DEBUG(_config_.verbosity, "Delayed preclusters on side 0: " << temp_preclusters.size());
    for(auto & prec : temp_preclusters)
    {
      _event_->add_precluster(prec, false, 0);
    }
    temp_preclusters = separate_hits(delayed_hits_side1);
    DT_LOG_DEBUG(_config_.verbosity, "Delayed preclusters on side 1: " << temp_preclusters.size());
    for(auto & prec : temp_preclusters)
    {
       _event_->add_precluster(prec, false, 1);
    }
    
    return;
  }
  
  void depth_first_search(int node, 
				   const std::vector<std::vector<int>> & adjacency, 
				   std::vector<bool> & visited, 
				   std::vector<int> & cluster_indices)
  {
    std::stack<int> stack;
    stack.push(node);
    visited[node] = true;

    while (!stack.empty())
    {
      int current = stack.top();
      stack.pop();
      cluster_indices.push_back(current);
      
      for(int neighbour : adjacency[current])
      {
        if(!visited[neighbour])
        {
          visited[neighbour] = true;
          stack.push(neighbour);
        }
      }
    }
  }

  // spatial clustering - separates a vector of tracker hits into spatialy distant subgroups (vectors)
  std::vector<std::vector<TrackerHitHdl>> Algos::separate_hits(const std::vector<TrackerHitHdl> & hits)
  {

    const double distance_treshold = _config_.clustering.max_distance;
    const size_t n = hits.size();

    // building adjacency graph
    // adjacency[i] = list of indices connected to hit i
    std::vector<std::vector<int>> adjacency(n); 
    
    for(size_t i = 0; i < n; ++i)
    {
      for(size_t j = i+1; j < n; ++j)
      {
        double dx = hits[i]->get_x() - hits[j]->get_x();
        double dy = hits[i]->get_y() - hits[j]->get_y();
        double distance = std::hypot(dx, dy);

        if (distance < distance_treshold)
        {
          adjacency[i].push_back(j);
          adjacency[j].push_back(i);
        }
    	}
    }

    // building clusters
    std::vector<bool> visited(n, false);
    std::vector<std::vector<TrackerHitHdl>> clusters;

    for(size_t i = 0; i < n; ++i)
    {
      if(visited[i]) continue;

      std::vector<int> cluster_indices;
      depth_first_search(i, adjacency, visited, cluster_indices);

      std::vector<TrackerHitHdl> cluster;
      for(int index : cluster_indices)
      {
        cluster.push_back(hits[index]);
      }

      clusters.push_back(std::move(cluster));
    }

    return clusters;
  }
  
  

//____________________________________________________________________________
//////////////////////////////////////////////////////////////////////////////
//// step 2: clustering //////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

  // master function of step 2
  void Algos::Legendre_transform_clustering()
  {
    for(auto & precluster : _event_->get_preclusters())
    {
      // clusterizes separatelly every prompt precluster
      if( precluster->is_delayed() ) continue;
     
      std::vector<TrackerHitHdl> & unclustered_hits = precluster->get_unclustered_tracker_hits(); 
      std::vector<ClusterHdl> & clusters = precluster->get_clusters();

      // recursively applies Legendre transform and spatial clustering to find linear clusters
      // creates clusters and removes its tracker hits from "unclustered_hits"
      clusterize_precluster(unclustered_hits, clusters);
    }
    return;
  }
  
  void Algos::clusterize_precluster(std::vector<TrackerHitHdl> & tracker_hits, std::vector<ClusterHdl> & clusters)
  {
    if( tracker_hits.size() < 3u ) return;
    
    const double association_distance = _config_.clustering.hit_association_distance;
    
    // Legendre transform on hits to find phi and r candidates
    double phi_estimate, r_estimate;
    find_cluster_Legendre( tracker_hits, phi_estimate, r_estimate );
    
    // identifying which tracker hits belong to phi, r line candidate
    std::vector<TrackerHitHdl> cluster_hits;
    separate_close_hits_to_line(tracker_hits, cluster_hits, phi_estimate, r_estimate, association_distance);
    
    // additional spatial separation to filter coincidental associations to the cluster          
    std::vector<std::vector<TrackerHitHdl>> sub_clusters = separate_hits(cluster_hits);
    
    // finding the biggest group
    int largest = 0;
    auto size = 0u;
    for(auto i = 0u; i < sub_clusters.size(); ++i)
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
    
    bool valid_Z = false;
    std::vector<TrackerHitHdl> invalid_Z_cluster;
    
    for(auto j = 0u; j < sub_clusters.size(); ++j)
    {
      if(j == largest && sub_clusters[j].size() > 2u)
      {        
        cluster_found = true;

        int no_valid_Z = std::count_if( sub_clusters[j].begin(), sub_clusters[j].end(), 
                                      [](const auto & hit){ return hit->has_valid_Z(); } );
                      
        if(no_valid_Z > 1)
        {
          valid_Z = true;
        
          // create and add new cluster to the precluster
          ClusterHdl cluster = std::make_shared<Cluster>( sub_clusters[j], phi_estimate, r_estimate );          
          clusters.push_back( cluster );
          DT_LOG_DEBUG(_config_.verbosity, "Cluster found with " << sub_clusters[j].size() << " hits and estimate (" << phi_estimate << " rad, " << r_estimate << " mm)");
        }
        else
        {
          invalid_Z_cluster = sub_clusters[j];
        }
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
      std::vector<std::vector<TrackerHitHdl>> sub_groups = separate_hits(tracker_hits); 
      tracker_hits.clear();
      for(auto & sub_group : sub_groups)
      {
        clusterize_precluster(sub_group, clusters);
        tracker_hits.insert(tracker_hits.end(), sub_group.begin(), sub_group.end());
      }
      if( not valid_Z )
      {
        tracker_hits.insert(tracker_hits.end(), invalid_Z_cluster.begin(), invalid_Z_cluster.end());
      }
    }
  }
  
  
  void Algos::separate_close_hits_to_line(std::vector<TrackerHitHdl> & hits,
                                          std::vector<TrackerHitHdl> & hits_separated,
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
  
  
  void Algos::find_cluster_Legendre(const std::vector<TrackerHitHdl> & hits, double & phi_estimate, double & r_estimate)
  {  
    // enclosing tracker hits in a smallest possible rectangle (min_x, max_y) x (min_y, max_y)
    std::pair<std::vector<TrackerHitHdl>::const_iterator,
              std::vector<TrackerHitHdl>::const_iterator> minmax_X, minmax_Y;
    minmax_X = std::minmax_element(hits.begin(), hits.end(), [](const auto & hit1, const auto & hit2)
    { 
      return hit1->get_x() <= hit2->get_x();
    });
    minmax_Y = std::minmax_element(hits.begin(), hits.end(), [](const auto & hit1, const auto & hit2)
    { 
      return hit1->get_y() <= hit2->get_y(); 
    });
    
    const double min_x = hits[minmax_X.first - hits.begin()]->get_x() - _geom_.tc_radius;
    const double max_x = hits[minmax_X.second - hits.begin()]->get_x() + _geom_.tc_radius;
    const double min_y = hits[minmax_Y.first - hits.begin()]->get_y() - _geom_.tc_radius;
    const double max_y = hits[minmax_Y.second - hits.begin()]->get_y() + _geom_.tc_radius;
    
    // center of the box enclosing the tracker hits
    const double center_X = (max_x + min_x) / 2.0;
    const double center_Y = (max_y + min_y) / 2.0;
    
    // size of region (delta_phi x delta_R) to be investigated
    double delta_phi = M_PI;
    double delta_r = std::hypot(max_x - min_x, max_y - min_y); 
    
    // setting up the sinogram
    prompt_sinogram_manager.reset_state();
    prompt_sinogram_manager.set_center( center_X, center_Y );    
    prompt_sinogram_manager.set_tracker_hits( hits );
    prompt_sinogram_manager.set_delta_r( delta_r );
    prompt_sinogram_manager.set_delta_phi( delta_phi );
    prompt_sinogram_manager.set_run( _event_->get_run_number() );
    prompt_sinogram_manager.set_event( _event_->get_event_number() );
    
    // processing the sinogram
    prompt_sinogram_manager.process();
    auto [peak_phi, peak_r, peak_value] = prompt_sinogram_manager.get_peak();
    phi_estimate = peak_phi;
    r_estimate = peak_r; 
    
    DT_LOG_DEBUG(_config_.verbosity, "Peak found: (" << phi_estimate << " rad, " << r_estimate << " mm)");
  }

//____________________________________________________________________________
//////////////////////////////////////////////////////////////////////////////
//// step 3: MLM line fitting + ambiguity checking and solving ///////////////
//////////////////////////////////////////////////////////////////////////////

  // master function of step 3
  void Algos::make_MLM_fits()
  {
    for(auto & precluster : _event_->get_preclusters())
    {
      // fits every prompt cluster
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
          //TODO: chi2 limiter could be implemented here to not propagate bad fits further
          // ,but would require a dedicated function for calculation (fits dont yet have chi2) 
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
    const double h = 0.000001; // has to be small but not too small - ideally roughly half the available decimal precision to maximize numerical stability and calculation precision
    for(int i = 0; i < 3; ++i)
    {
      double derivative = lik.log_likelihood_derivative(phi_0);
      phi_0 += -h * derivative / (lik.log_likelihood_derivative(phi_0 + h) - derivative);
    }
    const double sin_phi0 = std::sin(phi_0);
    const double cos_phi0 = std::cos(phi_0);
    const double tan_phi0 = std::tan(phi_0);

    double r_0 = (lik.Rr + lik.Rx * sin_phi0 - lik.Ry * cos_phi0) / lik.R;
    
    fit->phi = phi_0;
    fit->r = r_0;
    
   
    
    fit->a = tan_phi0;
    fit->b = -r_0 / cos_phi0;
    
    // calculating chi squared
    double min_likelihood = -lik.log_likelihood_value(phi_0);
    double no_of_measurements = lik.no_R + lik.no_Z;
    double chi_squared = min_likelihood / no_of_measurements;
    fit->chi_squared = chi_squared;
   
    // if at least 2 tracker hits have usable Z position, 3D ML fit is calculated
    if( lik.no_Z > 1 )
    {
    
      double denominator = (lik.Cov_Zyy * sin_phi0 * sin_phi0) 
                          + (2.0 * lik.Cov_Zxy * sin_phi0 * cos_phi0) 
                          + (lik.Cov_Zxx * cos_phi0 * cos_phi0);
      
      if( denominator != 0.0 )
      {
        double tan_theta = (lik.Cov_Zzy * sin_phi0 + lik.Cov_Zzx * cos_phi0) / denominator;
        double h_temp = (lik.Zz / lik.Z) - (lik.Zx * cos_phi0 + lik.Zy * sin_phi0) * tan_theta / lik.Z;
        
        fit->h = h_temp;
        fit->theta = std::atan(tan_theta);

        fit->c = tan_theta / cos_phi0;
        fit->d = h_temp - (r_0 * tan_phi0 * tan_theta);
      }
    }
    
    // in case of only 1 tracker hit Z position the fit goes horizontally at the height of the one tracker hit
    else if( lik.no_Z == 1 ) 
    {
      
      fit->h = lik.Zz / lik.Z;
      fit->theta = 0.0;
      
      fit->c = 0.0;
      fit->d = lik.Zz / lik.Z;
    } 
    
    DT_LOG_DEBUG(_config_.verbosity, "Linear fit created: (" << fit->phi << " rad, " 
                                                             << fit->r << " mm, " 
                                                             << fit->theta << " rad, "
                                                             << fit->h << " mm)" );
    DT_THROW_IF(lik.no_Z == 0, std::logic_error, "No Z infromation available for fitting!");
    
    return;
  }
  
  void Algos::detect_ambiguity_type(ClusterHdl cluster)
  {
    const std::vector<TrackerHitHdl> hits = cluster->get_tracker_hits();
    
    // detecting ambiguity type
    const double x0 = hits.front()->get_x();
    const double y0 = hits.front()->get_y();
    bool ambiguous;
    
    // type 1 == mirror image along line x = x0 
    ambiguous = true;
    for(auto i = 1u; i < hits.size(); ++i)
    {	
      if( x0 != hits[i]->get_x() )
      {
        ambiguous = false;
        break;
      }
    }
    if( ambiguous == true )
    {
      cluster->set_ambiguity_type( Ambiguity::Vertical );
      DT_LOG_DEBUG(_config_.verbosity, "Cluster ambiguity: Vertical");
      return;
    }
    
    // type 2 == mirror image along line y = y0 
    ambiguous = true;
    for(auto i = 1u; i < hits.size(); ++i)
    {	
      if( y0 != hits[i]->get_y() )
      {
        ambiguous = false;
        break;
      }
    }
    if( ambiguous == true )
    {
      cluster->set_ambiguity_type( Ambiguity::Horizontal );
      DT_LOG_DEBUG(_config_.verbosity, "Cluster ambiguity: Horizontal");
    	return;
    }
    
    // type 3 == mirror image along line y = x + (y0-x0) 
    ambiguous = true;
    for(auto i = 1u; i < hits.size(); ++i)
    {	
      if( y0 - hits[i]->get_y() != x0 - hits[i]->get_x() )
      {
        ambiguous = false;
        break;
      }
    }
    if( ambiguous == true )
    {
      cluster->set_ambiguity_type( Ambiguity::AntiDiagonal );
      DT_LOG_DEBUG(_config_.verbosity, "Cluster ambiguity: Anti-diagonal");
      return;
    }
    
    // type 4 == mirror image along line y = -x + (y0-x0) 
    ambiguous = true;
    for(auto i = 1u; i < hits.size(); ++i)
    {	
      if( y0 - hits[i]->get_y() != hits[i]->get_x() - x0 )
      {
        ambiguous = false;
       break;
      }
    }
    if( ambiguous == true )
    {
      cluster->set_ambiguity_type( Ambiguity::MainDiagonal );
      DT_LOG_DEBUG(_config_.verbosity, "Cluster ambiguity: Main-diagonal");
    	return;
    }

    cluster->set_ambiguity_type( Ambiguity::None );
    DT_LOG_DEBUG(_config_.verbosity, "Cluster ambiguity: None");
    return;
  }
  
  void Algos::create_mirror_fit(ClusterHdl cluster)
  {  
    if( cluster->get_ambiguity_type() == Ambiguity::None )
    {
      return;
    }
    if( cluster->get_linear_fits().size() != 1 ) 
    {
      return;
    }
    
    ConstLinearFitHdl fit = cluster->get_linear_fits().front();
    const std::vector<TrackerHitHdl> hits = cluster->get_tracker_hits();
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
    case Ambiguity::Vertical:
      mirror_a = -a0;
      mirror_b = 2.0 * a0 * x0 + b0;
      break;
    case Ambiguity::Horizontal:
      mirror_a = -a0;
      mirror_b = 2.0 * y0 - b0;
      break;
    case Ambiguity::AntiDiagonal:
      mirror_a = 1.0 / a0;
      mirror_b = y0 - x0 - (b0 / a0) + (y0 - x0) / a0;
      break;
    case Ambiguity::MainDiagonal:
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
    if( cluster->get_ambiguity_type() != Ambiguity::Vertical )
    {  
      mirror_fit->likelihood.Rr *= -1.0;
      mirror_fit->likelihood.Rrx *= -1.0;
      mirror_fit->likelihood.Rry *= -1.0;
      mirror_fit->likelihood.Cov_Rrx *= -1.0;
      mirror_fit->likelihood.Cov_Rry *= -1.0;
    }
    
    return;
  }
  
//____________________________________________________________________________
//////////////////////////////////////////////////////////////////////////////
//// step 4: Linear fits are associated to tracker hits //////////////////////
////         and combined into a precluster solutions   //////////////////////
//////////////////////////////////////////////////////////////////////////////
  
  // master function of step 4
  void Algos::combine_into_precluster_solutions()
  {
    // proccess is done separately for each precluster
    for(auto & precluster : _event_->get_preclusters())
    {
      // counting the number (N) of ambiguous clusters 
      int no_ambiguities = 0;
      for(auto & cluster : precluster->get_clusters())
      {
        if( cluster->get_ambiguity_type() != Ambiguity::None )
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
        precluster_solution->get_unclustered_tracker_hits() = precluster->get_unclustered_const_tracker_hits();
      }
      
      // filling the precluster solutions with different track combinations
      int N1 = 1;
      for(const auto & cluster : precluster->get_clusters())
      {
        // for each unambiguous cluster its single track is added to all solutions
        if(cluster->get_linear_fits().size() == 1u)
        {
          for(auto & precluster_solution : precluster_solutions)
          {  
            const auto hits = cluster->get_const_tracker_hits();
            TrackHdl track = std::make_shared<Track>(cluster->get_linear_fits().front(), hits);
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
          for(int i = 0; i < N1; ++i)
          {
            int N2 = no_precluster_solutions / (N1 * 2);
            for(int j = 0; j < N2; ++j)
            {  
              const auto hits = cluster->get_const_tracker_hits();
              TrackHdl track = std::make_shared<Track>(cluster->get_linear_fits()[0], hits);
              track->sort_associations();
              precluster_solutions[k]->add_track(track);
              ++k;
            }
            for(int j = 0; j < N2; ++j)
            {  
              const auto hits = cluster->get_const_tracker_hits();
              TrackHdl track = std::make_shared<Track>(cluster->get_linear_fits()[1], hits);
              track->sort_associations();
              precluster_solutions[k]->add_track(track);
              ++k;
            }
          }
          N1 *= 2;
        }
      }
    }
    return;
  }
  
//____________________________________________________________________________
//////////////////////////////////////////////////////////////////////////////
//// step 5: Kink finding and connecting into polyline trajectories ///////////
//////////////////////////////////////////////////////////////////////////////
  
  // master function of step 5 (straight track mode): Track -> non-kinked Trajectories
  void Algos::create_line_trajectories(TrajectoryType type)
  {    
    // proccess is done separately for each precluster
    for(auto & precluster : _event_->get_preclusters())
    {
      switch(type)
      {
        case electron:
          if(precluster->is_delayed()) continue;    
          break;
        case alpha:
          if(precluster->is_prompt()) continue;
          break;
        case both:
          break;
      }
      
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
  
  // master function of step 5 (polyline track mode): Track -> kinked Trajectories
  void Algos::create_polyline_trajectories()
  {
    // proccess is done separately for each precluster
    for(auto & precluster : _event_->get_preclusters())
    {
      if(precluster->is_delayed()) continue;
      
      // also separately for each precluster solution
      for(auto & precluster_solution : precluster->get_precluster_solutions())
      {
      
        std::vector<TrackHdl> unprocessed_tracks = precluster_solution->get_unprocessed_tracks(); 
        std::vector<std::vector<TrackHdl>> trajectory_candidates = find_polyline_candidates(unprocessed_tracks,
        																																				  precluster->get_side());
        
        for(auto & trajectory_candidate : trajectory_candidates)
        {
          // new Trajectory object for each found trajectory candidate
          TrajectoryHdl trajectory = std::make_shared<Trajectory>(trajectory_candidate);
          precluster_solution->get_trajectories().push_back(trajectory);
          
          // creating trajectory points
          if(trajectory->has_kink())
          {
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
              ++it;
            }
          }
        }       
      }
    }
    return;
  }
  
  std::vector<std::vector<TrackHdl>> Algos::find_polyline_candidates(std::vector<TrackHdl> & tracks, const int side) const
  {
    const double vertical_threshold           = _config_.polylines.max_vertical_distance;
    const double max_distance                 = _config_.polylines.max_tracker_hits_distance;
    const double min_distance_from_foil       = _config_.polylines.min_distance_from_foil;
    const double min_distance_from_main_walls = _config_.polylines.min_distance_from_main_walls;
    const double min_distance_from_X_walls    = _config_.polylines.min_distance_from_X_walls;
  
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
    
    // geo locators to access geometry information
    const snemo::geometry::gg_locator & ggLocator = _geom_.geo_loc->geigerLocator();    
    
    const snemo::geometry::calo_locator & caloLocator = _geom_.geo_loc->caloLocator();
    const double mainwall_x_cord = caloLocator.getXCoordOfWall(1);
    
    const snemo::geometry::xcalo_locator & xcaloLocator = _geom_.geo_loc->xcaloLocator();
    const double X_wall_y_cord = xcaloLocator.getYCoordOfWall(1, 1);
    
    // x coordinate of the outside border of tracker (side 1) 
    double tracker_x_max = ggLocator.getXCoordOfLayer(1, 8) + _geom_.tc_radius;
    tracker_x_max = std::min( tracker_x_max, mainwall_x_cord - min_distance_from_main_walls );
    
    // x coordinate of the inside border (source foil gap) of tracker (side 1)
    double tracker_x_min = ggLocator.getXCoordOfLayer(1, 0) - _geom_.tc_radius;
    tracker_x_min = std::max( min_distance_from_foil, tracker_x_min );

    // y coordinate of the border of tracker (side 1)
    double tracker_y_max = ggLocator.getYCoordOfRow(1, 112) + _geom_.tc_radius;
    tracker_y_max = std::min( tracker_y_max, X_wall_y_cord - min_distance_from_X_walls );
    
    // connection candidate (each pair (i, j) of tracks) are filtered based on several criteria 
    for(int i = 0; i < no_tracks; ++i)
    {
      TrackHdl track1 = tracks[i];
      for(int j = i+1; j < no_tracks; ++j)
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
        
        // detecting close hit on track1
        bool match1_found = false;
        for(auto & association : track1->get_associations())
        {
          if( distance_2D( kink_point_hdl, association.point) < max_distance )
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
          if( distance_2D( kink_point_hdl, association.point) < max_distance )
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
    for(int i = 0; i < no_tracks; ++i)
    {
      connection_counter[i] = std::accumulate(connections[i].begin(), connections[i].end(), 0);
    }
    
    // connecting the tracks = "trajectorizing"
    std::vector<bool> trajectorized(no_tracks, false);
    for(int i = 0; i < no_tracks; ++i)
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
        for(int count = 0; count < no_tracks-1; ++count)
        {
          for(int j = 0; j < no_tracks; ++j)
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
  
  
  // calculates the angle between the three point in radians 
  double calculate_angle(const Point & point1, const Point & point2, const Point & point3)
  {
    Point vec1 = {point1.x - point2.x,
                  point1.y - point2.y,
                  point1.z - point2.z };
    Point vec2 = {point3.x - point2.x,
                  point3.y - point2.y,
                  point3.z - point2.z };
    double angle = vec1.x * vec2.x + vec1.y * vec2.y + vec1.z * vec2.z;
    double norm1 = std::hypot(vec1.x, vec1.y, vec1.z);
    double norm2 = std::hypot(vec2.x, vec2.y, vec2.z);
    angle /= (norm1 * norm2);
    
    angle = std::max(-1.0, std::min(1.0, angle));
    return std::acos(angle);
  } 
  
  // takes a trajectory with multiple track segments and calculates the kink points and end points 
  void Algos::build_polyline_trajectory(PreclusterSolutionHdl & precluster_solution, TrajectoryHdl & trajectory)
  {
    DT_THROW_IF(trajectory->get_segments().size() < 2, std::logic_error, "No polyline trajectory for creating trajectory points");
    
    const double max_angle = _config_.polylines.max_kink_angle; 
    
    std::vector<PointHdl> & tr_points = trajectory->get_trajectory_points();
    std::vector<TrackHdl> & segments = trajectory->get_segments();

    // finding the kink points
    std::vector<PointHdl> kink_points;
    for(auto i = 0u; i < segments.size() - 1; ++i)
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
    for(auto i = 1u; i < segments.size() - 1; ++i)
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
        angle = M_PI - calculate_angle( *current_point, *kink_point, *next_kink_point );     
      }
      else 
      {
        angle = M_PI - calculate_angle( *current_point, *kink_point, *next_alternative_point );
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
    
    double angle = M_PI - calculate_angle( *current_point, *kink_point, *last_point );   
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
  
  
  
//____________________________________________________________________________
//////////////////////////////////////////////////////////////////////////////
//// step 6: Trajectory refinement ///////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
    
  // master function of step 6
  //  - clustering refinement, trajectory connecting, fit quality calculation
  //  - Potentially maximum likelihood refinement of the polyline fits...
  void Algos::refinements()
  { 
  	// refining electron tracks
    for(auto & precluster : _event_->get_preclusters())
    {
      if(precluster->is_prompt())
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
          
          // created connection changes the association points -> needs to be updated
          for(auto & trajectory : precluster_solution->get_trajectories())
          {
            if(trajectory->has_kink())
            {
              trajectory->update_segments();
            }
          }
        }
      }
    }  
    
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
    const double distance_threshold             = _config_.polylines.max_trajectories_middlepoint_distance;
    const double trajectory_connection_distance = _config_.polylines.max_trajectory_endpoints_distance;
    const double max_angle                      = _config_.polylines.max_trajectory_connection_angle;
    
    auto i = 0u;
    auto & trajectories = precluster_solution->get_trajectories(); 
    while(i < trajectories.size())
    {       
      bool connection_found = false;
      auto & trajectory_points1 = trajectories[i]->get_trajectory_points();
      for(auto j = i + 1; j < trajectories.size(); ++j)
      {
        auto & trajectory_points2 = trajectories[j]->get_trajectory_points();
        
        // case A: front - front
        PointHdl point1 = trajectory_points1.front();
        PointHdl point2 = trajectory_points2.front();
        double distance = distance_3D(point1, point2); 
        if( distance < trajectory_connection_distance )
        {
          Point middle_point = get_middle_point(*point1, *point2);
          double angle = M_PI - calculate_angle(*trajectory_points1[1],
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
          double angle = M_PI - calculate_angle(*trajectory_points1[1] ,
                                                 middle_point,
                                                 *trajectory_points2[trajectory_points2.size()-2]);
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
          double angle = M_PI - calculate_angle(*trajectory_points1[trajectory_points1.size()-2] ,
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
          
          double angle = M_PI - calculate_angle(*trajectory_points1[trajectory_points1.size()-2],
                                                  middle_point,
                                                 *trajectory_points2[trajectory_points2.size()-2]);
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
        ++i;
      }
    }
    return;
  }
  
  // checks whether all tracker hits are associated to correct segment inside a polyline trajectory
  // if not - associates it to the correct segment or removes it from the trajectory
  void Algos::remove_wrong_hits_associations(PreclusterSolutionHdl & precluster_solution, TrajectoryHdl & trajectory)
  {
    // tolerance primarily for comparing two doubles 
    const double numerical_tolerance = 0.01;
    
    std::vector<TrackHdl> & segments = trajectory->get_segments();
    std::vector<PointHdl> & traj_points = trajectory->get_trajectory_points();
    std::vector<ConstTrackerHitHdl> & unclustered_tracker_hits = precluster_solution->get_unclustered_tracker_hits();
    
    for(auto i = 0u; i < segments.size(); ++i)
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
        if( (t_start <= t_j + numerical_tolerance) && (t_j - numerical_tolerance <= t_end) ) 
        {
          ++j;
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
  
  // goes through unclustered tracker hits of a precluster solution and checks
  // if it can be associated to some segment of some trajectory
  void Algos::refine_clustering(PreclusterSolutionHdl & precluster_solution)
  {
    const double distance_threshold = _config_.clustering.hit_association_distance;
    const double max_extention_distance = _config_.polylines.max_extention_distance;
    
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

      for(auto i = 0u; i < segments.size(); ++i)
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
      const double hit_x = hit->get_x();
      const double hit_y = hit->get_y();
      const double hit_R = hit->get_R();
      bool clustered = false;
      for(auto j = 0u; j < trajectories.size(); ++j)
      {
        if( clustered ) break;
        
        auto & segments = trajectories[j]->get_segments();
        auto & traj_points = trajectories[j]->get_trajectory_points(); 
        for(auto k = 0u; k < segments.size(); ++k)
        {
          const double cos_phi = std::cos(segments[k]->get_phi());
          const double sin_phi = std::sin(segments[k]->get_phi());
          
      // 1. checking if the hit is close to the line (segment without bounds)
          //TODO possible to add vertical distance as well
          double distance = std::abs(segments[k]->get_r() - hit_x * sin_phi + hit_y * cos_phi) - hit_R;
          if( std::abs(distance) > distance_threshold )
          { 
                continue;
          }
                  
     // 2. checking if the hit is within bounds of the segment
          bool association_found = false;
          const double t = hit_x * cos_phi + hit_y * sin_phi;          
          const double t_start = t_starts[j][k];
          const double t_end = t_ends[j][k];
          
          auto position = 0u;
          std::vector<Association> & associations_k = segments[k]->get_associations();

          // hit might not be within the current segment ends, but might align well if we prolong the segment a little 
          // only the first and last segment can be prolonged! cases (k == 0) or (k == segments.size() - 1)
          bool trajectory_extention_needed = false;
          // new end_point in case of extention
          PointHdl end_point;             
       
        // A) tracker hit is in the middle of existing segment (no extentions)
          if(std::min(t_start, t_end) <= t &&
             std::max(t_start, t_end) >= t)
          {
            association_found = true;
            
            // finding where the hit belongs on the line (tracker_hits and tracker_hit_points are sorted vectors)
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
          
        // B) tracker hit can be added by extending the first segment 
          else if( k == 0 )
          {
            // checking if the potential extention is in the right direction 
            if( (t_end > t_start && t < t_start) ||
                (t_end < t_start && t >= t_start) )
            {
              const double distance = std::abs(t - t_start);
              if(distance < max_extention_distance)
              {
                association_found = true;
                trajectory_extention_needed = true;
                position = 0;
                end_point = traj_points.front();
              }
            }
          }
          
        // C) tracker hit can be added by extending the last segment 
          else if( k == segments.size() - 1 )
          {
            // checking if the potential extention is in the right direction
            if( (t_end > t_start && t >= t_end) ||
                (t_end < t_start && t <= t_end) )
            {
              const double distance = std::abs(t - t_end);
              if(distance < max_extention_distance)
              {
                association_found = true;
                trajectory_extention_needed = true;
                position = associations_k.size();
                end_point = traj_points.back();
              }
            }
          } 
          
          
      // 3. if good option found, the tracker hit is associated
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
  
//____________________________________________________________________________
//////////////////////////////////////////////////////////////////////////////
//// step 7: Trajectory refinement ///////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
    
  // master function of step 7
  //  - fit quality calculation
  void Algos::evaluate_trajectories()
  { 
  	// refining electron tracks
    for(auto & precluster : _event_->get_preclusters())
    {
      // evaluating electron and alpha tracks
      for(auto & precluster_solution : precluster->get_precluster_solutions())
      {
        // calculating all fit quality metrics
        for(auto & trajectory : precluster_solution->get_trajectories())
        {
          evaluate_trajectory(trajectory);
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
    
    // chi_squared is the sum of chi_squares of all segments divided by (number_of_measurements - DOF)
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
    

//____________________________________________________________________________
//////////////////////////////////////////////////////////////////////////////
//// step 8: Combining precluster solutions into all solutions ///////////////
//////////////////////////////////////////////////////////////////////////////
  
  bool compare_solutions(const SolutionHdl solution1, const SolutionHdl solution2);
  
  // master function of step 8
  void Algos::create_solutions()
  {
    // creating all unique combinations of all precluster solutions
    combine_precluster_solutions_into_solutions();
    
    // calculates some quality indicators of the solutions
    evaluate_solutions();
        
    // putting better solutions at the front (custom compare function) 
    auto & solutions = _event_->get_solutions();
    std::sort(solutions.begin(), solutions.end(), compare_solutions); 
    
    // removing suboptimal solutions from the collection
    if( _config_.keep_only_best_solutions )
    {
      remove_suboptimal_solutions();
    }
  }
  
  void Algos::combine_precluster_solutions_into_solutions()
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
  }
  
  bool compare_solutions(const SolutionHdl solution1, const SolutionHdl solution2)
  {
    //TODO not sure how exactly shoudl this work
    
    // 1. criteria: total number of identified linear segments
    // (more succesull fits -> better)
    unsigned int no_segments1 = solution1->get_num_of_segments();
    unsigned int no_segments2 = solution2->get_num_of_segments();
    if( no_segments1 != no_segments2 )
    {
      return (no_segments1 > no_segments2);
    }
    
    // 2. criteria: total number of trajectories 
    // (less trajectories -> more connected -> better)
    unsigned int no_trajectories1 = solution1->get_num_of_trajectories();
    unsigned int no_trajectories2 = solution2->get_num_of_trajectories();
    if( no_trajectories1 != no_trajectories2 )
    {
      return (no_trajectories1 < no_trajectories2);
    }
    
    // 3. criteria: total number of unclustered hits
    // (less unclustered hits -> better)
    unsigned int no_unclustered_hit1 = solution1->get_num_of_unclustered_hit();
    unsigned int no_unclustered_hit2 = solution2->get_num_of_unclustered_hit();
    if( no_unclustered_hit1 != no_unclustered_hit2 ) 
    {
      return (no_unclustered_hit1 < no_unclustered_hit2);
    }

    // 4. criteria: total summed chi2 of the fits
    // (lower total chi2 -> better)
    double chi2_summed1 = solution1->get_total_summed_chi2();
    double chi2_summed2 = solution2->get_total_summed_chi2();
    return (chi2_summed1 < chi2_summed2);
  
  }
  
  // calculating overall reconstruction qualities
  void Algos::evaluate_solutions()
  {
    auto & solutions = _event_->get_solutions();
    for(auto & solution : solutions)
    {
      unsigned int no_trajectories = 0;
      unsigned int no_segments = 0;
      unsigned int no_unclustered_hit = 0;
      double chi2_sum = 0.0;
      
      for(const auto & precluster_solution : solution->get_precluster_solutions())
      {
        const auto & trajectories = precluster_solution->get_trajectories();
        no_trajectories += trajectories.size();
        for(const auto & trajectory : trajectories)
        {
          no_segments += trajectory->get_segments().size();
          chi2_sum += trajectory->get_chi_squared();
        }
        no_unclustered_hit += precluster_solution->get_unclustered_tracker_hits().size();
      }
      
      solution->set_num_of_trajectories( no_trajectories );
      solution->set_num_of_segments( no_segments );
      solution->set_num_of_unclustered_hit( no_unclustered_hit );
      solution->set_total_summed_chi2( chi2_sum );
    }
  }
  
  // deletes solutions with not fully connected trajectories
  void Algos::remove_suboptimal_solutions()
  {
    auto & solutions = _event_->get_solutions();
    if( solutions.size() < 2 ) return;
    
    // solutions are order based on the number of segments (more -> better)
    unsigned int min_num_of_segments = solutions.front()->get_num_of_segments();
    auto first_suboptimal_sol = std::find_if(solutions.begin()+1, solutions.end(),
                        [min_num_of_segments](auto& sol){ 
                          return (sol->get_num_of_segments() < min_num_of_segments);
                        } );
    solutions.erase(first_suboptimal_sol, solutions.end());
    
   
    // solutions are order based on the number of trajectories (fewer the better)
    unsigned int max_num_of_trajectories = solutions.front()->get_num_of_trajectories();
    first_suboptimal_sol = std::find_if(solutions.begin()+1, solutions.end(),
                    [max_num_of_trajectories](auto& sol){ 
                      return (sol->get_num_of_trajectories() > max_num_of_trajectories);
                    } );
    solutions.erase(first_suboptimal_sol, solutions.end());
    
    // TODO: can be expanded to have more strict pruning
  }
  

//____________________________________________________________________________
//////////////////////////////////////////////////////////////////////////////
//// Alpha clustering  //////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

  // master function of alpha clustering 
  void Algos::alpha_clustering()
  {
    for(auto & precluster : _event_->get_preclusters())
    {
      // clusterizes separatelly every prompt precluster
      if( precluster->is_prompt() ) continue;
  
      // early exit for too few hits
      std::vector<TrackerHitHdl> & unclustered_hits = precluster->get_unclustered_tracker_hits(); 
      if(unclustered_hits.size() < 3)
      {
        DT_LOG_DEBUG(_config_.verbosity, "Delayed cluster not found - too few usable hits");
        continue;
      }
     
      // clustering - attempting to contraint phi range and remove misaligned hits
      std::vector<TrackerHitHdl> cluster_hits;
      auto [phi_min_ptr, phi_max_ptr] = find_delayed_cluster(unclustered_hits, cluster_hits); 
      
      // early exit if no large enough group is found
      if(cluster_hits.size() < 3) continue;
      
      // creating the delayed cluster
      DelayedClusterHdl delayed_cluster = std::make_shared<DelayedCluster>(cluster_hits, phi_min_ptr, phi_max_ptr);        
      precluster->get_clusters().push_back( delayed_cluster );

      // track estimation
      delayed_cluster->find_center();
      find_reference_time_bounds(delayed_cluster);
      estimate_delayed_track(delayed_cluster);
    }
  }

  
  // "Hough transform clustering" - Legendre transform clustering alternative for hits without drift radii
  // transforms the triggered cell (square) into space of all lines (phi, r) crossing this square
  // algo returns angular range [phi_min, phi_max] such that line (phi, r) from this range crosses the most cells geometrically possible   
  std::pair<double, double> Algos::find_delayed_cluster(std::vector<TrackerHitHdl> & hits, std::vector<TrackerHitHdl> & cluster_hits)
  {
    //number of different values of phi among which the cluster is being searched for
    const uint32_t bins_phi = _config_.alphas.clustering_resolution_phi;
    
    double phi_min_out = -M_PI / 2.0;
    double phi_max_out = M_PI / 2.0;
    if( hits.size() < 3 ) 
    {
      DT_LOG_DEBUG(_config_.verbosity, "cluster not found - too few usable hits");
      return {phi_min_out, phi_max_out};
    }
    
    const double tc_radius = _geom_.tc_radius;
        
    int max_count[bins_phi] = {0};
    double argmax_R[bins_phi];
    int global_max = 0;
    for(int i = 0; i < bins_phi; ++i)
    {
      std::vector<double> boundaries;
      std::vector<int> hit_count;
      
      double phi = i * M_PI / double(bins_phi);             
      double sin_phi = std::sin(phi);
      double cos_phi = std::cos(phi);
      double boundary_down, boundary_up;
              
      boundaries.push_back(-10000);
      hit_count.push_back(0);
      boundaries.push_back(10000);
              
      for(const auto & hit : hits)
      {
        double x = hit->get_x();
        double y = hit->get_y();
        if(i < bins_phi / 2)
        {
          boundary_down = (x - tc_radius) * sin_phi - (y + tc_radius) * cos_phi;
          boundary_up   = (x + tc_radius) * sin_phi - (y - tc_radius) * cos_phi;
        }
        else
        {         
          boundary_down = (x - tc_radius) * sin_phi - (y - tc_radius) * cos_phi;
          boundary_up   = (x + tc_radius) * sin_phi - (y + tc_radius) * cos_phi;  
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
              
      for(auto j = 0u; j < hit_count.size(); ++j)
      {
        if(hit_count.at(j) > max_count[i])
        {
          max_count[i] = hit_count.at(j);
          argmax_R[i] = (boundaries.at(j) + boundaries.at(j+1)) / 2.0;
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
      phi_bin_min = bins_phi - 1;
      while(max_count[phi_bin_min] == global_max)
      {
        phi_bin_min--;
      }

      phi_max = double(phi_bin_max) * M_PI / double(bins_phi);
      phi_min = double(phi_bin_min) * M_PI / double(bins_phi) - M_PI;
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

      phi_min = (phi_bin_min - 1.0) * M_PI / double(bins_phi);
      phi_max = (phi_bin_max) * M_PI / double(bins_phi);
    }

    double R_0 = argmax_R[phi_bin_max - 1];
    double phi_0 = double(phi_bin_max - 1) * M_PI / double(bins_phi);
        
    double sin_phi_0 = std::sin(phi_0);
    double cos_phi_0 = std::cos(phi_0);
    
    for(auto it = hits.begin(); it != hits.end(); )
    {
      auto hit = *it;
      
      double boundary_down, boundary_up;
      double x = hit->get_x();
      double y = hit->get_y();
      if(phi_bin_max-1 < bins_phi/2)
      {
        boundary_down = (x - tc_radius) * sin_phi_0 - (y + tc_radius) * cos_phi_0;
        boundary_up   = (x + tc_radius) * sin_phi_0 - (y - tc_radius) * cos_phi_0;
      }
      else
      {             
        boundary_down = (x - tc_radius) * sin_phi_0 - (y - tc_radius) * cos_phi_0;
        boundary_up   = (x + tc_radius) * sin_phi_0 - (y + tc_radius) * cos_phi_0;  
      }
      if(boundary_down <= R_0 && R_0 <= boundary_up)
      {
        cluster_hits.push_back(hit);
        hits.erase(it);
      }
      else
      {
        it++;
      }
    }
    
    int valid_Z = std::count_if(cluster_hits.begin(), cluster_hits.end(), 
                           [](const auto & hit){ return hit->has_valid_Z(); } );
    
    if(cluster_hits.size() >= 3 && valid_Z > 1)
    {
      phi_min_out = phi_min;
      phi_max_out = phi_max;
      DT_LOG_DEBUG(_config_.verbosity, "Delayed cluster found with " << cluster_hits.size() << 
                  " tracker hits and angular restriction: [" << phi_min_out << ", " << phi_max_out << "] rad");
    }
    else
    {
      hits.insert(hits.end(), cluster_hits.begin(), cluster_hits.end());
      cluster_hits.clear();
      DT_LOG_DEBUG(_config_.verbosity, "Delayed cluster not found - too few hits or not enough Z measurements" );
    }
    return {phi_min_out, phi_max_out};
  }
  
  // estimates the possible range for the missing reference time based on anodic times of provided hits
  void Algos::find_reference_time_bounds(DelayedClusterHdl delayed_cluster)
  {
    const double min_possible_drift_time = _config_.alphas.min_possible_drift_time; // in nanoseconds
    const double max_possible_drift_time = _config_.alphas.max_possible_drift_time; // in nanoseconds
    
    auto & tracker_hits = delayed_cluster->get_tracker_hits();
    
    std::sort(tracker_hits.begin(), tracker_hits.end(), 
              [](const auto& hit1, const auto& hit2){
                return hit1->get_delayed_time() < hit2->get_delayed_time(); });
   
    //std::for_each(tracker_hits.begin(), tracker_hits.end(), [](const auto& hit){ std::cout << hit->get_delayed_time() << std::endl; return;} );
    
    // optimal case
    double lower_bound = tracker_hits.back()->get_delayed_time() - max_possible_drift_time;           
    double upper_bound = tracker_hits.front()->get_delayed_time() - min_possible_drift_time;
    
    // TODO what to do to make this better??
    for(auto it = tracker_hits.rbegin()+1; it != tracker_hits.rend(); ++it)
    {
      if( lower_bound > upper_bound )
      {
        lower_bound = (*it)->get_delayed_time() - max_possible_drift_time; 
      }
      else
      {
        break;
      } 
    }
                  
    delayed_cluster->set_reference_time_min( lower_bound );
    delayed_cluster->set_reference_time_max( upper_bound );
    
    DT_LOG_DEBUG(_config_.verbosity, "Reference time bounds set: [" << delayed_cluster->get_reference_time_min()
                                       << " , " << delayed_cluster->get_reference_time_max() << "] ns");
    
    DT_THROW_IF(lower_bound > upper_bound, std::logic_error, "Invalid reference time bounds found!");
  }
  
  
  // find estimate of a linear fit for delayed cluster
  // finds optimal reference time and sets the drift radii - cluster can be further treated same as a prompt one
  void Algos::estimate_delayed_track(DelayedClusterHdl delayed_cluster)
  { 
    // constant same for all events
    const double initial_phi_step = _config_.alphas.phi_step;
    const double delta_r          = _config_.alphas.delta_r;
    const double time_step        = _config_.alphas.time_step;
  
    // event specific constants
    const double delta_phi = delayed_cluster->get_phi_max() - delayed_cluster->get_phi_min();
    const double peak_phi = (delayed_cluster->get_phi_min() + delayed_cluster->get_phi_max()) / 2.0;
    const double time_start = delayed_cluster->get_reference_time_min();
    const double time_end = delayed_cluster->get_reference_time_max();
    const double center_X = delayed_cluster->get_center_x();
    const double center_Y = delayed_cluster->get_center_y();
    
    auto & hits = delayed_cluster->get_tracker_hits(); 

    // setting up the sinograms - (constants with respect to one event)
    Sinogram::Settings & set = delayed_sinogram_manager.get_settings();
    set.resolution_phi = static_cast<uint32_t>(delta_phi / initial_phi_step);
    delayed_sinogram_manager.get_dual_space().resize(set.resolution_phi * set.resolution_r);
    delayed_sinogram_manager.set_center( center_X, center_Y );    
    delayed_sinogram_manager.set_tracker_hits( hits );
    delayed_sinogram_manager.set_run( _event_->get_run_number() );
    delayed_sinogram_manager.set_event( _event_->get_event_number() );
    
    int bins = (time_end - time_start) / time_step;
    
    //TH1F best_values = TH1F("best_values", "best_values", 
    //                   bins, time_start, time_end);
    
    // setting up the grid search through time range
    int time_bin = 0;
    double max = 0;
    double best_ref_time = time_start;
    double best_r = 0.0;
    double best_phi = 0.0;
    for(double time = time_end; time >= time_start; time -= time_step)
    {
      time_bin++;
      delayed_cluster->set_reference_time(time);
      delayed_cluster->update_drift_radii();
    
      // resetting the sinogram  
      delayed_sinogram_manager.reset_state( peak_phi, 0.0 );
      delayed_sinogram_manager.set_delta_phi( delta_phi );
      delayed_sinogram_manager.set_delta_r( delta_r );
      
      // processing the sinogram
      delayed_sinogram_manager.process();
      auto [phi_estimate, r_estimate, peak_value] = delayed_sinogram_manager.get_peak();
      
      if(peak_value > max)
      {
        max = peak_value;
        best_phi = phi_estimate;
        best_r = r_estimate;
        best_ref_time = time;
      }
      
      //best_values.AddBinContent(bins - time_bin + 1, peak_value);
    }
    
    /*
    {
      TCanvas canvas("alpha", "alpha", 1000, 800);
      canvas.cd();
      best_values.SetStats(0);
      best_values.Draw();
      const char* image_name = Form("Events_visu/alphas_Ev%d.png", _event_->get_event_number() );
      canvas.SaveAs(image_name);
    }
    */
    // saving the best found estimate 
    delayed_cluster->set_reference_time(best_ref_time);
    delayed_cluster->update_drift_radii();
    delayed_cluster->set_phi_estimate(best_phi);
    delayed_cluster->set_r_estimate(best_r);
    
    // reseting the sinogram after alpha track
    delayed_sinogram_manager.reset_sinogram_counter();
    
    DT_LOG_DEBUG(_config_.verbosity, "Delayed track estimate: phi = " << best_phi << " rad, r = " << best_r << " mm");
  
  }
  
  
  
  // different track estimation approach - more analytical (and experimental currently) - tries all 8 combinations and chooses
  std::vector<std::pair<double, double>> Algos::estimate_three_hit_tracks(const std::vector<TrackerHitHdl>& hits)
  {
    std::vector<std::pair<double, double>> output_estimates;

    if( hits.size() != 3 ) output_estimates;    
    
    const double log_likelihood_threshold = -10.0;
    const unsigned int phi_bins = 100;
    
    std::array<double, 8> maxima = {0};
    std::array<double, 8> maxima_phi = {0};
    std::array<double, 8> maxima_r = {0};
    for(int i = 0; i < 8; ++i)
    {
      std::vector<bool> signs = {bool(i & 1),
                                 bool(i & 2),
                                 bool(i & 4)};
      Likelihood likelihood(hits, signs);
      
      std::array<double, phi_bins> function_grid = {0};
      double phi_step = M_PI / static_cast<double>(phi_bins);
      double phi = -M_PI / 2.0;
      for(unsigned int j = 0; j < phi_bins; ++j)
      {
        function_grid[j] = likelihood.log_likelihood_value(phi);
        phi += phi_step;
      }
     
      auto it_max = std::max_element(function_grid.begin(), function_grid.end());
      maxima[i] = *it_max;
      
      double max_phi = (-M_PI / 2.0) + phi_step * static_cast<double>(it_max - function_grid.begin());
      maxima_phi[i] = max_phi;
      
      double max_r = (likelihood.Rx * std::sin(max_phi)) 
                   - (likelihood.Ry * std::cos(max_phi)) 
                   + likelihood.Rr;
      max_r /= likelihood.R; 
      maxima_r[i] = max_r;
    }
    
    for(int i = 0; i < 8; ++i)
    {
      if( maxima[i] > log_likelihood_threshold)
      {
        output_estimates.emplace_back(maxima_phi[i], maxima_r[i]);
        //std::cout << maxima[i] << ": (" << maxima_phi[i] << ", " << maxima_r[i] << ")" << std::endl; 
      }
    }
    
    return output_estimates;
  }
  
} //  end of namespace tkrec

