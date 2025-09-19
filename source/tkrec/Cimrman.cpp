// Interface from Falaise
#include "tkrec/Cimrman.h"

namespace tkrec {

  DPP_MODULE_REGISTRATION_IMPLEMENT(Cimrman, "Cimrman")

  // The working class which embeds private resources to
  // do the track reconstruction
  struct Cimrman::pimpl_type
  {
    Cimrman & tkrec; ///< Reference to the father recontruction module
    Geometry geom; ///< Geometry informations (ideally, should be extracted from the Falaise's geometry manager)
    Event event; ///< Working event to be reconstructed (working event data model)
    std::unique_ptr<Algos> palgo; ///< Reconstruction algorithms
    size_t eventCounter = 0;

    // Constructor
    pimpl_type(Cimrman & tkrec_);
    
  };

  Cimrman::pimpl_type::pimpl_type(Cimrman & tkrec_)
    : tkrec(tkrec_)
  {
    // Initialize the geometry informations from the father reconstruction module
    tkrec._init_geom_(geom);

    // Instantiate an 'algo' object:
    palgo = std::make_unique<Algos>(geom);
    
    // Initialize it from the configuration stored in the father reconstruction module
    palgo->initialize(tkrec._config_.recConfig);

    return;
  }

  const Cimrman::config_type & Cimrman::config() const
  {
    return _config_;
  }
  
  void Cimrman::_set_defaults_()
  {
    // Initialize the reference to the Falaise's geometry manager service
    _geoManager_ = snemo::service_handle<snemo::geometry_svc>{};
    return;
  }

  void Cimrman::_init_geom_(Geometry & geom_)
  {
    DT_LOG_DEBUG(_config_.verbosity, "Initializing geometry info...");
    const geomtools::mapping & geoMapping = _geoManager_->get_mapping();
    const geomtools::id_mgr & geoIdMgr = _geoManager_->get_id_mgr();
    std::string locator_plugin_name;
    const snemo::geometry::locator_plugin * geoLocator
      = snemo::geometry::getSNemoLocator(*_geoManager_.instance(), locator_plugin_name);
    const snemo::geometry::calo_locator & caloLocator = geoLocator->caloLocator();
    const snemo::geometry::xcalo_locator & xcaloLocator = geoLocator->xcaloLocator();
    const snemo::geometry::gveto_locator & gvetoLocator = geoLocator->gvetoLocator();
    const snemo::geometry::gg_locator & ggLocator = geoLocator->geigerLocator();

    geom_.has_Bi_source = false;
    std::vector<geomtools::geom_id> sourceCalibrationCarrierGids;
    uint32_t sourceCalibrationCarrierType = geomtools::geom_id::INVALID_TYPE; // (1110)
    if (geoIdMgr.has_category_info("source_calibration_carrier")) {
      geom_.has_Bi_source = true;
      sourceCalibrationCarrierType =
	geoIdMgr.get_category_info("source_calibration_carrier").get_type();
      geomtools::geom_id sourceCalibrarionCarrierGidPattern(sourceCalibrationCarrierType,
							    0, // module
							    geomtools::geom_id::ANY_ADDRESS,  // track
							    geomtools::geom_id::ANY_ADDRESS); // position
      geoMapping.compute_matching_geom_id(sourceCalibrarionCarrierGidPattern, sourceCalibrationCarrierGids);
      if (sourceCalibrationCarrierGids.size()) {
	      const geomtools::geom_info & srcCalibCarrierGinfo = geoMapping.get_geom_info(sourceCalibrationCarrierGids.front());
	      const geomtools::logical_volume & srcCalibCarrierLog = srcCalibCarrierGinfo.get_logical();
	      const geomtools::i_shape_3d & srcCalibCarrierShape = srcCalibCarrierLog.get_shape();    
	      const geomtools::box & srcCalibCarrierBox = dynamic_cast<const geomtools::box &>(srcCalibCarrierShape);
	      double x = srcCalibCarrierBox.get_z();
	      double y = srcCalibCarrierBox.get_y();
	      double z = srcCalibCarrierBox.get_x();
	      DT_LOG_DEBUG(_config_.verbosity, "x = " << x);
	      DT_LOG_DEBUG(_config_.verbosity, "y = " << y);
	      DT_LOG_DEBUG(_config_.verbosity, "z = " << z);
	      geom_.Bi_source_x = x;
	      geom_.Bi_source_y = y;
	      geom_.Bi_source_z = z;
      }
    }
    
    geom_.tc_radius = ggLocator.cellRadius() / CLHEP::mm;
    DT_LOG_DEBUG(_config_.verbosity, "tc_radius = " << geom_.tc_radius);

    geom_.mw_sizex = caloLocator.blockThickness() / CLHEP::mm; // 194.0 (31.0);
    geom_.mw_sizey = caloLocator.blockWidth() / CLHEP::mm; // 256.0;
    geom_.mw_sizez = caloLocator.blockHeight() / CLHEP::mm; // 256.0;

    geom_.gv_sizex = gvetoLocator.blockHeight() / CLHEP::mm; // 308.0;
    geom_.gv_sizey = gvetoLocator.blockWidth() / CLHEP::mm; // 310.0;
    geom_.gv_sizez = gvetoLocator.blockThickness() / CLHEP::mm; // 150.0;
 
    geom_.xw_sizex = xcaloLocator.blockWidth() / CLHEP::mm; // 200.0;
    geom_.xw_sizey = xcaloLocator.blockThickness() / CLHEP::mm; // 150.0;
    geom_.xw_sizez = xcaloLocator.blockHeight() / CLHEP::mm; // 208.5;

    DT_LOG_DEBUG(_config_.verbosity, "mw_sizex = " << geom_.mw_sizex);
    DT_LOG_DEBUG(_config_.verbosity, "mw_sizey = " << geom_.mw_sizey);
    DT_LOG_DEBUG(_config_.verbosity, "mw_sizez = " << geom_.mw_sizez);
    DT_LOG_DEBUG(_config_.verbosity, "gv_sizex = " << geom_.gv_sizex);
    DT_LOG_DEBUG(_config_.verbosity, "gv_sizey = " << geom_.gv_sizey);
    DT_LOG_DEBUG(_config_.verbosity, "gv_sizez = " << geom_.gv_sizez);
    DT_LOG_DEBUG(_config_.verbosity, "xw_sizex = " << geom_.xw_sizex);
    DT_LOG_DEBUG(_config_.verbosity, "xw_sizey = " << geom_.xw_sizey);
    DT_LOG_DEBUG(_config_.verbosity, "xw_sizez = " << geom_.xw_sizez);
    
    // TODO stored locators for cell and OM positions
    geom_.geo_loc = geoLocator;
    
    return;
  }

  Cimrman::Cimrman()
    : dpp::base_module() 
  {
    _set_defaults_();
    return;
  }

  Cimrman::~Cimrman()
  { 
    if (this->is_initialized())
      {
        this->reset();
      }
    return;
  }

  void Cimrman::read_config(const datatools::properties& config_)
  {
    if (config_.has_key("verbosity")) {
      std::string verbosityLabel = config_.fetch_string("verbosity");
      auto parsedVerb = datatools::logger::get_priority(verbosityLabel);
      DT_THROW_IF(parsedVerb == datatools::logger::PRIO_UNDEFINED,
		  std::logic_error,
		  "Undefined verbosity label " << std::quoted(verbosityLabel));
      _config_.verbosity = parsedVerb;
    }

    if (config_.has_key("CD_label")) {
      _config_.CD_label = config_.fetch_string("CD_label");
    }
 
    if (config_.has_key("TCD_label")) {
      _config_.TCD_label = config_.fetch_string("TCD_label");
    }
 
    if (config_.has_key("TTD_label")) {
      _config_.TTD_label = config_.fetch_string("TTD_label");
    }
 
    // Extract properties with prefix 'eventrec.' :
    datatools::properties recParamConfig;
    config_.export_and_rename_starting_with(recParamConfig, "eventrec.", "");
    _config_.recConfig.parse(recParamConfig);
    return;
  }

  void Cimrman::initialize(const datatools::properties & config_,
				 datatools::service_manager & services_,
				 dpp::module_handle_dict_type &  /*moduleDict*/
				 ) 
  {
    read_config(config_);
    
    // Geometry manager :
    _geoManager_ = snemo::service_handle<snemo::geometry_svc>{services_};
    // Instantiate the PIMPL working material:
    _work_ = std::make_unique<pimpl_type>(*this);
 
    this->_set_initialized(true);
    return;
  }

  void Cimrman::reset() 
  {   
    DT_THROW_IF(!is_initialized(), std::logic_error,
		"Module '" << get_name() << "' is not initialized !");
    this->_set_initialized(false);
    // Destroy the PIMPL working material:
    _work_.reset();
    _set_defaults_();
    return;
  }


  dpp::base_module::process_status Cimrman::process(datatools::things & workItem) 
  {
    _work_->eventCounter++;
    DT_LOG_DEBUG(_config_.verbosity, "============ New event #" << _work_->eventCounter);
    _populate_working_event_(workItem);
    _work_->palgo->process(_work_->event);	
	
    namespace snedm = snemo::datamodel;

    // Fill TKcluster data into TCD bank
    const auto & falaiseCDbank = workItem.get<snedm::calibrated_data>(_config_.CD_label);
	
    // Create or reset TCD bank
    auto & the_tracker_clustering_data
      = ::snedm::getOrAddToEvent<snedm::tracker_clustering_data>(_config_.TCD_label, workItem);
    the_tracker_clustering_data.clear();

    // Fill TKtrack data into TCD bank:
    _fill_TCD_bank_(falaiseCDbank, the_tracker_clustering_data);

    // Create or reset TTD bank
    auto & the_tracker_trajectory_data
      = ::snedm::getOrAddToEvent<snedm::tracker_trajectory_data>(_config_.TTD_label, workItem);
    the_tracker_trajectory_data.clear();

    // Fill TKtrack data into TTD bank
    _fill_TTD_bank_(the_tracker_clustering_data, the_tracker_trajectory_data);

    return falaise::processing::status::PROCESS_OK;
  }

  void Cimrman::_populate_working_event_(const datatools::things &workItem)
  {
    // Reset the working event
    _work_->event.reset();

    // Access to the event header
    const auto & header = workItem.get<snemo::datamodel::event_header>("EH");
    _work_->event.set_event_ids(header.get_id().get_run_number(),
				header.get_id().get_event_number());

    if(workItem.has(_config_.CD_label))
    {
	    DT_LOG_DEBUG(_config_.verbosity, "Has CD bank");
	    using namespace snemo::datamodel;

	    const auto & falaiseCDbank = workItem.get<calibrated_data>(_config_.CD_label);
	    DT_LOG_DEBUG(_config_.verbosity, "Nb calo hits = " << falaiseCDbank.calorimeter_hits().size());
	    for(const auto & calohit : falaiseCDbank.calorimeter_hits())
      {
        int SWCR[4] = {-1,-1,-1,-1};
        switch( calohit->get_geom_id().get_type() )
        {
          case 1302: 
	          SWCR[0] = calohit->get_geom_id().get(1);
	          SWCR[2] = calohit->get_geom_id().get(2);
	          SWCR[3] = calohit->get_geom_id().get(3);
	          break;
          case 1232:
	          SWCR[0] = calohit->get_geom_id().get(1);
	          SWCR[1] = calohit->get_geom_id().get(2);
	          SWCR[2] = calohit->get_geom_id().get(3);
	          SWCR[3] = calohit->get_geom_id().get(4);
	          break;
          case 1252:
	          SWCR[0] = calohit->get_geom_id().get(1);
	          SWCR[1] = calohit->get_geom_id().get(2);
	          SWCR[2] = calohit->get_geom_id().get(3);
	          break;
        }
        auto OMhitPtr = std::make_shared<OMHit>(SWCR);		
		    
        _work_->event.add_OM_hit( OMhitPtr );
      }

	    DT_LOG_DEBUG(_config_.verbosity, "Nb tracker hits = " << falaiseCDbank.tracker_hits().size());
	    
	    for (const auto & trhit : falaiseCDbank.tracker_hits() )
      {
        int SRL[3] = {trhit->get_side(), trhit->get_row(), trhit->get_layer()};
        auto hit = std::make_shared<TrackerHit>(SRL);
        hit->set_CDbank_tr_hit( trhit ); // TODO I cannot add as datatools::handle<const snemo::datamodel::calibrated_tracker_hit>
	      
        if(trhit->has_xy())
        {
          hit->set_x(trhit->get_x() / CLHEP::mm);
          hit->set_y(trhit->get_y() / CLHEP::mm);
        }
        
        if(trhit->is_prompt()) 
        {
          hit->set_as_prompt();
        }
        else if(trhit->is_delayed() && trhit->has_delayed_time())
        {
          hit->set_delayed_time(trhit->get_delayed_time());
        }
        
        double sigmaR = _config_.recConfig.default_sigma_r;
        if(not std::isnan(trhit->get_r()))
        {
          hit->set_R( trhit->get_r() / CLHEP::mm );
          hit->set_valid_R();
          if(not std::isnan(trhit->get_sigma_r()))
          {
            sigmaR = trhit->get_sigma_r() / CLHEP::mm;
          }
          // if (_config_.force_default_sigma_r) {
          //   sigmaR = _config_.recConfig.default_sigma_r;
          // }
          hit->set_sigma_R( sigmaR );
        }
        else
        {
          hit->set_R(datatools::invalid_real());        	
          hit->set_sigma_R(datatools::invalid_real());
        }
		    
        if(not std::isnan(trhit->get_z()))
        {
          hit->set_Z( trhit->get_z() / CLHEP::mm);
          hit->set_valid_Z();
          hit->set_sigma_Z( trhit->get_sigma_z() / CLHEP::mm );
        }
        else
        {
          hit->set_Z(datatools::invalid_real());
          hit->set_sigma_Z(datatools::invalid_real());
        }
        _work_->event.add_tracker_hit(hit);
      }

    }
    DT_LOG_DEBUG(_config_.verbosity, "Working event has been populated");
    
    return;
  }

  void Cimrman::_fill_TCD_bank_(const snemo::datamodel::calibrated_data & falaiseCDbank,
				        snemo::datamodel::tracker_clustering_data & the_tracker_clustering_data) const
  {
    namespace snedm = snemo::datamodel;
    // creating one clustering solution for each TK solution based on associated hits of individual trajectories
    // (1 falaise cluster = all associated tracker hits of 1 TK trajectory)
    std::vector<SolutionHdl> & solutions = _work_->event.get_solutions(); 
    for(auto i = 0u; i < solutions.size(); ++i)
    {
      SolutionHdl & solution = solutions[i]; 
      // creating empty clustering solution and adding to TCD bank
      auto htcs = datatools::make_handle<snedm::TrackerClusteringSolution>();
      the_tracker_clustering_data.append_solution(htcs, true);
      
      //snedm::tracker_clustering_solution & clustering_solution = the_tracker_clustering_data.get_default(); // TODO: get_default should be handled differently
      htcs->set_solution_id(the_tracker_clustering_data.size() - 1);
      
      auto & all_unclustered_hits = htcs->get_unclustered_hits();
      
      for(auto & precluster_solution : solution->get_precluster_solutions())
      {
      	// adding unclustered tracker hits into the solution
      	auto & unclustered_hits = precluster_solution->get_unclustered_tracker_hits();
      	for(auto & hit : unclustered_hits)
      	{
		    	all_unclustered_hits.push_back(hit->get_CDbank_tr_hit());
      	}
      	
        // creating one cluster for each TK trajectory
        std::vector<TrajectoryHdl> & trajectories = precluster_solution->get_trajectories();
        for(auto & trajectory : trajectories)
        {
          // creating empty cluster and adding to cluster collection
          snedm::TrackerClusterHdl cluster_handle = datatools::make_handle<snedm::tracker_cluster>();
          htcs->get_clusters().push_back(cluster_handle);
          cluster_handle->set_cluster_id(htcs->get_clusters().size() - 1);

          std::vector<TrackHdl> & traj_segments = trajectory->get_segments();
          for(auto & segment : traj_segments)
          {
            std::vector<Association> & associations = segment->get_associations();
            for(auto & association : associations)
            {
              cluster_handle->hits().push_back(association.tracker_hit->get_CDbank_tr_hit());
            }
          }
          cluster_handle->set_geom_id(geomtools::geom_id(1201, 0, cluster_handle->hits().front()->get_side()));
          
          // marking the cluster as prompt or delayed based on one of its hits 
          if(trajectory->get_segments().front()->get_associations().front().tracker_hit->is_prompt()) // TODO thats not nice
          {
		        cluster_handle->make_prompt();
          }
          else
          {
          	cluster_handle->make_delayed();
          }
        }
      }
    }
    
    if(not solutions.empty())
    {
    	the_tracker_clustering_data.set_default(0);
    }
    
    return;
  }

  void Cimrman::_fill_TTD_bank_(snemo::datamodel::tracker_clustering_data& the_tracker_clustering_data,
				      snemo::datamodel::tracker_trajectory_data& the_tracker_trajectory_data) const
  {	
    namespace snedm = snemo::datamodel;
    std::vector<SolutionHdl> & solutions = _work_->event.get_solutions(); 
    for(auto i = 0u; i < solutions.size(); ++i)
    {
      // clustering and tracking solutions are created with the same ordering as TK solutions
      
      // getting TK solution
      ConstSolutionHdl solution = solutions[i];
      
      // getting corresponding falaise clustering solution
      const snedm::TrackerClusteringSolutionHdl & cluster_solution = the_tracker_clustering_data.solutions()[i];

      // creating falaise trajectory solution 
      auto trajectory_solution = datatools::make_handle<snedm::tracker_trajectory_solution>();
      the_tracker_trajectory_data.add_solution(trajectory_solution); 
      trajectory_solution->set_solution_id(i);
      trajectory_solution->set_clustering_solution(cluster_solution);

      for(auto & precluster_solution : solution->get_precluster_solutions())
      {
        // creating one falaise trajectory for each TK trajectory
        std::vector<ConstTrajectoryHdl> trajectories = precluster_solution->get_trajectories();
        for(auto & trajectory : trajectories)
        {
          // creating new falaise tracker_trajectory and connecting to falaise cluster
          auto h_trajectory = datatools::make_handle<snedm::tracker_trajectory>();
          trajectory_solution->grab_trajectories().push_back(h_trajectory);
          
          int traj_ID = trajectory_solution->get_trajectories().size() - 1;
          h_trajectory->set_id(traj_ID);
          h_trajectory->set_cluster_handle(cluster_solution->get_clusters()[traj_ID]);

          //TODO: there must be a better way (like a link to Precluster from PreclusterSolution )
          int side = trajectory->get_segments().front()->get_associations().front().tracker_hit->get_SRL()[0]; 
          
          h_trajectory->set_geom_id(geomtools::geom_id(1201, 0, side));
          
          snedm::track_fit_infos & fit_infos = h_trajectory->get_fit_infos(); 
          fit_infos.set_chi2( trajectory->get_chi_squared() );
          int ndof = 1 + 3 * trajectory->get_segments().size();
          fit_infos.set_ndof( ndof ); 
                   
          //TODO there is no polyline option
          fit_infos.set_algo(snedm::track_fit_algo_type::TRACK_FIT_ALGO_LINE); 
          fit_infos.set_best(true);

          // creating and setting tracker_pattern as a line or polyline
          snedm::TrajectoryPatternHdl h_pattern;
          if(trajectory->has_kink())
          {
            auto polyline_pattern = new snedm::polyline_trajectory_pattern;
            h_pattern.reset(polyline_pattern);

            // creating the track points of polyline trajectory
            geomtools::polyline_3d& polyline = polyline_pattern->get_path();
            for (const auto& point : trajectory->get_trajectory_points())
            {
                geomtools::vector_3d geom_point = {
                  point->x * CLHEP::mm,
                  point->y * CLHEP::mm,
                  point->z * CLHEP::mm};
                polyline.add(geom_point);
            }
          }
          else
          {
            auto line_pattern = new snedm::line_trajectory_pattern;
            h_pattern.reset(line_pattern);

            // creating the track points of line trajectory
            ConstPointHdl start = trajectory->get_trajectory_points().front();
            ConstPointHdl end   = trajectory->get_trajectory_points().back();

            geomtools::line_3d & line_3d = line_pattern->get_segment();
            line_3d.set_first(
                  start->x * CLHEP::mm,
                  start->y * CLHEP::mm,
                  start->z * CLHEP::mm);
            line_3d.set_last(
                  end->x * CLHEP::mm,
                  end->y * CLHEP::mm,
                  end->z * CLHEP::mm);
          }
          h_trajectory->set_pattern_handle(h_pattern);
        }
      }
    }
    
    if(not solutions.empty())
    {
    	the_tracker_trajectory_data.set_default_solution(0);
    }
    
    return;
  }

} //  end of namespace tkrec
