// Cimrman headers
#include "geometry/Geometry.h"
#include "datamodel/Point.h"

// Standard headers
#include <cmath>

namespace cimrman {

  Geometry::Geometry(snemo::service_handle<snemo::geometry_svc> & geoManager, datatools::logger::priority _verbosity)
    : verbosity(_verbosity)
  {
    
    DT_LOG_DEBUG(verbosity, "Initializing geometry info...");
    
    const geomtools::mapping & geoMapping = geoManager->get_mapping();
    const geomtools::id_mgr & geoIdMgr = geoManager->get_id_mgr();
    std::string locator_plugin_name;
    const snemo::geometry::locator_plugin * geoLocator
      = snemo::geometry::getSNemoLocator(*geoManager.instance(), locator_plugin_name);
    const snemo::geometry::calo_locator & caloLocator = geoLocator->caloLocator();
    const snemo::geometry::xcalo_locator & xcaloLocator = geoLocator->xcaloLocator();
    const snemo::geometry::gveto_locator & gvetoLocator = geoLocator->gvetoLocator();
    const snemo::geometry::gg_locator & ggLocator = geoLocator->geigerLocator();

    this->has_Bi_source = false;
    std::vector<geomtools::geom_id> sourceCalibrationCarrierGids;
    uint32_t sourceCalibrationCarrierType = geomtools::geom_id::INVALID_TYPE; // (1110)
    if (geoIdMgr.has_category_info("source_calibration_carrier")) {
      this->has_Bi_source = true;
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
	      DT_LOG_DEBUG(verbosity, "x = " << x);
	      DT_LOG_DEBUG(verbosity, "y = " << y);
	      DT_LOG_DEBUG(verbosity, "z = " << z);
	      this->Bi_source_x = x;
	      this->Bi_source_y = y;
	      this->Bi_source_z = z;
      }
    }
    
    this->tc_radius = ggLocator.cellRadius() / CLHEP::mm;
    DT_LOG_DEBUG(verbosity, "tc_radius = " << this->tc_radius);

    this->mw_sizex = caloLocator.blockThickness() / CLHEP::mm; // 194.0 (31.0);
    this->mw_sizey = caloLocator.blockWidth() / CLHEP::mm; // 256.0;
    this->mw_sizez = caloLocator.blockHeight() / CLHEP::mm; // 256.0;

    this->gv_sizex = gvetoLocator.blockHeight() / CLHEP::mm; // 308.0;
    this->gv_sizey = gvetoLocator.blockWidth() / CLHEP::mm; // 310.0;
    this->gv_sizez = gvetoLocator.blockThickness() / CLHEP::mm; // 150.0;
 
    this->xw_sizex = xcaloLocator.blockWidth() / CLHEP::mm; // 200.0;
    this->xw_sizey = xcaloLocator.blockThickness() / CLHEP::mm; // 150.0;
    this->xw_sizez = xcaloLocator.blockHeight() / CLHEP::mm; // 208.5;

    DT_LOG_DEBUG(verbosity, "mw_sizex = " << this->mw_sizex);
    DT_LOG_DEBUG(verbosity, "mw_sizey = " << this->mw_sizey);
    DT_LOG_DEBUG(verbosity, "mw_sizez = " << this->mw_sizez);
    DT_LOG_DEBUG(verbosity, "gv_sizex = " << this->gv_sizex);
    DT_LOG_DEBUG(verbosity, "gv_sizey = " << this->gv_sizey);
    DT_LOG_DEBUG(verbosity, "gv_sizez = " << this->gv_sizez);
    DT_LOG_DEBUG(verbosity, "xw_sizex = " << this->xw_sizex);
    DT_LOG_DEBUG(verbosity, "xw_sizey = " << this->xw_sizey);
    DT_LOG_DEBUG(verbosity, "xw_sizez = " << this->xw_sizez);
    
    // TODO stored locators for cell and OM positions
    this->geo_loc = geoLocator;
  }
  
  double Geometry::distance_to_MW(const datamodel::Point & point /*, Side side*/) const
  {
    const snemo::geometry::calo_locator & caloLocator = geo_loc->caloLocator();
    const double mainwall_x_cord = caloLocator.getXCoordOfWall(1);
    return mainwall_x_cord - std::abs(point.x);
  }
  
  double Geometry::distance_to_XW(const datamodel::Point & point /*, Side side*/) const
  {
    const snemo::geometry::xcalo_locator & xcaloLocator = geo_loc->xcaloLocator();
    const double X_wall_y_cord = xcaloLocator.getYCoordOfWall(1, 1);
    return X_wall_y_cord - std::abs(point.y);
  }
  
  double Geometry::distance_to_SF(const datamodel::Point & point /*, Side side*/) const
  {
    return std::abs(point.x);
  }
  
  bool Geometry::is_inside_tracker(const datamodel::Point & point /*, Side side*/) const
  {
    // TODO there must be better way to check tracker volume
    const snemo::geometry::gg_locator & ggLocator = geo_loc->geigerLocator();
   
    // x coordinate of the outside border of tracker (side 1) 
    const double tracker_x_max = ggLocator.getXCoordOfLayer(1, 8) + tc_radius;    
    if( std::abs(point.x) > tracker_x_max ) return false;
    
    // x coordinate of the inside border (source foil gap) of tracker (side 1)
    const double tracker_x_min = ggLocator.getXCoordOfLayer(1, 0) - tc_radius;
    if( std::abs(point.x) < tracker_x_min ) return false;

    // y coordinate of the border of tracker (side 1)
    const double tracker_y_max = ggLocator.getYCoordOfRow(1, 112) + tc_radius;
    if( std::abs(point.y) > tracker_y_max ) return false;
    
    // TODO what about Z coordinate???
    
    // if all passes, returns true
    return true;
  }
  
  std::array<double, 2> Geometry::get_cell_position(const int SRL[3]) const
  {
    static const snemo::geometry::gg_locator & ggLocator = geo_loc->geigerLocator();
    std::array<double, 2> pos = {0};
    pos[0] = ggLocator.getXCoordOfLayer(SRL[0], SRL[2]);
    pos[1] = ggLocator.getYCoordOfRow(SRL[0], SRL[1]);
    return pos;
  }
  
  std::array<double, 3> Geometry::get_MW_OM_position(const int SWCR[4]) const
  {    
    static const snemo::geometry::calo_locator & caloLocator = geo_loc->caloLocator();
    std::array<double, 3> pos = {0};
    pos[0] = caloLocator.getXCoordOfWall(SWCR[0]);	
    pos[1] = caloLocator.getYCoordOfColumn(SWCR[0], SWCR[2]);
    pos[2] = caloLocator.getZCoordOfRow(SWCR[0], SWCR[3]);
    return pos;
  }

  std::array<double, 3> Geometry::get_XW_OM_position(const int SWCR[4]) const
  {    
    static const snemo::geometry::xcalo_locator & xcaloLocator = geo_loc->xcaloLocator();
    std::array<double, 3> pos = {0};
    pos[0] = xcaloLocator.getXCoordOfColumn(SWCR[0], SWCR[1], SWCR[2]);
    pos[1] = xcaloLocator.getYCoordOfWall(SWCR[0], SWCR[1]);
    pos[2] = xcaloLocator.getZCoordOfRow(SWCR[0], SWCR[1], SWCR[3]);
    return pos;
  }
  
  std::array<double, 3> Geometry::get_GV_OM_position(const int SWCR[4]) const
  {    
    static const snemo::geometry::gveto_locator & gvetoLocator = geo_loc->gvetoLocator();
    std::array<double, 3> pos = {0};
    pos[0] = gvetoLocator.getXCoordOfColumn(SWCR[0], SWCR[1], SWCR[2]);	
    pos[1] = gvetoLocator.getYCoordOfColumn(SWCR[0], SWCR[1], SWCR[2]);
    pos[2] = gvetoLocator.getZCoordOfWall(SWCR[0], SWCR[1]);
    return pos;
  }
  
} //  end of namespace cimrman

