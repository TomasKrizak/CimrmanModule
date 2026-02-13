// Cimrman headers
#include "visualization/Visu.h"
#include "geometry/Geometry.h"
#include "datamodel/Event.h"
#include "datamodel/TrackerHit.h"
#include "datamodel/OMHit.h"
#include "datamodel/Solution.h"
#include "datamodel/Point.h"
#include "datamodel/Association.h"
#include "datamodel/Trajectory.h"
#include "datamodel/Track.h"

// ROOT headers
#include "TCanvas.h"
#include "TColor.h"
#include "TEllipse.h"
#include "TH2D.h"
#include "TAttLine.h"
#include "TGLViewer.h"
#include "TGeoManager.h"
#include "TGeoVolume.h"
#include "TFile.h"
#include "TROOT.h"
#include "TPolyLine3D.h"
#include "TPolyLine.h"
#include "TBox.h"
#include "TLatex.h"
#include "TLine.h"
#include "TPoint.h"
#include "TStyle.h"
#include "TF1.h"
#include "TGraph.h"
#include "TPolyMarker3D.h"

// Standard headers
//#include <functional>
#include <memory>
#include <vector>
#include <array>


using namespace cimrman::datamodel;

namespace cimrman {

  Visu::Visu(const Geometry & geom_)
    : _geom_(geom_)
  {
    
    return;
  }

  bool Visu::has_event() const
  {
    return _event_ != nullptr;
  }

  void Visu::set_event(const Event & event_)
  {
    _event_ = &event_;
    return;
  }
 
  void Visu::make_top_projection(const uint32_t hits_option_,
                                 const uint32_t tracking_option_) const
  {
      for(auto i = 0u; i < _event_->get_solutions().size(); i++)
      {
          make_top_projection(i, hits_option_, tracking_option_);
      }
      return;
  }

  void Visu::make_top_projection(const uint32_t solution_id,
                                 const uint32_t hits_option_,
                                 const uint32_t tracking_option_) const
  {
    if (solution_id >= _event_->get_solutions().size()) return;
    ConstSolutionHdl solution = _event_->get_solutions()[solution_id];

    gROOT->SetBatch(true);
    TCanvas canvas("canvas","", 5800, 1600);     
    canvas.Range(-2900.0, -700.0, 2900.0, 900.0);
    TLatex title(-2600.0, 700.0,
                 Form("Run %d | Event %d | Solution %d",
                      _event_->get_run_number(),
                      _event_->get_event_number(),
                      solution_id));
    title.SetTextSize(0.07);
    title.Draw();

    // Drawing mainwall and mainwall calo hits
    std::vector<std::shared_ptr<TBox>> calorimeter_blocks;
    for(int om_side = 0; om_side < 2; om_side++) 
    {
      for(int om_column = 0; om_column < 20; om_column++) 
      {
        int SWCR[4] = {om_side, -1, om_column, 0};
        auto ohit = std::make_shared<OMHit>(SWCR);           
        std::array<double, 3> pos = _geom_.get_MW_OM_position( SWCR );
        
        ohit->set_x( pos[0] );	
        ohit->set_y( pos[1] );

        std::shared_ptr<TBox> calo = std::make_shared<TBox>(
                          ohit->get_y() - _geom_.mw_sizey / 2.0, 
                        -(ohit->get_x() + _geom_.mw_sizex / 2.0), 
                          ohit->get_y() + _geom_.mw_sizey / 2.0, 
                        -(ohit->get_x() - _geom_.mw_sizex / 2.0));

        calorimeter_blocks.push_back(calo);
                    
        bool is_hit = false;
        for(const auto & omHit : _event_->get_OM_hits())
        {
          if(om_side == omHit->get_SWCR()[0] && 
                  -1 == omHit->get_SWCR()[1] && 
           om_column == omHit->get_SWCR()[2] ) 
          {
            is_hit = true;
          }
        }
        if(is_hit)
        {
          calo->SetFillColor(kRed);       
        }
        calo->SetLineColor(kBlack); 
        calo->SetLineWidth(2);
        calo->Draw("same");
      }
    }
            
    // Drawing Xwall and Xwall calo hits
    for(int om_side = 0; om_side < 2; om_side++) 
    {
      for(int om_wall = 0; om_wall < 2; om_wall++) 
      {
        for(int om_column = 0; om_column < 2; om_column++) 
        {
          int SWCR[4] = {om_side, om_wall, om_column, 0};
          auto ohit = std::make_shared<OMHit>(SWCR); 
          std::array<double, 3> pos = _geom_.get_XW_OM_position( SWCR );
          
          ohit->set_x( pos[0] );	
          ohit->set_y( pos[1] );
                                                          
          std::shared_ptr<TBox> calo = std::make_shared<TBox>(
                            ohit->get_y() - _geom_.xw_sizey / 2.0, 
                          -(ohit->get_x() + _geom_.xw_sizex / 2.0), 
                            ohit->get_y() + _geom_.xw_sizey / 2.0, 
                          -(ohit->get_x() - _geom_.xw_sizex / 2.0));
                          
          calorimeter_blocks.push_back(calo);
                          
          bool is_hit = false;
          for(const auto & omHit : _event_->get_OM_hits())
          {
            if(om_side   == omHit->get_SWCR()[0] && 
               om_wall   == omHit->get_SWCR()[1] && 
               om_column == omHit->get_SWCR()[2] )  
            {
              is_hit = true;
            }
          }
          if(is_hit)
          {
            calo->SetFillColor(kRed);   
          }
          calo->SetLineColor(kBlack);     
          calo->SetLineWidth(2);
          calo->Draw("same");     
        } 
      }
    }
        
    // Drawing tracker and tracker hits
    std::vector<std::shared_ptr<TEllipse>> tracker_hits;
    for(int cell_num = 0; cell_num < 2034; cell_num++)
    {
      auto thit = std::make_shared<TrackerHit>(cell_num); 
      
      int* SRL = thit->get_SRL();
      std::array<double, 2> pos = _geom_.get_cell_position(SRL);
      thit->set_x( pos[0] );
      thit->set_y( pos[1] );
      
      double radius = _geom_.tc_radius;
      double sigma = std::numeric_limits<double>::quiet_NaN();
      bool is_hit = false;
      bool is_broken = false;
      bool is_associated = false;
      bool has_height = false;
      bool is_delayed = false;

      for(const auto & hit : _event_->get_tracker_hits())
      {
        if(cell_num == hit->get_cell_num())
        {
          is_hit = true;
          if( !hit->has_valid_R() )
          {
            is_broken = true;
          }
          else
          {
            radius = hit->get_R();
            sigma = hit->get_sigma_R();
          }     
          std::vector<ConstTrajectoryHdl> all_trajectories = solution->get_trajectories();
          for(auto & trajectory : all_trajectories)
          {
            const std::vector<Association> all_associations = trajectory->get_associations();
            for(const auto & association : all_associations)
            {
              if(association.tracker_hit == hit)
              {
                is_associated = true;
              }
            }
          }

          if( hit->has_valid_Z() )
          {
            has_height = true;
          }    
          
          if( !hit->is_prompt() )
          {
            is_delayed = true;
          }         
          break;
        }
      }

      // tracker hits
      if(is_hit)
      {
        if(is_broken)
        {
          auto tracker_cell = std::make_shared<TEllipse>(thit->get_y(),
                                        -thit->get_x(), radius, radius);
          tracker_hits.push_back(tracker_cell);
          tracker_cell->SetLineWidth(1);
          if(hits_option_ & show_hit_invalid)
          {
            tracker_cell->SetFillColor(kOrange);
            tracker_cell->Draw("same");
          }
        }
        else
        {
          auto tracker_cell = std::make_shared<TEllipse>(thit->get_y(),
                                      -thit->get_x(),
                                      radius + sigma,
                                      radius + sigma);
          tracker_hits.push_back(tracker_cell);
          tracker_cell->SetLineWidth(0);  
                                            
          if(hits_option_ & show_hit_used)
            tracker_cell->SetFillColor(kRed);       
              
          if(!has_height && (hits_option_ & show_hit_novertical))
            tracker_cell->SetFillColor(kMagenta);             
                     
          if(is_associated && (hits_option_ & show_hit_associated))
          {
            tracker_cell->SetFillColor(kGreen);         

            if(!has_height && (hits_option_ & show_hit_novertical))
              tracker_cell->SetFillColor(kTeal);                 
          }          
          
          if(is_delayed && (hits_option_ & show_hit_delayed))
          {
            tracker_cell->SetFillColor(kBlue);  
          }
          
          tracker_cell->Draw("same");
          if( radius - sigma > 0.0 )
          {
            auto tracker_cell_in = std::make_shared<TEllipse>(thit->get_y(),
                                                         -thit->get_x(),
                                                         radius - sigma,
                                                         radius - sigma);   
            tracker_hits.push_back(tracker_cell_in);
            tracker_cell_in->SetLineWidth(0);
            tracker_cell_in->Draw("same");
          }
        }
      }
      
      // non-triggered cells
      else
      {
        auto tracker_cell = std::make_shared<TEllipse>(thit->get_y(),
                                                    -thit->get_x(),
                                                    radius,
                                                    radius);
                                                    
        tracker_hits.push_back(tracker_cell);
        tracker_cell->SetLineWidth(1);
        tracker_cell->Draw("same");
      }
    }
        
    // Drawing Bi sources
    std::vector<std::shared_ptr<TBox>> sources;
    for(int column = 0; column < 6; column++)
    {
      auto Bi_source = std::make_shared<TBox>((column - 2.5) * _geom_.Bi_source_dist_y - _geom_.Bi_source_y / 2.0,
                                             -(-_geom_.Bi_source_x / 2.0),
                                             (column - 2.5) * _geom_.Bi_source_dist_y + _geom_.Bi_source_y / 2.0,
                                             -(_geom_.Bi_source_x / 2.0));
      sources.push_back(Bi_source);
              
      Bi_source->SetFillColor(kBlue); 
      Bi_source->SetLineWidth(2);
      Bi_source->Draw("same");
    }
        
    // Drawing trajectories
    std::vector<std::shared_ptr<TPolyLine>> polylines;
    if(tracking_option_ & show_tracking_trajectories)
    {
      std::vector<ConstTrajectoryHdl> trajectories = solution->get_trajectories();
      for (const auto & trajectory : trajectories)
      {
        std::shared_ptr<TPolyLine> traj = std::make_shared<TPolyLine>();
        polylines.push_back(traj);
        traj->SetLineColor(kOrange - 3);
        traj->SetLineWidth(4);
        std::vector<ConstPointHdl> trajectory_points = trajectory->get_trajectory_points();
        for(auto j = 0u; j < trajectory_points.size(); j++)
        {
          ConstPointHdl & point = trajectory_points[j];
          traj->SetPoint(j, point->y, -point->x);
        }
      traj->Draw("same");
      }
    }
        
    // Drawing avalanche origin points
    std::vector<std::shared_ptr<TGraph>> avalanche_origins;
    if(tracking_option_ & show_tracking_hit_avalanches)
    {
      std::vector<ConstTrajectoryHdl> trajectories = solution->get_trajectories();
      for (auto & trajectory : trajectories)
      {
        for(const auto & track : trajectory->get_segments())
        {
          std::vector<Association> associations = track->get_associations();
          if(associations.size() > 1)
          {
            std::shared_ptr<TGraph> graph = std::make_shared<TGraph>();
            avalanche_origins.push_back(graph);
            for (auto j = 0u; j < associations.size(); j++)
            {     
              ConstPointHdl point = associations[j].point;
              graph->SetPoint(j, point->y, -point->x);
            }
            graph->SetMarkerColor(kRed);
            graph->SetMarkerStyle(kFullCircle);
            graph->SetMarkerSize(1);
            graph->Draw("sameP");
          }
        }
      }
    }

    canvas.SaveAs(Form("./Events_visu/Run-%d_event-%d_solution-%d_2D.png",
                        _event_->get_run_number(),
                        _event_->get_event_number(),
                        solution_id));
        
    calorimeter_blocks.clear();
    tracker_hits.clear();
    sources.clear();
    //for(auto & track : lines) delete track;
    polylines.clear();
    avalanche_origins.clear();
    return;
  }
             
  void Visu::build_event(const uint32_t tracking_option_) const
  {
    for(auto i = 0u; i < _event_->get_solutions().size(); i++)
    {
        build_event(i, tracking_option_);
    }
    return;
  }

  void Visu::build_event(const uint32_t solution_id,
                         const uint32_t tracking_option_) const
  {     
    if(solution_id >= _event_->get_solutions().size()) return;
    
    ConstSolutionHdl solution = _event_->get_solutions()[solution_id];
    gROOT->SetBatch(true);
        
    TFile* file = new TFile(Form("./Events_visu/Run-%d_event-%d_solution-%d_3D.root",
                            _event_->get_run_number(),
                            _event_->get_event_number(),
                            solution_id), "RECREATE");

    TGeoManager *geom = new TGeoManager();
    TGeoMaterial *matVacuum = new TGeoMaterial("matVacuum", 0, 0, 0);    
    TGeoMedium *Vacuum = new TGeoMedium("Vacuum",1, matVacuum);
    TGeoVolume *top = gGeoManager->MakeBox("top", Vacuum, 1500, 1750, 2500);
    geom->SetTopVolume(top);

    int object_counter = 0;

    // Drawing calorimeter
    for(int omnum = 0; omnum < 712; omnum++) 
    {
      // Skipping calo hits
      bool is_hit = false;
      for(const auto & hit : _event_->get_OM_hits())
      {
        if(omnum == hit->get_OM_num())
        {
          is_hit = true;
        }
      }
      if(is_hit) continue;
              
      TGeoVolume * calo = nullptr;
      OMHit ohit = OMHit(omnum); 
      int* SWCR = ohit.get_SWCR();
      std::array<double, 3> pos;
         
      if(omnum < 520) // Mainwall
      {
        pos = _geom_.get_MW_OM_position( SWCR );
        calo = gGeoManager->MakeBox(Form("OM:%d.%d.%d.%d",
			      SWCR[0], SWCR[1], SWCR[2], SWCR[3]),
			    Vacuum,
			    _geom_.mw_sizex / 2.0,
			    _geom_.mw_sizey / 2.0,
			    _geom_.mw_sizez / 2.0);
      }
      else if(omnum < 648) // Xwall
      {
        pos = _geom_.get_XW_OM_position( SWCR );
        calo = gGeoManager->MakeBox(Form("OM:%d.%d.%d.%d",
			      SWCR[0], SWCR[1], SWCR[2], SWCR[3]),
			    Vacuum,
			    _geom_.xw_sizex / 2.0,
			    _geom_.xw_sizey / 2.0,
			    _geom_.xw_sizez / 2.0);
      }
      else // Gveto
      {
        pos = _geom_.get_GV_OM_position( SWCR );
        calo = gGeoManager->MakeBox(Form("OM:%d.%d.%d.%d",
			      SWCR[0], SWCR[1], SWCR[2], SWCR[3]),
		      Vacuum,
		      _geom_.gv_sizex / 2.0,
		      _geom_.gv_sizey / 2.0,
		      _geom_.gv_sizez / 2.0);
      }     
              
      TGeoHMatrix *trans = new TGeoHMatrix("Trans");
              
      trans->SetDx(pos[0]);
      trans->SetDy(pos[1]);
      trans->SetDz(pos[2]);
              
      calo->SetLineColor(kGray);      
              
      top->AddNode(calo, object_counter, trans);
      object_counter++;
    }

    // Drawing Bi calibration sources
    for(int row = 0; row < 7; row++)
    {    
      for(int column = 0; column < 6; column++)
      {
        TGeoVolume * Bi_source = gGeoManager->MakeBox("Bi_source",
					 Vacuum,
					 _geom_.Bi_source_x / 2.0,
					 _geom_.Bi_source_y / 2.0,
					 _geom_.Bi_source_z / 2.0);
            
        TGeoHMatrix *trans = new TGeoHMatrix("Trans");
        trans->SetDy( (double(column) - 2.5) * _geom_.Bi_source_dist_y );
        trans->SetDz( (double(row) - 3.0) * _geom_.Bi_source_dist_z );
                    
        Bi_source->SetLineColor(kCyan);
        Bi_source->SetLineWidth(3);
                    
        top->AddNode(Bi_source, object_counter, trans);
        object_counter++;
      }
    }

    // Adding calo hits
    for(const auto & hit : _event_->get_OM_hits())
    {
      TGeoVolume *box;
      const int* SWCR = hit->get_SWCR();
      std::array<double, 3> pos;
      if(hit->get_OM_num() >= 0 && hit->get_OM_num() < 520) // Mainwall
      {
        pos = _geom_.get_MW_OM_position( SWCR );
        box = gGeoManager->MakeBox(Form("OM:%d.%d.%d.%d",
				      SWCR[0], SWCR[1], SWCR[2], SWCR[3]),
		       Vacuum,
		       _geom_.mw_sizex / 2.0,
		       _geom_.mw_sizey / 2.0,
		       _geom_.mw_sizez / 2.0);
      }
      else if(hit->get_OM_num() < 648) // Xwall
      {
        pos = _geom_.get_XW_OM_position( SWCR );
        box = gGeoManager->MakeBox(Form("OM:%d.%d.%d.%d",
				      SWCR[0], SWCR[1], SWCR[2], SWCR[3]),
		       Vacuum,
		       _geom_.xw_sizex / 2.0,
		       _geom_.xw_sizey / 2.0,
		       _geom_.xw_sizez / 2.0);
      }
      else // Gveto
      {
        pos = _geom_.get_GV_OM_position( SWCR );
        box = gGeoManager->MakeBox(Form("OM:%d.%d.%d.%d",
				      SWCR[0], SWCR[1], SWCR[2], SWCR[3]),
		       Vacuum,
		       _geom_.gv_sizex / 2.0,
		       _geom_.gv_sizey / 2.0,
		       _geom_.gv_sizez / 2.0);
      }     
              
      TGeoHMatrix *trans = new TGeoHMatrix("Trans");
              
      trans->SetDx(pos[0]);
      trans->SetDy(pos[1]);
      trans->SetDz(pos[2]);
              
      box->SetLineColor(kRed);
      box->SetLineWidth(3);                                               

      top->AddNode(box, object_counter, trans);
      object_counter++;
    }

    // Adding tracker hits
    TPolyMarker3D* association_points = new TPolyMarker3D();
    int association_point_counter = 0;
    for(const auto & trhit : _event_->get_tracker_hits())
    {
      double radius = _geom_.tc_radius; // default value
      double sigma_R = trhit->get_default_sigma_R();
      double sigma_Z = trhit->get_default_sigma_R();
      if(trhit->has_valid_Z())
      {
        sigma_Z = trhit->get_sigma_Z();
      }
              
      if( trhit->has_valid_R() ) 
      {
        radius = trhit->get_R();
        sigma_R = trhit->get_sigma_R();
      }
            
      double radius_min = radius - sigma_R; 
      if( radius_min <= 0 )
      {
        radius_min = 0.0;
      }
    
      TGeoVolume *tracker_cell = geom->MakeTube(Form("cell:%d.%d.%d",
					       trhit->get_SRL()[0],
					       trhit->get_SRL()[1],
					       trhit->get_SRL()[2]),
					      Vacuum,
					      radius_min,
					      radius + sigma_R,
					      sigma_Z);
      TGeoHMatrix *trans = new TGeoHMatrix("Trans");
                      
      trans->SetDx( trhit->get_x() );
      trans->SetDy( trhit->get_y() );
      if( trhit->has_valid_Z() )
      {
        trans->SetDz( trhit->get_Z() );
      }
      else
      {
        trans->SetDz( 0.0 );
      }
              
      if( !trhit->has_valid_R() )
      {
        tracker_cell->SetLineColor(kOrange);
        tracker_cell->SetLineWidth(1);              
      }
      else if( !trhit->has_valid_Z() )
      {
        tracker_cell->SetLineColor(kMagenta);
        tracker_cell->SetLineWidth(1);
      }
      else
      {
        bool found = false;
        const auto& trajectories = solution->get_trajectories();
        for (const auto& trajectory : trajectories)
        {
          if(found) break;
          const auto& associations = trajectory->get_associations();
          for(const auto& association : associations)
          {
            if( association.tracker_hit == trhit )
            {
              found = true;
              const auto& association_point = association.point;
              association_points->SetPoint(association_point_counter++, 
                                           association_point->x, 
                                           association_point->y, 
                                           association_point->z);
              break;
            }
          }
        }
        if(found)
        {
          tracker_cell->SetLineColor(kGreen);
          tracker_cell->SetLineWidth(1);
        }
        else
        {
          tracker_cell->SetLineColor(kRed);
          tracker_cell->SetLineWidth(1);
        }
      }
              
      top->AddNode(tracker_cell, object_counter, trans);
      object_counter++;
    }
        
    // Close geometry and write to file
    geom->CloseGeometry(); 
    file->WriteObject(top, "demonstrator");

    association_points->SetMarkerStyle(20);    // marker style
    association_points->SetMarkerSize(1);    // size
    association_points->SetMarkerColor(kRed);
    file->WriteObject(association_points, "association_points");
    
    if(tracking_option_)
    {
      std::vector<ConstTrajectoryHdl> trajectories = solution->get_trajectories();
      for (auto i = 0u; i < trajectories.size(); ++i)
      {
        TPolyMarker3D* points = new TPolyMarker3D();

        std::vector<ConstPointHdl> trajectory_points = trajectories[i]->get_trajectory_points();
        TPolyLine3D* trajectory = new TPolyLine3D();        
        trajectory->SetLineColor(kOrange - 3);
        trajectory->SetLineWidth(4);
        for (auto j = 0u; j < trajectory_points.size(); ++j)
        {
          ConstPointHdl & point = trajectory_points[j];
          trajectory->SetPoint(j, point->x, point->y, point->z);
          
          points->SetPoint(j, point->x, point->y, point->z); 
        }
        points->SetMarkerStyle(20);    // marker style
        points->SetMarkerSize(1);    // size
        points->SetMarkerColor(kOrange);
        
        file->WriteObject( points, Form("trajectory_points-%d", i) );
        file->WriteObject( trajectory, Form("trajectory-%d", i) );
        delete trajectory;
        delete points;
      }     
    }

    // Close file and delete dynamically allocated objects
    file->Close();              
    delete file;        
    delete geom;        
    return;
  }

} // end of namespace cimrman
