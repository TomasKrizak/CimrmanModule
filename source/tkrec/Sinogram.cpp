// Cimrman headers
#include "tkrec/Sinogram.h"
#include "tkrec/TrackerHit.h"

// Standard headers
#include <iostream>

// Root headers
#include <TH2F.h>
#include <TCanvas.h>

// ClassImp(tkrec::Sinogram);

namespace tkrec {

  using namespace std;

    Sinogram::Settings & Sinogram::get_settings()
    {
      return settings;
    }
    
    // master function 
    void Sinogram::process()
    {
      while(iteration < settings.iterations)
      {
        // emptying the buffer
        reset_dual_space();
        
        // filling the sinogram one hit at a time
        for(auto & hit : tracker_hits)
        {
          if( hit->has_valid_R() ) 
          {
            add_hit_to_sinogram( hit );
          }
        }
        
        // creating png by converting the buffer into ROOT histogram
        if(settings.save_sinograms)
        {
          save_sinogram();
        }
        
        // finding new peak and updating the investigated range
        float* data = dual_space.data();
        int size = dual_space.size();
        uint32_t max_index = std::max_element( data, data + size ) - data;
        
        int peak_phi_index = max_index / settings.resolution_r;
        int peak_r_index = max_index % settings.resolution_r;
        
        peak_phi = index_to_phi( peak_phi_index );
        peak_r = index_to_r( peak_r_index );
        peak_value = dual_space[max_index];

        delta_phi = delta_phi / settings.zoom_factor; 
        delta_r = delta_r / settings.zoom_factor;      
        
        iteration++;
      }
      sinogram_ID++;
    }

    // converts phi index into the value of phi in the center of the phi column/row 
    inline double Sinogram::index_to_phi(uint32_t index_phi) const
    { 
      double min_edge = peak_phi - (delta_phi / 2.0);
      double step_size = delta_phi / static_cast<double>(settings.resolution_phi);
      return min_edge + ((static_cast<double>(index_phi) + 0.5) * step_size);
    }
    
    // converts r index into the value of r in the center of the r column/row 
    inline double Sinogram::index_to_r(uint32_t index_r) const
    {
      double min_edge = peak_r - (delta_r / 2.0);
      double step_size = delta_r / static_cast<double>(settings.resolution_r);
      return min_edge + ((static_cast<double>(index_r) + 0.5) * step_size);
    }
    
    // fast Padé approximation of exp(x) (order 5,5)
    // faster than std::exp with at least 5 digits of precision in range (-4.5, 0) = 3sigma range when used for Gauss 
    inline float fast_exp(const float x)
    {
      float nominator   = 1.0f + x*(0.5f + x*(1.0f/9.0f + x*(1.0f/72.0f + x*(1.0f/1008.0f + x*(1.0f/30240.0f)))));
      float denominator = 1.0f - x*(0.5f - x*(1.0f/9.0f - x*(1.0f/72.0f - x*(1.0f/1008.0f - x*(1.0f/30240.0f)))));
      return nominator / denominator;
    }
    
    // fast Padé approximation of exp(x) (order 3,3)
    // faster than std::exp, max error of 0.0085 in range (-4.5, 0) (3sigma range when used for Gauss)
    inline float faster_exp(const float x)
    {
      constexpr float c1 = 0.5f;
      constexpr float c2 = 1.0f / 10.0f;       
      constexpr float c3 = 1.0f / 120.0f;        

      float nominator   = 1.0f + x * (c1 + x * (c2 + x * c3));
      float denominator = 1.0f - x * (c1 - x * (c2 - x * c3));

      return nominator / denominator;
    }
    
    void Sinogram::add_hit_to_sinogram(TrackerHitHdl & tracker_hit)
    {
      double r_min = peak_r - (delta_r / 2.0);
      double r_max = peak_r + (delta_r / 2.0);
      double phi_min = peak_phi - (delta_phi / 2.0);
      double phi_max = peak_phi + (delta_phi / 2.0);
      
      double arr_sin[settings.resolution_phi];
      double arr_cos[settings.resolution_phi];
      for(int k = 0; k < settings.resolution_phi; ++k)
      {
        double phi = phi_min + ( (static_cast<double>(k) + 0.5) * delta_phi / static_cast<double>(settings.resolution_phi) );
        arr_sin[k] = std::sin(phi);
        arr_cos[k] = std::cos(phi);
      }
      
      for(int k = 0; k < settings.resolution_phi; ++k)
      {
        // r - legendre transform of the center of a circle (Hough transform)
        double r = (tracker_hit->get_x() - center_X) * arr_sin[k] - (tracker_hit->get_y() - center_Y) * arr_cos[k];
        double R_bin_width = delta_r / static_cast<double>(settings.resolution_r);
        
        bool overlap = false;
        for(int half = 0; half < 2; ++half)
        {    
          // mu - legendre transform of half circle (+R/-R)
          double mu = (r + (2.0 * half - 1.0) * tracker_hit->get_R());    

          // gauss is calculated only for -3 to 3 sigma region to cut time                            
          double r1 = mu - 3.0 * settings.sigma;
          double r2 = mu + 3.0 * settings.sigma;
          
          
          // if the 3sigma regions of the two halves of tracker hit overlap, we restrict the range to the middle (-+half bin for safety) 
          if(half == 0)
          {
            if(r2 > r) overlap = true;
      
            r2 = std::min(r2, r - 0.5 * R_bin_width); 
          }
          else
          {
            r1 = std::max(r1, r + 0.5 * R_bin_width);
          }

          // bin numbers coresponding to r1 and r2 values
          int bin1 = (static_cast<double>(settings.resolution_r) * (r1 - r_min) / delta_r);
          int bin2 = (static_cast<double>(settings.resolution_r) * (r2 - r_min) / delta_r) + 1;
          if(overlap && half == 0) bin2--;
            
          // if the 3 sigma borders (bin1 or bin2) are outside the investigate range of the histogtam, we restrict it to the border
          bin1 = std::max(0, bin1);
          bin2 = std::min(int(settings.resolution_r) - 1, bin2);

          // real values of r coresponding to each bin (lower and upper edge)
          double r_j1 = r_min + delta_r * static_cast<double>(bin1) / static_cast<double>(settings.resolution_r);
          double r_j2;
          
          // for large bins compared to the used sigma of gaussian bluring,
          // the function is integrated over the bin (in R direction)  
          if( R_bin_width > settings.sigma )
          {
            const double normalization = 1.0f / std::sqrt(2.0) * settings.sigma;
            for(int binj = bin1; binj < bin2 + 1; ++binj)
            {
              r_j2 = r_j1 + R_bin_width;
              
              // average probability density in a bin given by gauss distribution with mean in mu 
              float weight = ( std::erf( (r_j2 - mu) * normalization ) 
                              - std::erf( (r_j1 - mu) * normalization ) )
                            / (2.0 * R_bin_width);

              // result is 2D histogram of several sinusoid functions f(phi) in convolution with gauss in r                
              int globalBin = settings.resolution_r * k + binj;
              dual_space[globalBin] += weight;             
              
              r_j1 = r_j2;            
            }    
          }
          
          // for dense enough binning, the values are plotted without intergating
          // (saves A LOT of time - erf is expensive)
          else
          {
            for(int binj = bin1; binj < bin2 + 1; ++binj)
            {
              r_j2 = r_j1 + R_bin_width;
            
              // average probability density in a bin given by gauss distribution with mean in mu 
              double r_center = (r_j2 + r_j1) * 0.5f;
              float weight = (mu - r_center) / settings.sigma;
              weight = faster_exp( -0.5f * weight * weight ); // faster approximation of exp(x)
            
              // result is 2D histogram of several sinusoid functions f(phi) in convolution with gauss in r
              int globalBin = settings.resolution_r * k + binj;
              dual_space[globalBin] += weight;   
              
              r_j1 = r_j2;
            }                        
          }            
        }
      }
    }
    
    // returns the (peak_phi, peak_r) in global coordinates 
    // (r is return as a distance to the center of detector not the local sinogram center)
    std::tuple<double, double, double> Sinogram::get_peak()
    {    
      // result should be between -pi/2 and pi/2
      if( peak_phi > M_PI / 2.0 )
      {
        peak_phi -= M_PI;
        peak_r *= -1.0;
      }
      else if( peak_phi < -M_PI / 2.0 )
      {
        peak_phi += M_PI;
        peak_r *= -1.0;
      }
      
      return {peak_phi, peak_r + center_X * std::sin(peak_phi) - center_Y * std::cos(peak_phi), peak_value};
    }
    
    void Sinogram::set_run(int _run_number)
    {
      run_number = _run_number;
    }
    
    void Sinogram::set_event(int _event_number)
    {
      event_number = _event_number;
    }
    
    // saves a png of the current state of the sinogram
    void Sinogram::save_sinogram() const
    {
      double r_min = peak_r - (delta_r / 2.0);
      double r_max = peak_r + (delta_r / 2.0);
      double phi_min = peak_phi - (delta_phi / 2.0);
      double phi_max = peak_phi + (delta_phi / 2.0);
      
      TH2F sinograms("sinograms", "sinograms; phi; r",
         settings.resolution_phi, phi_min, phi_max,
         settings.resolution_r, r_min, r_max);
        
      float* sinograms_array = sinograms.GetArray();
      for(unsigned int i = 0; i < settings.resolution_phi; ++i)
      {
        for(unsigned int j = 0; j < settings.resolution_r; ++j)
        {
          float content = dual_space[(settings.resolution_r * i) + j];
          int globalBin = (settings.resolution_phi + 2) * (j + 1) + (i + 1);
          sinograms_array[globalBin] = content;          
        }
      }
             
      sinograms.SetEntries(settings.resolution_phi * settings.resolution_r);
      TCanvas canvas("sinograms", "sinograms", 1000, 800);
      canvas.cd();
      sinograms.SetStats(0);
      sinograms.SetContour(255);
      sinograms.Draw("COLZ");
      
      const char* image_name;
      if( prompt )
      {
        image_name = Form("Events_visu/prompt_sinogram-Rn%d_Ev%d_Tr%d_It%d.png",
                                    run_number, event_number, sinogram_ID, iteration );
      }
      else
      {
        image_name = Form("Events_visu/delayed_sinogram-Rn%d_Ev%d_time%d_It%d.png",
                                    run_number, event_number, sinogram_ID, iteration );
      }
      canvas.SaveAs(image_name);
      canvas.Close();    
    }
    
    // resets internal state 
    void Sinogram::reset_state(double _peak_phi, double _peak_r)
    {
      dual_space.resize(settings.resolution_phi * settings.resolution_r);
      iteration = 0u;
      peak_value = 0.0;
      peak_phi = _peak_phi;
      peak_r = _peak_r;
    }
    
    // resets internal buffer
    void Sinogram::reset_dual_space()
    {
      std::fill(dual_space.begin(), dual_space.end(), 0); 
    }
    
    void Sinogram::reset_sinogram_counter()
    {
      sinogram_ID = 0;
    }
    
    void Sinogram::set_delta_phi(double _delta_phi)
    {
      delta_phi = _delta_phi;
    }
    
    void Sinogram::set_delta_r(double _delta_r)
    {
      delta_r = _delta_r;
    }
    
    void Sinogram::set_center(double _center_X, double _center_Y)
    {
      center_X = _center_X;
      center_Y = _center_Y;
    }
    
    double Sinogram::get_delta_phi() const
    {
      return delta_phi;
    }
    
    double Sinogram::get_delta_r() const
    {
      return delta_r;
    }
    
    double Sinogram::get_center_X() const
    {
      return center_X;
    }
    
    double Sinogram::get_center_Y() const
    {
      return center_Y;
    }

    std::vector<float> & Sinogram::get_dual_space()
    {
      return dual_space;
    }
    
    const std::vector<float> & Sinogram::get_dual_space() const
    {
      return dual_space;
    }
    
    void Sinogram::set_tracker_hits(const std::vector<TrackerHitHdl> & _tracker_hits)
    {
      tracker_hits = _tracker_hits;
    }

} //  end of namespace tkrec
