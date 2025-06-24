# Cimrman Reconstruction Module

*"Algorithm so intuitive, you would think Jára Cimrman designed it."*

Cimrman is a track reconstruction module for SuperNEMO. It reconstructs charged particle trajectories using a combination of Legendre transform and maximum likelihood methods. Named in tribute to the great Czech thinker Jára Cimrman. The module takes **CD bank as an input**, reconstructs the trajectories within its own internal structure and **outputs both TCD and TTD data banks**, meaning it is responsible for clustering of tracker hits and for trajectory reconstruction at the same time. Currently it is able to produce **line and polyline trajectories**. It also detects ambiguous parts of tracker data and **provides all alternative possible solutions** (*tracker_trajectory_solution* and *tracker_clustering_solution*) in such cases. 

To obtain the final **PTD bank** you then need to apply **Charged Particle Tracking** module that extrapolates the verteces and creates the PTD bank.

You can find detailed description of the algorithm, its structure and the derivation of all methods in my [Master thesis](https://dspace.cvut.cz/handle/10467/123238).

## Module installation

In the future the module will be integrated in Falaise by default, but at the moment the module has to be installed independantly. For this reason, I temporarily provided the installation script install.sh

```
chmod 755 install.sh
./install.sh
```

## Module example usage

```
flreconstruct -i test_SD.brio -p reco.conf -o test_TTD.brio
```

## Module configuration

See the provided example configuration file *testing/test2/pipeline.conf*, which contains the default setting. 

1. **General configuration section:**
	* *visualization*: Creating and saving a png image and 3D object root file for each solution of each event.

     	**Warning:** intended for debugging purposes only, do not use on large datasets unless you want to generate enormous amount of files.
	* *save_sinograms*: saving the images of sinograms produced during Legendre trasnform based clustering.

     	**Warning:** intended for debugging purposes only, do not use on large datasets unless you want to generate enormous amount of files.
	* *force_default_sigma_r*: not using the drift radii uncertainties provided by Falaise and using default value instead
	* *default_sigma_r*: universal default uncertainty for drift radii
	* *chi_square_threshold*: currently not implemented. It is intended as a threshold for acceptable linear segment fits.
	
2. **Clustering section:**
	* *clustering_max_distance*: the maximum distance between two groups of tracker hits to be still considered as a part of one cluster 
		distance between tracker hits i and j = distance of their anodes:
		D(i,j) := sqrt( (x<sub>j</sub> - x<sub>i</sub>)<sup>2</sup> + (y<sub>j</sub> - y<sub>i</sub>)<sup>2</sup> )
		distance between two groups of tracker hits G<sub>1</sub> and G<sub>2</sub> = minimum distance between among all possible pairs of tracker cells from the two groups)
		D(G<sub>1</sub>, G<sub>2</sub>) := min{ D(i,j) | i from G<sub>1</sub>, j from G<sub>2</sub> }
	* *clustering_hit_association_distance*: clustering first finds a rough estimate for a line fit. Then associates tracker hits to it based on this distance threshold
	* *clustering_no_iterations*: number of iterations of the itterative zooming grid search algo used to find the maximum of the sinograms
	* *clustering_resolution_phi*: number of bins used for the grid search algo. starting range for phi is (0, phi)
	* *clustering_resolution_r*: number of bins used for the grid search algo. starting range for r is optimized to be minimal for each event.
	* *clustering_max_initial_precision_r*: If the initial range for r is low enough, number of used bins in r is reduced. In that case one bin = clustering_max_initial_precision_r
	* *clustering_zoom_factor*: each iteration of the search algo zooms in with this factor to achieve better precision
	* *clustering_uncertainty*: defines the blurring of the sinograms to compensate for errors and uncertainty caused by scattering
	
3. **Polyline reconstruction section:**
	* *polylines_max_extention_distance*: trajectory can be extended by associating additional tracker hit that satisfy the "clustering_hit_association_distance". This parameter defines the maximum possible elongation of the trajectory

	a) **Finding sharp kinks** - only vertical part of segments is changed
	* *polylines_max_vertical_distance*: maximum possible vertical distance of the two segments to be connected (computed in the point of their intersection in the horizontal projection)
	* *polylines_min_tracker_hits_distance*: tracker hits on both segments must be at least this close to the candidate kink to prevent the creation of distant "stranded" fake kinks
	* *polylines_max_kink_angle*: limits the 3D kink angle for the kink finder algo
	* *polylines_min_distance_from_foil*: prevents the kink to be created close to the foil
	
	b) **Finding tiny kinks** - all parameters of segments can be slightly changed			
	* *polylines_max_trajectories_middlepoint_distance*: maximum allowed horizontal shift of the segments to be connected into a polyline
   	* *polylines_max_trajectory_endpoints_distance*: maximum allowed 3D distance for two endpoints of two trajectories to be connected into polyline
	* *polylines_max_trajectory_connection_angle*: maximum alowed 3D angle for two trajectories to be connected into one polyline trajectory
	
	
