# research
In this repository, the main file to use is ray.py. This file takes in a data 
file and generates a grid of rays to create all-sky Aitoff projections of H I density, 
O VI density, and line-of-sight-velocity. It also creates a text file with the 
data from each ray. It also is able to make projections of H I density, O VI 
density, temperature, density, and metallicity. The file is set up to debug and time
the ray generation. It is also set to make plots of the ray data for a set fraction
of rays. 

The variables to set are:

pixels_per_dim      	  This determines the number of rays to generate in each direction 
		    	  (total number of rays = pixels_per_dim^2).

remove_first_N_kpc 	  This is the length of the ray to be removed, in order to exclude
		 	  data from the disk.


MakeProjections 	  Set to True to make projections of the dataset on each axis.

DeBug			  Set to True to debug the program while running it. 

WriteRays		  Set to True to save ray files. Needs to be True if ray files are
			  not already saved

fraction_ray_plots	  This determines the fraction of rays to make plots for. 

proj_scales 		  This determines the width of the projections being made.
			  
field_list 		  This sets the fields for which to extract data from the ray.
