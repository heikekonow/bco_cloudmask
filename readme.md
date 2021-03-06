# Cloudmask from radar observations at Barbados Cloud Observatory (BCO)

The [Barbados Cloud Observatory (BCO)](https://barbados.mpimet.mpg.de) is a meteorological observations site located at the eastern shore of Barbados. It is operated by the Max-Planck-Institute for Meteorology (MPI-M) in Hamburg, Germany since 2010. The main focus of observations at the BCO are trade wind cumuli and their surrounding atmosphere. To this end, a wide range of instruments are operated here. One of which is a cloud radar.

To make the classification of observed clouds by different paramters possible, this code has been developed to generate a cloud mask product from cloud radar observations. In the current state, only the two Ka-band radars are implemented. Adaption to the W-band radar is planned.

Processing is done in four steps:

   1. Concatenating data (functions/bco_cloudmask_concatData.m)
   2. Segmenting clouds (functions/bco_cloudmask_mask)
   3. Caculating cloud parameters (functions/bco_cloudmask_param)
   4. Saving data (functions/bco_cloudmask_save2netcdf)

Step 1 explicitly builds on the data structure of the BCO data at the MPI-M as it prepares available data for further processing. But the later steps can be applied to other kinds of radar data from different sources as well.

## 1. Concatenating data
Available measurement data from BCO radar observations are listed and analyzed. Since vertical resolutions changed during the BCO operations, only data with consistent vertical resolutions are concatenated. Every time, the vertical resolution changes, a new output product is generated.

To generate a data product, that has as little interruptions between files as possible, but still remains small enough to work with, yearly files are produced (as long as the vertical resolution remains constant).

To calculate cloud lengths from radar reflectivity measurements, wind observations are necessary. Therefore, wind measurements from surface meteorology data are read along with the radar data.

The concatenated radar data, together with surface meteorology, are saved in a temporary .mat file. The location for these temporary data files are specified in the beginning of functions/bco_cloudmask_concatData.m.


## 2. Segmenting clouds
Radar reflectivity data are converted to cloud masks where 0 denotes "no cloud" and 1 indicates "cloud". To avoid that short measurement interruptions are identified as beginning or ending of a cloud, this binary mask is cleaned up by applying morphological closing.

For closing of the cloud mask, a rectangular structuring element of 5x2 pixels is used. 5 vertical pixels correspond to 150 m. 2 horizontal pixels correspond to 20 s and with an estimated wind speed of 8 m/s, this are roughly 160 m cloud length. Although wind speed increases with height and thus 2 s measurement interruption would result in longer cloud gaps in higher altitudes. But this has not been incorporated in this code yet. It is assumed that this effect on the structuring is minor.

After morphological closing has been applied, cloud objects are segmented by applying connected components labeling ([Dillencourt et al., 1992](https://doi.org/10.1145/128749.128750)) to the cloud mask. For this, adjacent cloudy pixels are assigned the same label if they are connected via 8-connectivity. During this step, objects that comprise less than four pixels are discarded as clutter.


## 3. Calculating cloud parameters
In this step, parameters for each individual segmented cloud object (see previous step, 2. Segmenting clouds) are calculated. The calculated parameters are: cloud start time, cloud end time, cloud/chord length, cloud depth, cloud base height, cloud top height.

**Cloud start time** and **cloud end time** are the first and last time steps at which a cloud object passes over the radar.

Cloud length or **chord length** is calculated from duration (the time it took a cloud to pass over the radar) and wind speed. The wind speed *u(z)* at the height *z* is extrapolated from the 2 m reference height to cloud heights by applying the power law to the wind speed at reference height *u<sub>2m</sub>*:

![equation](https://latex.codecogs.com/gif.latex?u%28z%29%20%3D%20u_%7B2%20%5Cmathrm%7Bm%7D%7D%20%5Cleft%28%5Cfrac%7Bz%7D%7B2%5C%20%5Cmathrm%7Bm%7D%7D%5Cright%29%5E%7B0.11%7D).
<!-- \frac{u}{u_r} = \left(\frac{z}{z_r} \right )^{\alpha} -->
<!-- u(z) = u_{2 \mathrm{m}} \left(\frac{z}{2\ \mathrm{m}}\right)^{0.11} -->
<!-- $$u(z) = u_{2m} (\frac{z}{2}^{0.11}) $$ -->
The exponent is set to 0.11 to assume near-neutral stratification and almost undisturbed flow over the ocean ([Hsu et al., 1994](https://journals.ametsoc.org/view/journals/apme/33/6/1520-0450_1994_033_0757_dtplwp_2_0_co_2.xml?tab_body=pdf)).

**Cloud base height** is the height of the lowest cloudy pixel of each cloud object. It should be noted that in this processing setup, no distinction is made between cloud pixels and precipitation pixels. Therefore, the lowest pixel of a cloud object can also be the lowest pixel with rain. Accordingly, **cloud top height** is the height of the highest cloudy pixel.
**Cloud depth** is then calculated from the difference between cloud top height and cloud base height for each segmented cloud object.

## 4. Saving data
During the final step, the data are read from temporary files and saved into NetCDF files. Additionally, attributes for each variable, as well as global attributes are set. The resulting NetCDF file contains variables with the dimension "cloud number" and with the dimension(s) "time" and/or "height".

The cloud parameters that were described in the previous step (3. Calculating cloud parameters) were calculated for each cloud object and thus have the dimension "cloud number". Additionally, a numerical mask with the associated cloud object IDs for each pixel are given with the dimensions time and height. Flags for radar status and missing wind are given with the dimension time.

The NetCDF header for one example file is given below:

```
dimensions:
	N_clouds = 55773 ;
	height = 428 ;
	time = 2963520 ;
variables:
	double cloudBase(N_clouds) ;
		cloudBase:units = "m" ;
		cloudBase:long_name = "Cloud base height" ;
	double cloudDepth(N_clouds) ;
		cloudDepth:units = "m" ;
		cloudDepth:long_name = "Cloud depth" ;
	double cloudEndTime(N_clouds) ;
		cloudEndTime:units = "seconds since 1970-1-1 0:00:00 UTC" ;
		cloudEndTime:long_name = "End of cloud" ;
	double cloudLength(N_clouds) ;
		cloudLength:units = "m" ;
		cloudLength:long_name = "Cloud length" ;
	double cloudStartTime(N_clouds) ;
		cloudStartTime:units = "seconds since 1970-1-1 0:00:00 UTC" ;
		cloudStartTime:long_name = "Beginning of cloud" ;
	double cloudTop(N_clouds) ;
		cloudTop:units = "m" ;
		cloudTop:long_name = "Cloud top height" ;
	double height(height) ;
		height:units = "m" ;
		height:long_name = "Height" ;
	double numMask(time, height) ;
		numMask:units = " " ;
		numMask:long_name = "Cloud object ID" ;
	double N_clouds(N_clouds) ;
		N_clouds:units = " " ;
		N_clouds:long_name = "Total number of clouds" ;
	double time(time) ;
		time:units = "seconds since 1970-1-1 0:00:00 UTC" ;
		time:long_name = "Time" ;
	double status(time) ;
		status:units = " " ;
		status:long_name = "Radar On (1), radar Off (0), no data file (2), radar scanning (3)" ;
	double wind_missing(time) ;
		wind_missing:units = " " ;
		wind_missing:long_name = "Wind measured (0), wind missing (1)" ;

// global attributes:
		:information = "Cloud mask applied to BCO radar" ;
		:contact = "heike.konow@mpimet.mpg.de" ;
		:location = "The Barbados Cloud Observatory, Deebles Point, Barbados, West Indies" ;
		:institution = "Max Planck Institute for Meteorology, Hamburg" ;
		:instrument = "MBR2 cloud radar"
```

## Usage
If you want to run the code, start by renaming config.sample.m to config.m and adapt the paths in this file to you system. Decide on the parameters for the processing in the beginning of the file bco_run_cloudmask.m and run.

## References
Dillencourt, M. B., Samet, H., & Tamminen, M. (1992). A general approach to connected-component labeling for arbitrary image representations. Journal of the ACM, 39(2), 253–280. https://doi.org/10.1145/128749.128750

Hsu, S. A., Meindl, E. A., & Gilhousen, D. B. (1994). Determining the power-law wind-profile exponent under near-neutral stability conditions at sea. Journal of Applied …, 33, 757–765. https://doi.org/10.1175/1520-0450(1994)
