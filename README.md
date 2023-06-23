# yasha-RBR-tools
Development Matlab toolbox for processing of pressure sensor data from RBR or other instruments including wave calculations

Derived from Yasha-RBR_tools_v9 (latest working version)
 
Toolbox to read rbr pressure sensor data and calculate waves and sealevel variables
Uses ocn toolbox to get waves which allows for pressure attenuation correction.

** process_RBR_waves_sealevel_v9.m -- this is the main script used to do all the processing. Must carefully edit all the settings in this m-file
 
% This m-file reads raw .rsk file downloaded from RBR pressure sensors
% using rsk-tools (downloadad from rbr website)
% and saves to netcdf file. then it reads the netcdf file and calculates
% waves (optional), tidal analysis, lowpass filter, highpass filter (optional)
% and makes a number of plots that are saved to a directory along
% with the processed data as matfile (optional) or the hourly wave data as a csv text file
% **There are a lot of m-files required so you need to set the path to make
% sure all the m-files in the yasha_RBR_tools are available.
% 
%  Waves now correctly calculated with pressure response for depth using
%  Oceanlyz toolbox
%  http://zarmesh.com/wp-content/uploads/2016/05/Wind-wave-analysis-in-depth-limited-water-using-OCEANLYZ-A-MATLAB-toolbox.pdf
%  https://oceanlyz.readthedocs.io/en/latest/
% 
## contents
process_RBR_waves_sealevel_v9.m  - main driver for complete wave analysis and averaging of water level, filtering,tidal analysis, save outputs to hourly text files

%% process waves using chari's scripts made into a function 
get_chari_wave_analysis.m
autospec.m
opclevel.m
oppsd.m
oppsd2.m
opsmooth.m
periodo.m
get_tidespecft.m  - function to produce time-frequency plots from sea level data
wavenumber.m
wavepar.m
wavepar2_yh.m
wavepar_yh.m


interp1gap.m
interp_r.m
lanczosfilter.m
lp_filter.m
highpass.m
nanmean.m
newdatestr.m
fillgapYLH.m
get_beginfinish.m 		-interactively select start and finish times
get_finish.m			-interactively select only finish time (start at beginning)
get_outliers_yh.m
date2doy.m
rotateXLabels.m
savefig2_ylh.m

read_RBRyh.m    - old version  for reading rsk files. outdated,  but works on older instruments. replaced by v2
read_RBRyh_v2.m - read raw rsk file and save to netcdf


importDOTfile.m					- read departlent of transport 5min tide gauge file
process_DOT_5min_waterlevel.m   - process DOT tide gauge file
detrend_notes_yasha.m			- how to detrend sea level, etc.


##								    Toolboxes used
t_tide_v1/ - tidal analysis
subaxis/  - better multi axis plots
rbr-rsktools-3.5.3/ -  tools to read rsk data from rbr sensors
plot_bom_wind_arrows_yh/ - plot seabreeze style wind arrows
oceanlyz_2_0/  - main toolbox for calculating waves
altmany-export_fig-76c775b/ - optional way to save figures  and retain how it looks on screen. requires ghostscript to be installed


----NOT used ----
junk/
outdated_mfiles/
extra_outdated_mfiles/



##								    Fremantle Port - Stirling pressure sensor data analysis

subset_clean_Stirling_annual_to_monthly_matfiles.m
process_Stirling_waves_sealevel_v9.m
driver_get_stirling_waves_v9_daily_loop.m
driver_get_stirling_waves_v9_hourly_loop.m
driver_read_save_stirling_pressure_v9_daily.m
driver_read_save_stirling_pressure_v9_hourly.m
cat_Stirling_Outputs_months2year.m
importStirling2020.m
import_daily_BAD_pressure.m
import_daily_pressure.m
import_pressure_logfiles.m

------------------------------------------------------------------------------------------
					random minor use functions - not sure if needed
------------------------------------------------------------------------------------------
doc_datacursormode.m
myupdatefcn.m
