

Processing of stirling pressure data for waves and sea level is as follows


1) READ HOURLY OR DAILY TEXT FILES AND SAVE TO MATFILE [MONTH CHUNKS]

driver_read_save_stirling_pressure_v9_daily.m  					[read daily text files] 
driver_read_save_stirling_pressure_v9_hourly.m


		import_pressure_logfiles.m    							[import hourly logfiles]
		import_daily_pressure.m 	  							[import daily text files 2021+ from Fremantle port Authority for Stirling Channel]
		import_daily_BAD_pressure.m  						    [import some of the badly formatted daily text files]
		importStirling2020.m		  							[try to import 2020 csv data that has no date only time of day - not that useful]
		subset_clean_Stirling_annual_to_monthly_matfiles.m		[divide chari's annual 2hz data from stirling to monthly matfiles and clean outliers]



2) CALCULATE WAVES AND SEALEVEL, SAVE HOURLY OUTPUTS TO TEXT AND MATFILES [IN MONTH CHUNKS]

driver_get_stirling_waves_v9_hourly_loop.m 						[directly read hourly test fiels and calculate waves -- this has not been updated since process_Stirling_waves_sealevel_v9.m created, so better to read and save the data then use that file ]
driver_get_stirling_waves_v9_daily_loop.m  						[directly read daily text files and calculate waves - - allows for comparison between chari method and ocn toolbox]
**process_Stirling_waves_sealevel_v9.m 							[reads pressure sensor matfile created with  [driver_read_save_stirling_pressure_v9_daily.m]
  									    						 or [driver_read_save_stirling_pressure_v9_hourly.m] and calculates
 																 waves using ocn toolbox, tidal analysis, lowpass filter, highpass filter (optional)
 																 and makes a number of plots that are saved to a directory along
 																 with the processed data as matfile (optional) or the hourly wave data as a csv text file
 										 						**There are a lot of m-files required so you need to set the path to make
  																sure all the m-files in the yasha_RBR_tools_v9 are available.] ** this is preferred tested script


3) COMBINED THE MONTHLY OUTPUTS INTO ANNUAL [HOURLY] TIME SERIES, PLOT, SAVE AST MATFILE TABLES AND CSV TEXTFILES

cat_Stirling_Outputs_months2year.m  							[concatinate hourly outputs (monthly) from wave analysis created by process_Stirling_waves_sealevel_v9.m into year data]








NOTES:
DATA CLEANING USES THESE FILES:
interp1gap.m  [performs interpolation over small gaps in 1D data]
get_outliers_yh.m [Use moving window to get indexes of outliers > nd standard deviations from the window median. Window  can be balanced or offset around center]



Data shared with chari here: /Users/00068592/Dropbox/Public/Chari_Stirling_wave_Pressure_data/

Outputs saved here: /Users/00068592/Documents/RESEARCH/OTHER_PEOPLE/CHARI/Cockburn_pressure_sensor/processed_Stirling_waves_test/
