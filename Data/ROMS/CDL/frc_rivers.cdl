netcdf frc_rivers {

dimensions:
	xi_rho = 386 ;
	xi_u = 385 ;
	xi_v = 386 ;
	eta_rho = 130 ;
	eta_u = 130 ;
	eta_v = 129 ;
	s_rho = 30 ;
	river = 18 ;
	river_time = UNLIMITED ; // (0 currently)

variables:
	double river(river) ;
		river:long_name = "river runoff identification number" ;
	double river_Xposition(river) ;
		river_Xposition:long_name = "river XI-position at RHO-points" ;
		river_Xposition:valid_min = 1. ;
		river_Xposition:valid_max = 385. ;
	double river_Eposition(river) ;
		river_Eposition:long_name = "river ETA-position at RHO-points" ;
		river_Eposition:valid_min = 1. ;
		river_Eposition:valid_max = 129. ;
	double river_direction(river) ;
		river_direction:long_name = "river runoff direction" ;
	double river_Vshape(s_rho, river) ;
		river_Vshape:long_name = "river runoff mass transport vertical profile" ;
	double river_time(river_time) ;
		river_time:long_name = "river runoff time" ;
		river_time:units = "days since 1992-01-01 00:00:00" ;
		river_time:add_offset = 2448623. ;
	double river_transport(river_time river) ;
		river_transport:long_name = "river runoff vertically integrated mass transport" ;
		river_transport:units = "meter3 second-1" ;
		river_transport:time = "river_time" ;
	double river_temp(river_time s_rho, river) ;
		river_temp:long_name = "river runoff potential temperature" ;
		river_temp:units = "Celsius" ;
		river_temp:time = "river_time" ;
	double river_salt(river_time s_rho, river) ;
		river_salt:long_name = "river runoff salinity" ;
		river_salt:time = "river_time" ;

// global attributes:
		:type = "ROMS FORCING file" ;
		:title = "NENA River Forcing" ;
		:grd_file = "roms_nena_grid_3.nc" ;
		:rivers = "(1) Mississippi River at Vicksburg, MS, (2) Susquehanna River at Conowingo, MD, (3) Tombigbee River at Coffeeville L&D near Coffeeville, (4) Apalachicola River near Sumatra, FL, (5) Connecticut River at Hartford, CT, (6) Hudson River at Green Island NY, (7) Penobscot River at Eddington, ME, (8) Delaware River at Trenton NJ, (9) Altamaha River at Doctortown, GA, (10) Potomac River near Wash, DC Little Falls Pump, (11) Savannah River near Clyo, GA, (12) Roanoke River at Roanoke Rapids, NC, (13) Pee Dee River at Peedee, SC, (14) Pearl River near Bogalusa, LA, (15) Suwannee River near Wilcox, FL, (16) Kennebec River at North Sidney, ME,(17) Escambia River near Molino, FL, (18) James River near Richmond, VA" ;
}
