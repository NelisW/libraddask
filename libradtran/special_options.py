
"""--------------------------------------------------------------------
 * $Id: special_options.py 3112 2015-05-19 13:17:35Z robert.buras $
 * 
 * This file is part of libRadtran.
 * Copyright (c) 1997-2012 by Arve Kylling, Bernhard Mayer,
 *                            Claudia Emde, Robert Buras
 *
 * ######### Contact info: http://www.libradtran.org #########
 *
 * This program is free software; you can redistribute it and/or 
 * modify it under the terms of the GNU General Public License   
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.        
 * 
 * This program is distributed in the hope that it will be useful, 
 * but WITHOUT ANY WARRANTY; without even the implied warranty of  
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the   
 * GNU General Public License for more details.                    
 * 
 * You should have received a copy of the GNU General Public License          
 * along with this program; if not, write to the Free Software                
 * Foundation, Inc., 59 Temple Place - Suite 330, 
 * Boston, MA 02111-1307, USA.
 *--------------------------------------------------------------------"""


from . import option_definition
from . import GUI_definition
import io
class setup_special_group():

	group_name = 'Special'

	def __init__(self):
		documentation = get_special_documentation()

		brdf_non_parents=['albedo','albedo_file','albedo_library',
                                  'brdf_rossli', 'brdf_rossli_file', 'brdf_rossli_hotspot',
                                  'brdf_cam', 'brdf_cam_solar_wind', 'cox_and_munk_u10_map', 
                                  'brdf_rpv_file', 'brdf_rpv', 'brdf_rpv_library','brdf_rpv_type',
                                  'brdf_hapke', 'brdf_hapke_file',
                                  'mc_albedo_file', 'mc_rossli_file', 'mc_rpv_type', 'mc_rpv_file',
                                  'bpdf_tsang_u10']

		albedo_library = option_definition.option(
			name='albedo_library',
			group='special',
			helpstr='Spectral albedos of different surface types.',
			documentation=documentation['albedo_library'], 
			tokens=option_definition.addToken(name='Input.alb.library', datatype=str, valid_range=['IGBP',io.IOBase]),
			parents=['uvspec'],
			non_parents=brdf_non_parents,
			childs=['brdf_rpv_type'],
			showInGui = False,
                        developer = True
		)
	
		cloud_fraction_map = option_definition.option(
			name='cloud_fraction_map',
			group='special',
			documentation=documentation['cloud_fraction_map'],
			tokens=option_definition.addToken(name='Input.filename[FN_CLOUD_FRACTION_MAP]',datatype=io.IOBase),
			developer=True,
                )

		albedo_map = option_definition.option(
			name='albedo_map',
			group='special',
			documentation=documentation['albedo_map'],
			tokens=[option_definition.addSetting(name='Input.alb.source', setting='ALBEDO_FROM_ALBEDO_MAP') ,
				option_definition.addToken(name='Input.filename[FN_ALBEDO_MAP]',datatype=io.IOBase),
				option_definition.addToken(name='Input.alb.netCDF_alb_name',datatype=io.IOBase,optional=True),
				option_definition.addToken(name='Input.alb.scale_factor',datatype=float,optional=True) ],
                        developer=True
                )

		altitude_map = option_definition.option(
			name='altitude_map',
			group='special',
			documentation=documentation['altitude_map'],
			tokens=[option_definition.addSetting(name='Input.alt.source', setting='ALT_FROM_MAP') ,
				option_definition.addToken(name='Input.filename[FN_ALTITUDE_MAP]',datatype=io.IOBase),
				option_definition.addToken(name='Input.alt.netCDF_alt_name',datatype=io.IOBase,optional=True),
				option_definition.addToken(name='Input.alt.scale_factor',datatype=float,optional=True) ],
                        developer=True
                )

		cox_and_munk_u10_map = option_definition.option(
			name='cox_and_munk_u10_map',
			group='special',
			helpstr='Location of netCDF file with wind speed for ocean BRDF.', 
			documentation=documentation['cox_and_munk_u10_map'], 
			tokens = [ option_definition.addToken(name='Input.filename[FN_U10_MAP]', datatype=io.IOBase), 
				option_definition.addSetting(name='Input.disort2_brdf', setting='BRDF_CAM', default='BRDF_NONE') ], 
			parents=['uvspec'],
			non_parents=brdf_non_parents,
                        non_parent_exceptions=['brdf_cam', 'brdf_cam_solar_wind', 'cox_and_munk_u10_map', 'mc_rossli_file'],
			showInGui= False,
                        developer=True
		)

		cox_and_munk_pcl_map = option_definition.option(
			name='cox_and_munk_pcl_map',
			group='special',
			documentation=documentation['cox_and_munk_pcl_map'],
			tokens=[option_definition.addSetting(name='Input.disort2_brdf', setting='BRDF_CAM') ,
				option_definition.addToken(name='Input.filename[FN_PIGMENTS_MAP]',datatype=io.IOBase),
				option_definition.addToken(name='Input.cm.pcl_netCDF_name',datatype=io.IOBase,optional=True),
				option_definition.addToken(name='Input.cm.pcl_scale_factor',datatype=float,optional=True) ],
                        developer=True
		)

		cox_and_munk_sal_map = option_definition.option(
			name='cox_and_munk_sal_map',
			group='special',
			documentation=documentation['cox_and_munk_sal_map'],
			tokens=[option_definition.addSetting(name='Input.disort2_brdf', setting='BRDF_CAM') ,
				option_definition.addToken(name='Input.filename[FN_SALINITY_MAP]',datatype=io.IOBase),
				option_definition.addToken(name='Input.cm.sal_netCDF_name',datatype=io.IOBase,optional=True),
				option_definition.addToken(name='Input.cm.sal_scale_factor',datatype=float,optional=True) ],
                        developer=True
		)

		emissivity_map = option_definition.option(
			name='emissivity_map',
			group='special',
			documentation=documentation['emissivity_map'],
			tokens=[option_definition.addSetting(name='Input.alb.source', setting='ALBEDO_FROM_EMISSIVITY_MAP') ,
				option_definition.addToken(name='Input.filename[FN_ALBEDO_MAP]',datatype=io.IOBase),
				option_definition.addToken(name='Input.alb.netCDF_alb_name',datatype=io.IOBase,optional=True),
				option_definition.addToken(name='Input.alb.scale_factor',datatype=float,optional=True) ],
                        developer=True
		)

		surface_temperature_map = option_definition.option(
			name='surface_temperature_map',
			group='special',
			documentation=documentation['surface_temperature_map'],
			tokens=[option_definition.addToken(name='Input.filename[FN_SURFACE_TEMP_MAP]',datatype=io.IOBase),
				option_definition.addToken(name='Input.netCDF_name_surf_T',datatype=io.IOBase,optional=True)],
                        developer=True
		)

		surface_type_map = option_definition.option(
			name='surface_type_map',
			group='special',
			documentation=documentation['surface_type_map'],
			tokens=[option_definition.addSetting(name='Input.alb.surf_type_map',setting='TRUE'),
				option_definition.addToken(name='Input.filename[FN_SURFACE_TYPE_MAP]',datatype=io.IOBase),
				option_definition.addToken(name='Input.alb.netCDF_surf_name',datatype=io.IOBase,optional=True)],
                        developer=True
		)

		ECHAM_sza = option_definition.option(
			name='ECHAM_sza',
			group='special',
			documentation=documentation['ECHAM_sza'],
			tokens=[option_definition.addToken(name='Input.filename[FN_SZA]',datatype=io.IOBase),
				option_definition.addSetting(name='Input.atm.sza',setting='NOT_DEFINED_FLOAT'),
				option_definition.addSetting(name='Input.atm.sza_source',setting='SZA_ECHAM') ],
			developer=True,
		)

		ECMWF_ozone_climatology = option_definition.option(
			name='ECMWF_ozone_climatology',
			group='special',
			documentation=documentation['ECMWF_ozone_climatology'],
			tokens=option_definition.addSetting(name='Input.atm.ECMWF_ozone_climatology',setting='TRUE'),
                        developer=True
		)

		ECMWF_wc_file = option_definition.option(
			name='ECMWF_wc_file',
			group='special',
			documentation=documentation['ECMWF_wc_file'],
			tokens=[option_definition.addSetting(name='isp',setting=option_definition.CaothType('wc')),
				option_definition.addSetting(name='Input.caoth[isp].source',setting='CAOTH_FROM_ECMWF'),
				option_definition.addToken(name='Input.caoth[isp].filename',datatype=io.IOBase) ],
                        developer=True
		)

		ECMWF_ic_file = option_definition.option(
			name='ECMWF_ic_file',
			group='special',
			documentation=documentation['ECMWF_ic_file'],
			tokens=[option_definition.addSetting(name='isp',setting=option_definition.CaothType('ic')),
				option_definition.addSetting(name='Input.caoth[isp].source',setting='CAOTH_FROM_ECMWF'),
				option_definition.addToken(name='Input.caoth[isp].filename',datatype=io.IOBase) ],
                        developer=True
		)

		ECMWF_levels_only = option_definition.option(
                        name='ECMWF_levels_only',
                        group='special',
                        documentation=documentation['ECMWF_levels_only'],
                        tokens=option_definition.addSetting(name='Input.atm.rs_add_upper_levels',setting='FALSE'),
                        developer=True
                )

		ECMWF_wind_file = option_definition.option(
                        name='ECMWF_wind_file',
                        group='special',
                        documentation=documentation['ECMWF_wind_file'],
                        tokens=option_definition.addToken(name='Input.filename[FN_ECMWF_WIND_MAP]',datatype=io.IOBase),
                        developer=True
                )

		satellite_geometry = option_definition.option(
                        name='satellite_geometry',
                        group='special',
                        documentation=documentation['satellite_geometry'],
                        tokens=option_definition.addToken(name='Input.filename[FN_SATGEOMETRY]',datatype=io.IOBase),
                        developer=True
                )

		satellite_pixel = option_definition.option(
                        name='satellite_pixel',
                        group='special',
                        documentation=documentation['satellite_pixel'],
                        tokens=[option_definition.addToken(name='Input.sat_pixel_x',datatype=int),
				option_definition.addToken(name='Input.sat_pixel_y',datatype=int) ],
                        developer=True
                )

		self.options = [ ECHAM_sza,
				ECMWF_ozone_climatology, ECMWF_wc_file, ECMWF_ic_file, ECMWF_wind_file, ECMWF_levels_only,
				cox_and_munk_u10_map, cox_and_munk_pcl_map, cox_and_munk_sal_map,
				cloud_fraction_map, altitude_map, albedo_map, emissivity_map,
				surface_temperature_map, surface_type_map,
				albedo_library, satellite_geometry, satellite_pixel,
				]

		#Special group belongs to developers options; should not be in Gui for EsasLight version; 
		for opt in self.options:
			opt.showInGui = False

	def __iter__(self):
		return iter(self.options)

def get_special_documentation():
	return {
		'albedo_library' : r'''
	Albedo libraries are a collection of spectral albedos of different surface types. 
	This option must be used either with \code{bdrf_rpv_type} 
	or \code{surface_type_map}, in order to select the specific surface type.
	There are two possibilities for libraries: 
	the built-in IGBP library or a user defined albedo library.

	The built-in library of the International Geosphere Biosphere Programme is 
	selected with 
	\fcode{
	     albedo_library IGBP
	}
	The IGBP library contains 20 surface types which are set by \code{brdf_rpv_type}:
	\fcode{
	 1 evergreen_needle_forest\\
	 2 evergreen_broad_forest\\
	 3 deciduous_needle_forest\\
	 4 deciduous_broad_forest\\
	 5 mixed_forest\\
	 6 closed_shrub\\
	 7 open_shrubs\\
	 8 woody_savanna\\
	 9 savanna\\
	10 grassland\\
	11 wetland\\
	12 cropland\\
	13 urban\\
	14 crop_mosaic\\
	15 antarctic_snow\\
	16 desert\\
	17 ocean_water\\
	18 tundra\\
	19 fresh_snow\\
	20 sea_ice\\
	}
	Surface types 1 - 17 are defined by the International Geosphere Biosphere 
	Programme (IGBP); additionally there are tundra, fresh\_snow, and sea\_ice surface types.
	The spectral albedo of the ground is determined as a function of solar 
	zenith angle, precitable water, and clouds. The spectral resolution 
	equals the grid of the correlated-k Fu/Liou parameterisation.
	This library originates from the NASA CERES/SARB Surface Properties Project, 
	see \citet{Belward1996}.

	For creating your own albedo library use \code{albedo_library path}, where \code{path} is 
	the path of the directory where the albedo data is stored. The files are expected to have the 
	names \code{albedo_01.dat, albedo_02.dat, ...} If \code{brdf_rpv_type 1} is specified
	the albedo from \code{albedo_01.dat} will be used, and so on.  
	Each file is required to have two columns:
	Column 1 is the wavelength in nm, and column 2 the corresponding 
	Lambertian surface albedo. The wavelength grid may be freely set. The 
	albedo will be interpolated linearely to the wavelength grid used for the 
	radiation calculation. Comments start with \code{\#}. Empty lines are ignored.
	This option is similar to \code{albedo_file}, except that it offers an
	easy way to use the option \code{surface_type_map} in combinition with albedo files.
		''',

		'albedo_map' : r'''
	\emph{This option is preliminary and still subject to change (no wavelength dependency yet)!}
	It gives the possibility to specify a wavelength independent albedo with the help of 
	a \emph{netCDF} file, which is used in combinition 
	with the options \code{latitude}, \code{longitude}, and \code{time}.
	\fcode{
	albedo_map file [variable_name]
	}
	Here \code{file} is the location of the \emph{netCDF} file,
	The optional argument allows the name of the albedo variable in 
	the \emph{netCDF} file to be specified (the default name is \code{AL}).
	The albedo must be provided as function of \code{latitude} and \code{longitude} \code{AL(lat, lon)}, 
	and may also depend on time \code{AL(time, lat, lon)}. 
	The latitude, longitude, and time grids must be provided as doubles 
	\code{double lat(lat)}, \code{double lon(lon)}, and \code{double time(time)}.
	\code{uvspec} reads the value at the nearest pixel to the given \code{latitude} and 
	\code{longitude}. No spatial interpolation or averaging of the values are performed.
	If a time-dependent albedo is provided, the albedo data nearest to the specified time 
	will be selected (or linear interpolated if \code{time_interpolate} is switched on). 
		''',

		'altitude_map' : r'''
	Specifies an altitude map which is used in combinition with  
	\code{latitude}, \code{longitude} in order to  
	select the altitude for the simulation. 
	No interpolation is done between the pixels of the map. 
	The format of the call is: 
	\fcode{ 
	altitude_map file [variable_name] 
	} 
	where \code{file} is the location of the altitude map file. The map is expected to  
	be in \emph{netCDF} format. The file must contain \code{double lat(lat)},  
	\code{double lon(lon)}, and the altitude variable, where \code{variable_name} is the  
	name of the surface elevation variable in the \emph{netCDF} file. The default name is Z. 
	The altitude variable must be \code{altitude(lat, lon)}.  
	For format discribtion see also the example map included in {\sl libRadtran},  
	\file{data/altitude/ELEVATION_GTOPO_10min.cdf}. 
	To use this map in \code{uvspec}, you may also use \code{altitude_map GTOPO}. This map has a resolution  
	of 10 arc minutes and the unit of the altitude is meter. Please note that this resolution 
	might not ne adequate for your application.  
	If an altitude in the map is below the lowest level of the  
	\code{atmosphere_file}, the atmospheric profiles are extrapolated assuming a constant  
	gradient for temperature and mixing ratios. 
		''',

		'cloud_fraction_map' : r'''
	Undocumented option.
		''',

		'cox_and_munk_u10_map' : r'''
	Specify wind speed (in m/s) for the \citet{cox54a,cox54b}
	ocean BRDF with the help of an \emph{netCDF} file, which is used in combinition 
	with the options \code{latitude}, \code{longitude}, and \code{time}.
	\fcode{
	cox_and_munk_u10_map file
	}
	where \code{file} is the location of the \emph{netCDF} file. 
	{\sl libRadtran} reads the value at the nearest pixel to the given \code{latitude} and 
	\code{longitude}. No spatial interpolation or averaging of the values is done.

	The file must contain the elements of the wind vector \code{U10} and
	\code{V10}.  These must be specified as functions of \code{latitude}
	and \code{longitude} \code{U10(lat, lon)}, \code{V10(lat, lon)}, or
	additionally may also depent on time \code{U10(time, lat, lon)},
	\code{V10(time, lat, lon)}.  If the variable time is present in the
	file, the wind speed will be interpolated according to the option
	\code{time_interpolate}. All grids must be provided as \code{double
	lat(lat)}, \code{double lon(lon)}, and \code{double time(time)}.
		''',

		'cox_and_munk_pcl_map' : r'''
	A possibility to specify pigment concentration (in mg/m3) for the Cox and Munk 
	ocean BRDF with the help of an \emph{netCDF} file, which is used in combinition 
	with options \code{latitude}, \code{longitude}, and \code{time}.
	\fcode{
	cox_and_munk_pcl_map file [variable_name]
	}
	where \code{file} is the location of the \emph{netCDF} file. 
	{\sl libRadtran} reads the value at the nearest pixel to the given \code{latitude} and 
	\code{longitude}. No spatial interpolation or averaging of the values is done.
	
	The default name of the pigment concentration variable is
	\code{chlorophyll}, but can be changed with the optional argument
	\code{variable_name}.  The pigment concentration must be provided as
	function of \code{latitude} and \code{longitude},
	\code{chlorophyll(lat, lon)}, or additionally may also depend on time
	\code{chlorophyll(time, lat, lon)}.  If a time-dependent pigment
	concentration is specified, the pigment concentration will be
	interpolated according to the option \code{time_interpolate}. All
	grids must be provided in the file as \code{double lat(lat)},
	\code{double lon(lon)}, and \code{double time(time)}.
		''',

		'cox_and_munk_sal_map' : r'''
	Specify ocean salinity (in ppt) for the \citet{cox54a,cox54b}
	ocean BRDF with the help of an \emph{netCDF} file, which is used in combinition 
	with the options \code{latitude}, \code{longitude}, and \code{time}.
	\fcode{
	cox_and_munk_pcl_map file [variable_name]
	}
	where \code{file} is the location of the \emph{netCDF} file. 
	{\sl libRadtran} reads the value at the nearest pixel to the given \code{latitude} and 
	\code{longitude}. No spatial interpolation or averaging of the values is done.
	
	The expected name of the pigment concentration variable is per default
	\code{salinity}, but can be changed with the optional argument
	\code{variable_name}.  The pigment concentration must be provided as
	function of \code{latitude} and \code{longitude}, \code{salinity(lat,
	lon)}, or additionally may also depent on time \code{salinity(time,
	lat, lon)}.  If a time-dependent salinity is specified, the salinity
	will be interpolated according to the option
	\code{time_interpolate}. All grids must be provided as \code{double
	lat(lat)}, \code{double lon(lon)}, and \code{double time(time)}.
		''',	

		'emissivity_map' : r'''
	\emph{This option is preliminary and still subject to change (no wavelength dependency yet)!} 
	Specify a wavelength independent emissivity with the help of  
	an \emph{netCDF} file, which is used in combinition  
	with the options \code{latitude}, \code{longitude}, and \code{time}. 
	\fcode{ 
	emissivity_map file [variable_name] 
	} 
	where \code{file} is the location of the \emph{netCDF} file.  
	With the optional argument \code{variable_name} the name of the emissivity variable in  
	the \emph{netCDF} file can be specified. (By default the expected name is \code{EMIS}.) 
	The emissivity must be specified as function of \code{latitude} and \code{longitude} \code{EMIS(lat, lon)},  
	or additionally may also depent on time \code{EMIS(time, lat, lon)}.  
	All grids must be provided as  
	\code{double lat(lat)}, \code{double lon(lon)}, and \code{double time(time)}. 
	{\sl libRadtran} reads the value at the nearest pixel to the given \code{latitude} and  
	\code{longitude}. No spatial interpolation or averaging of the values is done. 
	If the variable time is present in the file, the emissivity data nearest to the specified time  
	will be selected (or interpolated if \code{time_interpolate} is switched on).  
		''',

		'ECHAM_sza' : r'''
	Undocumented option.
		''',

		'ECMWF_ozone_climatology' : '''
	The Integrated Forecast System (IFS) of the ECMWF uses a ozone 
	climatology for radiative transfer instead of the ozone simulated by 
	the IFS. If this option is activated the ozone profile of the 
	\code{atmosphere_file} or \code{ECMWF_atmosphere_file} is replaced by 
	the ozone climatology by Fortuin and Langematz (1995). 
	(If there is also a \code{mol_file} for ozone, it modifies the ozone 
	climatology profile.) 
		''',

		'ECMWF_ic_file' : r'''
	\fcode{ 
	ECMWF_ic_file file 
	} 
	For further information see \code{ECMWF_wc_file}.
		''',

		'ECMWF_levels_only' : r'''
	The atmosphere considered in the simulation has the same height range 
	as the data in the 
	\code{ECMWF_atmosphere_file}/\code{radiosonde}-file. No further levels 
	are added above those.  This option has only an effect in combination 
	with \code{ECMWF_atmosphere_file} or \code{radiosonde} 
	(this option is identical to \code{radiosonde_levels_only}). 
		''',

		'ECMWF_wc_file' : r'''
	Reads in combination with the options \code{latitude}, 
	\code{longitude}, and \code{time} (all mandatory) the pressure, 
	temperature, and cloud liqid water content (CLWC) and cloud cover (CC) 
	from an ECMWF \emph{netCDF} data file. 
	\fcode{ 
	ECMWF_wc_file file 
	} 
	No spatial interpolation of the values is done.  The data nearest to 
	the specified \code{time} will be selected (or linearly interpolated 
	if \code{time_interpolate} is switched on).  In order to use the ECMWF 
	data without cloud overlap assumption, use \code{cloud_overlap off}. 
		''',

		'ECMWF_wind_file' : r'''
	Reads in combination with the options \code{latitude}, 
	\code{longitude}, and \code{time} (all mandatory) the wind components 
	U, V, and W from an ECMWF \emph{netCDF} data file. 
	\fcode{ 
	ECMWF_wind_file file 
	} 
	The data nearest to the specified \code{time} will be selected (or 
	linearly interpolated, if \code{time_interpolate} is switched on). 
		''',

		'satellite_geometry' : r'''
	With this option the satellite geometry is determinded. The argument for this option 
	\fcode{ 
	satellite_geometry  netCDF_file  
	} 
	is the location of a \code{netCDF_file}, which must contain latitude and longitude position 
	as well as zenith and azimuth viewing angle for each pixel. 
 
	%ak20100306 What is below does not belong in the general Users Guide yet. 
	%For format specification see the example file data/satellite/MSG\_seviri/MSG\_seviri\_geometry.nc. 
	%This option has to be used in combination with \code{satellite_pixel}. 
	%You can also use the abbreviation 
	%\code{satellite_geometry MSG} 
	%for MSG simulation, \emph{but in order to use this, you have to copy the netCDF file} 
	%\fcode{ 
	%cp -r /data/A3/satellite_geometry/satellites libRadtran/data/ 
	%} 
	%\emph{which is not in cvs version, as it is too large.} 
		''',

		'satellite_pixel' : r'''
	This option specifies which pixel of the satellite image that should be simulated.  
	\fcode{ 
	satellite_pixel  pixel_x pixel_y 
	} 
	The arguments \code{pixel_nr_x} and \code{pixel_nr_y} specifies the pixel position in the native  
	system of the satellite, which is determinded by the option \code{satellite_geometry}. 
		''',

		'surface_temperature_map' : r'''
	Specify a surface\_temperature map with a \emph{netCDF} file which is used in combinition 
	with the options \code{latitude}, \code{longitude}, and \code{time}.
	\fcode{
	surface_temperature_map file [variable_name]
	}
	where \code{file} is the location of the \emph{netCDF} file. 
	\code{libRadtran} reads the value at the nearest pixel to the given \code{latitude} and 
	\code{longitude}. No spatial interpolation or averaging of the values is done.
		''',

		'surface_type_map' : r'''
	Specify a surface type map, which is used in combinition with 
	\code{albedo_library}, \code{latitude}, and \code{longitude} in order to 
	select the surface type relevant for the simulation.
	No pixel interpolation is done.
	The format of the call is:
	\fcode{
	surface_type_map file [variable_name]
	}
	where \code{file} is the location of the surface type map file. The map is expected to 
	be in \emph{netCDF} format. The file must contain the variables \code{double lat(nlat)}, 
	\code{double lon(nlon)}, and \code{byte brdf_rpv_type (nlat, nlon)}. If the name of the 
	surface type variable is different, the optional argument can be used in order to specify 
	the variable name. For format specification see also 
	\file{data/albedo/IGBP_map/SURFACE_TYPE_IGBP_10min.cdf}.

	For using the IGBP map, the call is \code{surface_type_map IGBP}. This map has a resolution of 10 
	minutes and contains the surface types 1 to 18 defined in the \code{albedo_library IGBP}. 
	Fresh snow and sea ice are not included, as their extent is too variable. Attention:  
	That implies e.g. that the Arctic is considered ocean\_water and not sea\_ice!

	Locations on the pixel boundaries are interpreted as the pixel northward and
	eastward respectively. E.g. location 0 N, 0 E is interpreted like the pixel ranging from
	0 to 10min North and from 0 to 10min East.
		''',
}
