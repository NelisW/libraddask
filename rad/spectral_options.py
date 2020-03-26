
"""--------------------------------------------------------------------
 * $Id: spectral_options.py 3145 2015-08-04 13:35:30Z josef.gasteiger $
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

class setup_spectral_group():

	group_name = 'Spectral/Source'

	def __init__(self):
		documentation = get_spectral_documentation()

		wavelength = option_definition.option(
			name='wavelength',
			group='spectral',
			helpstr='Set the wavelength range by specifying first and last wavelength in nm.',
			documentation=documentation['wavelength'],
			gui_inputs=(GUI_definition.FloatInput(name='Input.wl.start', default=None, valid_range=[0, 1000000.0]),
				    GUI_definition.FloatInput(name='Input.wl.end', default=None, valid_range=[0, 1000000.0], optional=True),),
			tokens = [ option_definition.addToken(name='Input.wl.start', datatype=float, default='NOT_DEFINED_FLOAT', valid_range=[0,1e6]),
				option_definition.addToken(name='Input.wl.end', datatype=float, default='NOT_DEFINED_FLOAT', valid_range=[0,1e6], optional=True) ], #TODO:argument 2 should be bigger that argument 1
			parents=['uvspec'],
			non_parents=['wavelength_index'],
		)

		wavelength_step = option_definition.option(
			name='wavelength_step',
			group='spectral',
			helpstr='Set the wavelength step (in nm) in conjunction with the wavelength range.',
			documentation=['wavelength_step'],
			gui_inputs=(),
			tokens=[],
			parents=['uvspec'],
			non_parents=['wavelength_index'],
		)

		wavelength_index = option_definition.option(
			name='wavelength_index',
			group='spectral',
			helpstr='Set the wavelengths to be selected.',
			documentation=documentation['wavelength_index'],
			gui_inputs=(GUI_definition.IntegerInput(name='Input.wl.start_index'),
				    GUI_definition.IntegerInput(name='Input.wl.end_index',
						 optional=True),),
			tokens = [ option_definition.addToken(name='Input.wl.start_index', datatype=int, default='NOT_DEFINED_INTEGER',),
				option_definition.addToken(name='Input.wl.end_index', datatype=int, default='NOT_DEFINED_INTEGER', optional=True) ],
			parents=['uvspec'],
			non_parents=['wavelength'],
		)

		wavelength_grid_file = option_definition.option(
			name='wavelength_grid_file',
			group='spectral',
			helpstr='Location of file with wavelength grid used for the internal transmittance calculations.',
			documentation=documentation['wavelength_grid_file'],
			gui_inputs=(GUI_definition.FileInput(name='Input.filename[FN_WLTRANS]'),),
			tokens=option_definition.addToken(name='Input.filename[FN_WLTRANS]', datatype=io.IOBase),
			parents=['uvspec'],
			non_parents=['thermal_bands_file'],
		)

		thermal_bands_file = option_definition.option(
			name='thermal_bands_file',
			group='spectral',
			helpstr='File with the center wavelengths and the wavelength band intervals for calculations in the thermal range.',
			documentation=documentation['thermal_bands_file'],
			gui_inputs=(GUI_definition.FileInput(name='Input.filename[FN_WLBANDS]'),),
			tokens=option_definition.addToken(name='Input.filename[FN_WLBANDS]', datatype=io.IOBase),
			parents=['uvspec'],
			non_parents=['filter_function_file','slit_function_file'],
		)

		thermal_bandwidth = option_definition.option(
			name='thermal_bandwidth',
			group='spectral',
			documentation=documentation['thermal_bandwidth'],
			gui_inputs=(GUI_definition.FloatInput(name='Input.bandwidth', valid_range=[0, 1000000.0]), GUI_definition.ListInput(name='Input.bandwidth_unit', valid_range=['','nm', 'cm-1'], optional=True),),
			tokens = [ option_definition.addToken(name='Input.bandwidth', datatype=float, valid_range=[0,1e6]),
				option_definition.addLogical(name='Input.bandwidth_unit', logicals=['nm','cm-1'], setting='UNIT_PER_', optional=True) ],
			parents=['uvspec'],
		)

		source = option_definition.option(
			name='source',
			group='spectral',
			helpstr='Source of radiation.', 
			documentation=documentation['source'], 
			gui_inputs=(GUI_definition.ListInput(name='Input.source', valid_range=['thermal', 'solar']),
				    GUI_definition.FileInput(name='Input.filename[FN_EXTRATERRESTRIAL]', optional=True),
				    GUI_definition.ListInput(name='Input.spectrum_unit', valid_range=['','per_nm', 'per_cm-1', 'per_band'], optional=True),),
			tokens= [option_definition.addLogical(name='Input.source', logicals=['thermal', 'solar'], setting='SRC_'),
				 option_definition.addToken(name='Input.filename[FN_EXTRATERRESTRIAL]', datatype=io.IOBase, optional=True),
				 option_definition.addLogical(name='Input.spectrum_unit', logicals=['per_nm', 'per_cm-1', 'per_band'], setting='UNIT_', optional=True ) ], 
			parents=['uvspec'], 
			plot = {'plot_type': '2D',
				'optional_args': {'column_names': (
						"wavelength",
						"extraterrestrial flux",)
						  }
				}
		)

		mc_sun_angular_size = option_definition.option(
			name='mc_sun_angular_size',
			group='spectral',
			documentation=documentation['mc_sun_angular_size'],
			tokens=option_definition.addToken(name='Input.rte.mc.sun_radius',datatype=float),
			threedmystic=True,
			showInGui=False,
		)

		mc_lidar = option_definition.option(
			name='mc_lidar',
			group='spectral',
			helpstr='Use local estimator to simulate a lidar.',
			documentation=documentation['mc_lidar'],
			tokens=[option_definition.addSetting(name='Input.source', setting='SRC_LIDAR'),
				option_definition.addSetting(name='Input.rte.mc.escape', setting=0),
				option_definition.addSetting(name='Input.rte.mc.locest', setting='MCLIDAR_SPACE'),
				option_definition.addLogical(name='Input.rte.mc.locest', logicals=['ascope','polarize','space','falcon','simple','pilot','moon','test'], setting='MCLIDAR_', optional=True) ],
			parents=['uvspec'],
			non_parents=['source'],
			showInGui=False,
			islidar=True
		)

		mc_lidar_file = option_definition.option(
			name='mc_lidar_file',
			group='spectral',
			helpstr='File containing positions, looking directions, and opening angles of lasers and detectors for lidar simulations in MYSTIC.',
			documentation=documentation['mc_lidar_file'],
			tokens=option_definition.addToken(name='Input.rte.mc.filename[FN_MC_LIDAR]',datatype=io.IOBase),
			parents=['uvspec'],
			non_parents=['source'],
			showInGui=False,
			islidar=True
		)

		mc_radar = option_definition.option(
			name='mc_radar',
			group='spectral',
			helpstr='Switch on radar.',
			documentation=documentation['mc_radar'],
			tokens=[option_definition.addSetting(name='Input.source', setting='SRC_LIDAR'),
				option_definition.addSetting(name='Input.rte.mc.escape', setting=0),
				option_definition.addSetting(name='Input.rte.mc.locest', setting='MCRADAR')],
			parents=['uvspec'],
			non_parents=['source'],
			showInGui=False,
			islidar=True
		)

		self.options = [wavelength, wavelength_index, wavelength_grid_file,
				thermal_bands_file, thermal_bandwidth,
				source, 
				mc_lidar, mc_lidar_file, mc_radar, 
				mc_sun_angular_size, ]

	def __iter__(self):
		return iter(self.options)

def get_spectral_documentation():
	return {
		'mc_sun_angular_size' : r'''
	At the moment only useful together with \code{mc_panorama_view}.
	Set the angular radius of the sun in degrees. If omitted the radius is calculated 
	(more or less correctly) via the earth-sun distance (not well tested).
	If no \code{mc_backward_sun_shape_file} is given a spectral sun shape according to
	Koepke2001 is used.
		''',

		'mc_lidar' : r'''
	Use local estimator to simulate a lidar. If \code{mc_lidar} is set,
	you need to provide a lidar file for \code{mc_lidar_file}. A detailed
	documentation is available on request from Robert Buras.
		''',

		'mc_lidar_file' : r'''
	File containing positions, looking directions, and opening angles of 
	lasers and detectors for lidar simulations in MYSTIC.  Only meaningful 
	with \code{mc_lidar}. 
	\fcode{ 
	mc_lidar_file file 
	}
		''',

		'mc_radar' : r'''
	Switch on radar, works exactly like lidar, so use with \code{mc_lidar_file}.
        Use \code{mc_ris} to get good statistics. Has not been tested with \code{mc_vroom},
        it's recommended to switch that off. Use with \code{write_output_as_netcdf} to get
        additional reflectivity factor output in a mmclx-file (netcdf format).
		''',

		'source' : r'''
	Set the radiation source type
	\fcode{
	source  type 
	}
	where \code{type} is either \code{solar} or \code{thermal}.
	Solar radiation is per default output in W/(m2 nm) if no \code{mol\_abs\_param} is specified 
        or \code{mol\_abs\_param} options \code{crs}, \code{lowtran}, or \code{reptran} are specified.
        For all other \code{mol\_abs\_param} options
	the output is integrated over the wavelength band.
	Thermal radiation is per default output in W/(m2 cm-1), if REPTRAN is used or the bandwidth 
	is equal to 1 cm-1 (default for \code{mol\_abs\_param lowtran}). 
	Otherwise the output is the integrated flux over the
	wavenumber interval specified by \code{thermal\_bandwidth}, \code{thermal\_bands\_file},
	or by the \code{mol\_abs\_param} option (\code{kato}, \code{kato2}, \code{kato2.96}, 
	\code{fu}, or \code{avhrr\_kratz}).
	\fcode{
	source type [file] [unit]
	}
	The second argument \code{file} specifies the location of file holding the extraterrestrial spectrum.
        In general, \code{file} is required for solar calculations if \code{mol\_abs\_param} is not used.
	\code{file} is ignored if \code{mol\_abs\_param} other than \code{lowtran} oder \code{reptran} is specified.

	The file must contain two columns. Column 1 is the wavelength in nm, and column 2 
	the corresponding extraterrestrial flux. The user may freely use any units
	he/she wants for the extraterrestrial flux. The wavelength  specified grid
	defines the wavelength resolution at which results are returned. However, 
	the wavelength range is determined by \code{wavelength}. \code{file} may be 
	omitted for thermal radiation calculations (\code{source thermal}) as well as 
	\code{output_quantity transmittance} and \code{output_quantity reflectivity} calculations. If omitted, the 
	output resolution equals the internal wavelength grid which the model chooses 
	for the radiative transfer calculation.
	Comments start with \code{\#}. Empty lines are ignored.

	For some purposes it is useful to tell libRadtran the units of the spectrum.
	This can be done with the optional third argument. 
	Possible choises for \code{unit} are \code{per\_nm}, \code{per\_cm-1} or \code{per\_band}.
	If \code{unit} is set to \code{per\_nm}
	libRadtran assumes that the unit of the spectrum is W/(m2 nm), if set to \code{per\_cm-1}
	it assumes W/(m2 cm-1).

		''',
		'thermal_bandwidth' : r'''
	Specify a constant bandwidth in cm-1 for thermal calculations. 
	\fcode{
	thermal\_bandwidth value
	}
	The default is 1 cm-1.
	This option is ignored if used together with \code{mol\_abs\_param kato/kato2/kato2.96/fu/avhrr\_kratz}.
			''',				

		'thermal_bands_file' : r'''
	File with the center wavelengths and the wavelength band intervals to be used for 
	calculations in the thermal range. 
	\fcode{
	thermal\_bands\_file file
	}
	The following three columns are expected:
	center (or reference) wavelength, lower wavelength limit, upper wavelength limit [nm].
	\code{thermal\_bands\_file} defines the wavelength grid for the radiative transfer 
	calculation. The RTE solver is called for each of the wavelengths in the first column. 
	The atmospheric (scattering, absorption, etc) properties are also evaluated at these 
	wavelengths. For thermal radiation calculations, the Planck function is integrated 
	over the wavelength bands defined in the second and third columns. The result will
	therefore be a band-integrated irradiance which does only make sense when the 
	\code{source solar file} grid equals the \code{thermal\_bands\_file} grid.
				''',

		'wavelength_grid_file' : r'''
	Location of single column file that sets the wavelength grid used for the 
	internal transmittance calculations. 
	\fcode{
	wavelength\_grid\_file file
	}
	The wavelengths must be in nm.
	Do not use this option unless you know what you are doing.
	Comments start with \code{\#}. Empty lines are ignored.
				''',

		'wavelength' : r'''
	Set the wavelength range by specifying first and last wavelength in nm. 
	\fcode{
	wavelength lambda\_0 lambda\_1
	}
	The default output wavelength grid is that defined in \code{source solar file}, 
	unless \code{spline} is specified. Note that the radiative transfer calculations 
	are done on an internal grid which can be influenced with \code{wavelength\_grid\_file}
	or \code{mol\_tau\_file abs file}
				''',

		'wavelength_index' : r'''
	Set the wavelengths to be selected. To be used together with predefined wavelength grids,
	such as \code{wavelength\_grid\_file}, \code{mol\_tau\_file abs file} and particularly 
	useful in combination with the \code{mol\_abs\_param} option where often only a 
	specified number of wavelength bands is required. E.g., in combination with
	\code{mol\_abs\_param AVHRR\_KRATZ}, \code{wavelength\_index 15 15} will select wavelength
	index 15 which corresponds to channel 4, or \code{wavelength\_index 10 14} will select 
	those bands required for channel 3. Indices start from 1.
				''',
		}
