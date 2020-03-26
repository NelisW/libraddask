"""--------------------------------------------------------------------
 * $Id: cloud_options.py 3180 2015-08-24 09:34:08Z Claudia.Emde $
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
class setup_cloud_group():

	group_name = 'Water and ice clouds'

	def __init__(self):
		documentation = get_cloud_documentation()
		#uvspec function to initialize new profile

		wc_file = option_definition.option( 	#TODO: option 3D only ifthreedmystic
			name='wc_file',
			group='cloud',
			helpstr='Location of file defining water cloud properties.',
			documentation=documentation['wc_file'],
			tokens = [ option_definition.addSetting(name='isp', setting=option_definition.CaothType('wc') ),
				option_definition.addLogical( name='Input.caoth[isp].source', logicals=option_definition.ProfileType().get_valid_range(), setting='CAOTH_FROM_'), 
				option_definition.addToken( name='Input.caoth[isp].filename', datatype=io.IOBase ) ],
			parents=['uvspec'],
			childs=['wc_modify',\
				 'interpret_as_level',\
				 'wc_properties',\
				 'cloudcover'],
			plot = {'plot_type': '2D',
				'optional_args': {'column_names': (
						"altitude",
						"liquid water content",
						"effective radius",)
						  }
				}
		)

		wc_properties = option_definition.option(
			name='wc_properties',
			group='cloud',
			helpstr='Define how liquid water content and effective droplet radius are translated to optical properties. ',
			documentation=documentation['wc_properties'],
			tokens = [ option_definition.addSetting(name='isp', setting=option_definition.CaothType('wc') ),
				option_definition.addLogical( name='Input.caoth[isp].properties', logicals=['hu', 'echam4', 'mie', io.IOBase], setting='PROP_', gui_name='properties'), 
				option_definition.addLogical( name='Input.caoth[isp].interpolate' , logicals=['interpolate'], optional=True, gui_name='interpolate') ],
			parents=['wc_file'],
		)

		wc_modify = option_definition.option(	#TODO: valid_ranges for GUI!!
			name='wc_modify',
			group='cloud',
			helpstr='Modify water cloud optical properties.',
			documentation=documentation['wc_modify'],
			tokens = [ option_definition.addSetting(name='isp', setting=option_definition.CaothType('wc') ),
				option_definition.addLogical( name='id1', logicals=['gg','ssa', 'tau', 'tau550'], setting='MODIFY_VAR_', gui_name='variable' ),
				option_definition.addLogical( name='id2', logicals=['set', 'scale' ], setting='MODIFY_TYPE_', gui_name='scale/set' ),
				option_definition.addToken( name='Input.caoth[isp].modify[id1][id2]', datatype=float, gui_name='value' ) ],
			parents=['wc_file'],
			non_unique=True,
		)

		ic_file = option_definition.option( 	#TODO: option 3D only ifthreedmystic
			name='ic_file',
			group='cloud',
			helpstr='Location of file defining ice cloud properties.',
			documentation=documentation['ic_file'],
			tokens = [ option_definition.addSetting(name='isp', setting=option_definition.CaothType('ic') ),
				option_definition.addLogical( name='Input.caoth[isp].source', logicals=option_definition.ProfileType().get_valid_range(), setting='CAOTH_FROM_' ),
				option_definition.addToken( name='Input.caoth[isp].filename', datatype=io.IOBase ) ],
			parents=['uvspec'],
			childs=['ic_modify','ic_habit',\
				 'ic_fu', 'interpret_as_level',\
				 'ic_properties',\
				 'cloudcover'],
			plot = {'plot_type': '2D',
				'optional_args': {'column_names': (
						"altitude",
						"ice water content",
						"effective radius",)
						  }
				}
		)

		ic_properties = option_definition.option(
			name='ic_properties',
			group='cloud', 
			helpstr='Define how ice water content and effective droplet radius are translated to optical properties. ',
			documentation=documentation['ic_properties'],
			tokens = [ option_definition.addSetting(name='isp', setting=option_definition.CaothType('ic') ),
                                   option_definition.addLogical( name='Input.caoth[isp].properties', logicals=['fu', 'echam4', 'yang', 'key', 'baum', 'baum_v36', 'hey', 'yang2013', io.IOBase], setting='PROP_', gui_name='properties' ), 
                        option_definition.addLogical( name='Input.caoth[isp].interpolate' , logicals=['interpolate'], optional=True, gui_name='interpolate') ],
			parents=['ic_file'],
			childs=['ic_habit', 'ic_habit_yang2013'],
			continious_update=True,
		)

		ic_modify = option_definition.option(	#TODO: valid_ranges for GUI!!
			name='ic_modify',
			group='cloud',
			helpstr='Modify ice cloud optical properties.',
			documentation=documentation['ic_modify'],
			tokens = [option_definition.addSetting(name='isp', setting=option_definition.CaothType('ic') ),
				option_definition.addLogical( name='id1', logicals=['gg','ssa', 'tau', 'tau550'], setting='MODIFY_VAR_', gui_name='variable' ),
				option_definition.addLogical( name='id2', logicals=['set', 'scale' ], setting='MODIFY_TYPE_', gui_name='scale/set' ),
				option_definition.addToken( name='Input.caoth[isp].modify[id1][id2]', datatype=float, gui_name = 'value' ) ],
			parents=['ic_file'],
			non_unique=True,
		)

		ic_habit = option_definition.option( 
			name='ic_habit',
			group='cloud',
			helpstr='Ice crystal habit for the \citet{Yang2000}, \citet{Key2002}, \code{baum_v36}, \code{hey} parameterizations.',
			documentation=documentation['ic_habit'],
			tokens = [ option_definition.addSetting(name='isp', setting=option_definition.CaothType('ic') ),
				option_definition.addLogical( name='Input.caoth[isp].habit',  
					logicals=['solid-column', 'hollow-column', 'rough-aggregate', 'rosette-4', 'rosette-6', 
					          'plate', 'droxtal', 'dendrite', 'ghm'], setting='IC_HABIT_' ) ],
			parents=['ic_file'],
			speaker='ic_properties',
			enable_values=('key','yang','hey','baum_v36'),
		)

		ic_habit_yang2013 = option_definition.option( 
			name='ic_habit_yang2013',
			group='cloud',
			helpstr='Ice crystal habit for the \citet{yang2013} parameterization.',
			documentation=documentation['ic_habit_yang2013'],
			tokens = [ option_definition.addSetting(name='isp', setting=option_definition.CaothType('ic') ),
                                   option_definition.addLogical(name='Input.caoth[isp].habit',  
                                              logicals=['column_8elements', 'droxtal', 'hollow_bullet_rosette', 'hollow_column', 'plate', 'plate_10elements',
                                                        'plate_5elements', 'solid_bullet_rosette', 'solid_column'], setting='IC_HABIT_' ),
                                   option_definition.addLogical(name='Input.caoth[isp].roughness',
                                              logicals=['smooth', 'moderate', 'severe'], setting='IC_ROUGHNESS_') ],
                        parents=['ic_file'],
			speaker='ic_properties',
			enable_values=('yang2013',),
		)

                
		ic_fu = option_definition.option(	
			name='ic_fu',
			group='cloud',
			documentation=documentation['ic_fu'],
			tokens = [ option_definition.addSetting( name='isp', setting=option_definition.CaothType('ic') ),
				option_definition.addLogical( name='id1', logicals=['reff_def', 'deltascaling'], setting='IC_FU_' ),
				option_definition.addLogical( name='id2', logicals=['on', 'off' ], setting='SWITCH_' ),
				option_definition.addSetting( name='Input.caoth[isp].ic_fu[id1]', setting='id2' ) ],
			parents=['ic_file'],
			non_unique=True,
		)

		cloud_fraction_file = option_definition.option(
			name='cloud_fraction_file',
			group='cloud',
			documentation=documentation['cloud_fraction_file'],
			tokens=[option_definition.addToken( name='Input.filename[FN_CLOUD_FRACTION]', datatype=io.IOBase ),],
		)

		cloud_overlap = option_definition.option(
			name='cloud_overlap',
			group='cloud',
			documentation=documentation['cloud_overlap'],
			tokens=[option_definition.addLogical( name='Input.cloud_overlap', logicals=['rand', 'max', 'maxrand', 'off' ], setting='CLOUD_OVERLAP_' )],
		)

		cloudcover = option_definition.option(
			name='cloudcover',
			group='cloud',
			helpstr='Fraction of the horizontal sky area which is covered by clouds.',
			documentation=documentation['cloudcover'],
			tokens = [ option_definition.addToken(name='isp', datatype=option_definition.CaothType, valid_range=['wc','ic'] ),
				option_definition.addToken( name='Input.caoth[isp].cloudcover', datatype=float, valid_range=[0,1] ) ,
				option_definition.addSetting( name='Input.caoth[isp].ipa', setting=1 ),
				option_definition.addSetting( name='Input.ipa', setting=1, default=0 ) ],
			parents=['ic_file','wc_file'],
			non_unique=True,
		)



		self.options = [wc_file, wc_properties, wc_modify, 
				ic_file, ic_properties, ic_modify, 
				ic_habit, ic_habit_yang2013, ic_fu,
				cloud_fraction_file, cloud_overlap, cloudcover]
		
	def __iter__(self):
		return iter(self.options)

def get_cloud_documentation():
	return {
		'wc_file'	: r'''
	Location of file defining water cloud properties. 
	\fcode{
	wc\_file type file
	}
	\code{type} defines the file type, which can be one of the
	following:
	\begin{description}

	\item[1D] Location of file defining one-dimensional profile.

	The file must contain three columns: Column 1 is the altitude in km,
	column 2 the liquid water content (LWC) in grams per cubic meter, and
	column 3 the effective droplet radius in micrometer. Empty lines are
	ignored. Comments start with \code{\#}.  Note that the definition of
	cloud altitudes in \code{wc\_file 1D} refers to sea level, not to altitude
	above ground. E.g., when \code{altitude} is set to 1.63km, and the
	first cloud level is defined at 3km, the cloud would start at 1.37km
	above ground.  An example of a water cloud is given in
	\file{examples/WC.DAT}.
	
	Per default the cloud properties are interpreted as layer
	properties. Before version 1.4 the default was level properties: The
	optical depth of a layer was calculated using information from the
	upper and lower levels defining the layer.  
	To switch to the old behaviour, use \code{interpret\_as\_level}. 
	See section~\ref{seq:wc} about water clouds for a
	realistic example how the contents of the \code{wc\_file 1D} are converted
	to optical properties.
	
	\ifthreedmystic{

	\item[3D] Define a MYSTIC 3D water cloud input file. 
	In case the \code{rte\_solver} is not \code{mystic},
	independent column calculations are performed for the 3D field.
	(See \file{examples/UVSPEC\_WC\_IPA.INP} for an example.)
	The expected format of the cloud file is: 
	\fcode{
	Nx  Ny  Nz flag \\
	dx  dy  z(1) z(2) ... z(n) \\
	ix  iy  iz  ext  g     ssa  (if flag == 1)\\
	ix  iy  iz  ext  reff       (if flag == 2)\\
	ix  iy  iz  LWC  reff       (if flag == 3)\\
	}
	where \code{Nx}, \code{Ny} and \code{Nz} are the number of grid boxes in 
	\code{x}, \code{y}, and \code{z}-direction.
	The parameter \code{flag} determines the format of the 3rd and 
	following lines. In the second line \code{dx} and \code{dy} are the sizes of the boxes in x- 
	and y-direction in km. 
	In the third and following lines the indices \code{ix}, \code{iy}, and \code{iz} specify cloudy pixels. 
	The optical properties of the cloud, are given by the other parameters in the line, 
	where \code{ext} is the extinction coefficient [1/km], \code{g} the asymmetry parameter,
	\code{reff} the effective radius [micrometer], and \code{ssa} the single scattering albedo.
	The conversion from microphysical to optical properties is defined by \code{wc\_properties} 
	(identical to the 1D case). Note that the first dimension (x)
	propagates east, and the second dimension (y) north. For more
	information see section~\ref{sec:mystic}.

%	\item[ipa] Define a 3D input file.
%	Independent column calculations are performed for the 3D field.
%	As argument a name of a 3D cloud file must be given. This file has to
%	be in the format as needed by MYSTIC, see \code{wc\_file 3D}.
%	(See \file{examples/UVSPEC\_WC\_IPA.INP} for an example.)
	}
	
	\item[ipa\_files] A two-column file, defining water cloud property files (see \code{wc\_file 1D}) in the first 
	column and the correspoding weights in the second column. 
	The radiative transfer calculation is performed independently for each
	cloud column and the result is the weighted average of all independent
	columns. If \code{ic\_file ipa\_files} and \code{wc\_file ipa\_files} are both
	defined, both must have the same columns in the same order, otherwise
	\code{uvspec} will complain. See
	\file{examples/UVSPEC\_WC\_IC\_IPA\_FILES.INP} for an example.
	
	\item[moments] A way to specify water cloud extinction coefficient, single
	scattering albedo, and scattering phase function for each layer.

	The file specified by 
	\code{wc\_file moments} has two columns where column 1 is the altitude in km. The second
	column is the name of a file which defines the optical properties of the layer 
	starting at the given altitude. The files specified in the second column must 
	have the following format:
	\begin{description}
	\item{Column 1: } 
	The wavelength in nm. These wavelengths may be different from those in \code{source solar filename}. 
	Optical properties are interpolated to the requested wavelengths.
	\item{Column 2:} 
	The extinction coefficient of the layer in units km-1. 
	\item{Column 3:} 
	The single scattering albedo of the layer.
	\item{Column 4-(nmom+4):} 
	The moments of the scattering phase function. 
	\end{description}
	Note that if using the \code{rte\_solver disort} or \code{rte\_solver
	fdisort2} it makes good sense to make the number of moments larger
	than \code{number\_of\_streams} because all moments are used in the calculation. For
	\code{rte\_solver fdisort1} and \code{rte\_solver polradtran} the number
	of moments included in the calculations will be \code{number\_of\_streams}+1. Higher
	order moments will be ignored for these solvers.  Please note that the
	uppermost line of the \code{wc\_file moments} denotes simply the top altitude
	of the uppermost layer. The optical properties of this line are
	consequently ignored. There are two options for this line: either an
	optical property file with zero optical thickness is specified
	or "NULL" instead.

	\end{description}
		''',

		'wc_properties'	: r'''
	Define how liquid water content and effective radius are translated to optical properties.
	\fcode{
	wc\_properties property [interpolate]
	}
	Possible choices for \code{property} are:
	\begin{description}
	\parameter{hu} 
	Parameterization by \citet{Hu1993}; 
	this is the default setting. Note that 
	the parameterization is somewhat different for \code{mol\_abs\_param FU} than for all other 
	cases because in the latter case the parameterization from the newer (March 2000) 
	Fu and Liou code is used while otherwise the data are taken from the original 
	paper by \citet{Hu1993}. Note that this parameterization has been
	developed to calculate irradiances, hence it is less suitable for radiances.
	This is due to the use of the Henyey-Greenstein phase function as an approximation
	of the real Mie phase function.
	\parameter{echam4} 
	Use the very simple two-band parameterization of the ECHAM4 climate model, described 
	in \cite{Roeckner1996}; this is probably only meaningful if you want to compare
	your results with ECHAM4, the two bands are 0.2 - 0.68 micrometer and 0.68 - 4.0 micrometer;
	within these bands, the optical properties are assumed constant.
	\parameter{mie} 
	Use pre-calculated Mie tables; useful for \code{mol\_abs\_param}; 
	the tables are expected in \code{data\_files\_path}\file{/correlated\_k/.../}. \\
	For spectral or pseudo-spectral (\code{mol\_abs\_param sbdart}) calculations, 
	a set of pre-calculated tables is also available.
	%CE: commented some parts. Should be rewritten when new Mie tables are well
	% optimized.
	% the wavelength grid points of these 
	%data has been carefully selected such that the extinction cross section, 
	%single scattering albedo, and the asymmetry parameter are accurate to 1\% 
	%(compared to the fully-resolved Mie calculation) for all wavelengths
	%between 250nm and 100 micrometer. 
	For spectral or pseudo-spectral
	calculations optional argument \code{interpolate} has to be defined explicitely to
	initiate the interpolation of the optical properties to the internal wavelength grid.
	%Please note that this option may be extremely memory-consuming because for each 
	%internal wavelength a full set of Legendre moments of the phase function is 
	%stored (up to several thousands). 
	The Mie tables are not part of the standard distribution 
	(because of their large size) but they are freely available from http://www.libradtran.org. 
	This is the correct option to calculate radiances, to be preferred over the 
	Henyey-Greenstein approach of \citet{Hu1993}.
	\parameter{filename} 
	Read optical properties from specified filename; file format is as produced 
	by the \code{mie}-tool of the {\sl libRadtran} package (see \code{output\_user netcdf}).
	\end{description}

	With the optional argument \code{interpolate} the water cloud optical properties 
	are interpolated over wavelength; useful for precalculated optical property files. 
	Please note that this option may be extremely memory-consuming because for each internal wavelength 
	a full set of Legendre moments of the phase function is stored (up to several thousands). 
		''',

		'wc_modify'	: r'''
	Modify water cloud optical properties.
	\fcode{
	wc\_modify variable scale/set value
	}
	\code{variable} can be one of the following parameter:
	\begin{description}

	\parameter{gg}
	Modify the water cloud asymmetry factor for all wavelengths and altitudes.
	\begin{description}
	\item[set]
	\code{value} can be a float between -1.0 and 1.0. 
	Please note that this option is only applied if a Henyey-Greenstein 
	phase function is used but not if an explicit phase function is 
	defined e.g. with \code{wc\_file moments}. It doesn't make sense to modify only 
	the first moment of an explicit phase function.  
	This option is useful only for monochromatic 
	calculations or in wavelength regions where the optical 
	properties of water clouds can be considered constant, 
	e.g. the ultraviolet range. 
	\item[scale]
	Scale the water cloud asymmetry factor for all wavelengths and altitudes
	with \code{value} between 0.0 and 1.0.
	\end{description}
	
	\parameter{ssa}
	Modify the water cloud single scattering albedo for all wavelengths
	and altitudes.
	\begin{description}
	\item[set] \code{value} can be a float between 0.0 and 1.0. 
	This option is useful only for monochromatic 
	calculations or in wavelength regions where the optical properties of water clouds 
	can be considered constant, e.g. the ultraviolet range.
	\item[scale]
	Scale the water cloud single scattering albedo for all wavelengths and altitudes
	with \code{value} between 0.0 and 1.0.
	\end{description}

	\parameter{tau}
	Modify the total water cloud optical thickness.
	\begin{description}
	\item[set] Set optical thickness to a constant value for all wavelengths.
	The optical thickness defined here is the integral from the surface at the 
	user-defined \code{altitude} to TOA (top of atmosphere). This option is useful only 
	for monochromatic calculations or in wavelength regions where the optical properties 
	of water clouds can be considered constant, e.g. the ultraviolet range.
	\item[scale] Scale the water cloud optical thickness for all wavelengths and altitudes
	with \code{value} between 0.0 and 1000000.0. Also works for 3d clouds.
	\end{description}

	\parameter{tau550}
	Set the water cloud optical thickness at 550nm. 
	The optical thickness defined here 
	is the integral from the surface at the user-defined \code{altitude} 
	to TOA (top of atmosphere). Other wavelengths are scaled accordingly.
	Note that this option requires for technical reasons that the
	wavelength interval defined by \code{wavelength} does contain 550nm.
	\end{description}
	\fcode{
	wc\_modify tau550 set value
	}
			''',

		'ic_file'	: r'''
	Location of file defining ice cloud properties. 
	\fcode{
	ic\_file type file
	}
	
	\code{type} defines the file type, which is identical to \code{wc\_file type}.
	See \code{wc\_file} for choices of \code{type} and a description on the file structures.
		''',

		'ic_properties'	: r'''
	Defines how ice water content and effective particle radius are translated 
	to optical properties. 
	\fcode{
	ic\_properties property [interpolate]
	}
	Possible choices for \code{property} are
	\begin{description}
	\parameter{fu}
	Parameterization by \citet{Fu1996,Fu1998}, see \code{ic\_file}; 
	this is the default setting. Note that this is a parameterization
	which has been created to calculate fluxes but not radiances. 
	Note also that the optical properties in the solar range provided by 
	\citet{Fu1996} are delta-scaled properties (that is, the forward peak of
	the phase function is truncated and optical thickness, asymmetry
	parameter, and single scattering albedo are reduced accordingly),
	whereas \code{uvspec} uses non delta-scaled properties unless the option
	\code{ic\_fu deltascaling on} is specified. By default the parameterization
	by \citet{Fu1996} is treated consistently with all other ice cloud
	parameterizations. 
	For wavelengths up to 4 micrometer \citet{Fu1996} is used while for wavelengths 
	larger than 4 micrometer \citet{Fu1998} is chosen. Please note that 
	\citet{Fu1996} is based on ray-tracing calculations while \citet{Fu1998}
	is a mixture of ray-tracing and Mie calculations (which is required for 
	the infrared wavelengths where the geometrical assumption does not hold).
	Hence, both parameterizations are not fully consistent. Rather, differences
	of some \% are to be expected in the wavelength region where both 
	parameterizations overlap. Also, the wavelength dependence in the solar 
	and infrared parts is treated differently: In the solar part \citep{Fu1996}
	the optical properties are defined for wavelength bands - hence they 
	are assumed constant within each band. In the infrared \citep{Fu1998} 
	they are defined at certain wavelengths and linearely interpolated 
	in between. If you use this option, please see also the
	discussion of \code{ic\_fu deltascaling} and \code{ic\_fu reff\_def}.
	The allowed range for the effective radius is from 9.315 -  65.120 micrometer. 
	\parameter{echam4}
	Use the simple two-band parameterization of the ECHAM4 climate model, described 
	in \cite{Roeckner1996}; this is probably only meaningful if you want to compare
	your results with ECHAM4, the two bands are 0.2 - 0.68 micrometer and 0.68 - 4.0 micrometer.
	Within the two ECHAM4 bands, the optical properties are assumed constant.
	\parameter{key}
	Parameterization by \citet{Key2002}. This parameterization can also
	be used to calculate radiances because it uses a
	double-Henyey-Greenstein phase function which better represents both
	forward and backward peaks. This parameterization covers the wavelength region  
	from 0.2 to 5.0 micrometer and is available for the following \code{habit}:
	\code{solid-column}, \code{hollow-column}, \code{rough-aggregate}, \code{rosette-4}, 
	\code{rosette-6}, and \code{plate}.
	\parameter{yang}
	Parameterization similar to \citet{Key2002} but based on more recent
	single scattering calculations. Below 3.4 micrometer it actually equals
	the \citet{Key2002} parameterization while from 3.4 - 100 micrometer new 
	coefficients have been calculated with much higher wavelength
	resolution and better accuracy. Hence, \code{yang} should give a reasonably
	consistent approximation from 0.2 - 100 micrometer, suitable for spectrally
	resolved calculations of radiance and irradiance.
	The covered range for the effective radius depends on the \code{ic\_habit}.
	(In micrometer: \code{solid-column} [5.96, 84.22], \code{hollow-column} [4.97, 70.24],
	 \code{rough-aggregate} [3.55, 108.10], \code{rosettes-4} [2.77, 45.30], 
	\code{rosettes-6} [2.85, 46.01], \code{plate} [4.87, 48.18], \code{dendrites}
	 [0.45, 1.88], \code{droxtal} [9.48, 293.32], 
	\code{spheroid}  [6.58, 203.39]). 
	\parameter{baum}
	Use ice cloud parameterization from \citet{baum05a:_bulk,baum05b:_bulk},
	\url{http://www.ssec.wisc.edu/\~baum/Cirrus/IceCloudModels.html}.
	In combination with the radiative transfer solvers \code{disort}, \ifmystic{\code{montecarlo},}
	 and \code{fdisort2}, accurate phase functions are used.
        \parameter{baum_v36} Use cloud parameterization from
	\citet{heymsfield2013,yang2013,baum2014} covering the 
	spectral range from 0.2 to 100~$\mu$m and effective radii from
	5 to 60 $\mu$m. The parameterization
	assumes severly roughened and randomly oriented ice particles
	and includes the full phase matrix. Three set of models are
	available, they can be selected using \code{ic\_habit}:
	\code{ghm} is based on a general habit mixture involving 9
	habits, \code{solid-column} assumes severely roughened solid
	columns, and \code{rough-aggregate} is based on severly
        roughened aggregates. The default is \code{solid-column}. 
	%\parameter{baum\_hufit }
	%Similar to the option \code{baum} but here the phase function
	%is parameterized by 128 Legendre coefficients, calculated with the
	%delta-fit method from 
	%\citet{Hu2000}. This parameterization covers the region  
	%from 0.4 to 2.2 micrometer. If high accuracy is needed e.g. in the vicinity of the halo, 
	%the forward peak, or the backscatter peak, \code{ic\_properties baum}
	%is recommended.
	\parameter{hey}
	Use pre-calculated ice cloud optical properties including full phase
	matrices.
	The parameterization
	is currently only available for the spectral region from 0.2 to 5
	micrometers. The single scattering
	properties have been been generated by Hong Gang using the
	models by \citet{Yang2000}. The parameterization is based on
	simple gamma distributions 
	\begin{equation}
	n(r) = n_0 r^{\alpha} \exp\left(-\frac{(\alpha+3)r}{r_e}\right),
	\end{equation}
	where $n_0$ is found by normalization and $\alpha$ is set to 1. In case of
	spherical particles the parameter $r_e$ would be the effective
	radius. For aspherical particles, the parameter $r_e$ is found
	iteratively so that the size distribution yields the required
	effective radius. The
	parameterization is availabe for the following habits: \code{solid-column},
	\code{hollow-column}, \code{rough-aggregate}, \code{rosette-6},
	\code{plate}, and \code{droxtal}. The
	default habit is \code{solid-column}. Furthermore a general
        habit mixture \code{ghm} similar to the one defined by
	\citet{baum05b:_bulk} may be selected.  For the HEY
	parameterization the ice crystals are assumed to be smooth, in
	contrast to the severely roughened particles assumed by
	\code{baum_v36}. The habit can be specified using the option
	\code{ic\_habit}.
        \parameter{yang2013} Pre-calculated ice optical properties
	including full phase matrices for the spectral region from 0.2
	to 99 $\mu$m. The parameterization uses single scattering data by
	\citet{yang2013} and assumes gamma size distributions as the
	\code{hey} parameterization. It is available for 9 habits and
	3 degrees of roughness, which can be selected using the option
	\code{ic\_habit\_yang2013}. 
	%\parameter{ic-mie}
	%Use pre-calculated Mie tables; useful for \code{mol\_abs\_param}; 
	%the tables are expected in \code{data\_files\_path}\file{/correlated\_k/.../}.
	%For spectral or pseudo-spectral (\code{mol\_abs\_param sbdart}) calculations, 
	%a set of pre-calculated tables is also available; the wavelength grid points of these 
	%data has been carefully selected such that the extinction cross section, 
	%single scattering albedo, and the asymmetry parameter are accurate to 1\% 
	%(compared to the fully-resolved Mie calculation) for all wavelengths
	%between 250nm and 100 micrometer. 
	%For spectral or pseudo-spectral
	%calculations the optional argument \code{interpolate} has to be defined explicitely to
	%initiate the interpolation of the optical properties to the internal wavelength grid.
	%Please note that this option may be extremely memory-consuming because for each 
	%internal wavelength a full set of Legendre moments of the phase function is 
	%stored (up to several thousands). The Mie tables are not part of the standard distribution 
	%(because of their large size) but they are freely available from http://www.libradtran.org. 
	%Note that a Mie calculation assumes spherical ice particles, the scattering function of 
	%which differs systematically from non-spherical particles. Hence, \code{ic\_properties mie}
	%is usually not representative of natural ice clouds.
	\parameter{filename}
	Read optical properties from specified filename; file format is as produced 
	by the \code{mie} tool of {\sl libRadtran} (see \code{output\_user netcdf}). 
	% or by Frank Evans' \code{cloudprp}.
	\end{description}
	The default property is \code{fu}.

	Please note also that, in contrast to spherical particles, there is no unique 
	definition of effective size for non-spherical particles. In particular, the
	above parameterizations use different definitions which, however, differ only by 
	a constant factor. 
	\citet{Yang2000}, \cite{Key2002}, and \citet{baum05a:_bulk,baum05b:_bulk}
	 use the general definition 
	\begin{equation}
	   r_{\rm eff} = {{3}\over{4}}{{\int V(h) n(h) dh}\over{\int A(h) n(h) dh}}
	\end{equation}
	  where $h$ is the maximum dimension of an ice crystal, $n(h)$ is the
	  number of particles with maximum dimension $h$ in the size distribution,
	  and $V$ and $A$ are the volume and mean
	  projected area of the particles, respectively. The volume and area are
	  based on the spherical diameter with equivalent volume and the
	  spherical diameter with equivalent projected area as defined by 
	  \citet{Yang2000}. On the other hand, \citet{Fu1996,Fu1998} use
	  hexagonal columns and use the following definition
	\begin{equation}
	   r_{\rm eff} =  {{\int D^2 L n(L) dL}\over{2 \int (D L + {\sqrt{3}\over{4}} D^2) n(L) dL}}
	\end{equation}
	  where $D$ is the width of the ice crystal (that is, the maximum diameter of the 
	  hexagonal area) and $L$ is the length. The integrand in the numerator is proportional 
	  to the volume while that in the denominator is proportional to the projected area. 
	  Evaluating these formulas one finds that, for the same hexagonal particle, the effective
	  radius would be $3 \sqrt{3} / 4 = 1.299$ times larger following the 
	  \citet{Yang2000}, \citet{Key2002} definition rather than the \citet{Fu1996,Fu1998} definition.
	  As an example, an effective radius of 20$\mu m$ with 
	  \code{ic\_properties fu} and \code{ic\_fu reff\_def on} and
	  1.299 $\cdot$ 20$\mu m$ = 26$\mu m$ with \code{ic\_properties yang} would give comparable results
	  for hexagonal columns. 
	To use the original definition of the effective radius by \citet{Fu1996,Fu1998} use
	\code{ic\_fu reff\_def on}!

	With the optional argument \code{interpolate} the ice cloud optical properties are interpolated 
	over wavelength; useful for precalculated optical property files defined with \code{ic\_properties}. 
	Please note that this option may be extremely memory-consuming because for each internal wavelength 
	a full set of Legendre moments of the phase function is stored (up to several thousands). 
		''',

		'ic_modify'	: r'''
	Modify ice cloud optical properties.
	\fcode{
	ic\_modify variable scale/set value
	}
	This option is identical to \code{wc\_modify}.
	Please refer to \code{wc\_modify} for a detailed description of \code{variable}.
	
	If you use this option in combination with the ice
	cloud properties by \citet{Fu1996}, please make sure that you understand
	the explanation of \code{ic\_fu}.
		''',

		'ic_habit'	: r'''
	Ice crystal habit for the \citet{Yang2000}, \citet{Key2002} and
	\code{hey} parameterizations, see also  \code{ic\_properties key/yang/hey}.
	\fcode{
	ic\_habit type
	} 
	For Key/Yang \code{type} may be one of \code{solid-column}, 
	\code{hollow-column}, \code{rough-aggregate}, \code{rosette-4}, 
	\code{rosette-6}, \code{plate}, \code{droxtal}, \code{dendrite}
	and \code{spheroid}. Please note that this parameterization is only valid for
	a restricted size range, depending on the habit (see table 1 in 
	\citet{Key2002}. Also, some of the habits are only available for
	wavelengths below 5 micrometer (\code{rosette-4}) while others are only available
	for wavelengths larger than 3 micrometer (\code{droxtal}, \code{spheroid}).

	For \code{hey} the following habits can be chosen: \code{solid-column}, 
	\code{hollow-column}, \code{rough-aggregate}, \code{rosette-6}, 
	\code{plate}, \code{droxtal}, and the general habit mixture \code{ghm} 
        which follows the ``recipe'' by \citet{baum05a:_bulk}. All
	habits and the habit mixture are available for effective radii from 5 to 90 micrometers in
	the wavelength region from 0.2 to 5 micrometers.

        The parameterization \code{baum\_v36} includes the general habit mixture
        \code{ghm} and the habits \code{solid-column} and \code{aggregate}, 
        all crystals modeled with severe roughness.          
		''',

		'ic_habit_yang2013'	: r'''
	Ice crystal habit for the \citet{yang2013} parameterization,
        \code{ic\_properties yang2013}.
	\fcode{
	ic\_habit type roughness
	} 
	The following habits are available: \code{column\_8elements}, 
        \code{droxtal}, \code{hollow\_bullet\_rosette}, \code{hollow\_column}, 
        \code{plate}, \code{plate\_10elements}, \code{plate\_5elements}, \code{solid\_bullet\_rosette}, 
        and \code{solid\_column}.

        For each habit three degrees of roughness needs to be specified, options are \code{smooth}, 
        \code{moderate} and \code{severe}.         
        
		''',        

		'ic_fu' : r'''
	\fcode{
	ic\_fu reff\_def on/off
	}
	Specify wich definition of the effective radius is used.
	
	If \code{on} the parameterization uses the original definition of the effective radius
	as specified in \citet{Fu1996,Fu1998}. 

	Default is \code{off}. The same definition of the effective radius is used as the \citet{Key2002},
	\citet{Yang2000} and \citet{baum05a:_bulk,baum05b:_bulk}
	parameterizations; see discussion of \code{ic\_properties}.

	\fcode{
	ic\_fu deltascaling on/off
	}
	Specify if the \citet{Fu1996} optical properties are delta-scaled or not. 

	If \code{on} delta-scaling is switched on.

	If  \code{off} delta-scaling is switched off. The default is without delta-scaling.
	Please note that this was changed on July 22, 2008: 
	Before, delta-scaling was switched on by default which might have
	caused some confusion, because irradiance calculations were not
	consistent with the other ice cloud parameterizations implemented in
	\code{uvspec}. 
	Using the \citet{Fu1996} parameterization in combination with one of 
	\code{ic\_modify} you now get
	consistent results with all other ice cloud parameterizations.
	
	%THIS IS THE OLD DOCUMENTATION (BEFORE CHANGE OF DELTA-SCALING
	%				      DEFAULT).
	%
	%It has been confirmed that the difference with and 
	%without delta-scaling is typically on the order of less than 1 - 2\% 
	%and therefore it was decided to switch delta-scaling off by default.
	%Please note that this has nothing to do with the internal delta-scaling
	%by the solvers: e.g. disort2 always applies it's own internal delta-scaling 
	%anyway. (THE FOLLOWING NEEDS TO BE MODIFIED).
	%If you define a cloud only by its microphysical properties (ice water
	%content, effective radius), delta-scaling should certainly be switched
	%on and you do not need to read further. 
	%If, however, you want to
	%use the Fu (1996) parameterization in combination with one of 
	%\code{ic_set_tau/tau550/gg/ssa} or \code{ic_scale_gg/ssa} it might be
	%reasonable to switch delta-scaling off and you should make sure that
	%you understand the following. Citing from Fu (1996): "For nonspherical
	%particles in cirrus clouds, it is found that a simple representation
	%of the scattering phase function through the asymmetry factor is
	%inadequate (Fu and Takano 1994). As demonstrated in appendix A, the
	%fraction of scattered energy residing in the forward peak, f, needs to
	%be removed from the scattering parameters to incorporate the strong
	%forward peak contribution in multiple scattering." Or in other words,
	%the sharp forward peak is truncated and added to the unscattered
	%direct radiation. The remaining phase function (excluding the sharp
	%forward peak) can be safely approximated by a Henyey-Greenstein
	%function. The scaling implies a reduction of the optical thickness,
	%the asymmetry parameter, and the single scattering albedo. This
	%reduction can be rather severe, e.g. a factor of about 3 for the
	%optical thickness in the visible spectral range. This implies
	%seemingly inconsistent optical properties: For idential IWC content
	%and effective radius, \code{ic_properties key/yang} would give an
	%(unscaled) optical thickness about three times higher than
	%\code{ic_properties fu}. The effect on the radiation field, however,
	%will be comparable, due the consistent scaling of optical thickness,
	%asymmetry parameter, and single scattering albedo. If you, however,
	%adjust the optical thickness using e.g. \code{ic_set_tau}, the effect
	%on the radiation field will be completely different because the
	%(unscaled) optical thickness by Key (2002) has a completely different
	%meaning as the (scaled) optical thickness by Fu (1996). In such cases
	%it might be reasonable to switch scaling off. This is a complicated
	%and confusing topic and it is suggested that you play around a bit
	%with the options, read the Fu (1996) paper, and make heavy use of the
	%\code{verbose} feature.
		''',

                'cloud_fraction_file' : r'''
        File containing a cloud fraction profile. 
        \fcode{
        cloud\_fraction\_file file
        }
        Two columns are expected:
        altitude [km] and cloud fraction, including ice and water clouds. If \code{cloud\_fraction\_file}
        is defined, effective cloud properties are calculated assuming either random overlap or maximum random overlap
        of the cloud layers (see also \code{cloud\_overlap}). An example is provided in \file{examples/CF.DAT}.
                ''',

                'cloud_overlap' : r'''
        Cloud overlap assumption. 
        \fcode{
        cloud\_overlap type
        }
        Following types are implemented:
        \begin{description}
        \item[rand] Random overlap of cloud layers 
        \item[maxrand]  Maximum random overlap scheme
        \item[max]      Maximum overlap scheme
        \item[off]      Turn off cloud overlap for ECMWF clouds                                        
        \end{description}
        Per default the \code{cloud\_overlap} scheme is switched off.
                ''',

                'cloudcover'    : r'''
        Set the fraction of the horizontal sky area which is covered by clouds.
        \fcode{
        cloudcover typename value
        }
        \code{typename} describes the name of the cloud type, which can be "wc" and "ic":

        \begin{description}
        \item[ic] Cloud cover of ice cloud, where the cloud properties are taken
        from \code{ic\_file}.
        \item[wc] Cloud cover of water cloud, where the cloud properties are taken 
        from \code{wc\_file}.
        \end{description}

        When a cloud cover is specified, the result will be calculated by
        the independent pixel approximation (IPA), that is, as weighted average 
        of cloudless sky and overcast sky. 
        Please note that, if both \code{cloudcover ic} and
        \code{cloudcover wc} are set, both must be equal.
        
        This option is ignored, if the option \code{cloud\_fraction\_file} is used.
                ''',

	}
