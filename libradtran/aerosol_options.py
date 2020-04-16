"""--------------------------------------------------------------------
 * $Id: aerosol_options.py 3139 2015-07-15 11:15:57Z josef.gasteiger $
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

# Modified by Derek Griffith to include aerosol_sizedist_file and aerosol_refrac_index


from . import option_definition
from . import GUI_definition
import io

class aerosol_mixed_options(option_definition.option):
    """
    Logic for haze, season, vulcan and visibility
    """

    def __init__(self, *args, **kwargs):
        option_definition.option.__init__(self, *args, **kwargs)
        self.dependencies.extend(['aerosol_haze', 'aerosol_vulcan',
                                  'aerosol_season', 'aerosol_visibility'])

    def isMandatory(self, is_set, get_value):
        cond = [is_set(opt) for opt in ['aerosol_haze',
                                       'aerosol_vulcan', 'aerosol_season',
                                       'aerosol_visibility']]
        if all(cond):
            return False
        elif any(cond):
            return True
        else:
            return False


class setup_aerosol_group():
    group_name = 'Aerosol'

    def __init__(self):
        documentation = get_aerosol_documentation()

        aerosol_default = option_definition.option(
            name='aerosol_default',
            group='aerosol',
            helpstr='Set up a default aerosol',
            documentation=documentation['aerosol_default'],
            tokens=[option_definition.addSetting(name='Input.aer.standard', setting=1)],
            childs=['aerosol_angstrom', 'aerosol_king_byrne'],
            continious_update=True
        )

        aerosol_file = option_definition.option(
            name='aerosol_file',
            group='aerosol',
            helpstr='Location of aerosol file',
            documentation=documentation['aerosol_file'],
            tokens=[option_definition.addLogical(name='id', logicals=['gg', 'ssa', 'tau', 'explicit', 'moments'], setting='FN_AER_'),
                    option_definition.addToken(name='Input.aer.filename[id]', datatype=io.IOBase),
                    option_definition.addSetting(name='Input.aer.spec', setting=1)],
            non_unique=True,
        )

        #	aerosol_shettle = option_definition.option(
        #		'aerosol_file',	#name
        #	        'aerosol',      #group
        #		'Specify aerosol properties in atmosphere',	#helpstr
        #		documentation['aerosol_shettle'],     #documentation
        #		[ option_definition.addLogical( 'id', ['vulcan', 'haze', 'season'], 'AER_SHETTLE_' ),
        #		option_definition.addToken('Input.aer[id]', int) ], #token
        #	        '',     #parents
        #	        '',     #non_parents
        #	        ''      #childs
        #	)
        #

        aerosol_profile_modtran = option_definition.option(
            name='aerosol_profile_modtran',
            group='aerosol',
            helpstr='Squeeze aerosol profile',
            documentation=documentation['aerosol_profile_modtran'],
            tokens=[option_definition.addSetting(name='Input.aer.profile_modtran', setting=True)],
        )

        aerosol_angstrom = option_definition.option(
            name='aerosol_angstrom',
            group='aerosol',
            helpstr='Scale the aerosol optical depth using the {\AA}ngstr{\"o}m formula.',
            documentation=documentation['aerosol_angstrom'],
            tokens=[option_definition.addToken(name='Input.aer.alpha', datatype=float),
                    option_definition.addToken(name='Input.aer.beta', datatype=float, default='NOT_DEFINED_FLOAT', valid_range=[0, 1e6]),
                    option_definition.addSetting(name='Input.aer.spec', setting=1)],
            parents=['aerosol_default'],
        )

        aerosol_king_byrne = option_definition.option(
            name='aerosol_king_byrne',
            group='aerosol',
            helpstr='Scale the aerosol optical depth using the King Byrne formula.',
            documentation=documentation['aerosol_king_byrne'],
            tokens=[option_definition.addToken(name='Input.aer.alpha_0', datatype=float, default='NOT_DEFINED_FLOAT',
                             valid_range=[-1e6, 1e6]),
                    # alph_0 is roughly equivalent with log(angstrom beta) and can thus will mostly be negative
                    option_definition.addToken(name='Input.aer.alpha_1', datatype=float, default='NOT_DEFINED_FLOAT',
                             valid_range=[-1e6, 1e6]),
                    option_definition.addToken(name='Input.aer.alpha_2', datatype=float, default='NOT_DEFINED_FLOAT',
                             valid_range=[-800, 800]),
                    option_definition.addSetting(name='Input.aer.spec', setting=1)],
            parents=['aerosol_default'],
        )

        aerosol_modify = option_definition.option(  # TODO: valid_ranges for GUI!!
                                  name='aerosol_modify',
                                  group='aerosol',
                                  helpstr='Modify aerosol optical properties.',
                                  documentation=documentation['aerosol_modify'],
                                  tokens=[option_definition.addLogical(name='id1', logicals=['gg', 'ssa', 'tau', 'tau550'],
                                                     setting='MODIFY_VAR_', gui_name='variable'),
                                          option_definition.addLogical(name='id2', logicals=['set', 'scale'], setting='MODIFY_TYPE_',
                                                     gui_name='scale/set'),
                                          option_definition.addToken(name='Input.aer.modify[id1][id2]', datatype=float, gui_name='value'),
                                          option_definition.addSetting(name='Input.aer.spec', setting=1, default=0)],
                                  non_unique=True,
                                  )

        aerosol_haze = option_definition.option(
            name='aerosol_haze',
            group='aerosol',
            helpstr='Specify the aerosol type in the lower 2 km of the atmosphere',
            documentation=documentation['aerosol_haze'],
            # gui_inputs=(IntegerGUI_definition.ListInput(name='Input.aer.haze', default=None, valid_range=[1, 4, 5, 6]),),
            tokens=[option_definition.addToken(name='Input.aer.haze', datatype=int, valid_range=[1, 6])],
            parents=[],  # Please review! ['aerosol_vulcan', 'aerosol_season', 'aerosol_visibility'],
            childs=['aerosol_option_specification'],
            extra_dependencies=['aerosol_haze', 'aerosol_vulcan', 'aerosol_visibility', 'aerosol_season'],
        )

        aerosol_set_tau_at_wvl = option_definition.option(
            name='aerosol_set_tau_at_wvl',
            group='aerosol',
            helpstr='Set the aerosol optical thickness at lambda (nm)',
            documentation=documentation['aerosol_set_tau_at_wvl'],
            tokens=[option_definition.addToken(name='Input.aer.tau_wvl_lambda', datatype=float, default='NOT_DEFINED_FLOAT',
                             valid_range=[0, 1e6]),
                    option_definition.addToken(name='Input.aer.tau_wvl_tau', datatype=float, default='NOT_DEFINED_FLOAT',
                             valid_range=[0, 1e6])],
        )

        aerosol_season = option_definition.option(
            name='aerosol_season',
            group='aerosol',
            helpstr='Specify season',
            documentation=documentation['aerosol_season'],
            tokens=[option_definition.addToken(name='Input.aer.seasn', datatype=int, valid_range=[1, 2]),
                    option_definition.addSetting(name='Input.aer.spec', setting=1)],
            childs=['aerosol_option_specification'],
            extra_dependencies=['aerosol_haze', 'aerosol_vulcan', 'aerosol_visibility', 'aerosol_season'],
        )

        aerosol_species_file = option_definition.option(  # TODO: undefined number of optional arguments are allowed
                                        name='aerosol_species_file',
                                        group='aerosol',
                                        helpstr='Specify mass density profiles of a mixture of aerosol types',
                                        documentation=documentation['aerosol_species_file'],
                                        gui_inputs=(GUI_definition.ListInput(name='Input.aer.mixture_name',
                                                              valid_range=['continental_clean', 'continental_average',
                                                                           'continental_polluted', 'urban',
                                                                           'maritime_clean', 'maritime_polluted',
                                                                           'maritime_tropical', 'desert', 'antarctic'],
                                                              optional=False),),
                                        tokens=[option_definition.addToken(name='Input.aer.mixture_name', datatype=str,
                                                         valid_range=['continental_clean', 'continental_average',
                                                                      'continental_polluted',
                                                                      'urban', 'maritime_clean', 'maritime_polluted',
                                                                      'maritime_tropical', 'desert', 'antarctic',
                                                                      'desert_spheriods']),
                                                option_definition.addSetting(name='Input.aer.n_species', setting=-1,
                                                           default='NOT_DEFINED_INTEGER'),
                                                option_definition.addSetting(name='Input.aer.spec', setting=1)],
                                        )

        aerosol_species_library = option_definition.option(
            name='aerosol_species_library',
            group='aerosol',
            helpstr='Location of optical property files',
            documentation=documentation['aerosol_species_library'],
            tokens=[option_definition.addToken(name='Input.aer.filename[FN_AER_SPECIES_LIB]', datatype=io.IOBase, valid_range=['OPAC', io.IOBase])],
        )

        aerosol_visibility = option_definition.option(
            name='aerosol_visibility',
            group='aerosol',
            helpstr='Horizontal visibility in km',
            documentation=documentation['aerosol_visibility'],
            tokens=[option_definition.addToken(name='Input.aer.visibility', datatype=float, default='NOT_DEFINED_FLOAT',
                             valid_range=[0, 1e6])],
            childs=['aerosol_option_specification'],
            extra_dependencies=['aerosol_haze', 'aerosol_vulcan', 'aerosol_visibility', 'aerosol_season'],
        )

        aerosol_vulcan = option_definition.option(
            name='aerosol_vulcan',
            group='aerosol',
            helpstr='Aerosol situation above 2 km',
            documentation=documentation['aerosol_vulcan'],
            # gui_inputs=(IntegerGUI_definition.ListInput(name='Input.aer.vulcan', default=None, valid_range=[1, 2, 3, 4], optional=False),),
            tokens=option_definition.addToken(name='Input.aer.vulcan', datatype=int, valid_range=[1, 4]),
            childs=['aerosol_option_specification'],
            extra_dependencies=['aerosol_haze', 'aerosol_vulcan', 'aerosol_visibility', 'aerosol_season'],
        )

        aerosol_sizedist_file = option_definition.option(
            name='aerosol_sizedist_file',
            group='aerosol',
            helpstr='Aerosol size distribution file.',
            documentation=documentation['aerosol_sizedist_file'],
            tokens=option_definition.addToken(name='Input.aer.filename[Id]', datatype=io.IOBase),
            childs=['aerosol_option_specification'],
            extra_dependencies=[],
        )
        aerosol_refrac_index = option_definition.option(
            name='aerosol_refrac_index',
            group='aerosol',
            helpstr='Set the aerosol particle refractive index.',
            documentation=documentation['aerosol_refrac_index'],
            tokens=[option_definition.addToken(name='Input.aer.nreal', datatype=float),
                    option_definition.addToken(name='Input.aer.nimag', datatype=float, default='NOT_DEFINED_FLOAT', valid_range=[0, 1e6]),
                    option_definition.addSetting(name='Input.aer.spec', setting=1)],
            parents=['aerosol_default'],
        )
        self.options = [aerosol_default,
                        aerosol_file, aerosol_species_library, aerosol_species_file,
                        aerosol_haze, aerosol_season, aerosol_vulcan, aerosol_visibility,
                        aerosol_profile_modtran,
                        aerosol_angstrom, aerosol_king_byrne,
                        aerosol_modify, aerosol_set_tau_at_wvl, aerosol_sizedist_file,
                        aerosol_refrac_index]

    def __iter__(self):
        return iter(self.options)


def get_aerosol_documentation():
    return {
        'aerosol_default': r'''
	Set up a default aerosol according to \citet{shettle89}. 
	The default properties are 
	a rural type aerosol in the boundary layer, background aerosol above 2km,
	spring-summer conditions and a visibility of 50km. These settings may
	be modified with \code{aerosol\_haze}, \code{aerosol\_vulcan}, \code{aerosol\_season},
	and \code{aerosol\_visibility}.
		''',

        'aerosol_file': r'''
	Location of file defining aerosol optical properties.
	\fcode{
	   aerosol\_file type file
	}
	\code{type} defines the file type, which can be one of the following:
	\begin{description}
	\parameter{gg} Location of aerosol asymmetry parameter file.

	The file must have two columns.
	Column 1 is the altitude in km. 
	Column 2 is the asymmetry parameter of each layer. 
	The asymmetry parameter defined with this option is constant with wavelength.
	If you require spectral dependence please use \code{aerosol\_file explicit}.
	Comments start with \code{\#}. Empty lines are ignored.
	\parameter{ssa} Location of aerosol single scattering albedo file. 

	The file must have two columns.
	Column 1 is the altitude in km. The altitude grid must be exactly equal to 
	the altitude grid specified in the file \code{atmosphere\_file}.
	Column 2 is the single scattering albedo of each layer. 
	The single scattering albedo defined with this option is constant with wavelength. If you require spectral dependence please use \code{aerosol\_file explicit}.
	Comments start with \code{\#}. Empty lines are ignored.
	\parameter{tau} Location of aerosol optical depth file. 

	The file must have two columns.
	Column 1 is the altitude in km. The altitude grid must be exactly equal to 
	the altitude grid specified in the file \code{atmosphere\_file}.
	Column 2 is the aerosol optical depth of each layer. 
	To allow wavelength-dependent aerosol optical thickness please
	use either \code{aerosol\_angstrom} or \code{aerosol\_file explicit}.
	Comments start with \code{\#}. Empty lines are ignored.
	\parameter{moments} Set the aerosol phase function moments to the values specified in 
	the aerosol moments file.

	The file contains one column with arbitrary number of 
	Legendre terms of the phase function. The phase function 
	 $p(\mu)$ is
	\begin{equation}
	   p (\mu) = \sum_{m=0}^{\infty} (2m+1) \cdot k_m \cdot P_m (\mu) 
	\end{equation}
	where $k_m$ is the m'th moment and $P_m (\mu)$ is the m'th Legendre
	polynomial.  If not specified, a Henyey-Greenstein phase function is
	assumed where the asymmetry parameter g is either a default value
	depending on the aerosol type or it may be specified using
	\code{aerosol\_modify set gg}.  The phase function will be the same for all
	altitudes and wavelengths.  See \code{aerosol\_file explicit} if more
	flexibility is wanted. May only be used together with the
	\code{disort} or \code{fdisort2} solver in combination with the option
	\code{disort_intcor moments}.
	\parameter{explicit} A way to specify aerosol optical depth, single scattering albedo, 
	and phase function moments for each layer. 

	The file must have two columns where column 1 is the altitude in km. The second
	column is a the name of a file which defines the optical properties of the layer
	starting at the given altitude. The files specified in the second column must have the 
	following format:
	\begin{description}
	\item[Column 1:]
	The wavelength in nm. These wavelengths may be different from those in \code{source solar file}. 
	Optical properties are interpolated to the requested wavelengths.
	\item[Column 2:]
	The extinction coefficient of the layer in units km-1. 
	\item[Column 3:]
	The aerosol single scattering albedo of the layer.
	\item[Column 4-(nmom+4):]
	The moments of the aerosol phase function. 
	\end{description}
	For some simple examples see the files \file{examples/AERO\_*.LAYER}. Note that 
	if using the \code{rte_solver disort} it makes good sense to make the 
	number of moments larger than \code{number_of_streams}. For \code{rte_solver fdisort1} and 
	\code{rte_solver polradtran} the number of moments included in the calculations 
	will be \code{number_of_streams}+1. Higher order moments will be ignored for these solvers.
	Please note that the uppermost line of the \code{aerosol_file explicit}
	denotes simply the top altitude of the uppermost layer. The optical
	properties of this line are consequently ignored. There are two
	options for this line: either an optical property file with zero
	optical thickness is specified or "NULL" us used.
	\end{description}
		''',

        'aerosol_species_file': r'''
	Specify mass density profiles of a mixture of aerosol types. 
	\fcode{
	aerosol\_species\_file profile [aero\_1 aero\_2 ... aero\_n]
	}
	where \code{aero\_1} to \code{aero\_n} are the aerosol species to be included. 
	For each of these species, the optical properties are read from the 
	\code{aerosol\_species\_library}, e.g. the OPAC data set provided with 
	libRadtran. The profile file needs to include vertical profiles for each
	of these species. This file can be either in \emph{netCDF}-format 
	(automatically recognized filename extension \code{.nc} or \code{.cdf}) 
	or in ASCII format. The format of the ASCII file is: 
	\fcode{
	z1    dens(aero\_1, z1) dens(aero\_2, z1)  ... dens(aero\_n, z1)\\
	z2    dens(aero\_1, z2) dens(aero\_2, z2)  ...\\
	 .        .           .\\
	 .        .           .\\
	}
	where \code{z} is the height in km, and \code{dens} are the aerosol mass
	densities in g/m3. Please make sure to include one column for each of 
	the species \code{aero\_1} to \code{aero\_n} listed after 
	\code{aerosol\_species\_file}. 
	For netCDF input it is also possible to specify the unit 'kg kg$^{-1}$'; the data 
	are then automatically converted to g/m$^3$.
	
	Some default aerosol mixtures are provided, corresponding to the definitions in 
	\citet{hess98:_optic_proper_aeros_cloud}.
	They can simply be invoked by
	\fcode{
	aerosol\_species\_file mixture\_name
	}
	where \code{mixture\_name} can be one of the following: 
	\fcode{\\
	continental\_clean\\
	continental\_average\\
	continental\_polluted\\
	urban\\
	maritime\_clean\\
	maritime\_polluted\\
	maritime\_tropical\\
	desert\\
	antarctic
	}
	A variation of the desert mixture containing nonspherical particles is 
	\fcode{
	desert\_spheroids
	}
		''',
        'aerosol_species_library': r'''
	With this option the \emph{directory} is specified where the optical property
	files for all aerosols species used in the \code{aerosol\_species\_file} are
	expected: For each species defined in \code{aerosol\_species\_file},
	\emph{netCDF}-file \emph{species\_name}\code{.nc},
	(e.g. \code{INSO.nc}), which contains the optical properties of the aerosol
	species, has to be provided. The netcdf format is the one produced by 
	the {\sl libRadtran} \code{mie} tool.
	 
	At the libRadtran webpage we provide the OPAC data set 
	\citep{hess98:_optic_proper_aeros_cloud}
	which can be directly used with \code{uvspec}: 
	\fcode{
	aerosol\_species\_library OPAC
	}
	OPAC contains following aerosol species: 
	\fcode{INSO   insoluble \\
	WASO   water\_soluble\\
	SOOT   soot\\
	SSAM   sea\_salt\_accumulation\_mode\\
	SSCM   sea\_salt\_coarse\_mode\\
	MINM   mineral\_nucleation\_mode\\
	MIAM   mineral\_accumulation\_mode\\
	MICM   mineral\_coarse_mode\\
	MITR   mineral\_transported\\
	SUSO   sulfate\_droplets\\
	}
	Variations of the mineral aerosol species containing nonspherical particles are:
	\fcode{MINM_SPHEROIDS   mineral\_nucleation\_mode\\
	MIAM_SPHEROIDS   mineral\_accumulation\_mode\\
	MICM_SPHEROIDS   mineral\_coarse\_mode
	}
	The nonspherical mineral components are described by \citet{koepke2015}.
		''',

        'aerosol_season': r'''
	Specify season to get appropriate aerosol profile.
	\fcode{
	   aerosol\_season season
	}
	where \code{season} is either \code{1} or \code{2}:
	\begin{description}
	\item[1] Spring-summer profile.
	\item[2] Fall-winter profile.
	\end{description}
		''',

        'aerosol_visibility': r'''
	Horizontal visibility in km. Affects the profile according to \citet{shettle89}
	and the optical thickness. 
	\fcode{
	   aerosol\_visibility value
	}
		''',

        'aerosol_vulcan': r'''
	Aerosol situation above 2 km as defined in \citet{shettle89}.
	\fcode{
	   aerosol\_vulcan value
	}
	where \code{value} is an integer choosing between the following models
	\begin{description}
	\item[1] Background aerosols.
	\item[2] Moderate volcanic aerosols.
	\item[3] High volcanic aerosols.
	\item[4] Extreme volcanic aerosols.
	\end{description}
		''',

        'aerosol_haze': r'''
	Specify the aerosol type in the lower 2 km of the atmosphere as
	\fcode{
	   aerosol\_haze type
	}
	where \code{type} is an integer identifying the following aerosol types:
	\begin{description}
	\item[1] Rural type aerosols.
	\item[4] Maritime type aerosols.
	\item[5] Urban type aerosols.
	\item[6] Tropospheric type aerosols.
	\end{description}
	For a description of the different aerosol types see \citet{shettle89}.
		''',

        'aerosol_profile_modtran': r'''
	Squeeze aerosol profile up to 6 km when altitude is non-zero as in MODTRAN.
	Per default the aerosol profile is shifted upwards
	and remains unchanged. 
		''',

        'aerosol_angstrom': r'''
	Scale the aerosol optical depth using the {\AA}ngstr{\"o}m formula:
	\begin{equation}
	\tau = \beta \lambda^{-\alpha}
	\end{equation}
	where $\lambda$ is in units of micrometer \citep{Angstrom1929}. Specify 
	the {\AA}ngstr{\"o}m alpha and beta coefficients by
	\fcode{
	   aerosol\_angstrom alpha beta
	}
	The optical thickness defined here  is the integral from the user-defined 
	\code{altitude} to TOA (top of atmosphere). 
		''',

        'aerosol_king_byrne': r'''
	Scale the aerosol optical depth using the King Byrne formula \citep{king1976}:
	\begin{equation}
	\tau = e^{\alpha_0} \cdot \lambda^{\alpha_1} \lambda^{-\alpha}
	\end{equation}
	where $\lambda$ is in units of micrometer \citep{Angstrom1929}. Specify 
	the King Byrne $alpha_0$ $alpha_1$ $alpha_2$ coefficients by
	\fcode{
	   aerosol\_king\_byrne alpha_0 alpha_1 alpha_2
	}
	The optical thickness defined here  is the integral from the user-defined 
	\code{altitude} to TOA (top of atmosphere). 
		''',

        'aerosol_modify': r'''
	Modify aerosol optical properties.
	\fcode{
	aerosol\_modify variable scale/set value
	}
	This option is identical to \code{wc\_modify}.
	Please refer to \code{wc\_modify} for a detailed description of \code{variable}.
		''',
        'aerosol_set_tau_at_wvl': r'''
	Set the aerosol optical thickness at wavelength lambda (nm). Other wavelengths are scaled accordingly.
	Note that this option requires for technical reasons that the
	wavelength interval defined by \code{wavelength} does contain \code{lambda}.
	The optical thickness defined here is the integral from the user-definded 
	\code{altitude} to TOA (top of atmosphere). 
	\fcode{
	   aerosol\_set\_tau\_at\_wvl lambda tau
	}
		''',

        'aerosol_sizedist_file': r'''
        Calculate optical properties from size distribution and index of refraction using Mie
        theory. Here is an exception from the rule that ALL values defined above are overwritten
        because the optical thickness profile is re-scaled so that the optical thickness
        at the first internal wavelength is unchanged. It is done that way to give the user an
        easy means of specifying the optical thickness at a given wavelength.
        ''',

        'aerosol_refrac_index': r'''
        Calculate optical properties from size distribution and index of refraction using Mie
        theory. Here is an exception from the rule that ALL values defined above are overwritten
        because the optical thickness profile is re-scaled so that the optical thickness
        at the first internal wavelength is unchanged. It is done that way to give the user an
        easy means of specifying the optical thickness at a given wavelength.
        ''',
    }
