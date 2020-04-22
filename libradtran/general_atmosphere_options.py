
"""--------------------------------------------------------------------
 * $Id: general_atmosphere_options.py 3441 2018-12-10 15:21:05Z Claudia.Emde $
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

#from . import option_definition
#from . import GUI_definition

from libraddask.libradtran import option_definition
from libraddask.libradtran import GUI_definition

class setup_general_atm_group():

    group_name = 'General atmosphere'

    def __init__(self):
        documentation = get_general_atmosphere_documentation()

        atmos_region = option_definition.option(
            name='atmos_region',
            group='general atmosphere',
            helpstr='Define atmospheric region.', 
            documentation=documentation['atmos_region'],
            tokens = [option_definition.addSetting(name="Input.rte.mc.backward.yes", setting=1),
                      option_definition.addToken(name="Input.rte.mc.ixmin", datatype=int, default='NOT_DEFINED_INTEGER'),
                      option_definition.addToken(name="Input.rte.mc.iymin", datatype=int, default='NOT_DEFINED_INTEGER'),
                      option_definition.addToken(name="Input.rte.mc.ixmax", datatype=int, default='NOT_DEFINED_INTEGER'),
                      option_definition.addToken(name="Input.rte.mc.iymax", datatype=int, default='NOT_DEFINED_INTEGER') ] ,
            parents=['uvspec'],
            showInGui =False,
            threedmystic =True,
            )

        no_absorption = option_definition.option(
            name='no_absorption',
            group='general atmosphere',
            helpstr='Switch off absorption.', 
            documentation=documentation['no_absorption'],
            gui_inputs=(GUI_definition.TextInput(name='no_absorption',optional=True),),
            tokens = [option_definition.addSetting(name='ispoff', setting='get_caothoff_index(&Input.caothoff,&Input.n_caothoff, "all")'),
                      option_definition.addToken(name='ispoff', datatype=option_definition.CaothoffType, optional=True),
                      option_definition.addSetting(name='Input.caothoff[ispoff].no_absorption', setting=1)],
            parents=['uvspec'],
            non_unique=True,
            )

        no_scattering = option_definition.option(
            name='no_scattering',
            group='general atmosphere',
            helpstr='Switch off scattering.', 
            documentation=documentation['no_scattering'],
            gui_inputs=(GUI_definition.TextInput(name='no_scattering',optional=True),),
            tokens = [option_definition.addSetting(name='ispoff', setting='get_caothoff_index(&Input.caothoff,&Input.n_caothoff, "all")'),
                      option_definition.addToken(name='ispoff', datatype=option_definition.CaothoffType, optional=True),
                      option_definition.addSetting(name='Input.caothoff[ispoff].no_scattering', setting=1)],
            parents=['uvspec'],
            non_unique=True,
            )

        interpret_as_level = option_definition.option(
            name='interpret_as_level',
            group='general atmosphere',
            helpstr='Interpret profile properties as level properties.',
            documentation=documentation['interpret_as_level'],
            gui_inputs=(GUI_definition.TextInput(name='Input.caothoff[isp].layer'),),
            tokens = [option_definition.addToken(name='isp', datatype=option_definition.CaothType),
                      option_definition.addSetting(name='Input.caoth[isp].layer', setting='FALSE')], 
            parents=['wc_file','ic_file','profile_file'],
            non_unique=True,
        )
        
        zout_interpolate = option_definition.option(
            name='zout_interpolate', 
            group='general_atmosphere',
            helpstr='Interpolate atmospheric profiles.', 
            documentation=documentation['zout_interpolate'], 
            tokens=option_definition.addSetting(name='Input.atm.zout_interpolate', setting='ZOUT_INTERPOLATE', default='NO_ZOUT_INTERPOLATE'),
            parents=['uvspec'],
            )
    
        z_interpolate = option_definition.not_yet_lex2py_option(
            name='z_interpolate',
            group='general_atmosphere',
            documentation='', #documentation['z_interpolate'], #TODO: Missing documentation
            gui_inputs=(GUI_definition.TextInput(name=''),),
            parents=['uvspec'],
            non_unique=True,
        )

        atm_z_grid = option_definition.option(
            name='atm_z_grid', 
            group='atmosphere', 
            documentation=documentation['atm_z_grid'],
            gui_inputs=(GUI_definition.TextInput(name='Input.atm.z_atm_forced_sea'),),
            tokens= [option_definition.addToken(name='Input.atm.z_atm_forced_sea', datatype=option_definition.SignedFloats), 
                     option_definition.addSetting(name='Input.atm.nz_atm_forced_sea', setting='ntokens')],
            parents=['uvspec'],
            )

        reverse_atmosphere = option_definition.option(
            name='reverse_atmosphere', 
            group='atmosphere', 
            helpstr='Atmosphere is turned on the head.',
            documentation=documentation['reverse_atmosphere'],
            tokens=option_definition.addSetting(name='Input.rte.reverse', setting=1, default=0), 
            parents=['uvspec'],
            )

        self.options = [atmos_region, no_absorption, no_scattering, 
                interpret_as_level, zout_interpolate, z_interpolate,
                atm_z_grid, reverse_atmosphere]

    def __iter__(self):
        return iter(self.options)

def get_documentation():
    return get_general_atmosphere_documentation()

def get_general_atmosphere_documentation():
    return {
        'atmos_region' : r'''
        \emph{Very experimental option, ask Robert Buras or Luca Bugliaro for usage.}

        When reading 2D or 3D atmospheric input files from netCDF, load only
        a certain xy-region from the files. The arguments are given as follows:
        \fcode{
        atmos_region xmin ymin xmax ymax
    }
        The values for \code{xmin}, \code{ymin}, \code{xmax} and \code{ymax} indicate
        pixel indices, are 0-based and inclusive (i.e. a region of \code{0 0 1 2} selects
        pixels \code{(0,0), (0,1), (0,2), (1,0), (1,1), (1,2)}).
        ''',

        'no_absorption' : r'''
    Switch off absorption. Please note that this option simply sets
    the absorption optical thickness to 0.
    \fcode{
    no\_absorption [name]
    }
    If \code{name} is not set, all absorption (molecular, aerosol, cloud, ice cloud, and any profile) is switched off.

    If used together with \code{xxx\_modify set tau} this might be a bit confusing but 
    probably the most logical way. E.g. when using \code{aerosol\_default} and 
    \code{aerosol\_modify set tau 1}, the aerosol optical thickness is set to 1,
    with 0.940539 scattering and 0.059461 absorption. If \code{no\_absorption} 
    is added, the absorption optical thickness is set to 0 while the scattering
    optical thickness is preserved at 0.940539 (even though 1 was specified by 
    the user). We find this the most logical solution of the problem because 
    by switching \code{no\_absorption} off and on one tests the effect of the 
    absorber in an isolated way, rather than mixing absorption and scattering. 
    The same is true for water and ice clouds. Note, that thermal emission of 
    molecules is also switched off. 

    Possible choises for the optional argument \code{name} are:
    \begin{description}
    \item[mol] Switch off molecular absorption.
    \end{description}
        ''',

        'no_scattering' : r'''
    Switch scattering off for 1D profiles.
    \fcode{
    no\_scattering [name]
    }
    If \code{name} is not set, all scattering (molecular, aerosol, cloud, ice cloud, and any profile) is switched off. 

    Possible choises for the optional argument \code{name} are:
    \begin{description}
    \item[mol] Switch off molecular scattering.
    \item[aer] Switch off scattering by aerosols.
    \item[wc]  Switch off scattering by water clouds.
    \item[ic]  Switch off scattering by ice clouds.
    \item[profile] Switch off scattering by any profile defined in \code{profile typename}.
    \end{description}
        ''',

        'interpret_as_level'    : r'''
    Interpret profile properties as level properties (this was the default
    behaviour before version 1.4). 
    \fcode{
    interpret\_as\_level profile
    }
    profile can be either \code{wc}, \code{ic} or any profile type specified in \code{profile\_file}.

    If \code{interpret\_as\_level wc} is defined, a \code{wc\_file} would be interpreted as
    follows:
    \begin{verbatim}
    #      z     LWC    R\_eff
    #     (km)  (g/m^3) (um)  
           5.000    0      0   
           4.000   0.2   12.0 
           3.000   0.1   10.0 
           2.000   0.1    8.0 
    \end{verbatim}
    The value 0.2 g/m$^3$ refers to altitude 4.0km, as e.g. in a
    radiosonde profile. The properties of each layer are calculated as
    average over the adjacent levels. E.g. the single scattering
    properties for the model layer between 3 and 4km are obtained by
    averaging over the two levels 3km and 4km. To allow easy definition of
    sharp cloud boundaries, clouds are only formed if both liquid water
    contents above and below the respective layer are larger than
    0. Hence, in the above example, the layers between 2 and 3 as well as
    between 3 and 4km are cloudy while those between 1 and 2km and between
    4 and 5km are not.

    Note that since version 1.4 the default is to interpret profile properties as 
    layer properties. For example wc properties are assumed to be constant over the layer. 
    The layer reaches from the level, where the properties are defined in the
    \code{wc\_file} to the level above that one.  The following lines
    \begin{verbatim}
    #      z     LWC    R_eff
    #     (km)  (g/m^3) (um) 
           4.000   0.0   0.0
           3.000   1.0  10.0
    \end{verbatim}
    define a cloud in the layer between 3 and 4 km with sharp boundaries.
        ''',

        'atm_z_grid' : r'''
    With this option the vertical resolution of the \code{atmosphere\_file} data is changed to 
    the levels (in km above sea surface) given as argument. This might be useful in oder to reduce the number 
    of levels (save computational time)\ifmystic{ or in order to easily adjust the atmosphere profile to the resolution 
    of a Monte Carlo cloud file \code{wc\_file 3D} or \code{ic\_file 3D}}.
    \fcode{
    atm\_z\_grid 0 2 4 6 8 10 20 30 ...
    }
        ''',

        'zout_interpolate' : r'''
    The z-grid of optical properties is determined by the \code{atmosphere\_file},
    and, if specified, by other profile files like \code{mol\_file}, \code{rh\_file}, 
    or \code{refractive\_index\_file}. 
    Additional levels might be introduced by the \code{zout} 
    option and the second argument of the \code{altitude} option. By default 
    (if \code{zout\_interpolate} is not specified) levels introduced 
    by the \code{zout} option will not affect the optical property 
    profiles, that is, the optical properties are constant within the layers specified by 
    the \code{atmosphere\_file} and profile files.
    If \code{zout\_interpolate} is specified, the atmospheric profiles (tracegases, temperature ...) 
    are interpolated to the levels introduced by \code{zout}, and optical 
    properties are determined from the interpolated atmospheric properties.
    If \code{heating\_rate}, \code{rte\_solver polradtran}, \code{rte\_solver rodents},
    \code{rte\_solver twostrebe}, \code{rte\_solver twomaxrnd}, \code{rte\_solver twomaxrnd3C}
    is specified, \code{zout\_interpolate} 
    will also be automatically activated. 
    \code{zout\_interpolate} generally causes smoother variation of the optical properties.
        ''',

        'reverse_atmosphere' : r'''
    Option for the strong and bold. Reverses the atmospheric input to the radiative
    transfer solvers. That is, the atmosphere is turned on the head.
    Yes, that is actually useful for some purposes. If you think
    you need this contact the author. Otherwise, do not use.
        ''',

    }
