"""--------------------------------------------------------------------
 * $Id: molecular_options.py 3403 2018-04-30 07:45:57Z Claudia.Emde $
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
import io

class setup_molecular_group():

    group_name = 'Molecular atmosphere'

    def __init__(self):
        documentation = get_molecular_documentation()

        atmosphere_file = option_definition.option(
            name='atmosphere_file',
            group='atmosphere',
            helpstr='Location of the atmosphere data file.', 
            documentation=documentation['atmosphere_file'], 
            tokens= option_definition.addToken(name='Input.atmosphere_filename', datatype=str, valid_range=['subarctic_winter', 'subarctic_summer', 'midlatitude_summer', 'midlatitude_winter', 'tropics', 'US-standard', io.IOBase]),
            parents=['uvspec'], 
            plot = {'plot_type': '2D',
                'optional_args': {'column_names': (
                        "altitude",
                        "pressure",
                        "temperature",
                        "air",
                        "ozone",
                        "oxygen",
                        "water vapour",
                        "CO2",
                        "NO2",)
                          }
                },
            mandatory=True,
        )

        atmosphere_file_3D = option_definition.option(
            name='atmosphere_file_3D',
            group='atmosphere',
            helpstr='Location of the 3D atmosphere data file.', 
            documentation=documentation['atmosphere_file_3D'], 
            tokens= [option_definition.addToken(name='Input.atmosphere3d_filename', datatype=str),
                     option_definition.addSetting(name='Input.atmosphere3d', setting='True'),
                     option_definition.addSetting(name='isp', setting='get_caoth_index(&Input.caoth,&Input.n_caoth,"molecular_3d",0)'),
                     option_definition.addSetting(name='Input.caoth[isp].source', setting='CAOTH_FROM_3D'),
                     option_definition.addSetting(name='Input.caoth[isp].filename', setting='"profile_mol3d_dummy.dat"'),
                     option_definition.addSetting(name='Input.caoth[isp].properties', setting='PROP_HU')],
            parents=['uvspec'],
            showInGui = False,
            developer = False,
            threedmystic =True
                )
                

        radiosonde = option_definition.not_yet_lex2py_option(
            name='radiosonde',
            group='atmosphere',
            helpstr='Specify density profiles.',
            documentation=documentation['radiosonde'],
            gui_inputs=(GUI_definition.TextInput(name=''),),
            parents=['uvspec'],
        )

        radiosonde_levels_only = option_definition.option(
            name='radiosonde_levels_only',
            group='atmosphere', 
            helpstr='', 
            documentation=documentation['radiosonde_levels_only'], 
            tokens=option_definition.addSetting(name='Input.atm.rs_add_upper_levels',setting='FALSE'), 
            parents=['uvspec'],
        )

        mol_file = option_definition.option(
            name='mol_file',
            group='molecular',
            documentation=documentation['mol_file'],
            gui_inputs=(GUI_definition.ListInput(name='mol_id',
                          valid_range=['O3', 'O2', 'H2O',
                               'CO2', 'NO2', 'BRO',
                               'OCLO', 'HCHO', 'O4',
                               'SO2', 'CH4', 'N2O',
                               'CO', 'N2'],
                          optional=False),
                       GUI_definition.FileInput(name='Input.atm.filename[mol_id]'),
                       GUI_definition.ListInput(name='Input.atm.unit_profile[mol_id]', valid_range=['','cm_3', 'm_3', 'MMR', 'VMR', 'RH'], optional=True),),
            tokens = [option_definition.addLogical(name='mol_id', logicals=['O3', 'O2', 'H2O', 'CO2', 'NO2', 'BRO', 'OCLO', 'HCHO', 'O4', 'SO2', 'CH4', 'N2O', 'CO', 'N2'], setting='MOL_' ) ,
                      option_definition.addToken(name='Input.atm.filename[mol_id]', datatype=io.IOBase),
                      option_definition.addLogical(name='Input.atm.unit_profile[mol_id]', logicals=['cm_3', 'm_3', 'MMR', 'VMR', 'RH'], optional=True)],
            parents=['uvspec'],
            non_unique=True,
        )

        pressure = option_definition.option(
            name='pressure',
            group='atmosphere', 
            helpstr='Surface pressure',
            documentation=documentation['pressure'], 
            gui_inputs=(GUI_definition.FloatInput(name='Input.pressure', default='NOT_DEFINED_FLOAT', valid_range=[0, 1000000.0]),),
            tokens=option_definition.addToken(name='Input.pressure', datatype=float, default='NOT_DEFINED_FLOAT', valid_range=[0,1e6]),
            parents=['uvspec'],
        )

        refractive_index_file = option_definition.option(
            name='refractive_index_file',
            group='atmosphere', 
            helpstr='', 
            documentation=documentation['refractive_index_file'], 
            gui_inputs=(GUI_definition.TextInput(name='Input.filename[FN_REFIND]'),),
            tokens=option_definition.addToken(name='Input.filename[FN_REFIND]', datatype=str), 
            parents=['uvspec'],
        )

        crs_model = option_definition.option(
            name='crs_model',
            group='molecular',
            helpstr='Specify cross section.',
            documentation=documentation['crs_model'],
            gui_inputs=(GUI_definition.ListInput(name='mol_id', valid_range=['no2', 'o3', 'o4', 'rayleigh'], optional=False), GUI_definition.ListInput(name='crs_model', valid_range=['Bass_and_Paur', 'Molina', 'Daumont', 'Serdyuchenko', 'Bogumil', 'Bodhaine', 'Bodhaine29', 'Nicolet', 'Penndorf', 'Burrows', 'Vandaele', 'Greenblatt', 'Thalman'], optional=False),),
            tokens= [option_definition.addLogical(name='mol_id', logicals=['no2', 'o3', 'o4', 'rayleigh'], setting='CRS_MOL_'),
                     option_definition.addLogical(name='Input.crs_model[mol_id]',  logicals=[ 'Bodhaine', 'Nicolet', 'Penndorf', 'Bodhaine29', 'Bass_and_Paur', 'Molina', 'Daumont', 'Serdyuchenko', 'Bogumil', 'Burrows', 'Vandaele', 'Greenblatt', 'Thalman'], setting='CRS_MODEL_')],
            parents=['uvspec'],
            non_unique=True,
        )

        crs_file = option_definition.option(
            name='crs_file',
            group='molecular',
            documentation=documentation['crs_file'],
            gui_inputs=(GUI_definition.ListInput(name='mol_id', valid_range=['O3', 'O2', 'H2O', 'CO2', 'NO2', 'BRO', 'OCLO', 'HCHO', 'O4', 'SO2', 'CH4', 'N2O', 'CO', 'N2'], optional=False), GUI_definition.FileInput(name='Output.crs.filename[mol_id]'),),
            tokens = [option_definition.addLogical(name='mol_id', logicals=['O3', 'O2', 'H2O', 'CO2', 'NO2', 'BRO', 'OCLO', 'HCHO', 'O4', 'SO2', 'CH4', 'N2O', 'CO', 'N2'], setting='MOL_' ) ,
                      option_definition.addToken(name='Output.crs.filename[mol_id]', datatype=io.IOBase)],
            parents=['uvspec'],
            non_unique=True,
        )

        rayleigh_depol = option_definition.option(
            name='rayleigh_depol',
            group='atmosphere',
            helpstr='Rayleigh depolarization factor.', 
            documentation=documentation['rayleigh_depol'], 
            gui_inputs=(GUI_definition.FloatInput(name='Input.rayleigh_depol'),),
            tokens=option_definition.addToken(name='Input.rayleigh_depol', datatype=float, default='NOT_DEFINED_FLOAT'), 
            parents=['uvspec'],
        )

        mol_abs_param = option_definition.option(
            name='mol_abs_param',
            group='spectral',
            helpstr='Set correlated_k scheme. ',
            documentation=documentation['mol_abs_param'],
            tokens = [option_definition.addLogical(name='Input.ck_scheme', logicals=['kato', 'kato2', 'kato2.96','kato2andwandji','fu','avhrr_kratz','sbdart','lowtran','reptran','reptran_channel','crs',io.IOBase], setting='CK_'),
                      option_definition.addToken(name='Input.ck_reptran_arg', datatype=str, optional=True)],
            parents=['uvspec'],
            childs= ['ck_lowtran_absorption','ck_fu_h2o_continuum'],
            continious_update=True,
        )
                
        reptran_file = option_definition.option(
            name='reptran_file',
            group='spectral',
            helpstr='File containing representative wavelengths.',
            documentation=documentation['reptran_file'],
            tokens=option_definition.addToken(name='Input.filename[FN_REPTRAN]', datatype=io.IOBase), 
            parents=['uvspec'],
            showInGui = False,
                )

        ck_lowtran_absorption = option_definition.option(
            name='ck_lowtran_absorption',
            group='molecular',
            helpstr='Switch off absorption by individual minor trace gases.', 
            documentation=documentation['ck_lowtran_absorption'],
            gui_inputs=(GUI_definition.ListInput(name='Input.absorption_gas', valid_range=['O4', 'N2', 'CO', 'SO2', 'NH3', 'NO', 'HNO3']),
            GUI_definition.ListInput(name='On/Off',valid_range=['on','off'], default='On'),),
            tokens = [option_definition.addLogical(name='mol_id', logicals=['O4','N2','CO','SO2','NH3','NO','HNO3'], setting='CK_ABS_' ) ,
                      option_definition.addLogical(name='Input.ck_abs[mol_id]', logicals=['on','off'], setting='SWITCH_')],
            parents=['mol_abs_param'],
            speaker="mol_abs_param",
            enable_values=("lowtran",),
            non_unique=True,
        )

        ck_fu_h2o_continuum = option_definition.option(
            name='ck_fu_h2o_continuum',
            group='molecular',
            helpstr='', #TODO
            documentation=documentation['ck_fu_h2o_continuum'],
            tokens=option_definition.addLogical(name='Input.ck_h2ocont',logicals=['on','off','v2.1','v2.4'],setting='CK_H2OCONT_'),
            parents=['uvspec'],
            speaker='mol_abs_param',
            enable_values=('fu',), #TODO:stimmt das?
            developer=True,
        )

        mol_tau_file = option_definition.option(
            name='mol_tau_file',
            group='molecular',
            helpstr='Location of Molecular optical depth file.', 
            documentation=documentation['mol_tau_file'],
            gui_inputs=(GUI_definition.ListInput(name='id', valid_range=['sca', 'abs'], optional=False), GUI_definition.FileInput(name='Input.filename[id]'),),
            tokens = [option_definition.addLogical(name='id', logicals=[ 'sca','abs'], setting='FN_MOL_TAU_'),
                      option_definition.addToken(name='Input.filename[id]', datatype=io.IOBase)],
            parents=['uvspec'],
            non_unique=True,
        )

        mol_modify = option_definition.option(
            name='mol_modify',
            group='molecular',
            helpstr='Modify column of molecular specie', 
            documentation=documentation['mol_modify'],
            gui_inputs = ( GUI_definition.ListInput(name='moltype', valid_range=[ 'O3','O2','H2O','CO2','NO2','BRO','OCLO','HCHO','O4','SO2','CH4','N2O','CO','N2'], optional=False),
                           GUI_definition.FloatInput(name='value',  valid_range=[0, 1000000.0]),
                           GUI_definition.ListInput(name='unit', valid_range=[ 'DU', 'CM_2', 'MM'], optional=False)),
            tokens = [option_definition.addLogical(name='id', logicals=[ 'O3','O2','H2O','CO2','NO2','BRO','OCLO','HCHO','O4','SO2','CH4','N2O','CO','N2'], setting='MOL_'),
                      option_definition.addToken(name='Input.atm.column[id]', datatype=float),
                      option_definition.addLogical(name='Input.atm.unit_column[id]', logicals=[ 'DU', 'CM_2', 'MM'], setting='MOL_UNIT_')],
            parents=['uvspec'],
            non_unique=True,
        )

        mixing_ratio = option_definition.option(
            name='mixing_ratio',
            group='molecular',
            helpstr='Mixing ratio of molecular specie', 
            documentation=documentation['mixing_ratio'],
            gui_inputs = (GUI_definition.ListInput(name='moltype', valid_range=[ 'O2','H2O','CO2','NO2','CH4','N2O','F11','F12','F22'], optional=False),
                          GUI_definition.FloatInput(name='value',  valid_range=[0, 1000000.0])),
            tokens = [option_definition.addLogical(name='id', logicals=[ 'O2','H2O','CO2','NO2','CH4','N2O','F11','F12','F22'], setting='MX_'),
                      option_definition.addToken(name='Input.mixing_ratio[id]', datatype=float, valid_range=[0,1e6])],
            parents=['uvspec'],
            non_unique=True,
        )


        self.options = [atmosphere_file, atmosphere_file_3D,
            radiosonde, radiosonde_levels_only,
            mol_file, mixing_ratio, mol_modify,
            pressure,
            refractive_index_file,
            crs_model,
            crs_file,
            rayleigh_depol,
            mol_abs_param,
            ck_lowtran_absorption,
            ck_fu_h2o_continuum,
            mol_tau_file,
                            reptran_file,
            ]

    def __iter__(self):
        return iter(self.options)


def get_documentation():
    return get_molecular_documentation()

def get_molecular_documentation():
    return {
        'ck_lowtran_absorption' : r'''
    Switch off absorption by individual minor trace gases which are currently only
    included when \code{mol\_abs\_param lowtran} is chosen. The syntax is
    \fcode{
      ck\_lowtran\_absorption species on/off
    }
    where species may be one of O4, N2, CO, SO2, NH3, NO, HNO3. By default all
    are switched on. 
    
    This option may also be used to turn on/off absorption by O4 in spectral 
    resolution. It is on by default.
        ''',
        'atmosphere_file' : r'''
    Location of the atmospheric data file. 
    \fcode{
    atmosphere\_file file
    }
    The file must have at least three columns containing the altitude,
    pressure, and temperature. Missing profiles are filled with 0 (e.g., if you did not specify 
    the ozone profile, there will be no ozone absorption!), with exception of the air density which 
    is calculated from pressure and temperature. Other trace gases may be set by \code{mol\_file}.  
    The columns are interpreted as follows:
    \begin{description}
    \item[1] Altitude above sea level in km
    \item[2] Pressure in hPa
    \item[3] Temperature in K
    \item[4] air density in cm$^{-3}$
    \item[5] Ozone density in cm$^{-3}$
    \item[6] Oxygen density in cm$^{-3}$
    \item[7] Water vapour density in cm$^{-3}$
    \item[8] CO2 density in cm$^{-3}$
    \item[9] NO2 density in cm$^{-3}$
    \end{description}
    The atmosphere is specified
    top-down, that is, the top level is the first line in the file, the bottom
    (surface) level the last line. All properties refer to model \emph{level} z,
    not to model \emph{layer}. It is important that the correct units are 
    used, otherwise unpredictable results are guaranteed.
    Comments start with \code{\#}. Empty lines are ignored. Please note that there 
    is some redundancy: For air as an ideal gas the density $\rho$, can be
    calculated from pressure and temperature, $\rho = p / kT$. \code{uvspec} will check 
    if this relation is fulfilled and will stop if it is not. 
    {\sl libRadtran} provides the six standard atmospheres by \citet{Anderson1986}:
    \begin{description}
    \item[afglt]   Tropical            (\code{tropics})                                           
           \item[afglms]  Midlatitude Summer  (\code{midlatitude\_summer})
    \item[afglmw]  Midlatitude Winter  (\code{midlatitude\_winter})
    \item[afglss]  Subarctic Summer    (\code{subarctic\_summer})
    \item[afglsw]  Subarctic Winter    (\code{subarctic\_winter})
    \item[afglus]  U.S. Standard       (\code{US-standard})                   
    \end{description}
    which may be chosen by for example
    \fcode{
    atmosphere\_file tropics
    }
    or by specifying the full file name. These atmosphere files are found in 
    \file{data/atmmod}.
    If no \code{atmosphere\_file} is defined, {\sl uvspec} will automatically select one. 
    If the information  \code{time}, \code{latitude} and \code{longitude} are provided in 
    the input file {\sl uvspec} will choose from the first 5 files, otherwise it takes 
    the U.S. Standard atmosphere.
        ''',
                'atmosphere_file_3D' : r'''
        Specify filename for a 3D molecular atmosphere. The file includes 3D fields of pressure, temperature and water vapor. Other species are not yet implemented in 3D and are specified as 1D altitude profiles using the option \code{atmosphere\_file}. The format of the 3D atmosphere file is as follows: 
        \fcode{Nx  Ny  Nz  flag \\
               dx  dy  z(1) z(2) z(3) ... z(n) \\
               ix  iy  iz  p  T  H2O \\
               ...
                }
        where \code{Nx}, \code{Ny} and \code{Nz} are the number of grid boxes in 
        \code{x}, \code{y}, and \code{z}-direction.
        The parameter \code{flag} is not yet used.
        In the second line \code{dx} and \code{dy} are the sizes of the boxes in x- 
        and y-direction in km. In the third and following lines the indices \code{ix}, \code{iy}, and \code{iz} specify atmosphere pixels. \code{p} is the pressure in hPa, \code{T} the temperature in K and \code{H2O} the water vapor concentration in kg/kg.
 
        See also the examples \code{examples/UVSPEC_MC_ABS3D.INP} and \code{examples/UVSPEC_MC_ABS3D_THERMAL.INP} which use the 3D atmosphere file \code{examples/UVSPEC_MC_ABS3D_AFGLUS3D.INP} as input. This atmosphere file includes the same atmospheric profile (US standard) in 3$\times$2 pixels and the examples check for solar and thermal radiation, whether the 3D input yields the same results as 1D calculations.

        Currently the implementation has some restrictions:
        \begin{itemize}       
        \item the altitude profiles of the 1D atmosphere file, the 3D atmosphere file and also other 3D profile files (including \code{wc_file 3D} and \code{ic_file 3D}) must include the same vertical grids
        \item the conversion of heating rates to K/d is only approximate, because it uses pressure and temperature of the first pixel, rather than the 3D fields
        \item ... may be more?
        \end{itemize}

        {\sl \code{atmosphere_file_3D} is an experimental option! Please check your results carefully and contact Claudia Emde in case you find any bugs or inconsistencies.}
            ''',    

        'radiosonde' : r'''
    This option allows to change the temperature and pressure profile, and optionally to 
    specify one or more density profiles. The entry in the input file looks like this:
    \fcode{
    radiosonde filename [gas\_species] [unit] ...
    }
    Currently the following gas\_species are included: ozone (O3), nitrogen dioxide (NO2), 
    water vapor (H2O), bromine oxide (BRO), chlorine dioxide (OCLO), formaldehyde (HCHO), 
    carbon dioxide (CO2), sulphur dioxide (SO2), and the oxygen dimer (O4). 
    Each gas species is identified by its abbrevations given in parentheses above.
    Unit is an optional argument to defines the unit of the density. The profiles can
    be given in particles per cm$^3$ (CM-3), in particles per m$^3$ (M-3), as volume 
    mixing ratio (VMR), as mass mixing ratio in kg/kg (MMR), or as relative humidity 
    (RH) (only for water). The default unit is RH for water vapour, 
    MMR for ozone, and CM3 for all other gases.
    The radiosonde file must have (2 + number of gases) columns:
    \begin{description}
    \item[1] pressure in hPa
    \item[2] temperature in Kelvin
    \item[3, 4, ...] density of trace gas in the specified unit
    \end{description}
    A new z-grid will be calculated, starting at \code{altitude} and assuming a linear temperature variation
    between levels. The air density will be recalculated according to the ideal gas law, and the density of 
    the well mixed gases O2 and CO2 will be scaled accordingly.
    The atmospheric data above the radiosonde data is taken from the \code{atmosphere\_file} level by level, starting 
    at the first pressure level above the radiosonde data. The z-grid of the \code{atmosphere\_file} in 
    this height region is shifted accordingly.
    Also if the density in the radiosonde file is specified as -1 at a level, 
    the value from the \code{atmosphere\_file} is used.
    Possible calls are
    \fcode{
    radiosonde ../examples/radiosonde.dat
    }
    just in order to change the temperature and pressure profile, or
    \fcode{
    radiosonde ../examples/radiosonde2.dat H2O RH O3 MMR NO2
    }
    where water vapour density will be given as relative humidity, ozone as mass mixing ratio, 
    and NO2 in cm$^{-3}$ (default).
        ''',

        'radiosonde_levels_only' : r'''
    The atmosphere considered in the simulation has the same height range as the data in 
    the \code{radiosonde}-file. No further levels are added above those.
    This option has only an effect in combination with \code{radiosonde}.
        ''',

        'mol_file' :    r'''
    Specify density profiles (or matrix, see below) of various trace gases to be included in the radiative 
    transfer calculation.
    \fcode{
    mol\_file gas\_species filename [unit]
    }
    At the moment following \code{gas\_species} are included: ozone (O3), nitrogen dioxide (NO2), water vapor (H2O),
    bromine oxide (BRO), chlorine dioxide (OCLO), formaldehyde (HCHO), carbon dioxide (CO2), 
    sulphur dioxide (SO2), oxygen (O2), the oxygen dimer (O4), methane (CH4), nitrous oxide (N20), 
        carbon monoxide (CO), and nitrogen (N2). 
    The gas species is identified 
    by their abbrevations given in the parenthesis above.

    The model expects a density file with two columns:
    \begin{description}
    \item[1] Altitude above sea level in km.
    \item[2] The density of trace gas [in the specified unit]
    \end{description}
    The altitude grid may be different from that in \code{atmosphere\_file}. All densities inside the range 
    of the \code{mol\_file} are replaced. For all other altitudes the values from the 
    \code{atmosphere\_file} are used. If the density is specified as -1 at a level, 
    the value from \code{atmosphere\_file} is used.
        Altitude ranges not covered by the \code{atmosphere\_file} are ignored.

    \code{unit} is an optional argument to define the unit of the density. The profiles can
    be given in particles per cm$^{3}$ (\code{cm\_3}), in particles per m$^{3}$ (\code{m\_3}), as volume mixing ratio (\code{vmr}), as mass mixing 
    ratio (\code{mmr}), or as relative humidity (\code{rh}) (only for water). The default for \code{unit} is cm$^{-3}$.
    
    To scale the profile to a total column value use \code{mol\_modify}.
    
    For airmass factor calculations it is for some species necessary to account for the 
    variation of the profile with sza. This may be accomplished by specifying a \code{mol\_file} 
    in the following format:
    \fcode{
    0.0       SZA1      SZA2 ...\\
    z(1)    dens(1,1)    ...\\
    z(2)      .           .\\
     .        .           .
    }
    where z(i) are the altitude levels above sea level in km, SZA is the solar zenith 
    angle in degrees, and dens is the density [in the specified unit] of the trace gases as 
    function of solar zenith angle and altitude. 
    The matrix may only be specified for one species. It may however be combined with profiles
    of other species. 
    A density matrix can only be used in connection with \code{rte\_solver sdisort}!
        ''', 

        'pressure' : r'''
    The surface pressure (at the user-defined \code{altitude}) in hPa. 
    \fcode{
    pressure value
    }
    The pressure profile as well as air, O2 and CO2 density profiles 
    are scaled accordingly.
        ''',

        'refractive_index_file' : r'',

        'crs_model' : r'''
    Choose between various cross sections.
    \fcode{
    crs\_model species crs
    }
    Following \code{species} are included:
    \begin{description}

    \parameter{rayleigh} Specify the Rayleigh cross section.
    Choose between the following Rayleigh scattering cross sections (\code{crs}):
    \begin{description}
    \item[Bodhaine]        \citet{Bodhaine1999} Rayleigh scattering cross section using their Eqs. 22-23.
    \item[Bodhaine29]    \citet{Bodhaine1999} Rayleigh scattering cross section using their Eq. 29.
    \item[Nicolet]         \citet{Nicolet1984} Rayleigh scattering cross section.
    \item[Penndorf]        \citet{Penndorf1957} Rayleigh scattering cross section.
    \end{description}
    \citet{Bodhaine1999} is default.

    \parameter{o3} Choose ozone cross section. 
    \code{crs} can be one of the following:
    \begin{description}
    \item[Bass\_and\_Paur]    \citet{Bass1985} ozone cross section.
    \item[Molina]         \citet{Molina1986} ozone cross section.
    \item[Daumont]        Ozone cross section by \citet{Daumont1992}, \citet{Malicet1995}.
    \item[Bogumil]        Ozone cross section from \citet{Bogumil2003}.
        \item[Serdyuchenko]     Ozone cross section from Serdyuchenko.        
    \end{description}
    \citet{Molina1986} is default.

    \parameter{no2} Choose between the various NO2 cross sections.
    \code{crs} is one of:
    \begin{description}
    \item[Burrows]        \citet{Burrows1998} NO2 cross section.
    \item[Bogumil]        NO2 cross section from \citet{Bogumil2003}.
        \item[Vandaele]         NO2 cross section from Vandaele et al.
    \end{description}
    \citet{Burrows1998} is default.

        \parameter{o4} Choose between the various O4 cross sections.
    \code{crs} is one of:
    \begin{description}        
        \item[Greenblatt]     O4 cross section by \citet{greenblatt1990}.
        \item[Thalman]       O4 cross section by \citet{thalman2013}.
        \end{description}

        \end{description}        
        ''',

        'crs_file' : r'''
    May be used to specify cross sections of O3, O2, H2O, CO2, NO2, BRO, OCLO, HCHO, 
        O4, SO2, CH4, N2O, CO, or N2 to be used instead of those supplied with 
    {\sl libRadtran}. No temperature dependence may be specified. Use as follows:
    \fcode{
    crs\_file NO2 ../examples/no2\_crs.dat
    }
    The species, e.g. \code{NO2}, must be specified to identify the
    species for which the cross section applies.
    The cross section file has two columns:
    \begin{description}
    \item[1] wavelength (nm)
    \item[2] cross section (cm$^2$)
    \end{description}
        ''',

        'rh_file' : r'''
    File that defines a profile of relative humidity. 
    \fcode{
    rh\_file file
    }
    If specified, the water vapour 
    profile in \code{atmosphere\_file} is over-written. If -1 is specified at a level, the value 
    from \code{atmosphere\_file} is used. 
        ''',

        'ck_fu_h2o_continuum' : r'''
    Undocumented option to switch the H2O continuum on or off or select a specific 
    version of the continuum.
        ''',

        'mixing_ratio' : r'''
    Mixing ratio in ppm.
    \fcode{
    mixing\_ratio species value
    }
    \code{species} can be one of the following:
    \begin{description}
    \item[O2] The mixing ratio of O2 in ppm. Scale the profile so that the mixing
        ratio at the user-defined \code{altitude} assumes the specified value.
    \item[H2O] The mixing ratio of H2O in ppm. Scale the profile so that the mixing
        ratio at the user-define \code{altitude} assumes the specified value.
        \item[CO2] The mixing ratio of CO2 in ppm. Scale the profile so that the mixing
        ratio at the user-defined \code{altitude} assumes the specified value.
    \item[NO2] The mixing ratio of NO2 in ppm. Scale the profile so that the mixing
        ratio at the user-defined \code{altitude} assumes the specified value.
        \item[CH4] The mixing ratio of CH4 in ppm (default: 1.6 ppm). 
    \item[N2O] The mixing ratio of N2O in ppm (default: 0.28 ppm).
    \item[F11] The mixing ratio of F11 in ppm (default: 0.000268 ppm).
    \item[F12] The mixing ratio of F12 in ppm (default: 0.000503 ppm).
    \item[F22] The mixing ratio of F22 in ppm (default: 0.000105 ppm).
    \end{description}
        The \code{mixing_ratio} of F11, F12, and F22 and the default values for CH4 and N2O are ignored in case of \code{mol_abs_param reptran}.

        ''',

        'mol_modify' : r'''
    Set the total column of a density profile. The column is integrated between the 
    user-defined \code{altitude} and TOA (top of atmosphere). The syntax is
    \fcode{
    mol\_modify species column unit
    }
    where \code{species} is one of O3, O2, H2O, CO2, NO2, BRO, OCLO, HCHO, O4, SO2, 
        CH4, N2O, CO, or N2, see also \code{mol\_file}.
    The second argument is the total column value, and the third argument is the unit, 
    in which the column is given. The unit can be DU (Dobson units), CM\_2 (molecules/cm$^2$) or MM.

    Please note that the unit MM is only valid for species H2O and specifies the precipitable water 
    in kg / m2 (which is approximately 1mm).The water vapor profile is scaled accordingly. The precipitable water 
        is integrated from the user-defined \code{altitude} to TOA (top of atmosphere).

    The default units are DU for O3, and CM\_2 for all other gases. It is possible to have
    several \code{mol\_modify} commands in the input file (maximum one per species). The following sets
    the NO$_2$ total column to 1.2 DU.
    \fcode{
    mol\_modify NO2 1.2 DU
    }
        ''',

        'rayleigh_depol' : r'''
    Rayleigh depolarization factor.
    \fcode{
    rayleigh\_depol value
    }
    The Rayleigh scattering phase function is
    $p(\mu) = a + b  \mu^2$ where $a = 1.5{(1+\texttt{depol})/(2+\texttt{depol})}$ and  
    $b = 1.5{(1-\texttt{depol})/(2+\texttt{depol})}$. By default the depolarization is calculated
    using the expressions from \citet{Bodhaine1999}.
        ''',

        'mol_abs_param' : r'''
    To calculate integrated shortwave or longwave irradiance, or to simulate 
    satellite instrument channels, use
    \fcode{
       mol\_abs\_param type
    }
    to choose between the following types of schemes:
    \begin{description}
        \item[reptran]
        Representative wavelengths parameterization adapted for spectral bands.
        This parameterization is used by default if no \code{mol\_abs\_param} option is given 
        in the {\sl uvspec} input file.
        Different band widths may be selected by
        \fcode{
          mol\_abs\_param reptran [fine|medium|coarse]
        } 
        (fine: 1cm$^{-1}$; medium: 5cm$^{-1}$; coarse: 15cm$^{-1}$; coarse is default).
        The data files for coarse resolution are included in the libRadtran package.
        The files required for fine and medium resolution can be downloaded from the libRadtran homepage.
        Absorption data is mainly based on HITRAN 2004. Absorption by H2O, CO2, O3, N2O, CO, CH4, O2, N2, and NO2 
        is considered, and absorption by all other gases is zero. 
        By default volume mixing ratios of N2O, CO, CH4, and N2 (those are not in the 
        \code{atmosphere\_file}) from the US standard atmosphere are applied. 
        Use \code{mol\_file} or \code{mol\_modify} to change the gas profiles.
        In case of radiative transfer problems with solar source, the extraterrestrial spectrum from
        Kurudz is applied by default. This parameterization is described in detail by \citet{gasteiger2014}.
        \item[reptran\_channel]
        Representative wavelengths parameterization for satellite channels. Usage 
        \fcode{
          mol\_abs\_param reptran\_channel channel\_name
        }
        Channel-integrated quantities are obtained using \code{output\_process per\_band}.
        The file \file{data/correlated\_k/reptran/channel\_list.txt} provides a list of available channels; 
        more information on the channels is provided in \file{data/filter/}. 
        See \citet{gasteiger2014} for details about the approach.
        \item[crs]
        Switch off spectral parameterizations. Only molecular absorption cross sections from 
        \code{crs_file} (including the default ones) are considered.
    \item[kato] 
    \citet{Kato1999b} correlated-k distribution, shortwave; based on HITRAN 96. Please note that the 
    bands above 2.5 micrometer are not very reliable which, however, this has only little impact
    on integrated shortwave radiation.
    \item[kato2] 
    \citet{Kato1999b}, shortwave; optimized version (Seiji Kato, personal communication, 2003);
    please note that \code{kato2} only has 148 subbands (that is, calls to the \code{rte\_solver}) 
    compared to 575 for \code{kato} which translates to a increase in computational speed by 
    up to a factor of 4 with only little increase in uncertainty. The absorption data are
    based on HITRAN 2000. Please note that the bands above 2.5 micrometer are not very reliable which, 
    however, this has only little impact on integrated shortwave radiation.
    \item[kato2andwandji] 
        Similar to \code{kato2}, but the UV bands \#3 and \#4 use the improved parameterization 
        by \citet{WandjiNyamsi2015}.
    \item[kato2.96]
    \citet{Kato1999b}, shortwave; optimized version (Seiji Kato, personal communication, 2003);
    similar to \code{kato2} but based on HITRAN96. Please note that the bands above 2.5 micrometer 
    are not very reliable which, however, has only little impact on integrated shortwave radiation.
    \item[fu]
    \citet{fu92,fu93}, shortwave and longwave; fast parameterization, developed for climate models.
    \item[avhrr\_kratz] 
    \citet{Kratz1995}, AVHRR instrument channels
    \item[lowtran]
    Gas absorption parameterization from LOWTRAN; code adopted from SBDART \citep{Ricchiazzi1998b}; 
    please see the section on "Spectral resolution".
    \item[sbdart]
    Identical to LOWTRAN.
    \end{description}
    If \code{mol\_abs\_param} kato/kato2/kato2.96/fu/avhrr\_kratz is specified, the extraterrestrial 
    flux is taken from
    internally defined files specific for each parameterization, not 
    from \code{source solar file}. The output is the integrated irradiance for 
    each band. To get e.g. integrated shortwave irradiance, simply add all 
    bands of the \citet{Kato1999b} or the \citet{fu92,fu93}
    parameterization. The five AVHRR channels are weighted sums of the
    libRadtran output. Examples how to integrate the output in the
    \code{avhrr\_kratz} case are included in the {\sl uvspec} self check
    which is initiated with  
    \code{make check}.
        ''',

                 'reptran_file' : r'''
        Location of the representative wavelengths file. 
        \fcode{
          reptran\_file file
        }
        This option is useful together with 'mol\_abs\_param reptran' and 'mol\_abs\_param reptran\_channel' 
        options, if you want to use your own representative wavelengths parameterization.
        ''',

        'mol_tau_file' : r'''
    Location of molecular scattering or absorption optical depth file.
    \fcode{
    mol\_tau\_file sca/abs filename
    }
    \begin{description}
    \parameter{sca}    Usually, the Rayleigh scattering
    cross section is calculated from the air pressure provided in \code{atmosphere\_file}
    (scaled with \code{pressure}). Use this parameter only if you really want to specify
    the optical depth directly (e.g. for a model intercomparison). The
    optical thickness profile may be either monochromatic or spectral. 
    \parameter{abs} Usually, molecular absorption
    is calculated from trace gas concentrations provided in \code{atmosphere\_file}
    (scaled with \code{mol\_modify O3}, etc.). Use this option only if you want to specify
    the optical depth directly (e.g. for a model intercomparison) or for a line-by-line
    calculation. If a spectral \code{mol\_tau\_file} is specified, the wavelength 
    grid defined there is used as the internal wavelength grid for the radiative transfer
    calculation, if not defined otherwise with \code{wavelength\_grid\_file}.
    \end{description}

    The file can be either of the following three formats:
    \begin{description}
    \parameter{Monochromatic} 
    Column 1 is the altitude in km 
    %AK: This is probably old stuff prior to redistribute
    %where the altitude grid must be exactly equal 
    %to the altitude grid specified in \code{atmosphere\_file}. 
    Column 2 is the absorption optical depth of each layer.
    \parameter{Spectral, ASCII}
    The first line contains the level altitudes in decreasing order; the following lines 
    contain the wavelength [nm] in the first column and then the absorption optical depths
    of each layer. 
    \parameter{Spectral,  netcdf}
    An example is available at the libRadtran homepage, 
    the file \file{UVSPEC.O2A.afglms.cdf} is a line-by-line spectrum of the oxygen A-Band 
    around 760nm, calculated for the mid-latitude summer 
    atmosphere. The advantage of
    netcdf compared to ASCII is that it is much faster to read, and that the file 
    is a self-contained, including data and a description of the variables and arrays.
    It is therefore particularly useful for line-by-line calculations where usually 
    many spectral data points are involved. 
    %netcdf is a common platform independent format; the description, a library to read and 
    %write netcdf including some tools to generate netcdf is available at 
    %http://www.unidata.ucar.edu/packages/netcdf/. A \code{mol\_tau\_file abs} must obey 
    %certain rules; 
    \end{description}
    Comments start with \code{\#}. Empty lines are ignored.
        ''',

}
