__author__ = 'DGriffith'
"""
 *--------------------------------------------------------------------
 *
 * This file is part of libraddask, orginally from morticia .
 * Copyright (c) 2015-2016 by Derek Griffith and Ari Ramkiloawn
 *
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

"""
.. module:: librad
    :platform: Unix only if running of libRadtran/uvspec is required. Also Windows for setting up uvspec cases.
    :synopsis: This module provides access to libRadtran through the uvspec utility. For information on libRadtran go to
    http://www.libtradtran.org
Some of the code within this module and code imported by the module is provided with libRadtran in the src_py
directory. In order to run libRadtran cases it is necessary to have libRadtran installed on your computer.
There is no Windows-native version of libRadtran, so that generally implies that you are running on a Unix/Linux,
however, reading and writing of uvspec input files, as well as reading of uvspec output files is possible on any
platform.

Important note concerning the uu output:
The uu field in the libRadtran/uvspec Case object contains the output radiances. In the most general case, this is a
5-dimensional numpy array, with axes in the following order:
1) umu - the cosine of the propagation zenith angles, where umu=1 is looking downward
2) phi - the propagation azimuth angle, where 180 deg is viewing northwards
3) spectral variable ('wavelength'/'wvl', 'wavenumber'/'wvn' and 'lambda' are alternative names)
4) altitudes of viewing point ('zout', 'zout_sea', 'pressure'/'p' and 'cth' are alternative level variables
5) stokes polarisation component (up to 4 stokes components, called I, Q, U and V)

By 'propagation' is meant the direction in which light is travelling, which is 180 degrees different from the
'viewing' direction.

Trailing singleton dimensions are dropped. That is, if there is only one stokes component, then uu has only 4
dimensions. If there is also only one level output (zout etc.), then uu has only 3 dimensions and so forth.

Leading and sandwiched singleton dimensions are not dropped in the current implementation. Dropping of singleton
dimensions can break existing code.

The relevant code here is taken from libRadtran version 2.0

"""

# TODO Reading of montecarlo (mystic/mc) output files, particularly:
# mc.flx.spc - spectral irradiance, actinic flux at levels
# mc.rad.spc - spectral radiance at levels
# TODO Dealing with exceptions, warnings. Keywords missing from the options library etc.
# TODO dealing with setting of units based on whatever is known about inputs/outputs, thermal is W/m^2/cm^-1

# from . import writeLex  # This imports all the libradtran option definitions
# from . import writeLex  # This imports all the libradtran option definitions
from libraddask.libradtran import writeLex

import os
import easygui  # For file open dialogs
import numpy as np
import xarray as xr
import re
import warnings

from libraddask.rad import moglo
from libraddask.rad import xd

import copy
from itertools import chain  #Used in RadEnv constructor
import matplotlib.pyplot as plt
import sys

_isfloatnum = '^[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?$'  # regular expression for matching tokens to floating point

# Definitions of basic Shettle/Fenn aerosols in libRadtran
shettle_aerosol_haze = {'Rural': 1, 'Maritime': 4, 'Urban': 5,
                    'Tropospheric': 6}  # Shettle/Fenn type aerosols.
# Seasonal modifications for Shettle/Fenn aerosools
shettle_aerosol_season = {'Spring-summer': 1, 'Fall-winter': 2}
# Shettle aerosol type above 2 km - largely stratospheric component
# driven by large volcanic eruptions
shettle_aerosol_vulcan = {'Background': 1, 'Moderate': 2, 'High': 3, 'Extreme': 4}

uvsOptions = writeLex.loadOptions()  # Load all options into a global dictionary of uvspec option specifications.

# The following dictionary provides units for some of the source solar files provided with libRadtran
sourceSolarUnits = {
    'apm_0_5nm': ['mW', 'm^2', 'nm'],
    'apm_1nm': ['mW', 'm^2', 'nm'],
    'atlas_plus_modtran': ['mW', 'm^2', 'nm'],
    'atlas_plus_modtran_ph': ['count/s', 'cm^2', 'nm'],  # count photons
    'atlas2': ['mW', 'm^2', 'nm'],
    'atlas3': ['mW', 'm^2', 'nm'],
    'fu': ['W', 'm^2', ''],  # in band
    'kato': ['W', 'm^2', ''],  # in band
    'kurudz_0.1nm.dat': ['mW', 'm^2', 'nm'],
    'kurudz_1.0nm.dat': ['mW', 'm^2', 'nm'],
    'NewGuey2003.dat': ['mW', 'm^2', 'nm'],
    'Thekaekara.dat': ['mW', 'm^2', 'nm']
}

# Units generally depend on the units in the solar flux file, but may be dictated by the correlated-k band model
uvspecOutVars = {
    'lambda': 'Wavelength [nm]',  # Cannot be used in Python because lambda is a keyword use 'wvl
    'wavelength': 'Wavelength [nm]', # scalar per line
    'wavenumber': 'Wavenumber [cm^1]',  # abbreviate to wvn
    'wvn': 'Wavenumber [cm^-1]',
    'wvl': 'Wavelength [nm]',
    'edir': 'Direct beam irradiance (same unit as extraterrestrial irradiance, e.g mW/m^2/nm if using the "atlas3" spectrum in the /data/solar_flux/ directory.)',
    'edn': 'Diffuse downwelling irradiance, i.e. total minus direct beam (same unit as edir).',
    'eup': 'Diffuse upwelling irradiance (same unit as edir).',
    'uavg': 'The mean intensity. Proportional to the actinic flux: To obtain the actinic flux, multiply the mean intensity by 4 pi (same unit as edir).',
    'uavgdir': 'Direct beam contribution to the mean intensity. (same unit as edir).',
    'uavgdn': 'Diffuse downward radiation contribution to the mean intensity. (same unit as edir).',
    'uavgup': 'Diffuse upward radiation contribution to the mean intensity. (same unit as edir).',
    'umu': 'Cosine of the zenith angles for sightline radiance (intensity) calculations.',  # This is a vector of values
    'u0u': 'The azimuthally averaged intensity at n_umu user specified angles umu. (units of e.g. mW/m^2/nm/sr if using the "atlas3" spectrum in the /data/solar_flux/ directory.)', # vector
    'uu': 'The radiance (intensity) at umu and phi user specified angles (unit e.g. mW/m^2/nm/sr if using the "atlas3" spectrum in the /data/solar_flux/ directory.)', # vector
    'uu_down': 'The downwelling radiances (intensity) at cmu and phi angles (unit e.g. mW/m^2/nm/sr if using the "atlas3" spectrum in the /data/solar_flux? directory.)', # vector
    'uu_up': 'The upwelling radiances (intensity) at cmu and phi angles (unit e.g. mW/m^2/nm/sr if using the "atlas3" spectrum in the /data/solar_flux/ directory.)', # vector
    'cmu': 'Computational polar angles from polradtran',
    'down_fluxI': 'The total (direct+diffuse) downward (down flux) irradiances (Stokes I component). Same units as extraterrestrial irradiance.',
    'up_fluxI': 'The total (direct+diffuse) upward (up flux) irradiances (Stokes I component). Same units as extraterrestrial irradiance.',
    'down_fluxQ': 'The total (direct+diffuse) downward (down flux) irradiances (Stokes Q component). Same units as extraterrestrial irradiance.',
    'up_fluxQ': 'The total (direct+diffuse) upward (up flux) irradiances (Stokes Q component). Same units as extraterrestrial irradiance.',
    'down_fluxU': 'The total (direct+diffuse) downward (down flux) irradiances (Stokes U component). Same units as extraterrestrial irradiance.',
    'up_fluxU': 'The total (direct+diffuse) upward (up flux) irradiances (Stokes U component). Same units as extraterrestrial irradiance',
    'down_fluxV': 'The total (direct+diffuse) downward (down flux) irradiances (Stokes V component). Same units as extraterrestrial irradiance.',
    'up_fluxV': 'The total (direct+diffuse) upward (up flux) irradiances (Stokes V component). Same units as extraterrestrial irradiance.',
    'eglo': 'Global downwelling irradiance',
    'enet': 'Global downwelling minus upwelling (net downward irradiance)',
    'esum': 'Global downwelling plus upwelling',
    'fdir': 'Direct actinic flux (scalar irradiance)',
    'fglo': 'Global actinic flux (scalar irradiance)',
    'fdn': 'Downwelling actinic flux (scalar irradiance)',
    'fup': 'Upwelling actinic flux (scalar irradiance)',
    'uavgglo': 'Total (global) mean diffuse intensity (radiance) = actinic flux/4pi',
    'f': 'Actinic flux (scalar irradiance)',
    'sza': 'Solar zenith angle [deg]',
    'n_air': 'Air refractive index',
    'zout': 'Altitude above ground [km]',
    'zout_sea': 'Altitude above sea level [km]',
    'sph_alb': 'Spherical albedo of the complete atmosphere',
    'albedo': 'Albedo',
    'heat': 'Heating rate in K/day',
    'n_xxx': 'Number density of gas xxx [cm^-3]',
    'rho_xxx': 'Mass density of gas xxx [kg/m^3]',
    'mmr_xxx': 'Mass mixing ratio of gas xxx [kg/kg]',
    'vmr_xxx': 'Volume mixing ratio of gas xxx [m^3/m^3]',
    'p': 'Pressure [hPa]',
    'T': 'Temperature [K]',
    'T_d': 'Dewpoint temperature [K]',
    'T_sur': 'Surface temperature [K]',
    'sur_temperature': 'Surface Temperature [K]',
    'theta': 'Potential temperature [K]',
    'theta_e': 'Equivalent potential temperature [K]',
    'rh': 'Relative humidity over water [%]',
    'rh_ice': 'Relative humidity over ice [%]',
    'c_p': 'Specific heat capacity of the air',
    'CLWC': 'Cloud liquid water content [kg/kg]',
    'CLWD': 'Cloud liquid water density [g/m^3]',
    'CIWC': 'Cloud ice water content [kg/kg]',
    'CIWD': 'Cloud ice water density [g/m^3]',
    'TCC': 'Total cloud cover [0-1]',
    'cth': 'Cloud top height [km]'
}
#  the following gas species can appear for _xxx in the above list of output variables.
gasSpecies = ['air', 'O3', 'O2', 'H2O', 'CO2', 'NO2', 'BRO', 'OCLO', 'HCHO', 'O4']

# Insert the additional items programatically
outVars = uvspecOutVars.copy()  # first make a copy
for thisVar in uvspecOutVars:
    if thisVar[-4:] == '_xxx':
        for gas in gasSpecies:
            outVars[thisVar.replace('xxx', gas)] = outVars[thisVar].replace('xxx', gas)  # insert in copy
        del outVars[thisVar]  # delete the template entry
uvspecOutVars = outVars  # copy the updated dict back to the original

# Some uvspec keywords can appear multiple times e.g.
# mol_modify can appear for many gas species individually.
# Keywords with multiple forms of this nature are called boson_keywords
boson_keywords = ['mol_modify', 'brdf_rossli', 'polradtran', 'brdf_hapke', 'sslidar', 'sdisort',
                  'ck_lowtran_absorption', 'mol_tau_file', 'aerosol_modify', 'aerosol_file', 'mixing_ratio',
                  'mol_file', 'brfdf_rpv', 'aerosol_file', 'brdf_cam', 'crs_file']

# For the output_user option, only certain variables are allowed as first and second index variables,
# These are wavelength (lambda, wvn), zout, zout_sea, p (pressure)
# The only two index variables that are allowed together is one height variable and one wavelength/number variable
# A height as well as a wavelength/number output should only be specified if there are both multiple
# wavelengths and multiple output heights. FOr the TZS solver, the cloud top height (cth) may be index.
indexVars = ['lambda', 'wavelength', 'wvl', 'wvn', 'zout', 'zout_sea', 'p', 'cth']

# Define a regexp for determining if a token of the zout keyword is a single level
# It is a single level if it is either a floating point number or the words boa, sur, cpt or toa (for surface,
# cold point tropopause and top of atmosphere)
re_isSingleOutputLevel = '(^[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$)|(^boa$)|(^sur$)|(^surface$)|(^cpt$)|(^toa$)'

def lookup_nearest_in_file(filename, values_and_offsets, column_number=0):
    """ Look up the nearest value in a free form text file of numeric data

    :param filename: File in which to look up the values
    :param column_number: Column number in which to search, starting from 0 - the default is 0
    :param values_and_offsets: An array with the values to search for (nominal values) in the first column and the
        minimum offsets to allow from the nominal values. Negative offsets will return the largest value which is
        less than the nominal value by at least the absolute value of the offset. Postitive offsets will return the
        value which is greater than the nominal value by at least the absolute offset value.
    :return: List of lookup values which are nearest to the nominal values by at least the offset values.

    In the tradition of libRadtran data files, lines starting with # are considered to be comments.
    The data in the column is assumed to be in monotonic, increasing order.
    Very large files should probably not be the subject of this function.

    The following excample will fetch the Thuillier solar spectrum wavelengths that span the range of
    385 nm to 955 nm with a margin of 2 nm on either side. This is useful when setting the uspec 'wavelength'
    keyword, which must give wavelengths that are actually listed in the source solar sepctrum file.
    >>> import libraddask.rad.librad as librad
    >>> wavelengths  = librad.lookup_nearest_in_file('data/Solar_irradiance_Thuillier_2002.txt', [[385.0, -2.0], [955, 2.0]])

    """
    values_and_offsets = np.asarray(values_and_offsets, dtype=np.float)
    # Read the whole file and select required column
    the_data = np.genfromtxt(filename)[:, column_number]
    nearest = np.zeros(values_and_offsets.shape[0])
    for ind_value, value in enumerate(values_and_offsets[:, 0]):
        ind_nearest = np.sum(the_data <= value + values_and_offsets[ind_value, 1])   # split the file
        if values_and_offsets[ind_value, 1] < 0.0:
            ind_nearest -= 1
        nearest[ind_value] = the_data[ind_nearest]
    return nearest

def angstrom_law(wavelength, alpha, beta):
    """ Calculation of aerosol optical thickness according to the Angstrom law.
    The Angstrom law is a simple power law that typifies aerosol optical thickness (aka optical depth) variation
    with wavelength. The Angstrom law is expressed as
    .. math::

        \\tau_{aer}=\\beta\\lambda^{-\\alpha}

    The Angstrom law can be used to set aerosol optical thickness in libRadtran/uvspec using the `aerosol_angstrom`
    keyword.

    :param wavelength: Wavelength(s) at which to compute the aerosol optical thickness. If any of the wavelengths
        are greater than 100, then wavelengths are assumed to be in nm, otherwise wavelengths are assumed to be in
        microns.
    :param alpha: The Angstrom alpha exponent.
    :param beta: The Angstrom beta (aerosol optical thickness at a wavelength of 1000 nm) parameter
    :return: Aerosol optical thickness at given wavelengths with Angstrom alpha and beta parameters as provided.

    .. seealso::
        librad.angstrom_law_fit, the libRadtran manual
    """
    if np.any(wavelength > 100.0):
        wavelength = wavelength / 1000.0
    return beta * wavelength ** (-alpha)

def angstrom_law_fit(wavelength, aot):
    """ Uses scipy.optimize to fit the Angstrom power law to an array of aerosol optical thickness values given
    at an array of wavelengths.

    :param wavelength: Wavelengths at which the aerosol optical thickness (aka optical depth) is provided in the aot
        input. If any of the wavelengths is larger than 100, then wavelengths are assumed to be in nm, otherwise
        wavelengths are assumed to be in microns.
    :param aot: Aerosol optical thickness at the given wavelengths
    :return: Angstrom alpha and beta parameters that best fit the input AOT data in the least squares sense.

    .. seealso::
        librad.angstrom_law
    """
    from scipy.optimize import curve_fit
    if np.any(wavelength > 100.0):
        wavelength = wavelength / 1000.0
    popt, pcov = curve_fit(angstrom_law, wavelength, aot)
    return popt[0], popt[1]

def king_byrne_formula(wavelength, alpha_0, alpha_1, alpha_2):
    """ The King Byrne formula for aerosol optical depth variation with wavelength.
    The King Byrne formula is

    .. math::

        \\tau_{aer}=e^{\\alpha_{0}}\\lambda^{\\alpha_{1}}\\lambda^{-\\alpha_{2}}

    :param wavelength: Wavelengths at which the aerosol optical thickness (aka optical depth) is provided in the aot
        input. If any of the wavelengths is larger than 100, then wavelengths are assumed to be in nm, otherwise
        wavelengths are assumed to be in microns.
    :param alpha_0: The :math:`\\alpha_0` parameter
    :param alpha_1: The :math:`\\alpha_1` parameter
    :param alpha_2: The :math:`\\alpha_2` parameter
    :return: Aerosol optical thickness at given wavelengths calculated with the King Byrne formula

    """
    if np.any(wavelength > 100.0):
        wavelength = wavelength / 1000.0
    return np.exp(alpha_0) * wavelength**alpha_1 * wavelength**(alpha_2 * np.log(wavelength))

def king_byrne_formula_fit(wavelength, aot):
    """ Uses scipy.optimize to fit the King Byrne formula to an array of aerosol optical thickness values
    provided at the given wavelengths

    :param wavelength: Wavelengths at which the aerosol optical thickness (aka optical depth) is provided in the aot
        input. If any of the wavelengths is larger than 100, then wavelengths are assumed to be in nm, otherwise
        wavelengths are assumed to be in microns.
    :param aot: Aerosol optical thickness at the given wavelengths
    :return: The :math:`\\alpha_0`, :math:`\\alpha_1` and :math:`\\alpha_2`
    """
    from scipy.optimize import curve_fit
    if np.any(wavelength > 100.0):
        wavelength = wavelength / 1000.0
    alpha_parameters, pcov = curve_fit(king_byrne_formula, wavelength, aot)
    return alpha_parameters[0], alpha_parameters[1], alpha_parameters[2]

def koschmieder_vis(ext_550=None, aot_550=None, scale_height=None, rayleigh_ext=0.01159):
    ''' Compute visibility in km using the Koschmieder relationship
     .. math::

        V_K = \\frac{\\ln 50}{\\epsilon_{aer} + \\epsilon_{ray}}

     Or,

     .. math::

        V_K = \\frac{\\ln 50}{\\tau_{aer}/ H + \\epsilon_{ray}}

     where :math:`V_K` is the visibility in km, :math:`\\epsilon_{aer}` is the aerosol extinction coefficient at
     550 nm in units of inverse km, :math:`\\tau_{aer}` is the vertical aerosol optical depth at 550 nm,
     and :math:`\\epsilon_{aer}` is the Rayleigh extinction coefficient, in units of inverse km.

     :math:`H` is the scale height or effective mixing layer height assuming that all aerosols are in the mixing
     layer. The boundary layer height is usually a good approximation.

    Provide either ``ext_550`` OR ``aot_550`` together with ``scale_height``.
    The Rayleigh extinction coefficient is optional.

    :param ext_550: aerosol extinction coefficient :math:`\\tau_{aer}` at 550nm in units of inverse km
    :param aot_550: vertical aerosol optical thickness
    :param scale_height: effective boundary layer height (ABL/mixing height if all aerosols in the mixing layer) in km
    :param rayleigh_ext: Rayleigh extinction coefficient, defaults to 0.01159 per km.
    :return: Visibility in km according to the Koschmieder relationship

    The relationship used here is taken from the MODTRAN manual.
    '''
    try:
        vis = np.log(50) / (ext_550 + rayleigh_ext)
        if (aot_550 is not None) or (scale_height is not None):
            warnings.warn('aot_550 and/or scale_height inputs igonored in librad.koschmieder_vis calculation.')
    except TypeError:
        try:
            vis = np.log(50) / (aot_550 / scale_height + rayleigh_ext)
        except TypeError:
            warnings.warn('Invalid input combination for librad.koschmieder_vis calculation provided.')
            return None
    return vis

class Case(object):
    """ Class which encapsulates a run case of libRadtran/uvspec.
    This class has methods to read libRadtran/uvspec input files, write uvspec input files, run uvspec in parallel on
    multiple compute nodes and read uvspec output files. An important use-case is that of reading a uvspec input
    file called the "base case", altering the parameters of particular option keywords and then running the case
    and reading the outputs. This class is also used by the RadEnv class which encapsulates a radiant environment.
    Construction of radiant environment maps typically requires running an array of librad.Case instances.
    """
    # Definitions of some of the possible uvspec output variables

    def __init__(self, casename='', filename=None, optionlist=None):
        """ Instantiate a libRadtran/uvspec case, typically by reading a uvspec .INP file.

        :param casename: A user-defined name for the libRadtran/uvspec case
        :param filename: An optional filename from which to read the libRadtran/uvspec input
        :param optionlist: A list of option keywords and parameteres (tokens). The keyword existence
            is verified. Besides that, no error checking is performed automatically.
        :return: None

        .. todo::
            Implement checking of uvspec input keyword tokens (parameters)

        >>> # Read a libRadtran/uvspec case from a .INP file and display the expanded input
        >>> import libraddask.rad.librad as librad
        >>> libRadCase = librad.Case(filename='./examples/UVSPEC_AEROSOL.INP')  # Read uvspec input and expand includes, if any
        >>> print(libRadCase)   # This prints the uvspec input file, compare to contents of UVSPEC_AEROSOL.INP
        atmosphere_file ../data/atmmod/afglus.dat
        source solar ../data/solar_flux/atlas_plus_modtran
        mol_modify O3 300. DU
        day_of_year 170
        albedo 0.2
        sza 32.0
        rte_solver disort
        number_of_streams 6
        wavelength 299.0 341.0
        slit_function_file ../examples/TRI_SLIT.DAT
        spline 300 340 1
        quiet
        aerosol_vulcan 1
        aerosol_haze 6
        aerosol_season 1
        aerosol_visibility 20.0
        aerosol_angstrom 1.1 0.2
        aerosol_modify ssa scale 0.85
        aerosol_modify gg set 0.70
        aerosol_file tau ../examples/AERO_TAU.DAT

        """
        self.name = casename
        self.error_txt = []
        self.stderr = ''  # Error output from the uvspec run may be read into this
        self.run_return_code = -1  # Will be set when uvspec run is executed
        self.purge = True  # This will purge uvspec input, output and error files after run, unless set False
        self.options = []  # options is a list [option_name (string), option_tokens (list of strings),
        self.tokens = []  # option keyword parameters (tokens)
        self.optionobj = []  # option object from writeLex
        self.filorigin = []
        self.solver = 'disort'  # default, modified b y the rte_solver option keyword
        self.fluxline = ['wvl', 'edir', 'edn', 'eup', 'uavgdir', 'uavgdn', 'uavgup']  # default output
        self.wvl = []  # wavelengths, wavenumbers and output levels difficult to ascertain to start with
        self.wvn = []
        self.sza = 0.0  # default solar zenith angle for this case
        self.zout = np.zeros(1)  # If no zout (or equivalent) is given, uvspec produces output at the surface
        self.zout_sea = np.zeros(1)  # altitudes above seal level, same as zout unless ground 'altitude' provided
        self.pressure_out = []  # Only populated if pressure_out keyword is given
        self.output_user = ''  # set with the output_user keyword
        self.output_quantity = 'radiance' # default is radiances and irradiances in units determined by input solar file
        self.source = ''  # 'solar' or 'thermal'
        self.source_file = ''  # TOA irradiance source file
        self.n_umu = 0  # number of zenith angles (actually cosine of zenith angle)
        self.umu = []  # zenith angles for radiance calculations
        self.uu = np.array([])  # This will be populated if there is radiance data output by the case
        self.n_phi = 0  # number of azimuth radiance angles
        self.phi = [] # np.zeros(1)  # azimuth angles for radiance calculations
        self.n_levels_out = 1  # Assume only one output level, unless zout/zout_sea/pressure_out keyword is used.
        self.levels_out = ['boa']  # Tokens provided on an output level keyword (zout, zout_sea, pressure_out), default
        self.levels_out_type = 'zout'  # Assume that type of levels out data is altitude above surface
        self.n_wvl = '?'  # number of wavelengths is difficult to ascertain to start with
        self.n_stokes = 1  # default number of stokes parameters for polradtran
        self.stokes = ['I']  # default is to compute intensity component only (all solvers, including polradtran)
        self.radND = []  # N-dimensional radiance data (umu, phi, wvn, zout, stokes)
        self.rad_units = ['mW', 'm^2', 'nm']  # radiance/irradiance units, could be K for brightness temperatures
        self.altitude = np.zeros(1)  # Default surface (ground) height above sea level
        self.mol_abs_param = 'reptran'  # REPTRAN is the default molecular absorption parametrization
        self.spectral_res = 'coarse'  # The coarse 15 cm^-1 spectral resolution is the default
        self.reptran_channel = ''  # Used with mol_abs_param keyword
        self.spectral_channels = []  #  These are the spectral channels, such as for 'kato', 'fu' band models.
        self.spectral_axis = 'wvl'  # By default RT calculations are spectral, with wavelength as the variable
        self.wavelength_index = []  # set by 'wavelength_index' keyword
        self.wavelength_index_range = []  # for use where range() would be applicable
        self.wavelength_grid_file = ''  # This is the name of the wavelength grid file if the option is used
        self.wavelength_grid = None  # Only gets set by self.set_wavelength_grid(wavelength_grid_points)
        self.wavelength_file = ''  # This associated with the keyword 'wavelength filename'
        self.output_process = 'none'  # Default is to output spectral data. See output_process in uvspec manual
        self.has_ice_clouds = False  # Default, unless discovered otherwise
        self.has_water_clouds = False  # ditto
        self.has_clouds = False  # either water or ice clouds
        if filename is not None:
            if not filename:
                # Open a dialog to get the filename
                filename = easygui.fileopenbox(msg='Please select a uvspec input file.', filetypes=["*.INP"])
            elif filename[-4:].lower() != '.inp':
                filename += '.INP'
            opdata, line_nos, path = Case.read(path=filename)
            self.infile = path
            self.outfile = self.infile[:-4] + '.OUT'
            self.errfile = self.infile[:-4] + '.ERR'
            # option_object, source_file_nos (filename, line_number)]
            # Process the results into lists of options
            for (optnumber, option) in enumerate(opdata):
                self.options.append(option[0])  # the option keyword (string)
                self.tokens.append(option[1:])  # The tokens following the keyword (list of strings)
                self.filorigin.append(line_nos[optnumber])  # The file origin of this keyword
                if self.options[optnumber] in uvsOptions:
                    self.optionobj.append(uvsOptions[self.options[optnumber]])  # The option object (dict lookup)
                else:
                    print('Warning, keyword option ' +  option[0] + ' not found in options library.')
                # Make any possible preparations for occurance of this keyword
                self.prepare_for_keyword(option[0],option[1:])
        if not self.name:  # set the name as the basname of the filename, excluding the extension
            self.name = os.path.basename(self.infile)[:-4]
        #TODO Build the case from the option list ?

    def prepare_for_solver(self, tokens):
        self.solver = tokens[0]  # First token gives the solver
        if any([self.solver == thesolver for thesolver in ['disort', 'disort2', 'sdisort',
                                                           'spsdisort', 'fdisort1', 'fdisort2']]):
            self.fluxline = ['wvl', 'edir', 'edn', 'eup', 'uavgdir', 'uavgdn', 'uavgup']  # default output
        elif any([self.solver == thesolver for thesolver in ['twostr','rodents']]):
            self.fluxline = ['wvl', 'edir', 'edn', 'eup', 'uavg']
        elif self.solver == 'polradtran':
            self.prepare_for_polradtran()
        elif self.solver == 'sslidar':
            self.fluxline = ['center_of_range', 'number_of_photons', 'lidar_ratio']
            self.rad_units = ['count', '', '']
        elif self.solver == 'montecarlo':
            self.solver = 'mystic'  # Same as mystic, use mystic

    def prepare_for_source(self, tokens):
        """ Prepare for source, particularly units of various kinds, depending on the source.
        libRadtran/uvspec source options are mainly 'solar' and 'thermal' with some additional options
        for units.

        :param tokens: uvspec 'source' keyword option parameters (tokens)
        :return: None
        """

        # TODO : More work required here - see libRadtran manual under source keyword
        # TODO : Also need to verify default radiance output units (mW/m^2/sr/nm or W/m^2/sr/nm)
        if tokens[0] == 'solar':  # dealing with solar spectrum
            self.source = 'solar'
            if len(tokens) > 1:  # solar source file provided
                self.source_file = os.path.basename(tokens[1])  # extract the filename only
                # try to establish the radiometric units based on the filename
                if self.source_file in sourceSolarUnits:
                    self.rad_units = sourceSolarUnits[self.source_file]
                elif len(tokens) == 3:
                    self.rad_units[2] = {'per_nm': 'nm', 'per_cm-1':'cm^-1',
                                      'per_band': ''}[tokens[2]]
        elif tokens[0] == 'thermal':
            self.source = 'thermal'
            if self.output_quantity == 'radiance':
                self.rad_units = ['W', 'm^2', 'cm^-1']

    def prepare_for_mol_abs_param(self, tokens):
        """ Make preparations for the desired molecular absorption parametrization. There are various
        molecular absorption options in libRadtran/uvspec and the options will certainly evolve.
        The default mode is to perform a spectral calculation (fine wavelength grid) with output
        of spectral radiances and irradiances. Some of the other modes, such as `crs` also output
        spectral data. The `crs` mode actually switches off molecular line absorption and considers
        only spectrally continuous scattering and absorption. This is really only good for the UV/blue
        spectrum. The `reptran_channel` mode and the correlated-k modes (`kato` variants, `fu`, `avhrr_kratz`
        and `lowtran`/`sbdart`) do not produce spectral radiances and irradiances. They produce band quantities
        which may even be summed using the `output_process sum` directive. A sub-range of correlated-k
        bins/channels can be selected using the `wavelength_index` directive.

        See the libRadtran manual for further information on the relevant options.
        :return:
        """
        tokens = [token.lower() for token in tokens]
        if tokens[0] == 'reptran':
            self.mol_abs_param = 'reptran'
            self.spectral_axis = 'wvl'  # RT calculations are spectral as opposed to band
            if len(tokens) > 1:
                if tokens[1] == 'coarse':
                    self.spectral_res = 'coarse'
                elif tokens[1] == 'medium':
                    self.spectral_res = 'medium'
                elif tokens[1] == 'fine':
                    self.spectral_res = 'fine'
                else:
                    warnings.warn('Invalid REPTRAN spectral resolution qualifier found in mol_abs_param directive.')
        elif tokens[0] == 'crs':  # This option actually switches off line-absorption and only considers continua
            self.mol_abs_param = 'crs'  # Only good for UV really, but will output in any spectral region
            self.spectral_res = ''
            self.spectral_axis = 'wvn'
        elif tokens[0] == 'reptran_channel':
            self.mol_abs_param = 'reptran_channel'
            self.reptran_channel = tokens[1]
            self.spectral_res = ''
            self.spectral_axis = 'chn'
            # Not sure if spectral or band-integrated quantities are provided
        elif tokens[0] == 'kato':
            self.mol_abs_param = 'kato'
            self.spectral_res = ''
            self.spectral_axis = 'chn'
            self.rad_units[2] = ''  # band integrated quantity out
        elif tokens[0] == 'kato2':
            self.mol_abs_param = 'kato2'
            self.spectral_res = ''
            self.spectral_axis = 'chn'
            self.rad_units[2] = ''  # band integrated quantity
        elif tokens[0] == 'kato2.96':
            self.mol_abs_param = 'kato2.96'
            self.spectral_res = ''
            self.spectral_axis = 'chn'
            self.rad_units[2] = ''  # band integrated quanitity
        elif tokens[0] == 'fu':
            self.mol_abs_param = 'fu'
            self.spectral_res = ''
            self.spectral_axis = 'chn'
            self.rad_units[2] = ''  # band integrated quantity
        elif tokens[0] == 'avhrr_kratz':
            self.mol_abs_param = 'avhrr_kratz'
            self.spectral_res = ''
            self.spectral_axis = 'chn'
            self.rad_units[2] = ''  # band integrated quantity
        elif tokens[0] == 'lowtran' or tokens[0] == 'sbdart':
            self.mol_abs_param = 'lowtran'
            self.spectral_res = ''
            self.spectral_axis = 'chn'  # These are pseudo-spectral with 20 cm^-1 resolution
            # Spectral quantities are provided
        else:
            warnings.warn('Unknown mol_abs_param type encountered.')

    def prepare_for_output_process(self, tokens):
        """ Prepare for effects of the `output_process` keyword in the libRadtran/uvspec input file.

        :param tokens: Keyword tokens of the uvspec `output_process` keyword.
        :return:
        """
        self.output_process = tokens[0]
        process = tokens[0].lower()
        if process == 'sum':  # the units become per band (?)
            self.output_process = 'sum'
            self.rad_units[2] = ''  # TODO : Check that this is correct somehow
        elif process == 'integrate':
            if not self.spectral_axis == 'wvl':
                warnings.warn('Option output_process integrate probably not valid with band quantities')
            self.output_process = 'integrate'
            self.rad_units[2] = ''  # The third element in the rad_units list gives the spectral variable
        elif process == 'per_nm':
            self.output_process = 'per_nm'
            self.rad_units[2] = 'nm'
        elif process == 'per_cm-1':
            self.output_process = 'per_cm-1'
            self.rad_units[2] = 'cm^-1'
        elif process == 'per_band':
            self.output_process = 'per_band'
            self.rad_units[2] = ''
        elif process == 'none':
            self.output_process = 'none'
            # TODO : Probably need to revert to other defaults here
        else:
            warnings.warn('Invalid output_process option encountered in uvspec input file.')

    def prepare_for_keyword(self, keyword, tokens):
        """ Make any possible preparations for occurrences of particular keywords
        :param keyword: The uvspec option keyword (string)
        :param tokens: The parameters (tokens) for the keyword as a list of strings
        :return:
        """
        # Prepare for different output formats, depending on the solver
        if keyword == 'rte_solver':
            self.prepare_for_solver(tokens)
        if keyword == 'source':
            self.prepare_for_source(tokens)
        # Prepare for radiances
        if keyword == 'umu':
            self.n_umu = len(tokens)  # The number of umu values
            self.umu = np.array(tokens).astype(np.float64)
            self.pza = xd_identity(np.arccos(self.umu), 'pza', 'rad')
            self.paz = xd_identity(np.deg2rad(self.phi), 'paz', 'rad')
        if keyword == 'phi':
            self.n_phi = len(tokens)
            self.phi = np.array(tokens).astype(np.float64)
            self.paz = xd_identity(np.deg2rad(self.phi), 'paz', 'rad')
        if keyword == 'output_user':
            self.output_user = [token.replace('lambda', 'wvl') for token in tokens]  # lambda is a keyword
            self.output_user = [token.replace('wavenumber', 'wvn') for token in self.output_user]  # abbreviate wavenumber
            self.output_user = [token.replace('wavelength', 'wvl') for token in self.output_user]  # abbreviate wavelength
            self.fluxline = []
        if keyword == 'polradtran' and tokens[0] == 'nstokes':
            self.n_stokes = int(tokens[1])
            if self.n_stokes == 1:
                self.stokes = ['I']
            elif self.n_stokes == 3:
                self.stokes = ['I', 'Q', 'U']
            elif self.n_stokes == 4:
                self.stokes = ['I', 'Q', 'U', 'V']
            else:
                raise ValueError('uvspec input polradtran nstokes must be 1, 3 or 4.')
            self.prepare_for_polradtran()
        if (keyword == 'zout' or keyword == 'zout_sea' or keyword == 'pressure_out' or
            keyword == 'tzs_cloud_top_height'):  # Determine number of output levels
            self.levels_out_type = keyword
            if all([re.match(re_isSingleOutputLevel, token.lower()) for token in tokens]):
                self.n_levels_out = len(tokens)
                self.levels_out = tokens
            else:  # Number of output levels is not known, will have to auto-detect
                self.n_levels_out = 0  # This indicates that the number of output levels is not known
        if keyword == 'output_quantity':  # Try to set output units
            self.output_quantity = tokens[0]
            if self.output_quantity == 'brightness':
                self.rad_units = ['K', '', '']
            elif self.output_quantity == 'reflectivity':
                self.rad_units = ['', '', '']
            elif self.output_quantity == 'transmittance':
                self.rad_units = ['', '', '']
        if keyword == 'heating_rate':
            self.output_user = ['zout', 'heat']
            self.fluxline = []  #TODO note that heating rate outputs for multiple wavelengths have a header line
        if keyword == 'print_disort_info':
            self.fluxline = '?'  # this output format is unknown or too complex to handle
        if keyword == 'altitude':
            self.altitude = np.float64(tokens[0])  # This is the ground altitude above sea level
        if keyword == 'mol_abs_param':
            self.prepare_for_mol_abs_param(tokens)
        if keyword == 'wavelength_index':
            self.wavelength_index = [int(tokens[0]), int(tokens[1])]
            self.wavelength_index_range = range(int(tokens[0]), int(tokens[1]) + 1)
        if keyword == 'output_process':
            self.prepare_for_output_process(tokens)
        if keyword == 'wc_file':
            self.has_water_clouds = True
            self.has_clouds = True
        if keyword == 'ic_file':
            self.has_ice_clouds = True
            self.has_clouds = True
        if keyword == 'profile_file':
            if tokens[0] == 'wc':
                self.has_water_clouds = True
                self.has_clouds = True
            elif tokens[0] == 'ic':
                self.has_ice_clouds = True
                self.has_clouds = True
        if keyword == 'sza':
            self.sza = np.float64(tokens[0])

    def prepare_for_polradtran(self):
        """ Prepare for output from the polradtran solver

        :return:
        """
        self.fluxline = ['wvl']
        for stokes in self.stokes:
            self.fluxline.extend(['down_flux' + stokes, 'up_flux' + stokes])
        warnings.warn('The polradtran solver does not produce direct solar irradiances and will only produce' +
                      ' output if the atmosphere file contains the altitudes specified by zout (see "zout" ' +
                      ' in the uvspec manual). ')

    def append_option(self, option, origin=('user', None)):
        """ Append a libRadtran/uvspec options to this uvspec case. It will be appended at the end of the file

        :param option: A list containing the keyword and keyword parameters (tokens). Boson keywords must be passed
            in "welded". e.g. x.append_option(['mol_modify O3', '270.0',  'DU'])
        :param origin: A 2-tuple giving the origin of the option and a "line number" reference. Default ('user', None)
            uvspec options.
        :return:
        """
        # May have a boson keyword passed in, so split at spaces
        option0split = option[0].split()
        if len(option0split) > 1:
            # Lookup the option in the available options
            the_option = uvsOptions[option0split[0]]  # this will raise ValueError if it does not exist
        else:
            the_option = uvsOptions[option[0]]
        # Check to see if it has logicals and if so, check first token against the logicals
        #TODO checking of the logicals
        # Check if this is a
        self.optionobj.append(the_option)  # The option object
        self.options.append(option[0])  # the option keyword (string)
        self.tokens.append(option[1:])  # The tokens following the keyword (list of strings)
        self.filorigin.append(origin)  # The origin of this keyword
        # Make any possible preparations for occurance of this keyword
        self.prepare_for_keyword(option0split[0], option0split[1:] + option[1:])

    def alter_option(self, option, origin=('user', None)):
        """ Alter the parameters of a uvspec input option. If the option is not found, the option is appended with
        append_option instead.

        :param option: List of keyword and tokens (parameters) to provide to the option keyword (list of strings).
        :param origin: A 2-tuple noting the "origin" of the change to this keyword. Default ('user', None)
        :return: None
        """
        # check if this is a boson_keyword
        if option[0] in boson_keywords:
            # Weld together the first tokens with a space in between
            option[0] = option[0] + ' ' + option[1]
            del option[1]
        try:
            ioption = self.options.index(option[0])
        except ValueError:  # The option is not currently being used in this case
            self.append_option(option, origin)  # just append the option if not found
        else:
            self.tokens[ioption] = option[1:]  # The tokens following the keyword (list of strings)
            self.filorigin[ioption] = origin  # The origin of this keyword
            option0split = option[0].split()
            self.prepare_for_keyword(option0split[0], option0split[1:] + option[1:])

    def set_option(self, *superlist):
        """ Set a uvspec input option with flexible input format

        :param superlist: All set_option arguments are collected into a list which are coverted to strings
            (if not already strings), splitting each resulting string at spaces.
        :return: None

        Example:

        >>> import libraddask.rad.librad as librad
        >>> libRadCase = librad.Case(casename='MyExample')  # Creates a blank libRadtran/uvspec case
        >>> libRadCase.set_option('source solar', '../data/solar_flux/atlas_plus_modtran')
        >>> libRadCase.set_option('wavelength', 300, 400)  # can provide literal numerics
        >>> print(libRadCase)

        """
        from itertools import chain
        optionlist = [str(element).split() for element in superlist]
        # now flatten the list
        option = list(chain(*optionlist))
        # and use alter_option to do the work
        self.alter_option(option)

    def del_option(self, option, all=True):
        """ Delete a uvspec input option matching the given option.

        :param option: Keyword of option to be deleted
        :param all: A flag indicating if all matching options must be deleted or only the first occurrence. The
        default is to delete all matching occurrences.
        :return: True if an option was deleted or False if not
        """
        #TODO consider providing warning if options does not exist
        deletedsomething = False
        while option in self.options:
            deletedsomething = True
            ioption = self.options.index(option)
            self.options.pop(ioption)
            self.tokens.pop(ioption)
            self.filorigin.pop(ioption)
            self.optionobj.pop(ioption)
        if deletedsomething:
            # Run through all options and reconstruct preparations
            for (ioption, option) in enumerate(self.options):
                self.prepare_for_keyword(self.options[ioption], self.tokens[ioption])
        return deletedsomething

    def set_wavelength_grid(self, wvl_grid):
        """ Set the internal RT solver wavelength grid for this librad.Case

        This method sets the internal grid for the RT calculations in libRadtran. Study the libRadtran manual
        before trying to use this option. Once this function has been used, the option wavelength_grid_file
        will be activated and the grid file will be *casename_wvl_grid.dat*, where the casename is the name of
        the librad.Case. When the case is Run, the wavelength grid file will be written out to the current
        directory before uvspec is executed. By default (unless purge is disabled), the file will also
        automatically be deleted once the tun terminates.

        :param wvl_grid: An np.array of floats that provide the wavelengths in nm at which to perform the
            internal radiative transfer calculation.
        :return: None

        .. seealso: unset_wavelength_grid
        """
        self.wavelength_grid_file = self.name + '_wvl_grid.dat'
        self.wavelength_grid = np.vstack(wvl_grid.flatten())  # should be a numpy array of floats, force to column
        self.set_option('wavelength_grid_file', self.wavelength_grid_file)

    def unset_wavelength_grid(self):
        """ Remove the wavelength grid file and associated grid data for this librad.Case.

        :return: None

        .. seealso: set_wavelength_grid()
        """
        self.wavelength_grid_file = ''
        self.wavelength_grid = None
        self.del_option('wavelength_grid_file')

    def set_view_geometry(self, sza=0.0, saa=0.0, oza=0.0, oaa=0.0):
        ''' Set the Case view geometry relative to the sun

        Sets the sun and observer zenith and azimuth angles. Note that libRadtran/uvspec sets the light
        propagation azimuth angle in the ``phi`` and ``phi0`` inputs. This is 180$^\circ$ different in
        azimuth to the direction of viewing. The ``phi0`` input is the direction of solar light propagation
        from the sun and is :math:`180^\\circ` from the saa input here. However, the ``phi`` input is the direction
        of light propagation from target to sensor and is therefore *the same* as the oaa input to this
        function.

        The ``umu`` input to libRadtran/uvspec is the cosine of the zenith angle of light propagation
        from target to sensor. Hence ``umu`` is positive for downward-looking (upward-propagating) cases.

        Also, from the libRadtran manual :
        phi = phi0 indicates that the sensor looks into the direction of the sun, while
        phi-phi0 = 180 deg means that the sun is behind the sensor.

        :param sza: solar zenith angle in degrees from the zenith (range 0 to 90 degrees)
        :param saa: solar azimuth angle in degrees from north through east (range zero to )
        :param oza: observation zenith angle in degrees from the zenith. Note that this is the zenith angle
            of a vector pointing from the target/scene to the observer (satellite/sensor).
        :param oaa: observation azimuth angle in degrees from north through east. Note again that this is the
            azimuth angle of a vector pointing from the target/scene to the observer.
        :return: None
        '''
        self.set_option('sza', sza)  # deg. This one is straightforward
        # Now when entering solar and observation zenith angles, it is necessary to provide the azimuth of light propagation
        # rather than the azimuth of the view direction, which is 180 deg different
        phi0 = saa + 180.0
        if phi0 > 360.0:
            phi0 -= 360.0  # Keep phi0 in the range from 0.0 to 360.0
        self.set_option('phi0', phi0)  # solar radiation propagation azimuth from north through east
        self.set_option('phi', oaa)  # This is the azimuth of the satellite as seen from the target - also azimuth of light propgation
        self.set_option('umu', np.cos(np.deg2rad(oza))) # For downward-looking (upward propagating) umu is positive

    @staticmethod
    def read(path, includes_seen=[]):
        """ Reads a libRadtran input file. This will construct the libRadtran case from the contents of the .INP file
        Adapted from code by libRadtran developers.

        :param path: File path from which to read the uvspec input
        :param includes_seen: List of files already included (for recursion purposes to avoid infinite include loops)
        :return: data, line_nos, path
          where data is the full data in the file with includes, line_nos shows the source of every line and
          path is the path to the main input file.
        """
        # Get the full path
        path = os.path.abspath(path)
        folder = os.path.dirname(path)
        # print(includes_seen)
        with open(path, 'rt') as INPfile:
            opdata = INPfile.readlines()  # read in the entire file and process in memory afterwards
        opdata = [line.split() for line in opdata]  # Split lines into keywords and option parameters (tokens)
        line_nos = [(path, i) for i in range(1, len(opdata)+1)]
        # Remove lines with comments and merge continuous lines
        line = 0
        while line < len(opdata):  # Remove empty lines, comments and merge options split over multiple file lines
            # Skip empty lines
            if not opdata[line]:
                opdata.pop(line)
                line_nos.pop(line)
                continue
            # Skip comments  #TODO save comments as well in the librad.Case
            if opdata[line][0].startswith("#"):
                opdata.pop(line)
                line_nos.pop(line)
                continue
            # Remove comments from the line  #TODO save comments into the librad.Case
            elif [True for word in opdata[line] if (word.find("#") != -1)]:
                tmp = []
                for word in opdata[line]:
                    pos = word.find("#")
                    if pos != -1:
                        if word != "#":
                            tmp.append(word[:pos])
                        break
                    else:
                        tmp.append(word)
                opdata[line] = tmp
            # continuous line
            elif opdata[line][-1].endswith("\\"):
                opdata[line][-1] = opdata[pos][-1][:-1] # remove the \
                    # if the \ was preceded and continued by whitespace
                if opdata[line][-1] == "":
                    opdata[line].pop()
                    line_nos.pop()
                opdata[line].extend(opdata[line + 1])
            else:
                line += 1
        # Get the includes and include them at the point where the include keyword appears
        all_opdata = []
        all_line_nos = []
        this_includes = includes_seen[:] # These are the include files seen up to this point in the recursion
        for line in range(len(opdata)):  # Run through all lines in the file and perform include substitutions
            if opdata[line][0] == 'include':  # Have found an include file, read it and append to all_opdata
                include_file = opdata[line][1]  # Obtain the filename
                include_file = os.path.normpath(os.path.join(folder, include_file))  # Extend to full path name
                if this_includes.count(include_file):  # This file has been included before
                    # print(this_includes, include_file)
                    raise IOError('Attempted to include the same uvspec input file more than once in parent.')
                this_includes.append(include_file)  # Add the filename to the list of included files already seen
                # Read the data from the included file
                inc_opdata, inc_line_nos, x = Case.read(include_file, this_includes)  # Recursion
                # Insert the data at this point in the file option and line number lists
                all_opdata.extend(inc_opdata)
                all_line_nos.extend(inc_line_nos)
            else:  # Not an include option, just append
                all_opdata.append(opdata[line])
                all_line_nos.append(line_nos[line])
        return all_opdata, all_line_nos, path

        # The following is the old code that inserted all include files at the end, which is not the same
        # behaviour as a C #INCLUDE statement, which inserts at the point where the #INCLUDE appears.
        # buff = 0
        # this_includes = []
        # for line in range(len(opdata)):
        #     if (opdata[line + buff][0] == "include" and len(opdata[line + buff]) == 2):
        #         include_file = opdata.pop(line + buff)[1]  # Obtain the filename of the included file and pop option
        #         include_file = os.path.normpath(os.path.join(folder, include_file))  # Extend to full path name
        #         this_includes.append(include_file)  # Add the filename to the list of included files at this level
        #         line_nos.pop(line + buff)  # Also pop the line numbers from the list of line numbers
        #         buff -= 1  # Take into account that a line has been removed
        # for include_path in this_includes:
        #     if not os.path.exists(include_path):
        #         msg = "Include file '%s' in '%s' does not exist." % (include_path, path)
        #         raise IOError(msg)
        #         #print " " * len(includes) + msg
        #     # If the file has been included before.
        #     elif (this_includes + includes_seen).count(include_path) != 1:  # Count number of times file is included
        #         msg = "File %s included more than once in %s. Please fix this." % (include_path, path)
        #         raise IOError(msg)
        #         #self.error_txt.append(msg)
        #         #print " " * len(includes) + msg
        #     else:
        #         include_data = Case.read(include_path, includes_seen + this_includes)
        #         # The include file might contain errors and return None.
        #         if include_data:
        #             opdata.extend(include_data[0])  #TODO surely an insert rather than extend ?
        #             line_nos.extend(include_data[1])
        # print(opdata)
        # print(line_nos)
        # return opdata, line_nos, path

    def __repr__(self):
        """ libRadtran/uvspec input data

        :return: The uvspec input data as it would appear in an input file.
        """
        uvsinp = []
        for (ioption, keyword) in enumerate(self.options):
            theTokens =  re.sub('[\[\],]', '', (' '.join(self.tokens[ioption])))  # Remove any square brackets or commas
            optionLine = (keyword + ' ' + theTokens).replace('  ', ' ')  # replace any double spaces with single spaces
            uvsinp.append(optionLine)  #TODO add comments
        return '\n'.join(uvsinp)

    def write(self, filename=''):
        """ Write libRadtran/uvspec input to a file (.INP extension by default.
        If the filename input is given as '', a file save dialog will be presented

        :param filename: Filename to which to write the uvspec input
        :return:
        """
        if not filename:
            filename = easygui.filesavebox(msg='Please save the uvspec input file.', filetypes=["*.INP"])
        if filename[-4:].lower() != '.inp':
            filename += '.INP'
        with open(filename, 'wt') as uvINP:
            alldata = self.__repr__()
            uvINP.write(alldata)

    def distribute_flux_data(self, fluxdata):
        """ Distribute flux/user data read from uvspec output file to various data fields.
        This method will look at `output_user` options and attempt to assign flux/user data in a sensible way.
        .. note::

            There are potentially uvspec output formats that are not possible to process or to assign correctly.
            These are typically cases in which it is not possible to determine from the .INP and/or .OUT file
            how this data should be assigned.

        :param fluxdata: Flux (irradiance) data read from uvspec output file
        :return: None
        """
        #TODO rename this function to just distribute_data ?
        # First split the data amongst output levels or output wavelengths/wavenumbers
        # We assume and handle only those cases of output_user where the primary variable is
        # zout, lambda (wvl) or wavenumber (wvn)
        if self.output_user:  # distribute to user-defined output variables
            fields = self.output_user
        else:
            fields = self.fluxline
        datashape = fluxdata.shape
        if len(datashape) == 1:  # Need some special cases to deal with single line output files
            linecount = 1
        else:
            linecount = datashape[1]
        # Deal with the shape of the data output and try to reshape, depending on the number of
        # wavelengths and output levels (zout, zout_sea or pressure)
        # if len(datashape) == 2:
        #    if datashape[1] != len(fields):  # number of fields in data does not match number of fields
                # print(datashape[1], ' not the same as ', len(fields)) #TODO try to deal with this
        if (fields[0] == 'zout' or fields[0] == 'zout_sea' or fields[0] == 'p' or
            fields[0] == 'cth'):  # output level or pressure level is the primary variable
            if self.n_levels_out == 0:  # Unknown number of output levels
                # Try just using the number of unique values in the first column
                self.n_levels_out = len(np.unique(fluxdata[:,0]))
            fluxdata = fluxdata.reshape((self.n_levels_out, -1, linecount), order='F')
            # Want wavelength to be first variable, so swap first and second axes if there are 2 or more axes
            if fluxdata.ndim >= 2:
                fluxdata = np.swapaxes(fluxdata, 0, 1)
        elif fields[0] == 'wvl' or fields[0] == 'wvn':  # wavelength/wavenumber is the primary variable
            if self.n_levels_out == 0:  # Don't know number of output levels
                self.n_wvl = len(np.unique(fluxdata[:,0]))  # Try to determine number of wavelengths/wavenumbers
                fluxdata = fluxdata.reshape((-1, self.n_wvl, linecount), order='F')
                # Want wavelength to be first variable, so swap first and second axes if there are 2 or more axes
                if fluxdata.ndim >= 2:
                    fluxdata = np.swapaxes(fluxdata, 0, 1)
            else:
                fluxdata = fluxdata.reshape((self.n_levels_out, -1, linecount), order='F')
                # Want wavelength to be first variable, so swap first and second axes if there are 2 or more axes
                if fluxdata.ndim >= 2:
                    fluxdata = np.swapaxes(fluxdata, 0, 1)
        else:  # Assume secondary variable is zout
            fluxdata = fluxdata.reshape((self.n_levels_out, -1, linecount), order='F')  #TODO provide warning or something
            # Want wavelength to be first variable, so swap first and second axes if there are 2 or more axes
            if fluxdata.ndim >= 2:
                fluxdata = np.swapaxes(fluxdata, 0, 1)

        self.fluxdata = fluxdata  # retain the flux data in the instance, reshaped as well as possible
        if linecount == 1:  # here the data is actually distributed, single line a special case
            # Some output fields, such as umu, uu, u0u, uu_down, uu_up, cmu(?) are vectors and therefore occupy
            # multiple columns, so keep track of columns and try to
            colstart = 0  # keep track of starting column
            for (ifield, field) in enumerate(fields):
                if field == 'uu' or field == 'uu_up' or field == 'uu_down':  # Do uu_up and uu_down actually exist ?
                    # Number of columns is n_umu * n_phi
                    ncols = max(self.n_umu, 1) * max(self.n_phi, 1)
                elif field == 'u0u':  #TODO does this exist and if so, what is the size ?
                    pass
                elif field == 'umu':
                    ncols = self.n_umu
                else:  # All other output variables are assumed to occupy only one column
                    ncols = 1
                # setattr(self, field, np.squeeze(fluxdata[colstart:(colstart + ncols)]))
                setattr(self, field, fluxdata[colstart:(colstart + ncols)])
                colstart = colstart + ncols
        else:  # There is more than 1 line of output flux/user data
            # Some output fields, such as umu, uu, u0u, uu_down, uu_up, cmu(?) are vectors and therefore occupy
            # multiple columns, so keep track of columns and try to distribute in a reasonable way
            colstart = 0  # keep track of starting column
            for (ifield, field) in enumerate(fields):
                if field == 'uu' or field == 'uu_up' or field == 'uu_down':  # Do uu_up and uu_down actually exist ?
                    # Number of columns is n_umu * n_phi
                    ncols = max(self.n_umu, 1) * max(self.n_phi, 1)
                elif field == 'u0u':  #TODO does this exist and if so, what is the size ?
                    pass
                elif field == 'umu':
                    ncols = self.n_umu
                else:  # All other output variables are assumed to occupy only one column
                    ncols = 1
                # setattr(self, field, np.squeeze(fluxdata[:, :, colstart:(colstart + ncols)]))
                setattr(self, field, fluxdata[:, :, colstart:(colstart + ncols)])
                colstart = colstart + ncols
        # Clean up zout and wavelength/wavenumber data
        # Convert levels to real numbers
        # TODO : What number should TOA translate to ? Maximum height defined in the atmosphere file ? Fixed number ? np.inf ?
        # TOA -> np.inf seems to cause problems for xarray
        level_dict = {'boa': self.altitude, 'BOA': self.altitude, 'sur': self.altitude,
                      'SUR': self.altitude, 'surface': self.altitude, 'SURFACE': self.altitude,
                      'toa': np.inf, 'TOA': np.inf, 'cpt': np.nan, 'CPT': np.nan}
        level_values = np.array([])
        for level in self.levels_out:
            try:
                level_values = np.append(level_values, np.float64(level))
            except ValueError:
                if level in level_dict:
                    level_values = np.append(level_values, level_dict[level])  # Try translating
                else:
                    warnings.warn('Invalid level token found in ' + self.levels_out_type + ' keyword.')
        #print('level_values ++')
        #print(level_values)
        self.level_values = level_values
        #self.zout = np.unique(self.zout)  #TODO check that this does not reorder the zout values (especially pressure)
        #if len(self.zout) > 0:  #  TODO : Problems here with pressure_out, zout_sea etc.
        if self.levels_out_type == 'zout':
            self.zout = np.unique(self.level_values)
            self.n_levels_out = len(self.zout)
            self.zout_sea = self.zout + self.altitude
        elif self.levels_out_type == 'zout_sea':
            self.zout_sea = np.unique(self.level_values)
            self.n_levels_out = len(self.zout_sea)
            self.zout = self.zout_sea - self.altitude
        elif self.levels_out_type == 'pressure_out':
            self.pressure_out = np.unique(self.level_values)[::-1]  # Reverse order
            self.n_levels_out = len(self.pressure_out)
            self.zout = np.nan
            self.zout_sea = np.nan

        #    self.n_levels_out = len(self.zout)
        self.wvn = np.unique(self.wvn)
        self.wvl = np.unique(self.wvl)
        if len(self.wvl) > 0 and len(self.wvn) == 0:  # Calculate wavenumbers if wavelengths available
            self.wvn = 1e7 / self.wvl
        if len(self.wvn) > 0 and len(self.wvl) == 0:  # Calculate wavelengths if wavenumbers available
            self.wvl = 1e7 / self.wvn
        self.n_wvl = len(self.wvl)

    def rad_units_str(self, latex=False):
        """ Provide a radiance units string e.g. W/sr/m^2/nm.
        If output_quantity is set to 'brightness' or 'reflectivity', rad_units will be 'K' or '' respectively.

        :return: Radiance units as a string
        """

        rad_units = self.rad_units[0]
        if rad_units and rad_units[-1] == 'W':
            rad_units += '/sr'
            if self.rad_units[1]:
                rad_units += '/' + self.rad_units[1]
            if self.rad_units[2]:
                rad_units += '/' + self.rad_units[2]
        if latex:
            rad_units = ' [$' + rad_units + '$]'
        return rad_units

    def irrad_units_str(self, latex=False):
        """ Provide an irradiance units string e.g. W/m^2/nm.
        If output_quantity is set to 'brightness' or 'reflectivity', irrad_units will be 'K' or '' respectively.

        :return: Irradiance units as a string
        """

        irrad_units = self.rad_units[0]
        if irrad_units and irrad_units[-1] == 'W':
            if self.rad_units[1]:
                irrad_units += '/' + self.rad_units[1]
            if self.rad_units[2]:
                irrad_units += '/' + self.rad_units[2]
        if latex:
            irrad_units = ' [$' + irrad_units + '$]'
        return irrad_units

    def readout(self, filename=None):
        """ Read uvspec output and assign to variables as intelligently as possible.

         The general process of reading is:
          1) If the user has specified output_user, just assume a flat file and read using
             np.loadtxt or np.genfromtxt.
          2) If not user_output and the solver has no radiance blocks, assume a flat file and
             read using np.loadtxt. The variables to be read should be contained in the self.fluxline attribute.
          3) Otherwise, if the output has radiance blocks, read those depending on the radiance block
             format for the specific solver. Keep reading flux and radiance blocks until the file is exhausted.

         Once the data has all been read, the data is split up between the number of output levels and number of
         wavelengths. For radiance data, the order of numpy dimensions is *umu*, *phi*, *wavelength*, *zout* and *stokes*. That
         is, if a case has multiple zenith angles, multiple azimuth angles, multiple wavelengths and multiple output
         levels, the radiance property uu will have 4 dimensions. In the case of polradtran, there will be 5
         dimensions to include the stokes parameters (if more than 1, which is the I = intensity parameter).

        Output from uvspec depends on the solver and a number of other inputs, including the directive ``output_user`.
        For the solvers ``disort``, ``sdisort``, ``spsdisort`` and presumably also ``disort2``, the irradiance (flux) outputs default
        to
        ::
            lambda edir edn eup uavgdir uavgdn uavgup

        If radiances (intensities) have been requested with the umu
        (cosine zenith angles input), each line of flux data is followed
        by a block of radiance data as follows:
        ::
         umu(0) u0u(umu(0))
         umu(1) u0u(umu(1))
         . . . .
         . . . .
         umu(n) u0u(umu(n))

        where u0u is the azimuthally averaged radiance for the requested zenith
        angles.

        If azimuth angles (phi) have also been specified, then the
        radiance block is extended as follows:
        ::
                                phi(0)        ...     phi(m)
         umu(0) u0u(umu(0)) uu(umu(0),phi(0)) ... uu(umu(0),phi(m))
         umu(1) u0u(umu(1)) uu(umu(1),phi(0)) ... uu(umu(1),phi(m))
         . . . .
         . . . .
         umu(n) u0u(umu(n)) uu(umu(n),phi(0)) ... uu(umu(n),phi(m))

        Radiance outputs are not affected by output_user options.

        For the polradtran solver, the flux block is as follows:
        ::
           lambda down_flux(1) up_flux(1) ... down_flux(iS) up_flux(iS)

        where iS is the number of Stokes parameters specified using the
        'polradtran nstokes' directive.
        If umu and phi are also specified, the radiance block is as
        follows:
        ::
                                 phi(0)      ...      phi(m)
         Stokes vector I
         umu(0) u0u(umu(0)) uu(umu(0),phi(0)) ... uu(umu(0),phi(m))
         umu(1) u0u(umu(1)) uu(umu(1),phi(0)) ... uu(umu(1),phi(m))
         . . . .
         . . . .
         umu(n) u0u(umu(n)) uu(umu(n),phi(0)) ... uu(umu(n),phi(m))
         Stokes vector Q
         . . .
         . . .

        The u0u (azimuthally averaged radiance) is always zero for
        polradtran.

        For the two-stream solver (twostr), the flux block is
        ::
           lambda edir edn eup uavg

        The directive keyword ``brightness`` can also change output. The documentation simply states that radiances and
        irradiances are just converted to brightness temperatures.

        The keyword directive ``zout`` and it's parameters will influence output format as well. In general the output
        is repeated for each given value of ``zout`` or ``zout_sea``.

        The keyword directive ``output`` and its parameters will also have a major effect.

        The ``output sum`` keyword sums output data over the wavelength dimension. This in contrast to ``output integrate``,
        which performs a spectral integral.

        The keyword directive ``header`` should not be used at all. This produces some header information in the output
        that will cause errors. A warning is issued of the ``header`` keyword is used in the input.

        :param filename: File from which to read the output. Defaults to name of input file, but with the .OUT
        extension.
        :return: None

        """
        #TODO check for use of header keyword
        if self.fluxline == '?':
            print('Unknown output format. Skipping file read.')
            return
        if filename is None:
            filename = self.outfile
        elif filename == '':
            filename = easygui.fileopenbox(msg='Please select the uvspec output file.', filetypes=["*.OUT"])
        fluxdata = []
        if not os.path.isfile(filename):
            print('Output file does not exist. Run uvspec.')
            return
        if self.output_user:
            fluxdata = np.loadtxt(filename, dtype=np.float64)
            self.distribute_flux_data(fluxdata)
        elif ((self.n_phi == 0 and self.n_umu == 0) or self.solver == 'sslidar' or
              self.solver == 'mystic'):  # There are no radiance blocks (sslidar). Mystic puts radiances in other files.
            fluxdata = np.loadtxt(filename, dtype=np.float64)
            self.distribute_flux_data(fluxdata)
        # elif self.n_phi == 0:   # Not sure about format for n_umu > 0, n_phi == 0
        #  Look at example UVSPEC_FILTER_SOLAR.INP, which indicates manual is not correct
        #     fluxdata = np.loadtext(filename)
        #     # Take away the radiance data
        #     raddata = fluxdata[len(self.fluxline):]
        #     fluxdata = fluxdata[:len(self.fluxline)]
        #     self.raddata = raddata
        #     self.fluxdata = fluxdata
        #     self.distribute_flux_data(fluxdata)
        else:  # read radiance blocks  #TODO polradtran has different format
            phicheck = []
            radND = []  # full 3D/4D radiance data is here (umu, phi and wavelength)
            with open(filename, 'rt') as uvOUT:
                # Read and append a line of flux data
                txtline = uvOUT.readline()
                while txtline:
                    fluxline = np.array(txtline.split()).astype(np.float64)
                    if len(fluxdata) > 0:
                        fluxdata = np.vstack((fluxdata, fluxline))
                    else:
                        fluxdata = fluxline
                    # Read the radiance block
                    if self.n_phi > 0:  # There is a line of phi angles
                        philine = uvOUT.readline()  #TODO check that phi angles are correct
                        phicheck = np.array(philine.split()).astype(np.float64)
                    # Read the lines for the umu values (radiances in radiance block)
                    raddata = []  # All polarisation blocks are vstacked
                    for polComp in self.stokes:  # read a radiances block for each polarisation component
                        if self.solver == 'polradtran':  # Read the polarisation block header
                            polBlockHeader = uvOUT.readline().strip()
                            # print(polBlockHeader + '}')
                            if not re.match('^Stokes vector ' + polComp + '$', polBlockHeader):
                                raise IOError('Stokes vector header not found for component ' + polComp)
                        # raddata = []
                        for i_umuline in range(self.n_umu):  #TODO this is possibly wrong if there is no phi specified - see manual
                            umuline = uvOUT.readline()
                            radline = np.array(umuline.split()).astype(np.float64)
                            if len(raddata) == 0:
                                raddata = radline
                            else:
                                raddata = np.vstack((raddata, radline))
                    if len(radND) == 0:
                        radND = raddata
                    else:
                        radND = np.dstack((radND, raddata))
                    txtline = uvOUT.readline()  # Read what should be the next line of flux data

            self.distribute_flux_data(fluxdata)  # distribute the flux data, which should also determine
                                                 # the number of wavelengths definitively
            self.radND = radND
            self.phi_check = phicheck

            # See if one size fits all
            # self.u0u = radND[:,1].reshape(self.n_umu, self.n_stokes, self.n_wvl, -1, order='F').squeeze()  # should actually all be zero
            #print(radND.shape)
            #print(radND)
            if radND.ndim == 1:
                self.u0u = radND[1]
                self.uu = radND[2:]
            else:
                self.u0u = radND[:,1].reshape(self.n_umu, self.n_stokes, -1, self.n_wvl, order='F')
                self.u0u = self.u0u.transpose([0, 3, 2, 1])
                self.uu = radND[:,2:]
            # There is actually some radiance data
            if self.uu.size:  # checks how many elements actually
                self.uu = self.uu.reshape((self.n_umu, self.n_stokes, max(self.n_phi, 1), -1, self.n_wvl), order='F')
                # transpose so that the nstokes axis is last and wavelength is 3rd (also 3rd last)
                self.uu = self.uu.transpose([0, 2, 4, 3, 1])  # transpose so that the nstokes axis is last
            else:
                self.uu = np.array([])
            # At this point, uu should be 5 dimensional (possibly with singleton dimensions) in order of
            # 'umu', 'phi', 'wvl', 'zout' (or equivalent) and stokes
            # if self.uu.shape[3] == 0:  # singleton dimension in zout
            #     print('Singleton in zout dimension, n_wvl = ')
            #     print(self.n_wvl)
            #     print('self.uu')
            #     print(self.uu)
            #     print('self.u0u'0
            #     print(self.u0u)
                # self.uu = self.uu

            # while self.uu.shape[-1] == 1:  # Remove any hanging singleton dimensions at the end
            #    self.uu = self.uu.squeeze(axis=self.uu.ndim-1)

            # if self.solver == 'polradtran':
            #     self.u0u = radND[:,1].reshape(self.n_umu, self.nstokes, self.n_wvl, -1, order='F').squeeze()  # should actually all be zero
            #     if radND.ndim == 3:  # There are multiple n_phi values
            #         self.uu = radND[:,2:,:]
            #         # Reshape to multiple nstokes values, multiple wavelengths and multiple zout levels
            #         self.uu = self.uu.reshape((self.n_umu, self.nstokes, self.n_phi, self.n_wvl, -1), order='F').squeeze()
            #
            #     else: # There are only different umu values, no different phi values
            #         self.uu = radND[:,2:]
            #         # reshape to multiple nstokes
            #         self.uu = self.uu.reshape((self.n_umu, self.nstokes, max(self.n_phi, 1), self.n_wvl, -1), order='F').squeeze()
            #
            # else:
            #     self.u0u = radND[:,1].reshape(self.n_umu, self.nstokes, self.n_wvl, -1, order='F').squeeze()
            #     #if self.u0u.shape[1] / self.n_wvl > 1:  #TODO problem here with single wavelength data
            #     #    self.u0u = self.u0u.reshape((self.n_umu, self.n_wvl, -1), order='F')
            #     if radND.ndim == 3:
            #         self.uu = radND[:,2:,:]
            #         # If there are multiple zout levels, must reshape self.uu
            #         if self.uu.shape[2] / self.n_wvl > 1:
            #             self.uu = self.uu.reshape((self.n_umu, self.n_phi, self.n_wvl, -1), order='F')
            #     else:
            #         self.uu = radND[:,2:]
            #         if self.uu.shape[1] / self.n_wvl > 1:  # if multiple output levels, reshape the radiance data appropriately
            #             self.uu = self.uu.reshape((self.n_umu, self.n_wvl, -1), order='F')
        if self.solver == 'mystic':
            # TODO : Read mystic fluxes and radiances from mc.flx.spc and mc.rad.spc
            warnings.warn('Reading of mc.flx.spc and mc.rad.spc skipped. To be implemented.')
        # Perform further processing of outputs, mainly production of xr.DataArray versions of outputs.
        self.process_outputs()


    def readerr(self, filename=None):
        """ Read any error output from the uvspec run and attach it to the self object in the self.stderr property

        It is important to check error output to see if the uvspec run was successful or what diagnostics were
        provided.

        WARNING : If the 'verbose' option is used in the uvspec input file, the error output file can be VERY large.
        The 'verbose' option prints a large amount of information relating to setup of the RT problem. Please refer to
        the libRadtran manual in respect of using 'verbose'. It is encouraged in the beginning stages of setting up
        an RT problem to verify correctness, but with (perhaps) a single wavelength to reduce the stderr output
        volume.

        :param filename: The name of te error file to read. If not provided, the default filename will be used. If
            provided as the empty string '', a file/open dialog will be presented.
        :return: True if a file was read and False if the file was not found.
        """
        if filename is None:
            filename = self.errfile
        elif filename == '':
            filename = easygui.fileopenbox(msg='Please select the uvspec error output file.', filetypes=["*.ERR"])
        if not os.path.isfile(filename):
            print('Error output file does not exist. Run uvspec with stderr redirected to an output file.')  ##TODO use an exception
            self.stderr = ''
            return False
        # Read the entire error file and attach as the stderr field
        with open(filename, 'rt') as sterrfid:
            self.stderr = sterrfid.readlines()
        return True

    def run(self, stderr_to_file=True, write_input=True, read_output=True, block=True, purge=True, check_output=False):
        """ Run the libRadtran/uvspec case.

        This will run the libRadtran/uvspec Case instance provided. Some control is provided regarding the handling of
        the standard error output from uvspec. This function only returns when uvspec terminates.

        :param stderr_to_file: Controls whether standard error output goes to the screen or is written to a file
            having the same name as the input/output files, except with the extension .ERR. Default is true.
            Any error output
        :param write_input: Controls whether the input file is written out before execution. Default is True.
        :param read_output: Controls whether the output file is read after execution. Default is True.
        :param block: By default, this method waits until uvspec terminates. If set False, the uvspec process
            is released to background and read_output is set to False (regardless of user input).
        :param purge: If set True, the input and output files from this run will be deleted after the run is complete
            and the outputs have been read. Will only be honoured if read_output is also True. Default is True.
            To prevent purging for a particular librad.Case, set the individual case property to False.
        :param check_output: If set True, the uvspec command is executed using the subprocess.check_output call, which
            will place the standard output from the run in self.check_output. This is useful for diagnostic
            purposes.
        :return: Returns self. This is important for running across networks.
        """
        # Write input file by default
        # Note that the location of the following imports is actually important, since this run code may be
        # executed on a foreign host with a different operating system to the machine calling for the
        # RT computations.
        import subprocess
        import os

        # prepare filenames
        inputname = self.name+'.INP'
        outputname = self.name+'.OUT'
        errname = self.name+'.ERR'
        swnprocess = 'uvspec'
        #print(sys.version)
        print(f'Current working directory: "{os.getcwd()}"')
        
        errstr = f'\n\n**** Unable to spawn "{swnprocess}" process on input file "{inputname}". \n**** Probably not installed system-wide on platform.\n'

        if write_input:
            self.write(filename=inputname)
        # Write the wavelength grid file if grid data is provided
        if self.wavelength_grid is not None:
            np.savetxt(self.wavelength_grid_file, self.wavelength_grid, fmt='%13.7f')
        # Spawn a sub-process using the subprocess module
        if not stderr_to_file:
            try:
                with open(inputname, 'rt') as stin, \
                     open(outputname, 'wt') as stout:
                    return_code = subprocess.call([swnprocess], stdin=stin, stdout=stout)
            except OSError:  # the uvspec command likely does not exist
                warnings.warn(errstr)
                return_code = 1
        else:  # redirect the standard error output to a file as well
            try:
                with open(inputname, 'rt') as stin, \
                     open(outputname, 'wt') as stout, \
                     open(errname, 'wt') as sterr:
                     return_code = subprocess.call([swnprocess], stdin=stin, stdout=stout, stderr=sterr)
            except OSError:  # the uvspec command likely does not exist
                warnings.warn(errstr)
                return_code = 1
        if not return_code and read_output:
            self.readout(filename=outputname)  # Read the output into the instance if the return code OK
        if stderr_to_file and read_output:
            self.readerr(filename=errname)  # Read and attach any error output
        if self.purge and purge and read_output:  # Delete the input and output files
            try:
                os.remove(inputname)
                os.remove(outputname)
                if stderr_to_file:
                    os.remove(errname)
                if self.wavelength_grid is not None:  # Remove the wavelength grid file
                    os.remove(self.wavelength_grid_file)
            except OSError:
                pass  # Just move on if file delete fails.
        self.run_return_code = return_code  # Add the return code to self
        return self

    def process_outputs(self):
        """ Process outputs from libRadtran into moglo.Scalar and xr.DataArray objects.
        Currently only radiance outputs are processed, along with a few typical flux outputs, such as `edir`.

        Note that this method probably does not cover all libRadtran/uvspec inputs and outputs and will most
        likely continue to evolve, perhaps breaking existing code.
        :return: None
        """
        # Determine output levels if possible
        # 'cpt', being cold-point tropopause is mapped to np.nan and toa maps to np.inf


        # The most important output for librad is radiances (self.uu), so those get processed
        # first. This is a 5 dimensional numpy array with axes 'umu'. 'phi', 'wvl', 'zout' (or equivalent)
        # and 'stokes'. The 'wvl' axis could also be 'chn' (spectral channel) or 'wvn' (spectral wavenumber)
        # Create each of the axes individually using xd_identity
        # TODO : Also include putting irradiance or other user output into xr.DataArray objects
        # TODO : Want to deal with azimuthally averaged radiances as well if possible (to xr.DataArray)
        # Azimuthally averaged radiances occur when phi is not specified
        # Set up sepctral axis, output level axis and stokes component axis
        if self.spectral_axis == 'wvl':
            spectral_axis = xd_identity(self.wvl, 'wvl','nm')  # The spectral axis is wavelength
        elif self.spectral_axis == 'wvn':
            spectral_axis = xd_identity(self.wvn, 'wvn', 'cm^-1')  # The spectral axis is wavenumber
        elif self.spectral_axis == 'chn':  # The spectral axis is channel number
            if np.any(self.wvl):
                spectral_axis = xd_identity(self.wvl, 'wvl','nm')
            elif np.any(self.wvn):
                spectral_axis = xd_identity(self.wvn, 'wvn', 'cm^-1')
        else:
            spectral_axis = xd_identity(np.nan, 'unknown', 'unknown')  # Don't know what is going on, so just put in something
        levels = xd_identity(self.level_values, self.levels_out_type)  # Presume units correct
        self.levels = levels
        stokes = xd_identity(range(self.n_stokes), 'stokes')
        self.stokes = stokes
        # Convert uu radiance output to xr.DataArray
        if self.uu.size:  # OK, there is some radiance data (not azimuthally averaged)
            umu = xd_identity(self.umu, 'umu', '')  # TODO : This must change to the propagation zenith (polar) angle
            paz = self.paz
            pza = self.pza
            phi = self.phi
            # Determine the units of uu
            uu_units = self.rad_units_str()
            # Build the xr.DataArray
            qty_name = {'radiance': 'specrad', 'transmittance': 'trnx', 'reflectivity': 'reflx',
                        'brightness': 'brightness'}[self.output_quantity]
            #print 'self.uu' , self.uu, ' ndim=', self.uu.ndim
            #print(spectral_axis)
            xd_uu = xr.DataArray(self.uu, [pza, paz, spectral_axis, levels, stokes],
                                    name=qty_name, attrs={'units': uu_units})
            self.xd_uu = xd_uu
        # Try to process fluxline data into xr.DataArray objects
        single_col_flux_fields = ['edir', 'edn', 'eup', 'uavgdir', 'uavgdn', 'uavgup', 'uavgglo', 'down_fluxI',
                                  'down_fluxQ', 'down_fluxU', 'down_fluxV', 'up_fluxI',
                                  'up_fluxQ', 'up_fluxU', 'up_fluxV','enet', 'esum']
        for flux_field in self.fluxline:  # This could actually be output_user data as well
            flux_units = self.irrad_units_str()  # Remember that "flux" means irradiance here
            if flux_field.startswith('uavg'):
                flux_units += '/sr'
            if flux_field in single_col_flux_fields:  # These are single column outputs
                flux_data = getattr(self, flux_field)
                if flux_data.shape[2] == 1:  # Remove trailing singleton dimension
                    setattr(self, flux_field, flux_data.squeeze(axis=2))
                else:
                    warnings.warn('Non-singleton third dimension encountered in scalar flux data.')
                try:
                    xd_flux = xr.DataArray(getattr(self, flux_field), [spectral_axis, levels], name=flux_field,
                                         attrs={'units': flux_units, 'long_name': long_name[flux_field]})
                    setattr(self, 'xd_' + flux_field, xd_flux)
                except:
                    warnings.warn('Unable to convert flux data to xarray.')
        # TODO : Would be preferable to put stokes parameters all into single xr.DataArray for polradtran
        # TODO : Processing of 'mystic' outputs to xr.DataArray objects
        # TODO : Mean intensity is the actinic flux divided by 4 pi, should perhaps be expressed /sr in units.

    def split_case_by_wavelength(self, n_sub_ranges, overlap):
        """ Split a librad.Case into a list of cases with wavelength sub-ranges
        This function can be useful to distribute a case over a compute cluster. Equal wavelength interval sub-cases
        should have similar runtimes. Use ipyparallel to run the list on a cluster. If the wavelength range
        is sufficiently large, consider setting n_sub_ranges to the number of available processors.

        This function only works with a librad.Case that uses the `wavelength` option keyword. Cases that make use of
        other methods to set the wavelength range will not work.

        :param n_sub_ranges: Number of cases in the list, which is also the number of wavelength sub-ranges.
        :param overlap: Amount of wavelength overlap between subrange cases.
        :return: List of librad.Case

        .. seealso:: merge_caselist_by_wavelength
        """
        # Obtain the wavelength range of the master case.
        if 'wavelength' in self.options:
            # Obtain the list index
            i_wvl_option = self.options.index('wavelength')
            # Obtain the wavelength range
            wvl_range_tokens = self.tokens[i_wvl_option]
            wvl_min = float(wvl_range_tokens[0])
            wvl_max = float(wvl_range_tokens[1])
            # Get the split points

            wvl_splits = np.linspace(wvl_min, wvl_max, int(n_sub_ranges) + 1)
            import copy  # need copy to make deep copies of self
            caselist = [copy.deepcopy(self) for icase in range(int(n_sub_ranges))]
            # Now run through the cases and set the wavelength sub-ranges
            for i_split, this_case in enumerate(caselist):
                this_case.set_option('wavelength', wvl_splits[i_split] - overlap * np.sign(np.float(i_split)),
                                     wvl_splits[i_split + 1])
                this_case.name += '{:04d}'.format(i_split)  # really important to change the name to avoid clashes at runtime
            return caselist
        else:
            return []  # Return an empty list if the 'wavelength' option does not occur

    @staticmethod
    def merge_caselist_by_wavelength(caselist, attr_name):
        """ Merges data in a particular attribute (property) of a list of librad.Case

        This function can be used to merge data from a list of runs created using the function
        split_case_by_wavelength.

        :param caselist: list of librad.Case after all cases have been run
        :param attr_name: name of attribute to merge e.g. 'edir'
        :return: merged wavelengths, requested attribute as numnpy arrays

        .. seealso:: split_case_by_wavelength

        """
        wvl_merged = caselist[0].wvl  # will put merged array of wavelengths in here
        data_merged = getattr(caselist[0], attr_name)  # merged data to go in here
        caselist = caselist[1:]
        for this_case in caselist:
            wvl_merged = np.concatenate((wvl_merged, this_case.wvl))
            data_merged = np.concatenate((data_merged, getattr(this_case, attr_name)))
        # remove duplicate wavelength data
        wvl_merged, merge_indices = np.unique(wvl_merged, return_index=True)
        data_merged = data_merged[merge_indices, ...]
        return np.vstack(wvl_merged), data_merged


class RadEnv(object):

    """ RadEnv is a class to encapsulate a large number of uvspec runs to cover a large number of sightlines over the
    whole sphere. A radiance map over the complete sphere is called a radiant environment map. The uvspec utility can
    only handle a limited number of sighlines per run, determined by the maximum number of polar and azimuthal angles
    specified in the file /libsrc_f/DISORT.MXD. If these values are changed, DISORT and uvspec must be recompiled. If
    the values are set too large, the memory requirements could easily exceed your computer's limit (there is
    currently no dynamic memory allocation in DISORT). The situation for the cdisort solver is less clear.

    .. note::
        The `polradtran` radiative transfer (RT) solver does not produce direct irradiance outputs (`edir`),
        This presents a problem for calculation of path transmission (optical depth) between output atmospheric
        levels (e.g. as specified by the uvspec `zout` keyword). This solver only produces total fluxes (irradiances)
        for each of desired stokes parameters. Hence for calculation of polarised radiant environment maps (REMs), the only
        feasible option for `librad` is to use the `mystic` solver with the `mc_polarisation` option. Currently,
        the librad.Case uvspec output file reading functions do not cater for `mystic`.

    """

    def __init__(self, base_case, n_pol, n_azi, mxumu=48, mxphi=19, hemi=False, n_sza=0):
        """ Create a set of uvspec runs covering the whole sphere to calculate a full radiant environment map.
        Where the base_case is the uvspec case on which to base the environmental map, Name is the name to give the
        environmental map and n_pol and n_azi are the number of polar and azimuthal sightline angles to generate. The
        mxumu and mxphi are the maximum number of polar and azimuth angles to calculate in a single run of uvspec.
        The default values are mxumu = 48, and mxphi = 19. These values are taken from the standard libRadtran
        distribution (/libsrc_f/DISORT.MXD) maximum parameter file. If using the polradtran solver, the corresponding
        file is /libsrc_f/POLRADTRAN.MXD. Other solvers may have different restrictions. A warning will be issued if
        the solver is not in the DISORT/POLRADTRAN family.

        :param base_case: librad.Case object providing the case on which the environment map is to be based. Note
            that not any basecase can be used. As a general guideline, the basecase should have standard irradiance
            outputs (i.e. should not use the `output_user` keyword). It should also not the use `output_process` or
            `output_quantity` keywords, which change the units and/or format of the libRadtran/uvspec output.
             Minimal validation of the basecase is performed. However, use with `mol_abs_param` such as `kato` and
             `fu` is important for `librad` and these are supported (k-distribution or `correlated-k`
             parametrizations). Use of `output_process per_nm` is appropriate for `source thermal` REMs to get
             radiance units per nanometre rather than per inverse cm.
        :param n_pol: Number of polar angles (view/propagation zenith angles)
        :param n_azi: Number of azimuthal angles.
        :param mxumu: Maximum number of polar angles per case.
        :param mxphi: Maximum number of azimuthal angles per case.
        :param hemi: If set True, will generate only a single hemisphere being on one side of
            the solar principle plane. Default is False i.e. the environment map covers the full sphere.
            Note that if hemi=True, the number of REM samples in azimuth becomes n_azi :math:`\\times` 2.
            This is the recommended mode (hemi=True) for librad purposes, since it reduces execution time.
        :param n_sza: The number of solar zenith angles (SZA) at which to perform transmittance and path radiance
            computations. Each SZA will result in another run of the base case (no radiances)


        The solver cdisort may have dynamic memory allocation, so the warning is still issued because the situation
        is less clear.

        A note about radiative propagation angles and viewing angles, which are 180 deg opposite to each other.
        For radiance calculations define the cosine of the viewing zenith angle
        `umu` and the sensor azimuth `phi` and don't forget to also specify the solar azimuth
        `phi0`. `umu` > 0 means sensor looking downward (e.g. a satellite), `umu` < 0 means looking
        upward. `phi` = `phi0` indicates that the sensor looks into the direction of the sun,
        `phi` - `phi0` = 180 means that the sun is in the back of the sensor.

        `phi` is the propagation azimuth angle `paz`, except that `paz` is in radians and `phi` is in degrees.

        `pza` is the propagation zenith angle in radians.

        `vaz` is the view azimuth angle and is 180 :math:`^\\circ` different from `paz`. `vaz` is expressed
        in degrees. `vza`, the view zenith angle is 180 :math:`^\\circ` different from `paz' and is expressed
        in degrees. In order to keep all azimuth angles in increasing order, 'vaz' is in the range [-180, 180],
        while `phi` is in the range [0, 360] and `vaz` = `phi` - 180.

        The value of `umu` is the cosine of the propagation zenith angle. The value of `phi` is the true azimuth of
        the propagation direction.

        For all one-dimensional solvers the absolute azimuth does not matter, but only the relative azimuth
        `phi` - `phi0`.

        For `librad` work, it is strongly recommended that the `hemi` flag be set True. This will automatically
        set the `phi0` keyword to zero in the RadEnv cases when running uvspec. This will essentially halve the
        execution time for radiant environment maps of the same effective spatial resolution.
        """
        # Require that n_pol be even
        # This is to ensure that there is no umu = 0 (horizontal direction, illegal in libRAdtran)
        if n_pol % 2:  # check comment below about Matlab polar angles
            warnings.warn('Input n_pol to librad.RadEnv must be even. Increased by 1.')
            n_pol += 1
        if n_azi % 2:
            warnings.warn('Input n_azi to librad.RadEnv must be even. Increased by 1.')
            n_azi += 1
        self.base_case = copy.deepcopy(base_case)  # Keep a copy of the uvspec base_case
        self.levels_out_type = self.base_case.levels_out_type
        self.n_levels_out = self.base_case.n_levels_out
        self.solver = self.base_case.solver  # Radiative transfer solver
        self.trans_base_case = copy.deepcopy(base_case)  # Keep a copy for transmittance computation purposes
        view_zen_angles = np.linspace(0.0, 180.0, int(n_pol))  # Viewing straight down is view zenith angle of 180 deg
        prop_zen_angles = np.linspace(np.pi, 0.0, int(n_pol))  # Radiation travelling straight up is propagation zenith angle of 0
        umu = np.cos(prop_zen_angles)  # Negative umu is upward-looking, downwards propagating
        if hemi:
            prop_azi_angles = np.linspace(0.0, 180.0, int(n_azi))
            view_azi_angles = np.linspace(-180.0, 0.0, int(n_azi))
        else:
            prop_azi_angles = np.linspace(0.0, 360.0, int(n_azi))
            view_azi_angles = np.linspace(-180.0, 180.0,  int(n_azi))
        phi = prop_azi_angles
        n_azi_batch = np.int(np.ceil(np.float64(len(prop_azi_angles))/mxphi))
        n_pol_batch = np.int(np.ceil(np.float64(len(prop_zen_angles))/mxumu))
        # Create an list of lists with all these batches of librad.Case
        self.cases = [[copy.deepcopy(base_case) for i_azi in range(n_azi_batch)] for j_pol in range(n_pol_batch)]

        # TODO : Take care of phi0 input in the case of hemi=True
        for iazi, iazi_start in enumerate(range(0, len(phi), mxphi)):
            batch_azi = phi[iazi_start:np.minimum(iazi_start+mxphi, len(phi))]
            for ipol, ipol_start in enumerate(range(0, len(umu), mxumu)):
                batch_pol = umu[ipol_start:np.minimum(ipol_start+mxumu, len(umu))]
                # Set the umu and phi keyword parameters
                self.cases[ipol][iazi].alter_option(['phi'] + [str(x) for x in batch_azi])
                self.cases[ipol][iazi].alter_option(['umu'] + [str(x) for x in batch_pol])
                self.cases[ipol][iazi].infile = (self.cases[ipol][iazi].infile[:-4] +
                                                 '_{:04d}_{:04d}.INP'.format(ipol, iazi))
                self.cases[ipol][iazi].outfile = (self.cases[ipol][iazi].outfile[:-4] +
                                                 '_{:04d}_{:04d}.OUT'.format(ipol, iazi))
                self.cases[ipol][iazi].name = (self.cases[ipol][iazi].name +
                                                 '_{:04d}_{:04d}'.format(ipol, iazi))
                if hemi:  # Doing only one hemisphere along solar principal plane
                    self.cases[ipol][iazi].alter_option(['phi0', '0.0'])  # Sun shining towards North
        # Create a flattened list view of the cases
        # Put all the cases into a single list
        self.casechain = list(chain(*self.cases))  # This creates a linear view of the cases
        self.hemi = hemi
        self.n_azi = n_azi
        self.n_pol = n_pol
        self.n_azi_batch = n_azi_batch
        self.n_pol_batch = n_pol_batch
        self.phi = xd_identity(prop_azi_angles, 'phi', 'deg')
        self.umu = xd_identity(umu, 'umu', '')
        self.pza = xd_identity(prop_zen_angles, 'pza', 'rad')
        self.paz = xd_identity(np.deg2rad(prop_azi_angles), 'paz', 'rad')
        self.vza = xd_identity(view_zen_angles, 'vza', 'deg')
        self.vaz = xd_identity(view_azi_angles, 'vaz', 'deg')
        self.has_water_clouds = self.base_case.has_water_clouds
        self.has_ice_clouds = self.base_case.has_ice_clouds
        self.has_clouds = self.base_case.has_clouds
        self.n_sza = n_sza
        if n_sza:  # Setup the transmission run cases
            self.setup_trans_cases(n_sza=n_sza)
            # self.trans_vza_up = []
        self.uu = np.array([])

    def setup_trans_cases(self, n_sza=32):
        """ Setup a list of cases for computing the transmission matrices between every level defined in the
        REM (and at every wavelength and stokes parameter combination). The computation of  transmittance between
        levels is accomplished in `librad` using libRadtran/uvspec by computing the direct solar irradiance
        transmittance for multiple zenith angles.

        Note that if there is an optically thick cloud layer between two levels in the REM, the transmittance
        will compute as zero or very small. This will also result in the incorrect/noisy computation of path
        radiances between the two levels. This situation is unavoidable. The recommended approach is that
        REMS be computed with altitude levels that all lie between cloud layers (i.e. no levels in the REM
        span a cloud layer). The user, or the code which uses the REMs should see to this. Essentially it must
        be recognised that optical surveillance is not possible through optically thick cloud layers.

        The REM is provided with a cloud flag that indicates if the base case incorporates clouds. The approach to
        computing inter-level transmittance (optical depth is the stored parameter, since this scales more
        closely in a linear fashion with distance) in the presence of cloud is to vary the 'cloudcover'
        keyword parameter (water and ice clouds independently). This is done for solar zenith angle of zero.
        The optical thickness from TOA to the level in question is then computed as a function of cloud cover
        fraction (CCF). The optical depth between levels (i.e the optical depth of a layer between two levels)
        is computed as the difference in optical depth to TOA of the lower level minus the optical depth to TOA
        of the upper level. If there is no difference in the layer optical depth when the CCF is varied from
        zero to some positive value (say 0.1, but not as high as 1), then the layer is free of cloud.

        A flag per layer is thus obtained which indicates if the layer contains cloud. If it does, the transmittance
        will compute as zero between the two levels in question. This means that the cloud base is not resolved
        to better than the level resolution in the REM. It is probably then quite important to ensure that
        cloud base altitude statistics are available in the theatre climatology.

        An upgrade to cloud handling could be to read the cloud profile files in order to obtain the exact vertical
        location of the cloud layers.

        Another inherent and unavoidable problem with computation of path optical depths and radiances using
        libRadtran/uvspec is that precisely horizontal paths cannot be dealt with using one-dimensional RT
        solvers. Therefore in this case, the maximum range that can be dealt with depends on the height difference
        between the REM levels and the maximum solar zenith angle (SZA) used for computation of optical depth.

        A further implication of the above point is that path transmittances and path radiances cannot be interpolated
        between the SZA nearest the horizon and the horizon proper. Some form of logarithmic extrapolation could be
        performed, but could result in large errors due to failure of Beer's Law and other problems.

        Execution and attribution of transmittance cases will provide each level with transmittance to the level
        above that altitude level. Therefore the topmost level will have transmittances to TOA, but the
        bottom level will not have transmittances (optical depths) to BOA, unless the bottom level *is* BOA.
        It is recommended that all librad base cases for REM include BOA as an output level.

        :param n_sza: The number of solar zenith angles at which to compute path optical depth and radiances.
             the SZA values are computed equi-spaced in the cosine of the solar zenith angle rather that the
             SZA itself. This is to help with the problem that the slant range between levels increases
             in linear relation to the secant of the view zenith angle. The optical depths are later
             interpolated to the same
        :return:
        """

        # Start by removing any cloud options in the transmission base_case, as well as any radiance options
        # to reduce runtime.
        new_option_list = []
        new_tokens_list = []
        # TODO : setup up transmission case series with AND without clouds
        # Remove radiance options to speed up transmission series computations
        options_to_remove = ['umu', 'phi', 'phi0']
        for (ioption, option) in enumerate(self.trans_base_case.options):
            if any([option == removal_option for removal_option in options_to_remove]):
                pass  # TODO : Look for other cloud options that may be important
            else:
                new_option_list.append(option)
                new_tokens_list.append(self.trans_base_case.tokens[ioption])
        self.trans_base_case.options = new_option_list
        self.trans_base_case.tokens = new_tokens_list
        # There are no radiance calculations now, so set n_phi and n_umu to zero
        self.trans_base_case.n_phi = 0
        self.trans_base_case.phi = []
        self.trans_base_case.n_umu = 0
        self.trans_base_case.umu = []
        self.trans_base_case.pza = []
        self.trans_base_case.paz = []
        # Change the transmission base case to output_quantity reflectivity. THis seems counter-intuitive,
        # but take a look at the libRadtran/uvspec manual to see why. Essentially, the reflectivity
        # option computes the ratio of the horizontal irradiance to the horizontal irradiance at TOA.
        self.trans_base_case.alter_option(['output_quantity', 'reflectivity'])
        self.trans_base_case.alter_option(['sza', '0.0'])
        # Calculate the view zenith angles equi-spaced in the secant of the VZA
        # First upward looking
        cosine_up = np.linspace(1.0, 0.0, int(n_sza) + 1)[:-1]  # Drop the last element because it is horizontal (illegal)
        # Calculate the view-zenith angles in radians
        vpa_up = np.arccos(cosine_up)  # view polar angle looking upwards
        # If less than the minimum depression angle for the REM, add a point
        rem_pza_limit = np.pi/2.0 - 0.99 * (np.abs(self.pza[1].data - self.pza[0]))/2.0
        if vpa_up[-1] < rem_pza_limit:
            vpa_up = np.hstack((vpa_up, rem_pza_limit))
        self.trans_vza_up = xd_identity(vpa_up, 'vza', 'rad')
        # Now build a list of uvspec runs, based on trans_base_case
        # The list is called trans_cases
        self.trans_cases = []
        for i_case in range(len(vpa_up)):
            self.trans_cases.append(copy.deepcopy(self.trans_base_case))
            # Set the solar zenith angle
            sza = self.trans_vza_up[i_case].data  # Solar zenith angle
            self.trans_cases[i_case].alter_option(['sza', str(sza)])
            # Set cloudcover to zero if the case has water or ice clouds
            if self.trans_base_case.has_water_clouds:
                self.trans_cases[i_case].alter_option(['cloudcover', 'wc', '0.0'])
            if self.trans_base_case.has_ice_clouds:
                self.trans_cases[i_case].alter_option(['cloudcover', 'ic', '0.0'])
            # Set pseudospherical option above sza of 75 degrees and solver is disort or twostr
            if sza > 75.0 and any([self.solver == thesolver for thesolver in ['disort', 'disort2', 'sdisort',
                                                           'spsdisort', 'fdisort1', 'fdisort2', 'twostr']]):
                self.trans_cases[i_case].alter_option(['pseudospherical', ''])
            # Change the name and input and output filenames, the _x_ is for transmission runs
            self.trans_cases[i_case].infile = (self.trans_cases[i_case].infile[:-4] +
                                             '_x_{:04d}.INP'.format(i_case))
            self.trans_cases[i_case].outfile = (self.trans_cases[i_case].outfile[:-4] +
                                             '_x_{:04d}.OUT'.format(i_case))
            self.trans_cases[i_case].name = (self.trans_cases[i_case].name +
                                             '_x_{:04d}'.format(i_case))
        # The transmission cases should be ready to run a this point.
        # The run_ipyparallel method will run these cases, but not in parallel with
        # the radiance cases.
        # Next, create the cloudcover series to detect cloud layers
        self.cloud_detect_cases = []
        # If there are clouds in this case, run three cases with cloud cover fraction of 0.0, 0.5 and 1.0
        if self.has_clouds:
            self.cloud_detect_cases.append(copy.deepcopy(self.trans_base_case))
            self.cloud_detect_cases[0].alter_option(['cloudcover', '0.0'])
            self.cloud_detect_cases.append(copy.deepcopy(self.trans_base_case))
            self.cloud_detect_cases[1].alter_option(['cloudcover', '0.5'])
            self.cloud_detect_cases.append(copy.deepcopy(self.trans_base_case))
            self.cloud_detect_cases[2].alter_option(['cloudcover', '1.0'])

    def run_ipyparallel(self, ipyparallel_view, stderr_to_file=False, purge=False):
        """ Run a complete set of radiant environment map cases of libRadtran/uvspec using the `ipyparallel`
        Python package, which provides parallel computation from Jupyter notebooks and other Python launch
        modes.
        Typical code for setting up the view:
            .. code-block:: python

               from ipyparallel import Client
               paraclient = Client(profile='mycluster', sshserver='me@mycluster.info', password='mypassword')
               paraclient[:].use_dill()  # Need dill as a pickle replacement for our purposes here
               ipyparallel_view = paraclient.load_balanced_view()
               ipyparallel_view.block = True  # Must wait for completion of all tasks on the cluster

        Note that if new ipengines are started, use_dill() must be executed again. The use_dill() call
        should be a routine before every function map to the cluster.


        :param ipyparallel_view: an ipyparallel view of a Python engine cluster (see ipyparallel documentation.)

        :param stderr_to_file: If set to True, standard error output will be sent to a file. use only for debugging
            purposes.
        :param purge: Boolean. If set True, the actual libRadtran cases that are executed to make up the REM are
            deleted in order to reduce the size of the object. If the object is purged, it is not possible to rerun
            the REM. Default is False - no purging (or minimal purging) is performed.
        :return:
        """
        # The following does work, but the list casechain is completely reassigned
        # instead of being assigned element for element
        # TODO : Consider passing in the client instead, to set blocking and use_dill() EVERY time.
        self.casechain = ipyparallel_view.map(Case.run, self.casechain)
        # Now recreate the list of lists view
        self.cases = [[self.casechain[i_pol * self.n_azi_batch + i_azi] for i_azi in range(self.n_azi_batch)]
                                                                        for i_pol in range(self.n_pol_batch)]
        # Compile the radiance results into one large array
        self.uu = np.hstack([np.vstack([self.cases[i_pol][i_azi].uu for i_pol in range(self.n_pol_batch)])
                                                                    for i_azi in range(self.n_azi_batch)])
        # if self.hemi:  # Double up in the azimuth direction, but flip as well
        #     self.uu = np.hstack([self.uu, self.uu[:,::-1,...]])
        #     # Perform doubling up in all azimuth variables
        #     self.pza = xr.concat((self.pza, self.pza + np.pi), dim='pza')
        #     self.phi = xr.concat((self.phi, self.phi + 180), dim='phi')
        #     self.vaz = xr.concat((self.vaz, self.vaz + 180), dim='vaz')  # view azimuth angle
        # Delete the individual results in an attempt to save memory
        for case in self.casechain:
            del case.uu
        # Concatenate the cases in umu and phi
        self.xd_uu = xr.concat([xr.concat([case_uu.xd_uu for case_uu in self.cases[jj]], dim='paz')
                                                 for jj in range(len(self.cases))], dim='pza')
        #self.xd_uu = xr.DataArray(self.uu, [self.pza, self.paz, self.spectral, self.levels, self.stokes])
        # Replace the values in the xr.DataArray with the exact original values in the zenith and azimuth
        # directions. Not doing this gave rise to a very subtle bug in spherical harmonic fitting
        self.xd_uu['pza'] = self.pza
        self.xd_uu['paz'] = self.paz
        # Also need to obtain the irradiances from one of the cases - they should actually all be the same
        fluxdata = self.casechain[0].fluxdata  # Would really want this as a xr.DataArray
        fluxline = self.casechain[0].fluxline
        irrad_units = self.casechain[0].irrad_units_str()
        # Promote irradiance data from
        for flux_component in fluxline:
            if hasattr(self.casechain[0], 'xd_' + flux_component):
                setattr(self, 'xd_' + flux_component, getattr(self.casechain[0], 'xd_' + flux_component))
        self.fluxdata = fluxdata
        self.fluxline = fluxline
        self.irrad_units = irrad_units
        # Run the transmittance sequences if there are any
        if self.n_sza:
            self.trans_cases = ipyparallel_view.map(Case.run, self.trans_cases)
            # If there are clouds in the radiant environment, run the could OD detection cases
            # These cases reveal if there are layers in the REM that include clouds
            if self.has_clouds:
                self.cloud_detect_cases = ipyparallel_view.map(Case.run, self.cloud_detect_cases)
            # Compile the transmittance data
            self.compute_path_transmittance()
            # Compile the path radiance data
            self.compute_path_radiance()
        if purge:
            del self.casechain
            del self.cases

    def run_parallel(self, n_nodes=4):
        """ Run the RadEnv in multiprocessing mode on the local host.

        This is not yet tested, but should work with the multiprocessing package on the local host to use all
        the available cores. Will only work if libRadtran is installed on the local host.

        In order to use dill instead of pickle, it is necessary to use the pathos multiprocessing module
        instead of the standard multiprocessing module

        :param n_nodes: Number of compute nodes to use. Default is 4. Preferably set to number of cores you have
            available on the local host.
        :return:
        """
        from pathos.multiprocessing import ProcessingPool
        worker_pool = ProcessingPool(nodes=n_nodes)
        self.casechain = worker_pool.map(Case.run, self.casechain)
        # Now recreate the list of lists view
        self.cases = [[self.casechain[i_pol * self.n_azi_batch + i_azi] for i_azi in range(self.n_azi_batch)]
                                                                        for i_pol in range(self.n_pol_batch)]
        # Compile the radiance results into one large array
        self.uu = np.hstack([np.vstack([self.cases[i_pol][i_azi].uu for i_pol in range(self.n_pol_batch)])
                                                                    for i_azi in range(self.n_azi_batch)])
        # Delete the individual results in an attempt to save memory
        for case in self.casechain:
            del case.uu

    def sph_harm_fit(self, degree, method='trapz'):
        """ Fit spherical harmonics to the radiant environment map (REM).
        One set of coefficients per wavelength or spectral channel will be fitted. The coefficients for each spectral
        bin/channel comprise one coefficient for order :math:`-m` to :math:`+m` for :math:`m = 0, 1, 2, ..., n`.
        The total number of coefficients is :math:`(n+1)^2`.

        The convention used for the spherical harmonics is that of Sloan with the Ramamoorthi and Hanrahan
        normalization. This is a real-valued basis defined as follows:

        .. math::
            y_{n}^{m}=\\begin{cases}
            (-1)^{m}\\sqrt{2}\\Re(Y_{n}^{m}) & m>0\\\\
            (-1)^{m}\\sqrt{2}\\Im(Y_{n}^{m}) & m<0\\\\
            Y_{n}^{0} & m=0
            \\end{cases}=\\begin{cases}
            (-1)^{m}\\sqrt{2}\\cos m\\phi\\,P_{n}^{m}(\\cos\\theta) & m>0\\\\
            (-1)^{m}\\sqrt{2}\\sin|m|\\phi\\,P_{n}^{|m|}(\\cos\\theta) & m<0\\\\
            K_{n}^{0}P_{n}^{0}(\\cos\\theta) & m=0
            \\end{cases},

        where the complex basis functions :math:`Y^{m}_{n}` are defined as:

        .. math::
            Y_{n}^{m}(\\theta,\\phi)=K_{n}^{m}e^{im\\phi}P_{n}^{|m|}(\\cos\\theta),\\,n\\in\\mathbf{N},\\,-n\\leq m\\leq n,

        having the normalisation factor of:

        .. math::
            K_{n}^{m}=\\sqrt{\\frac{(2n+1)(n-|m|)!}{4\\pi(n+|m|)!}}.

        The definition of the complex basis functions is consistent with the
        `scipy.special <http://docs.scipy.org/doc/scipy-0.16.0/reference/generated/scipy.special.sph_harm.html>`_
        definition of the spherical harmonics. Therefore, for fitting of the Sloan/Ramamoorthi/Hanrahan basis, the
        first definition is used above, that is :math:`y_{n}^{m}` can be calculated from the `scipy.special` function
        :math:`Y_{n}^{m}` as:

        .. math::
            y_{n}^{m}=\\begin{cases}
            (-1)^{m}\\sqrt{2}\\Re(Y_{n}^{m}) & m>0\\\\
            (-1)^{m}\\sqrt{2}\\Im(Y_{n}^{m}) & m<0\\\\
            Y_{n}^{0} & m=0
            \\end{cases}.

        The fitted coefficients of the spherical harmonics are computed by multiplying the REM by each of the
        harmonics and performing double numerical integration over zenith and azimuth angle as:

        .. math::
            L_{n}^{m}=\\int_{\\theta=0}^{\\pi}\\int_{\\phi=0}^{2\\pi}L(\\theta,\\phi)y_{n}^{m}(\\theta,\\phi)\\sin\\theta d\\theta d\\phi,

        where :math:`L_{n}^{m}` are the coefficients and :math:`L(\\theta,\\phi)` is the REM.

        :param degree: Spherical harmonics up to this degree :math:`n`, for all orders :math:`m` will be fitted.
        :param method: Integration method by which the coefficients are computed. 'trapz' for trapezoidal integration,
            'sum' for simple summation and 'simpson' for Simpson's Rule. The 'trapz' method seems to be
            considerably more accurate than 'sum' or 'simpson'. Therefore 'trapz' is the default.
        :return:
        """
        from scipy.special import sph_harm
        azi_angles = self.paz.data  # Propagation zenith angles in radians
        pol_angles = self.pza.data   # Propagation polar (zenith) angles in radians
        azi_ang_delta = azi_angles[1] - azi_angles[0]
        pol_ang_delta = abs(pol_angles[1] - pol_angles[0])
        # print( azi_ang_delta, pol_ang_delta)
        # Angles to be meshgridded - is this really necessary
        # TODO : Check if meshgridding is really necessary
        azi_angles, pol_angles = np.meshgrid(azi_angles, pol_angles)
        sin_pol_angles = xr.DataArray(np.sin(pol_angles), [self.pza, self.paz])
        # TODO : Symmetry checking still required
        # if self.hemi:  # Sun-symmetric REM
        sph_harm_coeff = []
        # If calculating over one hemisphere of a symmetrical REM with sun in the north-south line,
        # The sine (imaginary) coefficients are zero and the cos (real components) are doubled
        # SO should just run m = 0 to n, take twice the real part if hemi, otherwise full complex calc
        for n in range(degree + 1):
            sph_harm_coeff.append([])  # Add another list of coefficients for order m = 0 to n
            for m in range(0, n + 1):
                sph_harm_coeff[n].append(None)  # Add another coefficient for order m
                condon_short = (-1.0)**m * np.sqrt(2)
                # Compute the complex spherical harmonic basis
                if m == 0:
                    y_mn = sph_harm(0, n, azi_angles, pol_angles)
                else:
                    y_mn = condon_short * sph_harm(m, n, azi_angles, pol_angles)  # Ramamoorthi normalisation
                if self.hemi:
                    y_mn = y_mn.real * 2.0
                # Create xr.DataArray
                y_mn = xr.DataArray(y_mn, [self.pza, self.paz])
                # Compute the integrand
                inner_integrand = y_mn * self.xd_uu * sin_pol_angles
                # Compute the coefficients by integration by various methods
                if method == 'sum':  # Integrate by simple summation
                    # First integrate over the propagation zenith angle using summation
                    outer_integrand = inner_integrand.sum(dim='pza') * pol_ang_delta
                    # Then integrate over propagation azimuth
                    sph_harm_coeff[n][m] = outer_integrand.sum(dim='paz') * azi_ang_delta
                elif method == 'trapz':  # Integrate using the trapezoidal rule (seems best)
                    outer_integrand = (inner_integrand.isel(pza=0) +
                                            2.0 * inner_integrand.isel(pza=slice(1, -1)).sum(dim='pza') +
                                       inner_integrand.isel(pza=-1)) * pol_ang_delta / 2.0
                    # Then integrate over the propagation azimuth angle
                    sph_harm_coeff[n][m] = (outer_integrand.isel(paz=0) +
                                              2.0 * outer_integrand.isel(paz=slice(1, -1)).sum(dim='paz') +
                                            outer_integrand.isel(paz=-1)) * azi_ang_delta / 2.0
                elif method == 'simpson':  # The following is integration by Simpsons rule
                    outer_integrand = (inner_integrand.isel(pza=0) +
                                            4.0 * inner_integrand.isel(pza=slice(1, -1, 2)).sum(dim='pza') +
                                            2.0 * inner_integrand.isel(pza=slice(2, -2, 2)).sum(dim='pza') +
                                       inner_integrand.isel(pza=-1)) * pol_ang_delta / 3.0
                    # Then integrate over the propagation azimuth angle
                    sph_harm_coeff[n][m] = (outer_integrand.isel(paz=0) +
                                              4.0 * outer_integrand.isel(paz=slice(1, -1, 2)).sum(dim='paz') +
                                              2.0 * outer_integrand.isel(paz=slice(2, -2, 2)).sum(dim='paz') +
                                            outer_integrand.isel(paz=-1)) * azi_ang_delta / 3.0
                else:
                    warnings.warn('Unknown integration method ' + method + ' encountered in sph_harm_fit.')
        self.sph_harm_coeff = sph_harm_coeff
        return sph_harm_coeff

    def sph_harm_fat(self, degree):
        """ This code was used for debugging purposes - ignore
        :param degree:
        :return:
        """
        from scipy.special import sph_harm
        azi_angles = self.paz.data  # Propagation zenith angles in radians
        pol_angles = self.pza.data   # Propagation polar (zenith) angles in radians
        azi_ang_delta = azi_angles[1] - azi_angles[0]
        pol_ang_delta = np.abs(pol_angles[1] - pol_angles[0])
        # print(azi_ang_delta, pol_ang_delta)
        # Angles to be meshgridded - is this really necessary
        # TODO : Check if meshgridding is really necessary
        azi_angles, pol_angles = np.meshgrid(azi_angles, pol_angles)
        sin_pol_angles = xr.DataArray(np.sin(pol_angles), [self.pza, self.paz])
        # TODO : Symmetry checking still required
        # if self.hemi:  # Sun-symmetric REM
        sph_harm_coeff_cos = []
        sph_harm_coeff_sin = []
        # If calculating over one hemisphere of a symmetrical REM with sun in the north-south line,
        # The sine (imaginary) coefficients are zero and the cos (real components) are doubled
        # SO should just run m = 0 to n, take twice the real part if hemi, otherwise full complex calc
        for n in range(degree + 1):
            sph_harm_coeff_cos.append([])
            sph_harm_coeff_sin.append([])# Add another list of coefficients for order m = 0 to n
            for m in range(0, n + 1):
                sph_harm_coeff_cos[n].append(None)# Add another coefficient for order m
                sph_harm_coeff_sin[n].append(None)
                # Compute the complex spherical harmonic basis
                if m == 0:
                    y_mn = sph_harm(0, n, azi_angles, pol_angles)  # Small imaginary components may arise - drop
                else:
                    y_mn = (-1.0)**m * np.sqrt(2.0) * sph_harm(m, n, azi_angles, pol_angles)  # Ramamoorthi normalisation
                if self.hemi:  # integrate over hemisphere, therefore need a factor of 2
                    y_mn = y_mn * 2.0  # equivalent of taking only the cosine components times 2
                # Create xr.DataArray
                y_mn_cos = xr.DataArray(y_mn.real, [self.pza, self.paz])
                y_mn_sin = xr.DataArray(y_mn.imag, [self.pza, self.paz])
                # Compute the integrand
                inner_integrand_cos = self.xd_uu * y_mn_cos * sin_pol_angles
                # Compute the coefficients
                # First integrate over the propagation zenith angle using summation (very slight overestimate)
                outer_integrand_cos = inner_integrand_cos.sum(dim='pza') * pol_ang_delta
                # Then integrate over the propagation azimuth angle
                sph_harm_coeff_cos[n][m] = outer_integrand_cos.sum(dim='paz') * azi_ang_delta

                # Compute the sin integrand
                inner_integrand_sin = self.xd_uu * y_mn_sin * sin_pol_angles
                # Compute the coefficients
                # First integrate over the propagation zenith angle using summation (very slight overestimate)
                outer_integrand_sin = inner_integrand_sin.sum(dim='pza') * pol_ang_delta
                # Then integrate over the propagation azimuth angle
                sph_harm_coeff_sin[n][m] = outer_integrand_sin.sum(dim='paz') * azi_ang_delta
        self.sph_harm_coeff_cos = sph_harm_coeff_cos
        self.sph_harm_coeff_cos = sph_harm_coeff_cos
        return sph_harm_coeff_cos, sph_harm_coeff_sin

    def compute_path_transmittance(self):
        """ Compute path transmittances from the set of libRadtran/uvspec runs executed for solar zenith angles
        of 0 to near 90 degrees.

        This method computes optical depth (-log(transmittance)) of all paths from a particular level, both
        upward and downward. Paths that lie in the horizontal "blind zone" are assigned OD of 0.0. These should
        actually be assigned OD of np.nan or perhaps np.inf.

        .. seealso::
            RadEnv.setup_trans_cases()

        :return: None
        """
        # Run through all the cases in the self.trans_cases and compile the direct solar irradiance data
        # This is actually transmittance data (edir is not a flux/irradiance with output_quantity reflectivity)
        for trans_case in self.trans_cases:
            # Add a propagation zenith angle for each sza. The pza is pi - sza (in radians)
            trans_case.xd_edir = trans_case.xd_edir.assign_coords(pza=np.deg2rad(trans_case.sza))
        # Concatenate results from all transmission runs
        self.xd_edir_trans = xr.concat([this_case.xd_edir for this_case in self.trans_cases], dim='pza')
        # Interpolate transmission results onto the pza grid for the RadEnv
        # This is not a "harmonisation" interpolation. The transmission grid is being interpolated
        # onto another grid in pza (propagation zenith angle)
        # TODO : Check out reliability of using quadratic/cubic interpolation for transmittance
        self.xd_trans_toa = xd_interp_axis_to(self.xd_edir_trans, self.xd_uu, axis='pza', interp_method='linear',
                                              fill_value=1.0, assume_sorted=False)
        self.xd_opt_depth = -np.log(self.xd_trans_toa)  # Compute the optical depths from a level to TOA
        # Subtract the optical depth of the level above it.
        # This provides the optical depth from any level to the level above it
        xd_opt_depth = -self.xd_opt_depth.diff(self.levels_out_type, label='lower')
        #  TODO : Set long_name and units of optical depth
        # xd_opt_depth.attrs['long_name'] = long_name['od']
        # Implicitly, the optical depth from the top level is known to TOA so concat those values
        levels_axis_num = xd_opt_depth.get_axis_num(self.levels_out_type)
        levels_out_type = self.levels_out_type  # just to shorten it
        self.xd_opt_depth = xr.concat((xd_opt_depth, self.xd_opt_depth[{levels_out_type: -1}]),
                                        dim=levels_out_type)
        # Now confront the problem of computing optical depths to the level below
        # For the lowest level, the optical depths to the level below are undefined, perhaps an appropriate
        # value is np.nan. Otherwise, the optical depth looking to the next lower level is found by
        # flipping the OD data from the lower level along the pza axis and summing to upper level
        for ilevel in range(self.n_levels_out-1, 0, -1):  # Start at the top and work down
            # Create a copy of the optical depths from the level beneath ilevel
            opt_depth_from_beneath = copy.deepcopy(self.xd_opt_depth[{levels_out_type: ilevel-1}])
            # Flip optical depth along the pza axis
            opt_depth_from_beneath.data = opt_depth_from_beneath[{'pza': slice(None, None, -1)}].data
            self.xd_opt_depth[{levels_out_type: ilevel}] += opt_depth_from_beneath
            # Replace optical depths of zero with nan
            # self.xd_opt_depth[{levels_out_type: ilevel}]
        # Now compute the transmittance between levels as np.exp(optical_depth)
        self.xd_trans = np.exp(-self.xd_opt_depth)
        # TODO : Put in correct long_name and and units (unitless actually)
        # Now run through the cases in the cloud detection sequence to find layers affected by clouds
        if self.has_clouds:
            for librad_case in self.cloud_detect_cases:
                pass  #  TODO : Cloud detection


    def compute_path_radiance(self):
        """ Compute path radiances for path segments between all altitudes in the REM.
        The path transmittances (optical depth) as well as the total radiances at each altitude are required to
        calculate path radiances. If :math:`L_{pi}^{\\downarrow}` is the downwelling path radiance at level :math:`i`
        originating between level :math:`i` and level :math:`i+1` and :math:`L_{pi}^{\\uparrow}` is the upwelling
        path radiance at level :math:`i` originating between level :math:`i` and level :math:`i-1`, then

        .. math::

            L_{pi}^{\\downarrow}=L_{i}^{\\downarrow}-L_{i+1}^{\\downarrow}\\tau_{i+1}^{\\downarrow},

        and

        .. math::

            L_{pi}^{\\uparrow}=L_{i}^{\\uparrow}-L_{i-1}^{\\uparrow}\\tau_{i-1}^{\\uparrow}.

        :math:`L_{i}^{\\uparrow}` is the lower hemisphere (upwelling hemisphere) of the REM and
        :math:`\\tau_{i}^{\\uparrow}` is the transmission between level :math:`i` and level :math:`i+1`.


        .. seealso:: RadEnv.compute_path_transmittance()

        :return: None
        """
        # First compute the product of radiance and transmittance at every level
        self.xd_uu_times_tau = self.xd_trans * self.xd_uu
        # The upper and lower hemispheres of the path radiance REMs have to be computed separately as shown in
        # the docstring equations. Typical up/downwelling indexing is xd_uu[xd_uu['pza'] < np.pi/2].
        # Run through the levels first from bottom level to top, computing the downwelling path radiance
        # to the level above.
        # Initialise path radiance to something having same units etc.
        self.xd_path_radiance = copy.deepcopy(self.xd_uu)  # Initialise to total radiance
        for ilevel in range(0, self.n_levels_out - 1):
            # Set up hemispherical indexing for this level
            hemi_index_level = {self.levels_out_type: ilevel, 'pza': self.xd_path_radiance['pza'] > np.pi/2}
            # Setup hemisperhical indexing for next level up
            hemi_index_above = {self.levels_out_type: ilevel + 1, 'pza': self.xd_path_radiance['pza'] > np.pi/2}
            self.xd_path_radiance[hemi_index_level] = (self.xd_uu[hemi_index_level] -
                                                       self.xd_uu_times_tau[hemi_index_above])
        # Now run from top level down to bottom, computing upwelling path radiance
        for ilevel in range(self.n_levels_out - 1, 1, -1):
            # Set up hemispherical indexing for this level
            hemi_index_level = {self.levels_out_type: ilevel, 'pza': self.xd_path_radiance['pza'] < np.pi/2}
            # Setup hemisperhical indexing for next level down
            hemi_index_below = {self.levels_out_type: ilevel - 1, 'pza': self.xd_path_radiance['pza'] < np.pi/2}
            self.xd_path_radiance[hemi_index_level] = (self.xd_uu[hemi_index_level] -
                                                       self.xd_uu_times_tau[hemi_index_below])



    def write_openexr(self, filename, chan_names=None, chan_per_exr=3, normalise=False, half=False, repeat_azi=1,
                      use_mitsuba_wvl=False):
        """Write a radiant environment as an OpenEXR file or set of OpenEXR files.

        Use of this method requires that the OpenEXR package be installed on your platform. This can be a problem
        for Windows. However, it should be possible to find a Python wheel (.whl) file for Windows which will
        enable installation if pip cannot do the job.

        :param filename: Name of the file to which the OpenEXR environment map should be written. If multiple exr
            files are written, this will be the filename prefix. No default.
        :param chan_names: String or list of strings giving the names to be used for the channels. For example
            give 'RGB' if channels 'R', 'G' and 'B' if RGB files are to be written. This is the same as giving
            ['R', 'G', 'B']. If there are insufficient channel names for all the wavelengths, the names are
            reused.
        :param chan_per_exr: Number of spectral channels to write per exr file. Default is 3.
        :param normalise: Boolean. Normalise all channels to the maximum value (in any one file). Default False.
            Normalisation will make ti possible to display the .exr file in IrfanView or mtsgui (Mitsuba GUI), but
            will destroy the radiometric correctness.
        :param half: Boolean. Use float16 instead of default float32. Halves data size at cost of radiometric
            resolution.
        :param repeat_azi: Integer. Repeat the data in the REM azimuthal direction. This is mostly a convenience for
            display purposes when the REM has no azimuthal variation (e.g. in the thermal bands).
        :param use_mitsuba_wvl: Boolean. If set True, write channels using the Mitsuba wavelength assignments i.e.
            remap wavelengths to the 360 nm to 830 nm range. Setting this flag True will also scale the radiance
            values to W/m^2/sr/nm, which are the canonical units for Mitsuba rendering in the context of MORTIlibradCIA.

        Notes
        -----
        Polarization is not dealt with. Only the first Stokes component (I) is written to the EXR files.
        The number of channel names should be equal to the number of channels per EXR file, but this is not
        enforced.
        Irfanview can display three-component, normalised EXR files only. The Mitsuba GUI (mtsgui) can display
        EXR files with any number of channels, but it is necessary to step through the channels (using [ and ])
        and they are displayed in grayscale. Even mtsgui will clip EXR files having radiance values exceeding 1.0.


        :return:
        """
        import OpenEXR
        import Imath
        wvl = self.xd_uu.wvl.data  # Might not exist, but assume it does
        n_wvl = np.float(wvl.size)  # Number of wavelengths
        wvl_units = self.xd_uu.wvl.units
        rad_units = self.xd_uu.units
        zout = self.xd_uu.zout.data  # output levels
        n_zout = zout.size  # Number of output levels
        n_exr_files = np.int(np.ceil(n_wvl / chan_per_exr))
        # Currently cannot handle polarization, stick to output of the first stokes component (I)
        xd_uu = self.xd_uu[:,:,:,:,0].data  # Zero in last index is stokes I
        if self.hemi:  # Need to double up by reflection left right
            # Note that this puts azimuth
            xd_uu = np.concatenate((xd_uu, xd_uu[:, -2:0:-1, ...]), axis=1)
        datashape = xd_uu[:,:,0,0].shape
        if half:  # Use 16-bit floating point values
            np_type = np.float16
            imath_type = Imath.PixelType.HALF
        else:  # Use 32 bit floating point values
            np_type = np.float32
            imath_type = Imath.PixelType.FLOAT
        imath_chan = Imath.Channel(Imath.PixelType(imath_type))
        # Generate channel names
        if use_mitsuba_wvl:  # Map wavelength channels to equally space channels in the 360 nmto 830 nm range
            channel_names=[]
            channel_interval = (830.0 - 360.0) / (chan_per_exr)
            chan_start = np.arange(360.0, 830.0, channel_interval)
            chan_stop = chan_start + channel_interval
            for i_chan, chan_start_wvl in enumerate(chan_start):
                channel_names.append('{:.2f}'.format(chan_start[i_chan]) + '-' + '{:.2f}'.format(chan_stop[i_chan]) +
                                     wvl_units)
            if rad_units[0:2] == 'mW':
                xd_uu = xd_uu/1000.0  # Convert to W/m^2 ....
                rad_units = rad_units[1:]

        else:  # Use provided channel names or default channel names
            if chan_names is None:
                channel_names = []
                if chan_per_exr == 3:
                    channel_names = 'RGB'
                else:
                    for i_chan, wvl_chan in enumerate(wvl):
                        channel_names.append('{:.3f}'.format(wvl_chan) + '-' + '{:.3f}'.format(wvl_chan) + wvl_units)
            else:
                channel_names = chan_names
        n_channel_names = len(channel_names)
        # Number of channel names should be equal to the number of channels per exr
        if n_channel_names != chan_per_exr:
            warnings.warn('Number of channel names should equal number of channels per EXR file when writing'
                          '.exr files from REMs.')
        # Run through EXR files and write them out
        for i_zout, zout_value in enumerate(zout):
            for i_exr in range(n_exr_files):
                # Generate the file header
                theheader = OpenEXR.Header(datashape[1] * repeat_azi, datashape[0])
                # Add some metadata
                theheader['generatedBy'] = 'libraddask.rad.librad'
                theheader['zout'] = str(zout_value)
                theheader['zout_units'] = 'km'
                theheader['rad_units'] = rad_units
                theheader['wvl_units'] = wvl_units
                theheader['wvl'] = 'np.' + repr(wvl)
                channels = []
                i_channels = []
                chan_data = {}
                chan_max = 0.0
                for i_chan in range(i_exr * chan_per_exr, (i_exr + 1) * chan_per_exr):
                    if i_chan > n_wvl - 1:
                        break
                    # Append the channel names
                    this_channel_name = channel_names[int(np.mod(i_chan, n_channel_names))]
                    channels.append(this_channel_name)
                    i_channels.append(i_chan)
                    # Build the channel data
                    single_chan_data = xd_uu[:,:,i_chan,i_zout].astype(np_type)
                    if repeat_azi > 1:
                        single_chan_data = np.tile(single_chan_data, (1, repeat_azi))
                    if normalise:
                        single_channel_max = np.max(single_chan_data)
                        chan_max = max(chan_max, single_channel_max)
                    chan_data[this_channel_name] = single_chan_data
                if normalise:
                    for channel_key in chan_data:
                        chan_data[channel_key] /= chan_max
                theheader['channels'] = dict([(c, imath_chan) for c in channels])
                theheader['i_channels'] = str(i_channels)
                exr_filename = filename + '{:04d}'.format(min(i_channels)) + '_' + '{:04d}'.format(max(i_channels)) + \
                                'z' + str(zout_value) + '.exr'
                #print(exr_filename)
                exr = OpenEXR.OutputFile(exr_filename, theheader)
                exr.writePixels(dict([this_channel, chan_data[this_channel].tostring()] for this_channel in channels))
                exr.close()


class HyperRadEnv(RadEnv):
    """A higher dimensional form of radiant environment map REM. This class inherits from RadEnv.
     A RadEnv instance has mandatory radiance data dimensions of propagation zenith angle (pza), propaation azimuth
     angle (paz), wavelength, height above surface (zout) and stokes parameter. These are the outputs that can all
     be obtained from a single run of libRadtran. The HyperRadEnv allows for increasing the dimensionality of the
     REM dataset using other (typically scalar) inputs to libRadtran. Most significant amongst these are solar
     zenith angle (sza), surface albedo (albedo) and surface temperature (sur_temperature). Solar zenith angle and
     albedo are only applicable in `source solar` runs, while sur_temperature is only applicable in `source thermal`
     REMs.

     Other possible higher dimensions that can be added include:
     Aerosol loading as specified by visibility or optical depth



    """

    def __init__(self, base_case, n_pol, n_azi, hyper_axes, mxumu=48, mxphi=19, hemi=False, n_sza=0):
        """ Create a set of uvspec runs covering the whole sphere to calculate a full radiant environment map.
        Where the base_case is the uvspec case on which to base the environmental map, Name is the name to give the
        environmental map and n_pol and n_azi are the number of polar and azimuthal sightline angles to generate. The
        mxumu and mxphi are the maximum number of polar and azimuth angles to calculate in a single run of uvspec.
        The default values are mxumu = 48, and mxphi = 19. These values are taken from the standard libRadtran
        distribution (/libsrc_f/DISORT.MXD) maximum parameter file. If using the polradtran solver, the corresponding
        file is /libsrc_f/POLRADTRAN.MXD. Other solvers may have different restrictions. A warning will be issued if
        the solver is not in the DISORT/POLRADTRAN family.

        :param base_case: librad.Case object providing the case on which the environment map is to be based. Note
            that not any basecase can be used. As a general guideline, the basecase should have standard irradiance
            outputs (i.e. should not use the `output_user` keyword). It should also not the use `output_process` or
            `output_quantity` keywords, which change the units and/or format of the libRadtran/uvspec output.
             Minimal validation of the basecase is performed. However, use with `mol_abs_param` such as `kato` and
             `fu` is important for `librad` and these are supported (k-distribution or `correlated-k`
             parametrizations). Use of `output_process per_nm` is appropriate for `source thermal` REMs to get
             radiance units per nanometre rather than per inverse cm.
        :param n_pol: Number of polar angles (view/propagation zenith angles)
        :param n_azi: Number of azimuthal angles.
        :param hyper_axes: Dictionary of libRadtran keywords, with a list (or numpy array) of values to assign  to
            that keyword to create multiple copies of the REM e.g. {'sza': [0.0 30.0 45.0]} will create 3 runs of the
            REM with these solar zenith angles. Note that the number of REMs goes up very rapidly. If 3 surface
            albedo values as well as 3 solar zenith angles, the compute overhead goes up by a factor of 9.
            In general, the keywords should have a scalar and numeric input (such as `albedo`, `sza` and
            `sur_temperature`). This input parameter is mandatory and there is no default.
        :param mxumu: Maximum number of polar angles per case.
        :param mxphi: Maximum number of azimuthal angles per case.
        :param hemi: If set True, will generate only a single hemisphere being on one side of
            the solar principle plane. Default is False i.e. the environment map covers the full sphere.
            Note that if hemi=True, the number of REM samples in azimuth becomes n_azi :math:`\\times` 2.
            This is the recommended mode (hemi=True) for librad purposes, since it reduces execution time.
        :param n_sza: The number of solar zenith angles (SZA) at which to perform transmittance and path radiance
            computations. Each SZA will result in another run of the base case (no radiances)


        """
