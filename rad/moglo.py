__author__ = 'DGriffith'

import numpy as np
from pint import UnitRegistry
ureg = UnitRegistry()
Q_ = ureg.Quantity

import warnings

def U_(units):
    """ Returns a pint.Quantity of magnitude one with the given units

    :param units: A string providing a unit of measure in the pint unit registry
    :return: A pint.Quantity object with a magnitude of 1 and the specified units.
    """
    return Q_(1.0, units)

#  globals and utility functions

# The following are the long names that will be provided for plotting purposes and netCDF conventions
long_name = {
    'wvl': 'Wavelength',
    'wvn': 'Wavenumber',
    'trn': 'Transmission',
    'rad': 'Radiance',
    'irr': 'Irradiance',
    'specrad': 'Spectral Radiance',
    'specirr': 'Spectral Irradiance',
    'spf': 'Linear Spatial Frequency',
    'spfa': 'Angular spatial Frequency',
    'mtf': 'Modulation Transfer Function',
    'otf': 'Optical Transfer Function',
    'ptf': 'Phase Transfer Function',
    'fno': 'Focal Ratio',
    'efl': 'Effective Focal Length',
    'lum': 'Luminance',
    'ill': 'Illuminance',
    'flux': 'Optical Flux',
    'fldy': 'Field Position in x',
    'fldx': 'Field Position in y',
    'fldr': 'Radial Field Position',
    'fldz': 'Defocus',
    'fldo': 'Field Orientation',  # Typically horizontal/vertical, across/along track, sagittal/tangential
    'obs': 'Obscuration Ratio',
    'pitchx': 'Pixel Pitch in x',
    'pixapx': 'Pixel Aperture in x',
    'pixapy': 'Pixel Aperture in y',
    'pitchy': 'Pixel Pitch in y',
    'asr': 'Absolute Spectral Response',
    'sqe': 'Spectral Quantum Efficiency',
    'chn': 'Spectral Channel Number',
    'rsr': 'Relative Spectral Response',
    'srf': 'Spectral Response Function',
    'phe': 'Photoelectrons',
    'dn': 'Digital Numbers',
    'wellcap': 'Pixel Well Capacity',
    'readnoise': 'RMS Readout Noise',
    'darkcurr': 'Pixel Dark Current',
    'psnu': 'Photo-Response Non-Uniformity',
    'dsnu': 'Dark Signal Non-Uniformity',
    'tref': 'Reference Temperature',
    'temp': 'Temperature',
    'deltat': 'Temperature Delta',
    'nyqx': 'Nyquist Frequency in x',
    'nyqy': 'Nyquist Frequency in y',
    'bitdepth': 'Bit Depth',
    'dgain': 'Digital Gain',
    'doffset': 'Digital Offset',
    'dnoise': 'RMS Digital-Equivalent Noise',
    'spfcut': 'Diffraction Cutoff Spatial Frequency',
    'vza': 'View Zenith Angle',  # degrees
    'vpa': 'View Polar Angle',  # same as vza, but in radians
    'vea': 'View Elevation Angle',  # complement of the VZA
    'vaz': 'View Azimuth Angle',
    'svaz': 'Solar-Relative View Azimuth Angle',
    'pza': 'Propagation Zenith Angle',  # in radians
    'umu': 'Cosine of the Propagation Zenith Angle',
    'phi': 'Propagation Azimuth Angle',  # in unit of degrees
    'paz': 'Propagation Azimuth Angle',  # but in units of radians
    'phi0': 'Solar Propagation Azimuth Angle',
    'stokes': 'Stokes Parameter',  # I=0, Q=1, U=2, V=3
    'zout': 'Height Above Surface',
    'zout_sea': 'Height Above Sea Level',
    'p': 'Pressure',
    'sph_alb': 'Atmospheric Spherical Albedo',
    'reflx': 'Reflectivity',
    'trnx': 'Transmissivity',
    'edir': 'Direct Solar Horizontal Irradiance',
    'edn': 'Downwelling Diffuse Irradiance',
    'eup': 'Upwelling Diffuse Irradiance',
    'eglo': 'Total Downwelling Irradiance',
    'enet': 'Net Downward Irradiance',
    'uavgdir': 'Direct Solar Actinic Sterance',
    'uavgdn': 'Diffuse Downward Actinic Sterance',
    'uavgup': 'Diffuse Upward Actinic Sterance',
    'uavgglo': 'Global Downward Actinic Sterance',
    'down_fluxI': 'Global Downwelling Irradiance Stokes I',
    'down_fluxQ': 'Global Downwelling Irradiance Stokes Q',
    'down_fluxU': 'Global Downwelling Irradiance Stokes U',
    'down_fluxV': 'Global Downwelling Irradiance Stokes V',
    'up_fluxI': 'Global Upwelling Irradiance Stokes I',
    'up_fluxQ': 'Global Upwelling Irradiance Stokes Q',
    'up_fluxU': 'Global Upwelling Irradiance Stokes U',
    'up_fluxV': 'Global Upwelling Irradiance Stokes V',
    'brightness': 'Brightness Temperature',
    'od': 'Optical Depth'

    # TODO : Include all libRadtran definitions as well from rad.librad.py
}

# The following are intended to be the standard names as per the Climate and Forecast (CF) convention
standard_names = {

}

# Default units are the units that are preferred and to which librad will convert if other units are given
default_units = {
    'wvl': 'nm',
    'wvn': 'cm^-1',
    'spf': '1/mm',
    'efl': 'mm',
    'rad': 'W/m^2/sr',
    'irr': 'W/m^2',
    'specrad': 'W/m^2/sr/nm',
    'specirr': 'W/m^/nm',
    'fldx': 'mm',
    'fldy': 'mm',
    'fldr': 'mm',
    'fldz': 'mm',
    'flux': 'W',
    'asr': 'A/W',
    'sqe': '',
    'lum': 'cd/m^2',
    'ill': 'lux',
    'fldo': 'deg',
    'phe': 'e',
    'dn': 'count',
    'pitchx': 'mm',
    'pitchy': 'mm',
    'pixapx': 'mm',
    'pixapy': 'mm',
    'wellcap': 'e',
    'readnoise': 'e',
    'darkcurr': 'e/s',  # per pixel
    'prnu': '%',  # Must actually be in percentage  This is a special case and needs special handling
    'dsnu': '%',  # Must actually be in percentage  This is a special case and needs special handling
    'tref': 'degC',  # Reference temperature
    'temp': 'degC',  # Operating temperautre
    'deltat': 'delta_degC',  # Change in temperature
    'nyqx': '1/mm',  # Nyquist x
    'nyqy': '1/mm',  # Nyquist y
    'bitdepth': 'bit',
    'dgain': 'e/count',
    'doffset': 'count',
    'dnoise': 'count',
    'mtf': '',
    'spfcut': '1/mm',
    'spfa': '1/mrad',
    'vza': 'deg',  # View zenith angle
    'vpa': 'rad',  # Same as vza, but in radians
    'vaz': 'deg',  # View azimuth angle
    'vea' : 'deg',  # View elevation angle (complement of the view zenith angle)
    'pza': 'rad',  # Propagation zenith angle
    'svaz': 'deg',  # Solar-relative view azimuth angle
    'umu': '',  # Cosine of the propagation zenith
    'phi': 'deg',  # Propagation azimuth angle, same as paz, except for units
    'paz' : 'rad',  # same as phi except for units
    'phi0': 'deg',  # Solar propagation azimuth angle
    'srf': '',  # spectral response function
    'stokes': '',  # Stokes parameter number (not the actual stokes parameter)
    'zout': 'km',  # Altitude above ground level
    'zout_sea': 'km',  # Altitude above sea level
    'p': 'hPa',  # Atmospheric pressure preferred units : hectopascals
    'reflx': '',  # reflectivity
    'trnx': '',  # transmissivity
    'brightness': 'K'
}

class Scalar(object):
    """ The Scalar class is for representation of scalar numeric values, together with units of measure, a
    mnemonic, being one of those listed in the librad long_name vocabulary.

    A feature of the Scalar class is that the data attribute can only be set using a.data = [value, units]

    """

    @property
    def data(self):
        return self.__data

    @data.setter
    def data(self, data_and_units):
        value = Q_(np.asarray(data_and_units[0], dtype=np.float64), data_and_units[1]) # Will blow up if units unknown to pint
        self.__units = data_and_units[1]
        if self.mnemonic in default_units:
            # Convert units
            value = value.to(default_units[self.mnemonic])  # Will blow up if units cannot be converted to default
            self.__units = default_units[self.mnemonic]
        self.__data = value.magnitude

    @property
    def units(self):
        return self.__units

    @units.setter
    def units(self, value):
        warnings.warn('Attempt to alter scalar units ignored.')

    @property
    def mnemonic(self):
        return self.__mnemonic

    @mnemonic.setter
    def mnemonic(self, value):
        warnings.warn('Attempt to alter scalar mnemonic ignored.')

    def __init__(self, mnemonic, data, units, attrs=None):
        """ Construct a scalar value using a librad mnemonic (from the long_name vocabulary), the scalar value
        and the units of measure.

        :param mnemonic: A string such as 'efl' for Effective Focal Length, which is defined in the librad
            long_name vocabulary
        :param data: The scalar numeric value to be stored
        :param units: The units string of the data. This string must be known to the Python Pint package.
        :param attrs: A user-defined dictionary of any other attributes for the scalar value. The mnemonic will
            be looked up in moglo.long_name and the long_name will be added as an attribute.
        :return:
        """
        self.__mnemonic = mnemonic
        self.__units = units
        if attrs is None:
            self.attrs = {}
        else:
            self.attrs = attrs
        try:
            self.attrs['long_name'] = long_name[mnemonic]  # Will blow up if the mnemonic is not known
        except KeyError:
            raise KeyError('Unknown scalar mnemonic ' + mnemonic + ' provided in Scalar instantiation.')
        value = Q_(np.asarray(data, dtype=np.float64), units)  # Will blow up if units unknown to pint
        if mnemonic in default_units:
            # Convert units
            value = value.to(default_units[mnemonic])  # Will blow up if units cannot be converted to default
            self.__units = default_units[mnemonic]
        self.__data = value.magnitude

    def __repr__(self):  #  TODO : Should we be using __repr__ or __str__
        return self.attrs['long_name'] + ' : ' + str(self.__data) + ' ' + self.__units

    # TODO : Consider doing all arithmetical magic methods for Scalar class.
    # TODO : Problem is that of catering for multiplcation of 2 scalar values, using isinstance considered unpythonic

    def __pos__(self):
        return +self.__data

    def __neg__(self):
        return -self.__data

    def __abs__(self):
        return abs(self.__data)

    def __mul__(self, other):
        if isinstance(other, Scalar):
            return self.__data * other.__data
        else:
            return self.__data * other  # Hope that other defines multiplication with a simple scalar numeric

    def __add__(self, other):
        if isinstance(other, Scalar):
            return self.__data + other.__data
        else:
            return self.__data + other

    def __radd__(self, other):
        if isinstance(other, Scalar):
            return other.__data + self.__data
        else:
            return other + self.__data

    def __sub__(self, other):
        if isinstance(other, Scalar):
            return self.__data * other.__data
        else:
            return self.__data * other

    def __div__(self, other):
        if isinstance(other. Scalar):
            return self.__data / other.__data
        else:
            return self.__data / other

    def __rdiv__(self, other):
        if isinstance(other, Scalar):
            return other.__data / self.__data
        else:
            return other / self.__data



# Define a global exception for unit mismatch
class UnitMismatch(Exception):
    pass

# Define a global exception for missing units
class MissingUnits(Exception):
    pass

class MissingDataArrayAxis(Exception):
    pass