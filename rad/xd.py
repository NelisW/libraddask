__author__ = 'DGriffith'

# Functions related to interpolation of xarray.DataArray and other utilities

# Some sort of global data dictionary (CF compliant ?), including short names and synonyms
# Could import some or all CF definitions from XML file.

import numpy as np
import xarray
from scipy.interpolate import RegularGridInterpolator, interp1d
import warnings
from operator import mul, add
import copy
import functools

from libraddask.rad import moglo #units and vocabulary

def xd_identity(np_vector, axis_name, units=None, attrs=None):
    """ Create an identity xarray.DataArray. That is, a DataArray vector in which both the values and axis
        coordinates are identical.

    :param np_vector: Vector of numeric data
    :param axis_name: Name for the axis - must be in vocabulary defined in moglo.py
    :param units: The units of the np_vector data. If the units are not in the default units, the Python pint package
        is used to make a conversion. If units are not given, it is assumed that the data is already in the
        default units for the quantity named in axis_name. It is better to provide the units if in any way
        unsure.
    :param attrs: Dictionary of additional attributes to attach to the DataArray
    :return:
    """
    if axis_name in long_name:
        the_long_name = long_name[axis_name]
    else:
        warnings.warn('Unknown axis name ' + axis_name + ' encountered in xd_identity creation.')
    if axis_name in default_units:
        the_units = default_units[axis_name]
    else:
        the_units = ''  # Assumed to be unitless quantity
    if units is None:
        units = the_units
    values = moglo.Q_(np_vector, units)  # Create a pint quantity with the given units
    values = values.to(the_units)  # actually the default units
    np_vector = values.magnitude
    if attrs is not None:
        the_attrs = attrs
    else:
        the_attrs = {}
    the_attrs.update({'long_name': the_long_name})
    the_attrs.update({'units': the_units})
    return xarray.DataArray(np_vector, [(axis_name, np_vector)], name=axis_name, attrs=the_attrs)


def xd_harmonise_interp(xd_list):
    """ Perform linear interpolation on merged set of axis points for two or more xarray DataArray objects.
        This function can be used to prepare (harmonise) multiple xarray.DataArray objects for multiplication or addition
        on a common set of coordinate axis points by linearly interpolating all DataArray objects onto the same
        set of points, obtained by merging and sorting the points from all input DataArray objects.

    The DataArry objects provided.    The scipy linear grid interpolator is used for this purpose. See:
    scipy.interpolate.RegularGridInterpolator
    :param xd_list:
    :return: Tuple of xarray.DataArray objects with merged and linearly interpolated values in all axes.
    Only unique values in the interpolation axis are used.

    """
    # TODO : enforce compatible attributes or not ? What attributes in returned object ?
    # TODO : Ignore axes that have non-numeric coordinates e.g. xdarray['axisname'].dtype.char in 'SUa',
    # TODO : which detects dtypes that are string, or xdarray['axisname'].dtype.kind in 'fc' (float or complex)
    # TODO : alternatively require that axis coordinates are always numeric, say with a list of labels as attrs
    # TODO : What about interpolation on times axes
    # TODO : Need to expand on extrapolation (and possibly also single-axis interpolation) schemes
    # Accumulate the index values from each of the given arrays, for each of the axes in the first array
    index_vals = {}  # dictionary of index coordinates for each axis
    index_float = {}  # determine if the index kind is a floating point type (complex included)
    #metadata = {}

    for xd_arr in xd_list:
        for axis in xd_arr.dims:
            # accumulate dictionary for all dimensions in the entire collection of DataArrays
            if not axis in index_vals:
                index_vals[axis] = xd_arr[axis]
            else:
                index_vals[axis] = np.hstack((index_vals[axis], xd_arr[axis]))
            # also accumulate the attributes (metadata)
            # metadata.update(xd_arr.attrs)
    # get the unique values in increasing numerical order using np.unique for each axis found in the whole set
    for axis in index_vals:
        index_vals[axis] = np.unique(index_vals[axis])
        index_float[axis] = index_vals[axis].dtype.kind in 'fc'
    # interpolate each of the DataArray objects onto the new grid (for whatever axes it does have)
    xd_return_list = []
    for xd_arr in xd_list:
        # Create the linear interpolator
        interpolator = RegularGridInterpolator([xd_arr[axis].values for axis in xd_arr.dims],
                                               xd_arr.values,
                                               method='linear', bounds_error=False, fill_value=0.0)
        merged_coordinates = np.meshgrid(*[index_vals[axis] for axis in xd_arr.dims],
                                               indexing='ij')
        interp_vals = interpolator(tuple(merged_coordinates))
        # reconstruct the xarray.DataArray with interpolated data
        xd_arr_interp = xarray.DataArray(interp_vals, [(axis, index_vals[axis]) for axis in xd_arr.dims],
                                       name=xd_arr.name, attrs=xd_arr.attrs)
        xd_arr_interp.attrs = xd_arr.attrs  # transfer the attributes verbatim
        xd_return_list.append(xd_arr_interp)

        # There may be axes not present in a specific DataArray. These are omitted for that DataArray and
        # simply allowed to broadcast when performing operations with other DataArrays
    return xd_return_list

def xd_interp_axis_to(from_xd, to_xd, axis, interp_method='linear', bounds_error=False, fill_value=0.0,
                      assume_sorted=True):
    """ Interpolate a single xarray.DataArray axis from one set of coordinates to another. Since interpolation
    occurs along a single axis, there is more flexibility in the method of interpolation that can be used.
    The `scipy.interpoalte.interp1d` class is used to perform the interpolation.

    :param from_xd: The xarray>DataArray object with originating data.
    :param to_xd: The xarray>DataArray object that will provide the new coordinates to which the interpolation will
        be carried out.
    :param axis: The name of the axis along which to perform the interpolation.
    :param interp_method: Is the kind of interpolation to perform. Options are as for sipy.interpolate.interp1d,
        namely 'linear', 'nearest', 'zero', 'slinear', 'quadratic' and 'cubic', where 'slinear', 'quadratic' and
        'cubic' produce spline interpolation of first, second or third order respectively. The default is 'linear'.
    :return: New xarray.DataArray with xd_from interpolated along given axis to coordinates provided by xd_to in
        the given axis.
    """
    from_dims = from_xd.dims
    from_axes = [copy.deepcopy(from_xd[this_axis]) for this_axis in from_dims]
    interp_func = interp1d(from_xd[axis].data, from_xd.data, kind=interp_method, axis=from_xd.get_axis_num(axis),
                            copy=False, bounds_error=bounds_error, fill_value=fill_value, assume_sorted=assume_sorted)
    new_data = interp_func(to_xd[axis].data)  # Interpolate along the named axis
    # Now reconstruct the xd_from DataArray
    from_axes[from_xd.get_axis_num(axis)] = to_xd[axis]  # Grab the new axis from the xd_to DataArray
    new_from_xd = xarray.DataArray(new_data, from_axes, attrs=from_xd.attrs)  # Use attributes from original
    return new_from_xd


def xd_harmonised_product(xd_list,raise_exception=True):
    """ Compute the harmonised product of a number of N-dimensional data arrays.
        The DataArrays are interpolated onto a common set of coordinates and then the product of the DataArrays
        is computed, returning a single DataArray with merged attributes. Unit mismatches are flagged with warnings.

    :param xd_list: List/tuple of xarray.DataArray objects to be multiplied
    :return: Product of xarray.DataArray objects with merged attributes
    :except UnitMismatch, MissingUnits:
    """
    # TODO : This function to be checked to correct "var_units" mistake

    #main_attrs = {}  # Will accumulate all main attributes here - not sure what to do with units ?
    unit_dict = {}  #  Dictionary of units
    axis_attrs = {}  # Dictionary of axis attribute dictionaries
    # Check units and merge metadata
    # have to merge attributes for main data and all axes individually
    for xd_arr in xd_list:
        #main_attrs.update(xd_arr.attrs)
        for axis in xd_arr.dims:
            if axis in axis_attrs:
                axis_attrs[axis].update(xd_arr[axis].attrs)  # Accumulate the attributes for each axis
            else:
                axis_attrs[axis] = xd_arr[axis].attrs
            if not axis in unit_dict:
                if 'units' in xd_arr[axis].attrs:
                    unit_dict[axis] = xd_arr[axis].attrs['units']
                else:
                    if raise_exception:
                        raise MissingUnits('Units not found for ' + xd_arr.name + ' on axis ' + axis)
                    else:
                        print(MissingUnits('Units not found for ' + xd_arr.name + ' on axis ' + axis))
            elif ('units' in xd_arr[axis].attrs['units']) and (unit_dict[axis] != xd_arr[axis].attrs['units']):
                # TODO : Consider throwing a unit mismatch error, or converting to desired units with pint
                warnings.warn('Unit mismatch found when taking xarray.DataArray harmonised product.')
                if raise_exception:
                    raise UnitMismatch('Unit mismatch encountered for ' + xd_arr.name + ' on axis ' + axis)
                else:
                    print(UnitMismatch('Unit mismatch encountered for ' + xd_arr.name + ' on axis ' + axis))
            else:
                if raise_exception:
                    raise MissingUnits('Units not found for ' + xd_arr.name + ' on axis ' + axis)
                else:
                    print(MissingUnits('Units not found for ' + xd_arr.name + ' on axis ' + axis))


    xd_factors = xd_harmonise_interp(xd_list)
    xd_product = functools.reduce(mul, xd_factors)  # take the product by reducing the list using the mul operator
    for axis in xd_product.dims:  # Put the merged attributes into each of the axes
        xd_product[axis].attrs = axis_attrs[axis]
    return xd_product


def check_convert_units(value_with_units, preferred_units):
    """ Check the units of a quantity and convert to preferred units using Python `pint`

    :param value_with_units: A list with a numeric value or numpy array in the first position and a string
        providing units in the second position. The unit string must be recognisable by the Python `pint` package.
    :param preferred_units: A string expressing the units to which `pint` should convert the scalar
    :return: Value expressed in the preferred units
    """

    # Use pint to convert
    try:
        value = moglo.Q_(np.asarray(value_with_units[0], dtype=np.float64), value_with_units[1])  # Will blow up if units not recognised
    except TypeError:
        warnings.warn('A scalar value was supplied without units. Example of correct scalar input is [40.0, "degC"]')
        return None
    value = value.to(preferred_units)
    return value.magnitude


def xd_check_convert_units(xd, axis_name, preferred_units):
    """ Check and convert units for one or more axes of an `xarray.DataArray`

    :param xd: An xarray.DataArray object having an axis called `axis_name` and a value in the `attrs` dictionary
    :param axis_name: Name of the axis or data to check/convert
    :param preferred_units: A string providing the preferred units that can be passed to `pint`
    :return: A xarray.DataArray, in which the values in the named axis have been converted to the preferred units
        The `axis_name_units` field is also updated.
    """

    # Create a pint.Quantity object using the data from the named array and use that to convert to
    # preferred units
    if axis_name == xd.name:  # The primary data must be converted
        Q_values = moglo.Q_(xd.data, xd.units)  # Can fetch units this way, but not set them
        Q_values = moglo.Q_values.to(preferred_units)
        xd.data = moglo.Q_values.magnitude
        xd.attrs['units'] = preferred_units
    else:  # Convert units of the named axis
        Q_values = moglo.Q_(xd[axis_name].data, xd[axis_name].units)
        Q_values = moglo.Q_values.to(preferred_units)
        xd[axis_name].data = moglo.Q_values.magnitude
        xd[axis_name].attrs['units'] = preferred_units


def xd_attrs_update(xd_list):
    """ Update long_name and units attributes of all axes in an xarray.DataArray
    :param xd: Input xarray.DataArray
    :return:
    """

    for xd_arr in xd_list:
        xd_arr.attrs['long_name'] = long_name[xd_arr.name]
        for axis in xd_arr.dims:
            xd_arr[axis].attrs['long_name'] = long_name[axis]  # Will blow up if axis mnemonic name not found
            xd_arr[axis].attrs['units'] = default_units[axis]  # Likewaise if units not found for this axis name


