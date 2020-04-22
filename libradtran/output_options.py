"""--------------------------------------------------------------------
 * $Id: output_options.py 3523 2019-12-20 13:14:51Z Claudia.Emde $
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

class setup_output_group():

    group_name = 'Output and misc.'

    def __init__(self):
        documentation = get_output_documentation()

        include = option_definition.not_yet_lex2py_option(
            name='include',
            group='output',
            helpstr='Include a file into uvspec input.',
            documentation=documentation['include'],
            gui_inputs=(GUI_definition.FileInput(name='filename'),),
            parents=['uvspec'],
            showInGui=False,
        )

        quiet = option_definition.option(
            name='quiet',
            group='output',
            helpstr='If specified, informative messages are turned off. See also \code{verbose}.',
            documentation=documentation['quiet'],
            tokens=option_definition.addSetting(name='Input.quiet', setting=1, default=0),
            parents=['uvspec'],
            non_parents=['verbose'],
        )

        verbose = option_definition.option(
            name='verbose',
            group='output',
            helpstr='If specified abundances of informative messages are output to stderr.',
            documentation=documentation['verbose'],
            tokens=option_definition.addSetting(name='Input.verbose', setting=1, default=0),
            parents=['uvspec'],
            non_parents=['quiet'],
        )

        write_ext_to_file = option_definition.option(
            name='write_ext_to_file',
            group='output',
            helpstr='If specified abundances of informative messages are output to stderr.',
            documentation=documentation['write_ext_to_file'],
            tokens=option_definition.addSetting(name='Input.write_ext_to_file', setting=1, default=0),
            parents=['uvspec'],
            islidar=True,
        )

        test_optical_properties = option_definition.option(
            name='test_optical_properties',
            group='output',
            helpstr='Write optical properties to file without solving RTE.',
            documentation=documentation['test_optical_properties'],
            tokens=option_definition.addSetting(name='Input.test_optical_properties', setting=1, default=0),
            parents=['uvspec'],
            developer=True,
        )

        print_disort_info = option_definition.option(
            name='print_disort_info',
            group='output',
            helpstr='Print various disort output. See disort.doc', 
            documentation=documentation['print_disort_info'],
            tokens=[option_definition.addToken(name="Input.rte.prndis", datatype=option_definition.Integers),
                option_definition.addSetting(name='Input.rte.nprndis',setting='ntokens')],
        )

        data_files_path = option_definition.option(
            name='data_files_path', 
            group='output', 
            documentation=documentation['data_files_path'],    
            gui_inputs=(GUI_definition.FileInput(name='Input.filename[FN_PATH]'),),
            tokens = option_definition.addToken(name='Input.filename[FN_PATH]' , datatype=io.IOBase),    
        )

        output_process = option_definition.option(
            name='output_process',
            group='output',
            helpstr='Decide how the output from \code{uvspec} is processed.',
            documentation=documentation['output_process'],
            gui_inputs=(GUI_definition.ListInput(name='Input.output_unit_processing', valid_range=['integrate', 'sum', 'rgbnorm', 'rgb_norm', 'rgb', 'per_nm', 'per_cm-1', 'per_band'], optional=False),),
            tokens=option_definition.addToken(name='Input.output_unit_processing', datatype=str, valid_range=['integrate','sum','rgbnorm','rgb_norm','rgb',
                                           'per_nm','per_cm-1','per_band']),
            parents=['uvspec'],
        )

        output_file = option_definition.option(
            name='output_file',
            group='output', 
            helpstr='Location of the output file', 
            documentation=documentation['output_file'], 
            gui_inputs=(GUI_definition.TextInput(name='Input.filename[FN_OUTFILE]'),),
            tokens=option_definition.addToken(name='Input.filename[FN_OUTFILE]', datatype=str),
            parents=['uvspec'],
        )

        output_format = option_definition.option(
            name='output_format', 
            group='output',
            helpstr='Output format',
            documentation=documentation['output_format'],
            tokens=option_definition.addLogical(name='Input.output_format', logicals=['ascii','flexstor','netCDF','sat_picture'], setting='OUTPUT_FORMAT_'), 
        )

        output_user = option_definition.not_yet_lex2py_option(
            name='output_user',
            group='output',
            helpstr='User defined output', 
            documentation=documentation['output_user'],
            gui_inputs=(GUI_definition.TextInput(name=''),),
            tokens=option_definition.addToken(name="", datatype=str),
            parents=['uvspec'],
        )

        output_quantity = option_definition.option(
            name='output_quantity',
            group='output',
            helpstr='Define output quantity to calculate instead of absolute quantities.',
            documentation=documentation['output_quantity'],
            tokens=option_definition.addLogical(name='Input.calibration', logicals=['reflectivity','transmittance','brightness'], setting='OUTCAL_'),
            parents=['uvspec'],
        )

        heating_rate = option_definition.option(
            name='heating_rate',
            group='output',    
            helpstr='Calculation of heating rates.',
            documentation=documentation['heating_rate'],
            gui_inputs=(GUI_definition.ListInput(name='Input.heating', valid_range=['none', 'local', 'layer_fd', 'layer_cd'], optional=True),),
            tokens = [option_definition.addSetting(name='Input.heating', setting='HEAT_LAYER_CD'),
                      option_definition.addLogical(name='Input.heating', logicals=['none','local','layer_fd','layer_cd'], setting='HEAT_', optional=True) ] ,
            parents=['uvspec'],
        )

        write_output_as_netcdf = option_definition.option(
            name='write_output_as_netcdf',
            group='output',
            helpstr='Output is written into netcdf file.',
            documentation=documentation['write_output_as_netcdf'],
            tokens=option_definition.addSetting(name='Input.write_output_as_netcdf', setting=1),
            parents=['uvspec'],
            islidar=True,
            showInGui=False,
        )

        slit_function_file = option_definition.option(
            name='slit_function_file',
            group='output',
            helpstr='If specified, the calculated spectrum is convolved with the function found in the slit_function_file.',
            documentation=documentation['slit_function_file'],
            gui_inputs=(GUI_definition.FileInput(name='Input.filename[FN_SLITFUNCTION]'),),
            tokens = [option_definition.addToken(name='Input.filename[FN_SLITFUNCTION]', datatype=io.IOBase),
                      option_definition.addSetting(name='Input.convolve', setting=1, default=0)],
            parents=['uvspec'],
            non_parents=['filter_function_file','thermal_bands_file'],
            plot = {'plot_type': '2D',
                'optional_args': {'column_names': (
                        "wavelength",
                        "slit function",)
                          }
                }
        )

        spline = option_definition.option(
            name='spline',
            group='output',
            helpstr='Spline interpolate the calculated spectrum.',
            documentation=documentation['spline'],
            gui_inputs=(GUI_definition.FloatInput(name='Input.spline_lambda_0', valid_range=[0, 1000000.0]), GUI_definition.FloatInput(name='Input.spline_lambda_1', valid_range=[0, 1000000.0]), GUI_definition.FloatInput(name='Input.spline_lambda_step', valid_range=[0, 1000000.0]),),
            tokens = [option_definition.addToken(name='Input.spline_lambda_0', datatype=float, valid_range=[0,1e6]), #TODO:Valid_range from option wavelength.
                      option_definition.addToken(name='Input.spline_lambda_1', datatype=float, valid_range=[0,1e6]),
                      option_definition.addToken(name='Input.spline_lambda_step', datatype=float, valid_range=[0,1e6]),
                      option_definition.addSetting(name='Input.spline', setting=1, default=0)],
            parents=['uvspec'],
            non_parents=['spline_file'],
        )

        spline_file = option_definition.option(
            name='spline_file',
            group='output',
            helpstr='Spline interpolate to arbitrary wavelengths, in nm, given as a single column in file.',
            documentation=documentation['spline_file'],
            gui_inputs=(GUI_definition.FileInput(name='Input.filename[FN_SPLINE]'),),
            tokens = [option_definition.addToken(name='Input.filename[FN_SPLINE]', datatype=io.IOBase),
                      option_definition.addSetting(name='Input.spline', setting=1, default=0) ] ,
            parents=['uvspec'],
            non_parents=['spline'],
        )

        filter_function_file = option_definition.option(
            name='filter_function_file',
            group='output',
            helpstr='If specified, the calculated spectrum is multiplied with a filter function defined in file.',
            documentation=documentation['filter_function_file'],
            tokens=[option_definition.addToken(name='Input.filename[FN_FILTERFUNCTION]',datatype=io.IOBase),
                    option_definition.addLogical(name='Input.filter_function_normalize', logicals=['normalize'], optional=True)],
            parents=['uvspec'],
            non_parents=['slit_function_file','thermal_bands_file'],
            plot = {'plot_type': '2D',
                'optional_args': {'column_names': (
                        "wavelength",
                        "filter function",)
                          }
                }
        )

        pressure_out = option_definition.not_yet_lex2py_option(
            name='pressure_out',
            group='output',
            helpstr='Specify the output levels in pressure coordinates.',
            documentation=documentation['pressure_out'],
            gui_inputs=(GUI_definition.TextInput(name=''),),
            parents=['uvspec'],
            non_parents=['zout','zout_sea'],
        )

        zout = option_definition.not_yet_lex2py_option(
            name='zout',
            group='output',    
            helpstr='Output altitude(s).', 
            documentation=documentation['zout'],
            gui_inputs=(GUI_definition.TextInput(name=''),),
            tokens=option_definition.addToken(name="", datatype=str),
            parents=['uvspec'],
            non_parents=['zout_sea','pressure_out'],
        )

        zout_sea = option_definition.not_yet_lex2py_option(
            name='zout_sea',
            group='output',    
            helpstr='Output altitude(s) above sea surface.', 
            documentation=documentation['zout_sea'],
            gui_inputs=(GUI_definition.TextInput(name=''),),
            tokens=option_definition.addToken(name="", datatype=str),
            parents=['uvspec'],
            non_parents=['zout','pressure_out'],
        )

        mc_backward_output = option_definition.option(
            name='mc_backward_output', 
            group='output',
            helpstr='Specify quantity to be calculated using backward Monte Carlo.', 
            documentation=documentation['mc_backward_output'], 
            tokens=[option_definition.addLogical(name='Input.rte.mc.backward.output', logicals=['edir','edn','eup','fdir','fdn','fup','act','abs','emis','heat','exp','exn','eyp','eyn','ednpv'], setting='MCBACKWARD_', gui_name='output'),    
                    option_definition.addLogical(name='Input.rte.mc.abs_unit', logicals=['W_per_m2_and_dz', 'W_per_m3', 'K_per_day'], setting='MCABS_UNIT_' , optional = True, gui_name='unit')], 
            parents=['mc_backward'], 
            speaker='rte_solver',
            enable_values=("mystic","montecarlo"),
            mystic =True
        )

        mc_forward_output = option_definition.option(
            name='mc_forward_output',
            group='output',
            helpstr='Specify quantity to be calculated using forward Monte Carlo.',
            documentation=documentation['mc_forward_output'],
            tokens = [option_definition.addLogical(name='Input.rte.mc.absorption', logicals=['absorption','actinic','emission','heating'], setting='MCFORWARD_ABS_', gui_name='output'), 
                      option_definition.addLogical(name='Input.rte.mc.abs_unit', logicals=['W_per_m2_and_dz', 'W_per_m3', 'K_per_day'] , setting='MCABS_UNIT_',  optional = True, gui_name='unit')],
            parents=['uvspec'], 
            non_parents=['mc_backward'],
            speaker='rte_solver',
            enable_values=("mystic","montecarlo"),
            mystic =True
        )

        mc_basename = option_definition.option(
            name='mc_basename', 
            group='output', 
            helpstr='Filename for MYSTIC output.',
            documentation=documentation['mc_basename'],
            tokens=option_definition.addToken(name='Input.rte.mc.filename[FN_MC_BASENAME]',datatype=io.IOBase),
            parents=['uvspec'], 
            speaker='rte_solver',
            enable_values=("mystic","montecarlo"),
            mystic =True
        )

        mc_backward_writeallpixels = option_definition.option(
            name='mc_backward_writeallpixels', 
            group='output', 
            helpstr='Write all pixels to the output files',
            documentation=documentation['mc_backward_writeallpixels'],    
            tokens=option_definition.addSetting(name='Input.rte.mc.backward.writeallpixels', setting=1),    
            parents=['mc_backward'], 
            speaker='rte_solver',
            enable_values=("mystic","montecarlo"),
            threedmystic =True
        )

        mc_std = option_definition.option(
            name='mc_std',
            group='output',
            helpstr='Calculate MYSTIC standard deviation.', 
            documentation=documentation['mc_std'],
            tokens=option_definition.addSetting(name='Input.rte.mc.std', setting=1, default=0), 
            parents=['uvspec'], 
            speaker='rte_solver',
            enable_values=("mystic","montecarlo"),
            threedmystic =True
        )

        mc_surfaceparallel = option_definition.option(
            name='mc_surfaceparallel',
            group='output',
            helpstr='Calculate irradiance parallel to the surface.',
            documentation=documentation['mc_surfaceparallel'],
            tokens=option_definition.addSetting(name='Input.rte.mc.surfaceparallel', setting=1),
            parents=['uvspec'],
            speaker='rte_solver',
            enable_values=("mystic","montecarlo"),
            threedmystic =True,
            showInGui=False,
        )

        mc_jacobian = option_definition.option(
            name='mc_jacobian',
            group='output',
            helpstr='Calculate jacobi matrix.',
            documentation=documentation['mc_jacobian'],
            tokens=option_definition.addSetting(name='Input.rte.mc.jacobian', setting=1),
            parents=['uvspec'],
            childs=['mc_jacobian_std'],
            speaker='rte_solver',
            enable_values=("mystic","montecarlo"),
            threedmystic =False,
            islidar=False,
            showInGui=False,
                        developer=True
        )

        mc_jacobian_std = option_definition.option(
            name='mc_jacobian_std',
            group='output',
            helpstr='Calculate jacobi matrix standard deviation.',
            documentation=documentation['mc_jacobian_std'],
            tokens=option_definition.addSetting(name='Input.rte.mc.jacobian_std', setting=1),
            parents=['mc_jacobian'],
            speaker='rte_solver',
            enable_values=("mystic","montecarlo"),
            threedmystic =True,
            islidar=True,
            showInGui=False,
                        developer=True
        )

        self.options = [include,
                    quiet, verbose, write_ext_to_file, test_optical_properties,
                print_disort_info,
                data_files_path,
                output_process,
                output_file, output_format, output_user,
                output_quantity,
                heating_rate,
                write_output_as_netcdf,
                slit_function_file, spline, spline_file,
                filter_function_file,  
                pressure_out,
                zout, zout_sea,  
                mc_backward_output,
                mc_forward_output, 
                mc_basename, 
                mc_backward_writeallpixels,
                mc_std,
                mc_surfaceparallel, mc_jacobian, mc_jacobian_std,
                ]

    def __iter__(self):
        return iter(self.options)


def get_documentation():
    return get_output_documentation()

def get_output_documentation():
    return {
        'data_files_path' : r'''
    The path to the directory where all \code{uvspec} internal data files live, e.g.
    the files that are in the subdirectories of the \file{data} directory
    that comes with the \code{uvspec} distribution. 
    \fcode{
    data\_files\_path path
    }
    The default for \code{path} is \code{../data/}.
        ''',

        'include' : r'''
    Include a file into the \code{uvspec} input.  
    \fcode{ 
    include file 
    } 
    Works exactly like the C \code{\#include} or  
    the Fortran \code{INCLUDE} statements. 
        ''',

        'mc_backward_writeallpixels' : r'''
    If set, write all pixels to the output files; otherwise, only those are written 
    which are actually calculated.
        ''',

        'mc_backward_output' : r'''
    Specify quantity to be calculated using backward Monte Carlo. 
    \fcode{
    mc\_backward\_output output [unit]
    }
    So far the
    following \code{output} options have been implemented: 
    \begin{description}
    \parameter{edir} 
      direct horizontal irradiance
    \parameter{edn} 
      diffuse downward irradiance (default)
    \parameter{eup} 
      diffuse upward irradiance
        \parameter{exp} 
      diffuse irradiance in positive x direction for grid box above specified altitude at the left (lower x) face; output is written to eup.
        \parameter{exn} 
      diffuse irradiance in negative x direction for grid box above specified altitude at the left (lower x) face; output is written to edn.
        \parameter{eyp} 
      diffuse irradiance in positive y direction for grid box above specified altitude at the front (lower y) face; output is written to eup.
        \parameter{eyn} 
      diffuse irradiance in negative y direction for grid box above specified altitude at the front (lower y) face; output is written to edn.
    \parameter{act} 
      actinic flux
    \parameter{abs} 
      absorption
        \parameter{emis} 
          emission
        \parameter{heat} 
          heating rates, that is absorption + emission
    \end{description}
    For \code{abs}, \code{emis}, \code{heat} an optional argument \code{W\_per\_m2\_and\_dz} (default), 
    \code{W\_per\_m3}, or \code{K\_per\_day} may be specified which converts the result e.g. to heating rates. 
        ''',

        'mc_forward_output' : r'''
    Specify quantity to be calculated using forward Monte Carlo.
    Forward output is calculated for each grid box; may need a considerable amount of memory,
        depending on the 3D cloud grid. Only meaningful with \code{rte\_solver mystic}.
    This option is only available for forward calculations. 
    For backward please use \code{mc\_backward\_output}.
    \fcode{
    mc\_forward\_output quantity [unit]
    }
    \code{quantity} can be one of the following:
    \begin{description}
    \parameter{absorption}
    Calculate MYSTIC absorption and write it to file \code{mc\_basename}.abs.spc 
        in the following format:
        \fcode{
        lambda ix iy iz absorption
        }
    For backward please use \code{mc\_backward\_output abs}.
    \parameter{actinic}
    Calculate MYSTIC actinic flux by dividing the absorbed energy by the 
        absorption coeffcient; this method is much better than the traditional photon counting 
        which usually comprises spikes (because in the latter method each photon is weighted with 
        $1/\cos(\theta)$ which may be a very large number);  
        For backward please use \code{mc\_backward\_output act}. 
    \parameter{emission}
    Calculate MYSTIC emission. Note that emission is calculated directly without tracing 
        any photons which makes this option very fast. Changing \code{mc\_photons} will therefore 
        not affect the result. Only meaningful with \code{source thermal}. 
    For backward please use \code{mc\_backward\_output emis}.
    \parameter{heating}
    Calculate MYSTIC heating rates. This is identical to \code{absorption} for \code{source solar}. 
    For \code{source thermal}, however, the emission of a photon in a grid box is counted as cooling. 
    For backward please use \code{mc\_backward\_output heat}.
    \end{description}
        The optional argument \code{unit} may be one of \code{W\_per\_m2\_and\_dz} (default), \code{W\_per\_m3}, or \code{K\_per\_day}.
    \code{unit} may be specified to convert the result e.g. to heating rates. 
        ''',

        'mc_basename' : r'''
    Filename for MYSTIC 3D output (default: mc).
    \fcode{
    mc\_basename basename
    }
        ''',

    'mc_jacobian' : r'''
    Calculate Jacobians. i.e. the derivative of radiance/irradiance with respect to absorption/scattering optical thickness. The derivatives are calculated for each vertical layer. 

        The option is implemented only for \code{rte\_solver montecarlo} in backward tracing mode \code{mc\_backward}.

        The derivatives are calculated separately for all species (i.e. molecules, aerosols and clouds). So far the option is only included for monochromatic simulations. 

        The result is provided in the output file \code{basename.jac}, which includes the following columns: \\

        x \hspace{1ex}  y \hspace{1ex}  z$_i$ \hspace{1ex}  $\frac{\partial{I}}{\partial{\tau_{s,mol,i}}}$ \hspace{1ex} $\frac{\partial{I}}{\partial{\tau_{a,mol,i}}}$ \hspace{1ex}  $\frac{\partial{I}}{\partial{\tau_{s,aer,i}}}$ \hspace{1ex}  $\frac{\partial{I}}{\partial{\tau_{a,aer,i}}}$
        
        ''',
    'mc_jacobian_std' : r'''
    Calculate jacobi matrix standard deviation. Only meaningful with 
    \code{mc\_jacobian}. 
        ''',

        'mc_std' : r'''
    Calculate standard deviation of the average.
        ''',

        'mc_surfaceparallel' : r'''
    Calculate irradiance parallel to the surface instead of horizonal irradiance. This
    option is obviously only interesting for topograpy and only for calculations at
    the surface. For other levels the option is ignored. 
        ''',

        'filter_function_file' : r'''
    If specified, the calculated spectrum is multiplied with a filter function 
    defined in file.
    \fcode{
    filter\_function\_file file [normalize]
    }
    The file must contain two columns. 
    Column 1 is the wavelength, in nm. Column 2 is the corresponding filter 
    function value. Comments start with \code{\#}. Empty lines are ignored.
    In combination with \code{output\_process sum} or \code{output\_process integrate} this option 
    is useful e.g. to calculate weighted irradiances or actinic fluxes or 
    to simulate broadband or satellite observations.

    If the optional second argument \code{normalize} is specified, 
    the integral of the filter function over wavelength is normalized such
    that \code{output\_process integrate} gives radiative properties per wavelength, averaged over
    the filter function.
            ''',

        'slit_function_file' : r'''
    If specified, the calculated spectrum is convolved with the function found
    in the \file{slit\_function\_file}. 
    \fcode{
    slit\_function\_file file 
    }
    The file must contain two columns. Column 1 is the wavelength, in nm, and relative to the
    center wavelength. Column 2 is the corresponding slit function value. It must
    be unity at the maximum. The wavelength steps in the slit function file must
    be equidistant. Comments start with \code{\#}. Empty lines are ignored. Please 
    note that prior to convolution the spectrum is interpolated to the wavelength
    steps of the slit function. For this reason, make sure that the resolution
    of the slit function is high enough even if the slit function is e.g. a 
    simple triangle which could in principle be described with 3 grid points. For an 
    example see \file{examples/TRI\_SLIT.DAT} and the \code{make\_slitfunction} tool.
        ''',

        'spline' : r'''
    \fcode{
    spline lambda\_0 lambda\_1 lambda\_step
    }
    Spline interpolate the calculated spectrum between wavelengths \code{lambda\_0} 
    and \code{lambda\_1} in steps of \code{lambda\_step},
    in nm. Specified as e.g.
    \fcode{
    spline 290. 365. 0.5
    }
    Here, the calculated spectrum is interpolated to wavelengths 290.0, 290.5, 291.0, 
    ..., 364.5, 365.0. For interpolation to arbitrary wavelengths use \code{spline\_file}.
    The specified wavelength interval must be within the one specified by \code{wavelength}.
        ''',

        'spline_file' : r'''
    Spline interpolate to arbitrary wavelengths, in nm, given as a single column in file 
    \file{spline\_file}. 
    \fcode{
    spline\_file file
    }
    The specified wavelengths must be within the range specified 
    by \code{wavelength}. Comments start with \code{\#}. Empty lines are ignored.
        ''',

        'test_optical_properties' : r'''
    Write optical properties to file. This is made for test suite. 
    rte_solver will be switched to null_solver.
        ''',

        'output_quantity' : r'''
    Convert radiances / irradiances to equivalent output quantity.
    \fcode{
    output\_quantity quantity
    }
    where \code{quantity} can be one of the following:
    \begin{description}
    \parameter{brightness}
    Convert radiances / irradiances to equivalent brightness temperatures.
    \parameter{reflectivity}
    Calculate transmission / reflectivity instead of absolute quantities. 
    For irradiances / actinic fluxes the transmission T is defined as 
    \begin{equation}
       T = {{E}\over{E_0 \cos{\theta}}}
    \end{equation}
      where $E$ is the irradiance / actinic flux, $E_0$ is the extraterrestrial flux,
      and $\theta$ is the solar zenith angle.
    The reflectivity R is defined as 
    \begin{equation}
       R = {{\pi \cdot L}\over{E_0 \cos{\theta}}}
    \end{equation}
     where $L$ is the radiance, $E_0$ is the extraterrestrial flux,
     and $\theta$ is the solar zenith angle.
    Obviously, reflectivities do not depend on Sun-Earth distance. Please 
    note the difference to \code{transmittance}.
    \parameter{transmittance}
    Calculate transmittance / reflectance instead of absolute quantities. 
    That is, set the extraterrestrial irradiance to 1 and do not correct for Sun-Earth distance:
    \begin{equation}
       T = \frac{E}{E_0}
    \end{equation}
      where $E$ is the irradiance / actinic flux / radiance and $E_0$ is the extraterrestrial flux.
    Please note the difference to \code{reflectivity}.
        \end{description}
        ''',

        'output_user' : r'''
    User defined output. 
    %This option is case sensitive.
    Here the user may specify the columns desired for output. 
    \fcode{
    output\_user format
    }
    where \code{format} is one or more of the following.
    \begin{description}
    \parameter{lambda} 
    Wavelength in nm.
    \parameter{wavenumber} 
    Wave number in cm$^{-1}$.
    \parameter{sza} 
    solar zenith angle
    \parameter{zout} 
    Output altitude in km.
    \parameter{edir, eglo, edn, eup, enet, esum}
    The direct, global, diffuse downward, and diffuse upward irradiance.
    Net is global - upward, sum is global + upward. 
    \parameter{uu} 
    Radiances:  uu(umu(0),phi(0)) ... uu(umu(0),phi(m)) ... uu(umu(n),phi(0)) ... uu(umu(n),phi(m))
    \parameter{fdir, fglo, fdn, fup, f}
    The direct, global, diffuse downward, diffuse upward, and total actinic flux.
    \parameter{uavgdir, uavgglo, uavgdn, uavgup, uavg} 
    The Direct, global, diffuse downward, diffuse upward, and total diffuse 
    mean intensity (= actinic flux / 4$\pi$).
    \parameter{spher\_alb}
    Spherical albedo of the complete atmosphere.
    \parameter{albedo} 
    Albedo.
    \parameter{heat} 
    Heating rate in K/day.
    \end{description}

    It is also possible to gain some information about the atmosphere and the clouds:
    \begin{description}
    \parameter{p} 
    pressure [hPa], ,
    \parameter{T, T\_d}
    temperature [K], dewpoint temperature [K]
    \parameter{T\_sur} 
    surface temperature [K]
    \parameter{theta } 
    potential temperature [K]
    \parameter{theta\_e} 
    equivalent potential temperature [K]
    \parameter{n\_xxx} 
    number density of the gas xxx [cm$^{-3}$]
    \parameter{rho\_xxx} 
    mass density of the gas xxx [kg/m$^3$]
    \parameter{mmr\_xxx} 
    mass mixing ratio of the gas xxx [kg/kg]
    \parameter{vmr\_xxx} 
    volume mixing ratio of the gas xxx [m$^3$/m$^3$]
    \parameter{rh} 
    relative humidity over water [percent]
    \parameter{rh\_ice} 
    relative humidity over ice   [percent]
    \parameter{c\_p} 
    specific heat capacity of the air (humidity and temperature dependent)
    \parameter{CLWC} 
    cloud liquid water content   [kg/kg]
    \parameter{CLWD} 
    cloud liquid water density   [g/m$^3$]
    \parameter{CIWC} 
    cloud ice water content      [kg/kg]
    \parameter{CIWD} 
    cloud ice water density      [g/m$^3$]
    \parameter{TCC} 
    total cloud cover            [0-1]
    \end{description}
    where \code{xxx} is one of 
    AIR, O3, O2, H2O, CO2, NO2, BRO, OCLO, HCHO, or O4.

    Default output is  
    \fcode{
    output\_user lambda lambda, edir, edn, eup, uavgdir, uavgdn, uavgup 
    }
    for \code{fdisort1}, \code{sdisort}, and \code{spsdisort}, whereas the default for 
    \code{twostr} is 
    \fcode{
    output\_user lambda, edir, edn, eup, uavg.
    }
    The lines containing radiances and the output of \code{rte\_solver polradtran} are not affected. 
        ''',

        'pressure_out' : r''' 
    Specify the output levels in pressure coordinates. The syntax is   
    \fcode{ 
    pressure\_out p1 p2 ... 
    }
    where '\code{p1 p2 ...}' are the output levels in hPa.  
    The pressure output levels must be sorted in decreasing order.  
    Output pressure levels must be within the range defined in the 
    \code{atmosphere\_file}. You can also use \code{toa} for  
    \code{top of atmosphere} and \code{sur} for \code{surface altitude} and \code{cpt}  
    for \code{cold point tropopause}. 
        ''',

        'print_disort_info' : r'''
    Specify one or more integers between 1 and 7.
    \fcode{
    print\_disort\_info  value
    }
    Print various disort input and output in disorts own format. See 
    \file{libsrc\_f/DISORT2.doc} for more information.
    \strong{Warning:} Produces a lot of output. 
        ''',

        'heating_rate' : r'''
    Calculation of heating rates. Output is only provided at altitudes specified by \code{zout}.
    To get heating rate profiles a number of altitudes must thus be specified by \code{zout}.
    Heating rates is the change of temperature with time in units of K/day.
    For spectral calculations the default output is a matrix: 
    \begin{Verbatim}
        0.0        zout1          zout2 ...
        lambda1    heating\_rates  ...
        lambda2      .
           .         .
           .         .
    \end{Verbatim}
    For integrated calculations (\code{output\_process sum} or \code{output\_process integrate}) the default output 
    is in two columns with column 1 being the altitude and column 2 the heating rates.
    The output of \code{heating\_rate} can also be specified with the \code{output\_user} option.
    Note that heating rates are only well-behaved up to altitudes
    for which the respective correlated-k options are valid. E.g. about 60 km for 
    \code{fu} and about 80 km for \code{kato}, \code{kato2}, \code{kato2.96}, and \code{lowtran}.
    Attention: For spectral calculations, the extraterrestrial spectrum is assumed to be in 
    mW/(m2 nm).

    Three different methods are implemented to calculate the heating rate, which can be selected
    with an optional keyword:
    \fcode{
    heating\_rate [method]
    }
    where \code{method} may be either \code{layer\_cd} (heating rates are derived from centered 
    differences of the flux (the default method), \code{local} (heating rates are derived 
    from the actinic flux), or \code{layer\_fd} (heating rates are derived from forward differences of the flux 
        over one layer. Attention: \code{heating\_rate local} introduces new levels into the profile which  
    slightly affects the model output.
    With \code{layer\_fd}, the output is 
    not representative for a \emph{level}, but for the \emph{layer} from the z-level of the line in the output file, 
    where it is written, up to next \emph{output level} above!
        ''',

        'output_process' : r'''
    Decide how the output from \code{uvspec} is processed:
    \fcode{
    output\_process type
    }
    where type is one of
    \begin{description}
    \parameter{sum} 
    Sum output over wavelength. Useful in combination with the 
    \code{mol\_abs\_param} option (\code{kato}, \code{kato2}, \code{kato2.96}, 
    \code{Fu}, \code{avhrr\_kratz}).
        In case of \code{mol\_abs\_param reptran}, the units are automatically converted to \code{per\_band} before summation. 
    \parameter{integrate} 
    Integrate output over wavelength for solar and over wavenumber for thermal simulations. 
    Useful for spectral calculations, \code{mol\_abs\_param lowtran} and \code{reptran}.
    \parameter{per\_nm } 
    Output is given in W/(m$^2$ nm) or mW/(m$^2$ nm) (W or mW is determined by the extraterrestrial spectrum.)
    \parameter{per\_cm-1 } 
    Output is given in W/(m$^2$ cm$^{-1})$ or mW/(m$^2$ cm$^{-1}$).
    \parameter{per\_band } 
    Output is given in W/m$^2$ or mW/m$^2$ per correlated-k or reptran band. (This option can not be used for 
    spectral calculations and \code{mol\_abs\_param} LOWTRAN in the solar range.)
    \parameter{none} 
    No processing - output spectral information (default).
    \end{description}
        ''',

        'quiet' : r'''
    If specified, informative messages are turned off. See also \code{verbose}.
        ''',

        'verbose' : r'''
    If specified abundances of informative messages are output to stderr. To make
    use of this information, you may want to write the standard \code{uvspec} output to 
    one file and the diagnostic messages to another. To do so, try
    \code{(./uvspec < uvspec.inp > uvspec.out) >\& verbose.txt} (depending on your shell you 
    might need a slightly different syntax). The irradiances 
    and radiances will be written to \file{uvspec.out} while all 
    diagnostic messages go into \file{verbose.txt}. See also \code{quiet}.
        ''',

        'output_file'    : r'''
    \fcode {             
        output\_file filename
    }
    uvspec output will be written to \code{filename}. The format can be set by the option 
    \code{output\_format} (default is ascii).  
        ''',

        'write_output_as_netcdf' : r'''
     Option only works with 
    \code{mc\_lidar} and \code{mc\_radar}. Output is written into netcdf file rather than dozens 
    of ASCII files. Saves disk space. 
        ''',

        'zout'    : r'''
    This option is used to specify the output altitudes in km \emph{above surface altitude}. 
    One or more altitudes may be specified in increasing magnitude. 
    \fcode{
    zout 0 1 2 3 4 5 ...
    }
    Output altitudes must be within the range defined in the
    \code{atmosphere\_file}. Note that \code{zout} does not restructure the
    atmosphere model. Hence, if you specify \code{zout 0.730} and have your
    atmosphere model in \code{atmosphere\_file} go all the way down to 
    sea level, i.e. 0.0km., output is presented at 0.730km and calculations
    performed with an atmosphere between 0.0 and 0.730 km (and above of course). 
    If you want calculations done for e.g. an elevated site you have to restructure
    the atmosphere model and make sure it stops at the appropriate altitude.
    This you may either due by editing the atmosphere file or by using 
    \code{altitude}. Note that for \code{rte\_solver polradtran} the atmosphere file must 
    contain the altitudes specified by \code{zout}. You can also use \code{toa} for 
    \code{top of atmosphere} and \code{sur} for \code{surface altitude} and \code{cpt} 
    for \code{cold point tropopause}.

    Instead of specifying the altitudes in km, it is also possible to use keywords as 
    argument for this option. Possible keywords are \code{atm\_levels}, 
    \code{all\_levels}, \code{model\_levels}, \code{model\_layers}, and \code{model\_levels\_and\_layers}. 
    For \code{atm\_levels}, all levels from the \code{atmosphere\_file} are used as output levels.  
    For \code{all\_levels}, all levels (including levels from \code{atmosphere\_file}, 
    \code{mol\_file}, cloud files, \code{altitude} options) are used as output levels. 
    For \code{model\_levels}, \code{model\_layers}, \code{model\_levels\_and\_layers} 
    the levels, layers, or both from the \code{ECMWF\_atmosphere\_file} are used as output level.
    Usage e.g.:
    \fcode{
    zout model\_levels [nlev\_max]
    }
    With the optional argument \code{nlev\_max} the user may specify the number of zout layers 
    from the ground.
        ''',

        'zout_sea'    : r'''
    like zout, but \emph{above sea surface}
        ''',

        'output_format'    : r'''
    \fcode {
        output\_format format
    }
    where \code{format} is either \code{ascii} (default) or \code{flexstor}. 
    Note that \code{flexstor} does not currently work when \code{umu} and/or \code{phi} is specified.
    There is also the
    possibility to write uvspec simulation results to an existing netCDF file. For
    that \code{format} must be \code{netCDF} and the option \code{output\_file} must be given
    and point to a file that contains a lat/lon/time grid.
    If \code{format} is set to \code{sat picture} then \code{output\_file} must be given and point
    to a NetCDF-File that contains a pixel x/pixel y/time grid.
        ''',

        'write_ext_to_file'    : r'''
    Option needed by ARLEM when calling uvspec. Writes extinction
    information to file instead of to output.
        ''',

    }
