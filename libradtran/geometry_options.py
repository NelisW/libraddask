"""--------------------------------------------------------------------
 * $Id: geometry_options.py 3519 2019-12-10 10:53:35Z Claudia.Emde $
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


class setup_geometry_group():

    group_name = 'Geometry'

    def __init__(self):
        documentation = get_geometry_documentation()

        sza = option_definition.option(
            name='sza',
            group='geometry',
            helpstr='Solar zenith angle',
            documentation=documentation['sza'],
            tokens = [option_definition.addToken(name='Input.atm.sza', datatype=float, default='NOT_DEFINED_FLOAT', valid_range=[0.,180.]), 
                      option_definition.addSetting(name='Input.atm.sza_source', setting='SZA_DIRECT_INPUT', default='SZA_BY_TIME_AND_LOCATION')], 
            parents=['uvspec'],
            non_parents=['sza_file','ECHAM_sza'],
        )

        sza_file = option_definition.option(
            name='sza_file', 
            group='geometry',
            helpstr='Location of solar zenith angle file.', 
            documentation=documentation['sza_file'], 
            tokens = [option_definition.addToken(name='Input.filename[FN_SZA]', datatype=io.IOBase),
                      option_definition.addSetting(name='Input.atm.sza_source', setting='SZA_FROM_SZA_FILE', default='SZA_BY_TIME_AND_LOCATION')], 
            parents=['uvspec'],
            non_parents=['sza', 'ECHAM_sza'],
        )

        phi0 = option_definition.option(
            name='phi0',
            group='geometry',
            helpstr='Solar azimuth angle',
            documentation=documentation['phi0'], 
            tokens=option_definition.addToken(name='Input.atm.phi0', datatype=float, valid_range=[-360.,360.]), 
            parents=['uvspec'],
        )

        phi = option_definition.option(
            name='phi',
            group='geometry',
            helpstr='User output azimuth angles.',
            documentation=documentation['phi'],
            tokens = [option_definition.addToken(name='Input.rte.phi', datatype=option_definition.SignedFloats),
                      option_definition.addSetting(name='Input.rte.nphi', setting='ntokens'),
                      option_definition.addSetting(name='Input.rte.maxphi', setting='ntokens')],
            parents=['uvspec'],
        )

        umu = option_definition.option(
            name='umu',
            group='geometry',
            helpstr='Cosine of user output polar angles.',
            documentation=documentation['umu'],
            tokens = [option_definition.addToken(name='Input.rte.umu', datatype=option_definition.SignedFloats),
                      option_definition.addSetting(name='Input.rte.numu', setting='ntokens'),
                      option_definition.addSetting(name='Input.rte.maxumu', setting='ntokens')],
            parents=['uvspec'],
            non_parents=[],
            childs=[]
        )

        mc_bw_umu_file = option_definition.option( 
            name='mc_bw_umu_file',
            group='geometry',
            helpstr='Define a different umu and phi for each pixel.', 
            documentation=documentation['mc_bw_umu_file'],
            gui_inputs=(GUI_definition.FileInput(name='filename'),),
            tokens=[option_definition.addToken(name="Input.rte.mc.filename[FN_MC_UMU]", datatype=io.IOBase),
                    option_definition.addSetting(name='Input.rte.mc.allocate_umu_and_phi',setting=True)],
            parents=['mc_backward'],
            speaker='rte_solver',
            enable_values=("mystic","montecarlo"),
            threedmystic =True
        )

        earth_radius = option_definition.option(
            name='earth_radius',
            group='geometry',
            helpstr='Radius of the Earth in km',
            documentation=documentation['earth_radius'], 
            tokens=option_definition.addToken(name='Input.r_earth', datatype=float, default=6370.0),
            parents=['uvspec'], 
        )

        latitude = option_definition.option(
            name='latitude',
            group='geometry',
            helpstr='Specify the latitude of the location to simulate.',
            documentation=documentation['latitude'],
            gui_inputs=( GUI_definition.ListInput(name='hemisphere', valid_range=['n','s']), GUI_definition.FloatInput(name='degrees',valid_range=[0,180]), 
                GUI_definition.FloatInput(name='minutes',valid_range=[0,60],optional=True), GUI_definition.FloatInput(name='seconds',valid_range=[0,60], optional=True),),
            tokens = [option_definition.addLogical(name='Input.lat_signum',logicals=['n','s'],setting='LATITUDE_SIGNUM_'),
                      option_definition.addToken(name='Input.lat_degrees', datatype=float),
                      option_definition.addToken(name='Input.lat_minutes', datatype=float, optional=True),
                      option_definition.addToken(name='Input.lat_seconds', datatype=float, optional=True)],
            parents=['uvspec'], 
            extra_dependencies=['latitude','longitude'],
            showInGui=False
        )

        longitude = option_definition.option(
            name='longitude',
            group='geometry',
            helpstr='Specify the longitude of the location to simulate.',
            documentation=documentation['longitude'],
            gui_inputs=( GUI_definition.ListInput(name='hemisphere', valid_range=['e','w']), GUI_definition.FloatInput(name='degrees',valid_range=[0,180]), 
                GUI_definition.FloatInput(name='minutes',valid_range=[0,60],optional=True), GUI_definition.FloatInput(name='seconds',valid_range=[0,60], optional=True),),
            tokens = [option_definition.addLogical(name='Input.lon_signum',logicals=['e','w'],setting='LONGITUDE_SIGNUM_'),
                      option_definition.addToken(name='Input.lon_degrees', datatype=float),
                      option_definition.addToken(name='Input.lon_minutes', datatype=float, optional=True),
                      option_definition.addToken(name='Input.lon_seconds', datatype=float, optional=True)],
            parents=['uvspec'],
            extra_dependencies=['latitude','longitude'],
            showInGui=False,
        )

        day_of_year = option_definition.option(
            name='day_of_year',
            group='geometry',
            helpstr='Correction fo sun-earth distance',
            documentation=documentation['day_of_year'], 
            tokens=option_definition.addToken(name='Input.UTC.tm_yday', datatype=int, default='NOT_DEFINED_INTEGER'),
            parents=['uvspec'], 
        )

        mc_bcond = option_definition.option(
            name='mc_bcond',
            group='geometry',
            helpstr='Define MYSTIC boundary conditions.', 
            documentation=documentation['mc_bcond'], 
            tokens=option_definition.addLogical(name='Input.rte.mc.bcond', logicals=['periodic','mirror','absorb'], setting='MCBCOND_'), 
            parents=['uvspec'], 
            speaker='rte_solver',
            enable_values=("mystic","montecarlo"),
            mystic =True,
            developer=True,
            showInGui=False,
        )

        mc_sample_grid = option_definition.option(
            name='mc_sample_grid',
            group='geometry',
            helpstr='Sample grid size (Nx Ny [dx dy]) for MYSTIC.', 
            documentation=documentation['mc_sample_grid'],
            tokens= [option_definition.addToken(name="Input.rte.mc.Nx_sample", datatype=int),
                     option_definition.addToken(name="Input.rte.mc.Ny_sample", datatype=int),
                     option_definition.addToken(name="Input.rte.mc.dx_sample", datatype=float, optional=True),
                     option_definition.addToken(name="Input.rte.mc.dy_sample", datatype=float, optional=True)],
            parents=['uvspec'],
            speaker='rte_solver',
            enable_values=("mystic","montecarlo"),
            threedmystic =True,
        )

        mc_spherical = option_definition.option(
            name='mc_spherical',
            group='geometry',
            helpstr='Use MYSTIC spherical geometry.',
            documentation=documentation['mc_spherical'],
            tokens=[option_definition.addLogical(name='id', logicals=['1D','3D'], setting='DIM_'),
                    option_definition.addSetting(name='Input.rte.mc.spherical[id]', setting=1)],
            gui_inputs=(GUI_definition.ListInput(name='dimension', valid_range=['1d']),),
            parents='mc_backward',
            childs=['mc_refraction'],
            speaker='rte_solver',
            enable_values=("mystic","montecarlo"),
            mystic =True,
            #showInGui=False,
        )

        mc_spherical3D_scene = option_definition.option(
            name='mc_spherical3D_scene',
            group='geometry',
            helpstr='Define scene on globe.', 
            documentation=documentation['mc_spherical3D_scene'],
            tokens=[option_definition.addSetting(name='Input.rte.mc.spherical3D_scene', setting=1),
                    option_definition.addToken(name='Input.rte.mc.spherical3D_scene_lon_min', datatype=float, gui_name='LON_W'),
                    option_definition.addToken(name='Input.rte.mc.spherical3D_scene_lon_max', datatype=float, gui_name='LON_E'),
                    option_definition.addToken(name='Input.rte.mc.spherical3D_scene_lat_min', datatype=float, gui_name='LAT_S'),
                    option_definition.addToken(name='Input.rte.mc.spherical3D_scene_lat_max', datatype=float, gui_name='LAT_N')],
            speaker='mc_spherical',
            enable_values=("3d",),
            islidar =True,
            showInGui=False,
        )

        mc_satellite_view = option_definition.option(
            name='mc_satellite_view',
            group='geometry',
            helpstr='Simulate a satellite image.', 
            documentation=documentation['mc_satellite_view'],
            tokens=[option_definition.addSetting(name='Input.rte.mc.pan_angles', setting='calloc (4, sizeof(double))'),
                    option_definition.addToken(name='Input.rte.mc.pan_angles[0]', datatype=float),
                    option_definition.addToken(name='Input.rte.mc.pan_angles[1]', datatype=float),
                    option_definition.addToken(name='Input.rte.mc.pan_angles[2]', datatype=float),
                    option_definition.addToken(name='Input.rte.mc.pan_angles[3]', datatype=float),
                    option_definition.addSetting(name='Input.rte.mc.panorama', setting='PAN_MODE_SATELLITE'),
                    option_definition.addSetting(name='Input.rte.mc.allocate_umu_and_phi', setting=True)],
            parents=['mc_backward'],
            non_parents=['mc_panorama_view'],
            speaker='rte_solver',
            enable_values=("mystic","montecarlo"),
            threedmystic =True,
            islidar=True,
            showInGui=False,
        )

        mc_satellite_position = option_definition.option(
            name='mc_satellite_position',
            group='geometry',
            helpstr='Define the position of the satellite.', 
            documentation=documentation['mc_satellite_position'],
            tokens=[option_definition.addSetting(name='Input.rte.mc.sensorposition', setting='MC_SENSORPOSITION_SPHERICAL'),
                option_definition.addToken(name='Input.rte.mc.sensorposition_lon', datatype=float),
                option_definition.addToken(name='Input.rte.mc.sensorposition_lat', datatype=float),
                option_definition.addToken(name='Input.rte.mc.sensorposition_alt', datatype=float)],
            non_parents=['mc_sensor_position'],
            islidar=True,
            showInGui=False,
        )

        mc_sun_position = option_definition.option(
            name='mc_sun_position',
            group='geometry',
            helpstr='Sun position', 
            documentation=documentation['mc_sun_position'],
            tokens=[option_definition.addSetting(name='Input.atm.sza_source', setting='SZA_DIRECT_INPUT'),
                    option_definition.addToken(name='Input.atm.phi0_spher', datatype=float),
                    option_definition.addToken(name='Input.atm.sza_spher', datatype=float)],
            islidar=True,
            speaker='mc_sperical',
            enable_values=('3D',),
            showInGui=False,
        )

        mc_sensordirection = option_definition.option(
            name='mc_sensordirection',
            group='geometry',
            helpstr='Define viewing direction of an irradiance sensor.', 
            documentation=documentation['mc_sensordirection'],
            tokens=[option_definition.addToken(name="Input.rte.mc.sensordirection_dx", datatype=float, gui_name='dx'),
                    option_definition.addToken(name="Input.rte.mc.sensordirection_dy", datatype=float, gui_name='dy'),
                    option_definition.addToken(name="Input.rte.mc.sensordirection_dz", datatype=float, gui_name='dz'),
                    option_definition.addSetting(name="Input.rte.mc.sensordirection", setting=1)],
            parents=['mc_backward'],
            speaker='rte_solver',
            enable_values=("mystic","montecarlo"),
            mystic = True,
        )

        mc_sensorposition = option_definition.option(
            name='mc_sensorposition',
            group='geometry',
            helpstr='Define the position of a sensor.', 
            documentation=documentation['mc_sensorposition'],
            tokens=[option_definition.addToken(name="Input.rte.mc.sensorposition_x", datatype=float, gui_name='x'),
                    option_definition.addToken(name="Input.rte.mc.sensorposition_y", datatype=float, gui_name='y'),
                    option_definition.addToken(name="Input.rte.mc.sensorposition_z", datatype=float, gui_name='z'),
                    option_definition.addSetting(name="Input.rte.mc.sensorposition", setting='MC_SENSORPOSITION_CARTESIAN')],
            parents=['mc_backward'],
            non_parents=['mc_satellite_position'],
            speaker='rte_solver',
            enable_values=("mystic","montecarlo"),
            threedmystic =True,
        )

        mc_reference_to_nn = option_definition.option(
            name='mc_reference_to_nn', 
            group='geometry',
            helpstr='The sampled pixels correspond to he surface pixels.',
            documentation=documentation['mc_reference_to_nn'],
            tokens=[option_definition.addSetting(name='Input.rte.mc.reference_to_NN', setting=1, default=0)], 
            parents=['mc_backward'],  
            speaker='rte_solver',
            enable_values=("mystic","montecarlo"),
            threedmystic =True
        )

        mc_panorama_forward = option_definition.option(
            name='mc_panorama_forward',
            group='geometry',
            helpstr='Simulate a 4pi panorama with forward MYSTIC.', 
            documentation=documentation['mc_panorama_forward'],
            tokens=[option_definition.addSetting(name='Input.rte.mc.pan_angles', setting='calloc (4, sizeof(double))'),
                    option_definition.addSetting(name='Input.rte.mc.pan_angles[0]', setting=0.0),
                    option_definition.addSetting(name='Input.rte.mc.pan_angles[1]', setting=360.0),
                    option_definition.addSetting(name='Input.rte.mc.pan_angles[2]', setting=0.0),
                    option_definition.addSetting(name='Input.rte.mc.pan_angles[3]', setting=180.0),
                    option_definition.addSetting(name='Input.rte.mc.panorama_forward', setting=1),
                    option_definition.addSetting(name='Input.rte.mc.allocate_umu_and_phi', setting=True)],
            non_parents=['mc_satellite_view','mc_backward','mc_escape'],
            speaker='rte_solver',
            enable_values=("mystic","montecarlo"),
            threedmystic =True
        )

        mc_panorama_view = option_definition.option(
            name='mc_panorama_view',
            group='geometry',
            helpstr='Simulate a panorama.', 
            documentation=documentation['mc_panorama_view'],
            tokens=[option_definition.addSetting(name='Input.rte.mc.pan_angles', setting='calloc (4, sizeof(double))'),
                    option_definition.addToken(name='Input.rte.mc.pan_angles[0]', datatype=float),
                    option_definition.addToken(name='Input.rte.mc.pan_angles[1]', datatype=float),
                    option_definition.addToken(name='Input.rte.mc.pan_angles[2]', datatype=float),
                    option_definition.addToken(name='Input.rte.mc.pan_angles[3]', datatype=float),
                    option_definition.addSetting(name='Input.rte.mc.panorama', setting='PAN_MODE_CAMERA'),
                    option_definition.addSetting(name='Input.rte.mc.allocate_umu_and_phi', setting=True)],
            parents=['mc_backward'],
            non_parents=['mc_satellite_view'],
            speaker='rte_solver',
            enable_values=("mystic","montecarlo"),
            threedmystic =True,
            childs=['mc_panorama_alignment','mc_panorama']
        )

        mc_panorama_alignment = option_definition.option(
            name='mc_panorama_alignment',
            group='geometry',
            helpstr='Panorama alignment.', 
            documentation=documentation['mc_panorama_alignment'],
            tokens=option_definition.addLogical(name='Input.rte.mc.pan_alignment',logicals=['mu','sun','zenith'],setting='PAN_ALIGNMENT_'),
            threedmystic=True,
            parents=['mc_panorama_view'],
        )           

        mc_panorama = option_definition.option(
            name='mc_panorama',
            group='geometry',
            helpstr='Panorama settings.', 
            documentation=documentation['mc_panorama'],
            tokens=[option_definition.addLogical(name='id',logicals=['distr_photons_over_pixel','no_pixel','quicklook','with_direct_rad','weight_with_cos','circumsolar_var_red'],setting='PAN_FLAG_'),
                    option_definition.addSetting(name='Input.rte.mc.pan[id]',setting=1)],
            parents=['mc_panorama_view'],
            threedmystic=True,
            non_unique=True,
        )

        mc_blitz_position = option_definition.option(
            name='mc_blitz_position',
            group='geometry',
            helpstr='Define lightning as source.', 
            documentation=documentation['mc_blitz_position'],
            tokens=[option_definition.addSetting(name='Input.rte.mc.blitz_position', setting='calloc (6, sizeof(double))'),
                    option_definition.addToken(name='Input.rte.mc.blitz_position[0]', datatype=float),
                    option_definition.addToken(name='Input.rte.mc.blitz_position[1]', datatype=float),
                    option_definition.addToken(name='Input.rte.mc.blitz_position[2]', datatype=float),
                    option_definition.addToken(name='Input.rte.mc.blitz_position[3]', datatype=float),
                    option_definition.addToken(name='Input.rte.mc.blitz_position[4]', datatype=float),
                    option_definition.addToken(name='Input.rte.mc.blitz_position[5]', datatype=float),
                    option_definition.addSetting(name='Input.source', setting='SRC_BLITZ'),
                    option_definition.addSetting(name='Input.rte.mc.allocate_umu_and_phi', setting=True)],
            parents=['mc_forward'],
            speaker='rte_solver',
            enable_values=("mystic","montecarlo"),
            threedmystic =True,
            showInGui=False,
        )

        self.options = [sza, sza_file, phi0, phi, umu, mc_bw_umu_file,
                earth_radius,
                latitude,longitude,
# missing time, 
                day_of_year,
                mc_bcond, mc_sample_grid, mc_spherical, mc_spherical3D_scene, 
                mc_satellite_view, mc_satellite_position, mc_sun_position,
                mc_sensordirection, mc_sensorposition, 
                mc_reference_to_nn, mc_panorama_view, 
                mc_panorama_alignment, mc_panorama, mc_panorama_forward,
                mc_blitz_position,
                ]

    def __iter__(self):
        return iter(self.options)

def get_documentation():
    return get_geometry_documentation()

def get_geometry_documentation():
    return {
        'day_of_year' : r'''
    Integer, to correct the calculated radiation quantities for the 
    Sun-Earth distance for the specified Julian day (1-365). 
    \fcode{
       day_of_year value
    }
    If not specified, the Earth-Sun distance is 1 AU (i.e. equinox
    distance), that is, no correction is applied to the extraterrestrial
    irradiance \code{source solar file}. Alternatively \code{time} may be used
    for that purpose.
        ''',
        'earth_radius' : r'''
    Specify the earth radius in km. 
    \fcode{
        earth_radius value
    }
    This is needed by all solvers in spherical geometry\ifmystic{, e.g. \code{mystic}
    in combination with option \code{mc_spherical}}. The default value is 6370 km.
        ''',
        'latitude' : r'''
    This option can be used to specify the latitude of the location to simulate.
    (This option only has an effects, if \code{longitude} is specified, too.)
    \fcode{
    latitude N/S deg [min] [sec]
    }
    where \code{N} and \code{S} stands for the northern and southern hemisphere, respectively.
    \code{deg min sec} is the position in degrees, arc minutes, and arc seconds.
    \code{deg} might also be a float number. \code{min} and \code{sec} may be obmitted.
    The \code{latitude} information will be used for the following: 
    
    %\code{latitude} in combination with \code{longitude}, \code{time}, and any \code{map}-option
    %is used to select the location where to read the input data.
    
    \code{latitude} in combination with \code{longitude} and \code{time} is used to calculate 
    the solar zenith angle, if no \code{sza} is specified (see also \code{time_interval}).
    
    \code{latitude} in combination with \code{longitude} and \code{time} is used to choose a
    suitable default atmosphere file, if no \code{atmosphere_file} is specified. 
        ''',

        'longitude' : r'''
    This option can be used to specify the longitude of the location to simulate.
    (This option only has an effects, if \code{latitude} is specified, too.)
    \fcode{
    longitude E/W deg [min] [sec]
    }
    where \code{E} and \code{W} stand for the eastern and western hemisphere, respectively.
    \code{deg min sec} is the position in degrees, arc minutes, and arc seconds.
    \code{deg} might also be a float number. \code{min} and \code{sec} may be obmitted.
    For possible usage of the \code{longitude} information, see \code{latitude}. 
        ''',

        'mc_bcond' : r'''
    Undocumented option to define MYSTIC boundary conditions: periodic (default),
    absorbing (\code{absorb}), mirror (\code{mirror}) -- make sure that you understand 
    what you are doing! While \code{mirror} gives reasonable results in the thermal spectral range,
    the direct sun is messed up because it changes direction (by mirroring) without being scattered.
        ''',

        'mc_blitz_position' : r'''
    Define lightning as source. The position of the the lightning light
    source is defined as a simple line, going from [ X\_1 Y\_1 Z\_1 ] to [
    X\_2 Y\_2 Z\_2 ] (in units of m). Usage:
    \fcode{
    mc_blitz_position X_1 Y_1 Z_1 X_2 Y_2 Z_2
    }
    Only works in forward mode, for radiances.
        ''',

        'mc_bw_umu_file'  : r'''
    Define a different umu and phi for each pixel. You specify the file
    containing a two-dimensional array with umu and phi for each sample
    pixel. Only works in connection with \code{mc_backward}.
    \fcode{
    mc_bw_umu_file file
    }
        ''',

        'mc_coherent_backscatter' : r'''
    Switches on coherent backscattering for lidar and radar. Works only 
    with \code{mc_lidar} or \code{mc_radar}, use only with \code{mc_polarisation}.
        ''',

        'mc_panorama' : r'''
    Use this option only together with \code{mc_panorama_view}.

    \fcode{
    mc_panorama circumsolar_var_red
    }
    Variance reduction method for the calculation of circumsolar radiation.
    Re-distributes the initial backward photon directions according to DDIS phase\_max.
    Use only if you calculate the average radiance for one large FOV (at least two times of the angular radius of the sun disk, e.g. 
    if your sensor convers 0$^\circ$ to 1$^\circ$ from the sun center or 0$^\circ$ to 5$^\circ$ from the sun
    center). For smaller FOVs (e.g 0$^\circ$ to 0.35$^\circ$) it may slow down convergence as this method does not yet consider
    the extent of the sun disk. Results should not be wrong but you may need more photons to achieve the same 
    standard deviation. Consideration of the sun's extent in this method should be easy to implement, however.
    For the simulations of radiance profiles (so called sunshapes) this option should probably not be used as it will
    only increase computing time without any benefit. For further information see Reinhardt (2013) Sect. 3.4.2 
    (Dissertation, LMU, urn:nbn:de:bvb:19-164380, http://edoc.ub.uni-muenchen.de/16438/).

    \fcode{
    mc_panorama distr_photons_over_pixel
    }
    With only \code{mc_panorama_view} all photons for the same pixel
    are emitted into exactly the same direction. The finite solid angle that
    is covered by one pixel is not accounted for. With 
    \code{mc_panorama distr_photons_over_pixel} the photons are distributed
    over the FOV covered by the pixel. 
    
    \fcode{
    mc_panorama no_pixel
    }
    Useful for calculation of radiance distributions which shall be exactly
    charted. Originally \code{mc_panorama_view} was meant to calculate images
    like produced by a CCD camera. You set the edges of your desired image
    and the resolution of the sensor. MYSTIC will then calculate the
    radiances at the center of each pixel. This implies that values for
    the edges of the image as given by the input file are never
    computed. The border angles for which radiances are computed are
    shifted half a pixel to the inside of the domain. With
    \code{mc_panorama no_pixel} the edging values as given in the input
    file are directly hit. A short example: Assume you only have three
    pixels. You choose your edges to be at $\phi=0$ and $\phi=180$. With
    \code{mc_panorama no_pixel} MYSTIC will calculate radiances at
    $\phi=0$, $\phi=90$ and $\phi=180$. Without this option you will get
    values for $\phi=30$, $\phi=90$ and $\phi=150$.
    
    \fcode{
    mc_panorama quicklook
    }
    This option allows to do a fast calculation in order to see whether you
    have set your parameters right. All panorama pixels are run in a fast
    mode: only one photon is started per pixel, it is moved till to a
    fixed optical depth (0.1), and then a local estimate is performed. In
    case the local estimate is done scattering (i.e.~the photon has not
    hit a surface before reaching tau=0.1), isotropic scattering is
    assumed. Then the photon is killed. We suggest you turn off Rayleigh
    scattering and aerosols while using this option. The result lets you
    see roughly where your panorama is looking at surface, sunlit clouds,
    shaded clouds, or slanted clouds.
    
    \fcode{
    mc_panorama with_direct_rad
    } 
    Output is direct plus diffuse radiance.
    Use together with \code{mc_panorama distr_photons_over_pixel}.

    \fcode{
    mc_panorama weight_with_cos
    } 
    Output is weighted with the cosine of the angle between photon and zenith of camera.
    Use together with \code{mc_panorama distr_photons_over_pixel}.
        ''',

        'mc_panorama_alignment' : r'''
    Use this option only together with \code{mc_panorama_view}. This
    option has been introduced to facilitate e.g.~the calculation of
    circum solar radiation. Three sub-options are available, \code{sun},
    \code{mu}, and \code{zenith}. When set \code{sun}, the camera is
    aligned to sun (i.e.~using the options \code{sza} and
    \code{phi0}). $\theta$ and $\mu$ as given by the option
    \code{mc_panorama_view} now relates to the sun's position. $\mu$ sets
    cosine of the angular distance from the center of the sun. $\mu=$-1
    implies looking directly into the sun. With $\phi$ we select the
    direction around the sun on a circle with angular distance ${\rm
    acos}(\mu)$ from the center of the sun. $\phi=0$ implies that the
    camera pixel is looking in a direction exactly below the sun, at
    $\phi=180$ the camera pixel will be looking exactly above the
    sun. Increasing $\phi$ means the direction of looking moves clock-wise
    around the sun-camera-line (seen from the sun). If the sun is in the
    zenith, the definition of $\phi$ is identical to that in \code{phi}.
    
    The sub-option \code{mu} aligns the camera geometry according to the
    options \code{umu} and \code{phi}. Be aware that \code{phi0} and
    \code{phi} are defined differently!
    
    The sub-option \code{zenith} makes only sense when using
    \code{mc_spherical 3D}. Then the camera is aligned to the zenith with
    respect to the lon-lat position of the camera.
        ''',

        'mc_panorama_view'  : r'''
    Simulate a panorama. Only works with \code{mc_sensorposition},
    \code{mc_sample_grid} and \code{mc_backward}. Can best be explained by
    following example:

    \fcode{\\
    mc_panorama_view PHI_1 PHI_2 THETA_1 THETA_2\\
    mc_sensorposition X Y Z \\
    mc_sample_grid NPHI NTHETA [] []\\
    mc_backward NPHI1 NTHETA1 NPHI2 NTHETA2\\
    }

    or, with numbers:

    \fcode{\\
    mc_panorama_view 45 135 80 180\\
    mc_sensorposition 3000 3000 26 \\
    mc_sample_grid 90 100 [] []\\
    mc_backward 10 0 10 99\\
    }

    In this example, \code{mc_panorama_view} defines the camera to cover
    the area between 10 degrees below the horizon (80 degrees from nadir,
    looking into the ground) and the zenith (180 degrees from nadir), and
    looking into the directions between south-west (45) and
    north-west (135). The camera is positioned
    at x=3km, y=3km and z=26m, see \code{mc_sensorposition}. The number of
    pixels that the camera has is defined by \code{mc_sample_grid}, hence
    the camera horizontal resolution is 1 degree in the
    horizon: (135-45)/90. The camera vertical resolution is also 1
    degree: (180-80)/100. This run only calculates one vertical line at
    phi=55.5 (see \code{mc_backward}). The square brackets in the line
    \code{mc_sample_grid} are usually dx and dy, they are set
    automatically in the presence of 3d clouds. If no clouds are present,
    you need to set them, if you use \code{mc_elevation_file}, you need to
    set them consistent with the elevation file. NOTE: If you have
    formerly used \code{mc_panorama}, this is what has changed: The order
    of the numbers has changed: Formerly, we wrote \code{mc_panorama MU_1
    MU_2 PHI_1 PHI_2}, where MU\_i=cos(THETA\_i). In the options
    \code{mc_sample_grid} and \code{mc_backward}, X and Y have been
    swapped. Although these changes are far from being nice, they are
    necessary to introduce a conform nomenclatur throughout
    libRadtran. Sorry for that. (CAUTION! The definition of phi for
    panoramas was formerly different than that of \code{phi}!  This has
    been changed on 19.11.2010!)
        ''',

        'mc_panorama_forward'  : r'''
    Simulate a panorama by forward tracing. Only works with
    \code{mc_sample_grid}. Can best be explained by
    following example:

    \fcode{\\
    mc_panorama_forward\\
    mc_sample_grid NPHI NTHETA [] []\\
    }

    In this example, \code{mc_panorama_forward} calculated an all-sky (4 pi) panorama with NTHETA * NPHI grid points.
    Please note that the results are averaged over the whole domain. \code{panorama_forward} therefore doesn't make
    much sense with 3D cloud or surface geometry.
        ''',

        'mc_reference_to_nn': r'''
    Only works with \code{mc_backward}. The sampled pixels correspond to
    the surface pixels. In other words, the photons are still counted at
    zout, but the sample pixel which they are counted in corresponds to
    the pixel at sea level. In other words, the sample grid is projected
    from nn to zout using umu and phi.
        ''',

        'mc_sample_grid'  : r'''
    Sample grid size (Nx Ny [dx dy]) for MYSTIC.
    \fcode{
    mc_sample_grid Nx Ny [dx dy]
    }
    Only meaningful with \code{rte_solver montecarlo}.
        ''',

        'mc_satellite_position' : r'''
    Define the position of the satellite. See \code{mc_satellite_view} for
    description.
    \fcode{
    mc_satellite_position longitude latitude alitude
    }
        ''',
        'mc_satellite_view' : r'''
    Simulate a satellite image. Only works with
    \code{mc_satellite_position}/\code{mc_sensorposition},
    \code{mc_sample_grid}, \code{mc_backward} and
    \code{mc_spherical 3D}. Can best be explained by following example:
    
    \fcode{\\
    mc_spherical 3D \\
    mc_satellite_view PHI_W PHI_E THETA_S THETA_N \\
    mc_sample_grid N_PHI N_THETA \\
    mc_satellite_position SAT_LON SAT_LAT SAT_ALT \\
    mc_sun_position SUN_LON SUN_LAT \\
    mc_backward \\
    }
    
    In this example, the earth is a sphere (not plane parallel). The
    satellite is in zenith above longitude SAT\_LON degrees, latitude
    SAT\_LAT degrees, and at an altitude of SAT\_ALT
    meters. (Alternatively, you can use \code{mc_sensorposition S_X S_Y
    S_Z}, with [S\_X,0,0] being toward MSG, [0,S\_Y,0] being above the
    Maledives, and [0,0,S\_Z] being above the north pole.) The image taken
    has a resolution of N\_THETA times N\_PHI, and the picture (taken with
    a camera with spherical picture geometry) spans from PHI\_W to PHI\_E
    degrees (positive being east of nadir) and from THETA\_S to THETA\_N
    degrees (positive being north of nadir). In this example, the sun is
    in zenith at the latitude SUN\_LAT degrees (positive meaning north) and
    longitude SUN\_LON degrees (positive meaning east). See also
    \code{mc_panorama_view}, which is a similar option.  Example for MSG:
    
    \fcode{\\
    mc_spherical 3D \\
    mc_satellite_view -10 10 -10 10 \\
    mc_sample_grid 4000 4000 \\
    mc_satellite_position 0 0 3.6e7 \\
    mc_sun_position -10 23\# about 1PM UTC on mid summer \\
    mc_backward \\
    }
    
    Only meaningful with \code{rte_solver montecarlo}.
        ''',

        'mc_sensordirection'  : r'''
    Define viewing direction of an irradiance sensor in Monte Carlo backward mode. 
    \fcode{
    mc_sensordirection x-value y-value z-value
    }
    Has been introduced for 
    irradiance calculations in topography and might not properly work with all options. 
    For radiance use the usual \code{umu} and \code{phi}. 
        ''',

        'mc_sensorposition'  : r'''
    Define the position of a sensor. 
    \fcode{
    mc_sensorposition x-value y-value z-value
    }
    Has been introduced for irradiance calculations in topography
    and might not properly work with all options. 
        ''',

        'mc_spherical' : r'''
    Spherical geometry in MYSTIC. 
    \fcode{
    mc_spherical 1D
    }
    Works only in "1D" - \code{wc_file 3D} 
    and \code{ic_file 3D} are not yet considered. If \code{mc_spherical} is selected
    \code{mc_backward} is switched on automatically. Viewing
    direction (\code{umu}, \code{phi}) and sun position
    (\code{sza}, \code{phi0}) are defined with respect to the
    sensor position specified by \code{zout}. For details about
    the implementation of spherical geometry please refer to
    \citet{emde2007}. 

    \iflidar{
    \fcode{
    mc_spherical 3D
    }
    3D spherical geometry in MYSTIC. Only works with
    \code{mc_satellite_position}/\code{mc_sensorposition}, and
    \code{mc_satellite_view}, and only for solar backward. You have to set
    the sun position using \code{mc_sun_position}. Still very
    experimental. Does not work with \code{mc_elevation_file}.
    }
        ''',

        'mc_spherical3D_scene' : r'''
    Define scene on globe. Only works with \code{mc_spherical 3D}. Without
    this option, 3d clouds and 2d albedo maps are folded around the full
    globe when using \code{mc_spherical 3D}. This option defines the limits
    of the clouds and albedo maps in terms of longitude and latitude. Example:
    
    \fcode{\\
    mc_spherical3D_scene LON_W LON_E LAT_S LAT_N \\
    }
    
    The scene spans from LON\_W degrees in the west to LON\_E degrees in
    the east (positive values are east, negative are west of Greenwich)
    and from LAT\_S degrees in the south to LAT\_N degrees in the
    north (positive are northern hemisphere, negative are southern
    hemisphere).
        ''',

        'mc_sun_position' : r'''
    Sun position when using \code{mc_spherical 3D}. See \code{mc_satellite_view}
    for an example and explanation.
    \fcode{
    mc_sun_position longitude latitude
    }
        ''',

        'phi' : r'''
    Azimuth output angles (in degrees) in increasing order. 
    \fcode{
    phi values
    }
    The radiance is output at \code{phi} and \code{umu}.
    \begin{itemize}
    \item Sensor in the North (looking South):  0 deg 
    \item Sensor in the East  (looking West):   90 deg 
    \item Sensor in the South (looking North):  180 deg 
    \item Sensor in the West  (looking East):   270 deg 
    \end{itemize}
    For all one-dimensional solvers the absolute azimuth does not matter, but only the relative azimuth \code{phi}-\code{phi0}.
        ''',

        'phi0' : r'''
    Azimuth angle of the sun (0 to 360 degrees). 
    \fcode{phi0 value}
    \begin{itemize}
    \item Sun in the South:    0 degrees
    \item Sun in the West:    90 degrees
    \item Sun in the North: 180 degrees
    \item Sun in the East: 270 degrees
    \end{itemize}
    For all one-dimensional solvers the absolute azimuth does not matter, but only the relative azimuth \code{phi}-\code{phi0}.
        ''',

        'sza' : r'''
    The solar zenith angle (degrees). 
    \fcode{
    sza value
    }
    The default solar zenith angle is 0.
        ''',

        'sza_file' : r'''
    Location of solar zenith angle file for wavelength-dependent solar
    zenith angle. 
    \fcode{
    sza_file file
    }
    This option is useful if you want to simulate an instrument 
    which scans so slowly that the solar zenith angle may change significantly during 
    the wavelength scan. 
    The file must have two or three columns. Column 1 is the wavelength, 
    in nm, and column 2 the corresponding  solar zenith angle. Optionally
    the third column may contain the corresponding solar azimuth angle. The solar 
    azimuth angle is only needed when calculating radiances. The wavelength 
    grid may be freely set. The solar zenith and azimuth angle will be interpolated to 
    the wavelength grid used for the radiation calculation.
    Comments start with \code{\#}. Empty lines are ignored.
        ''',

        'umu' : r'''
    Cosine of output polar angles in increasing order, starting with 
    negative (downwelling radiance, looking upward) values (if any) and on through 
    positive (upwelling radiance, looking downward) values. 
    Must not be zero.
    \fcode{
    umu values
    }
        ''',
    }
