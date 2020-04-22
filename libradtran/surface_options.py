
"""--------------------------------------------------------------------
 * $Id: surface_options.py 3500 2019-09-30 12:31:35Z Claudia.Emde $
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

class setup_surface_group():

    group_name = 'Surface'

    def __init__(self):
        documentation = get_surface_documentation()

        brdf_non_parents=['albedo','albedo_file','albedo_library',
                            'brdf_rossli', 'brdf_rossli_file', 'brdf_rossli_hotspot',
                            'brdf_ambrals', 'brdf_ambrals_file', 'brdf_ambrals_hotspot',
                            'brdf_cam', 'brdf_cam_solar_wind', 'cox_and_munk_u10_map', 
                            'brdf_rpv_file', 'brdf_rpv', 'brdf_rpv_library','brdf_rpv_type',
                            'brdf_hapke', 'brdf_hapke_file',
                            'mc_albedo_file', 'mc_albedo_spectral_file', 'mc_albedo_type',
                            'mc_rossli_file', 'mc_ambrals_file',
                            'mc_ambrals_type',
                            'mc_ambrals_spectral_file', 'mc_rpv_type',
                            'mc_rpv_spectral_file', 'mc_rpv_file',
                            'mc_triangular_surface_file',
                            'bpdf_tsang_u10']

        altitude = option_definition.option(
            name='altitude',
            group='surface',
            helpstr='Set the bottom level atmosphere provided in atmosphere_file to be at the given altitude.',
            documentation=documentation['altitude'], 
            gui_inputs=(GUI_definition.FloatInput(name='Input.alt.altitude', default='NOT_DEFINED_FLOAT', valid_range=[-1000000.0, 1000000.0]), GUI_definition.FloatInput(name='Input.alt.altitude_dz', valid_range=[0, 1000000.0], optional=True),),
            tokens = [option_definition.addSetting(name='Input.alt.source',setting='ALT_DIRECT_INPUT', default='ALT_NOT_DEFINED'), 
                      option_definition.addToken(name='Input.alt.altitude',datatype=float, default='NOT_DEFINED_FLOAT', valid_range=[-1e6,1e6]),
                      option_definition.addToken(name='Input.alt.altitude_dz', datatype=float, default=0.0, valid_range=[0,1e6], optional=True)],
            parents=['uvspec'], 
        )

        mc_elevation_file = option_definition.option(
            name='mc_elevation_file', 
            group='surface',
            helpstr='Define a MYSTIC 2D elevation input file.', 
            documentation=documentation['mc_elevation_file'], 
            gui_inputs=(GUI_definition.FileInput(name='Input.rte.mc.filename[FN_MC_ELEVATION]'),),
            tokens=option_definition.addToken(name='Input.rte.mc.filename[FN_MC_ELEVATION]',datatype=io.IOBase),
            parents=['uvspec'], 
            speaker='rte_solver',
            enable_values=("mystic","montecarlo"),
            threedmystic =True
        )

        albedo = option_definition.option(
            name='albedo',
            group='surface',
            helpstr='Lambertian surface albedo',
            documentation=documentation['albedo'], 
            gui_inputs=(GUI_definition.FloatInput(name='Input.alb.albedo', default=0.0, valid_range=[0, 1]),),
            tokens = [option_definition.addToken(name='Input.alb.albedo', datatype=float, default=0.0, valid_range=[0,1]),
                      option_definition.addSetting(name='Input.alb.source', setting='ALBEDO_CONSTANT', default='NOT_DEFINED_INTEGER')], 
            parents=['uvspec'], 
            non_parents=brdf_non_parents,
        )

        albedo_file = option_definition.option(
            name='albedo_file',
            group='surface',
            helpstr='Location of surface albedo file for wavelength dependent surface albedo.',
            documentation=documentation['albedo_file'], 
            gui_inputs=(GUI_definition.FileInput(name='Input.filename[FN_ALBEDO]'),),
            tokens = [option_definition.addToken(name='Input.filename[FN_ALBEDO]', datatype=io.IOBase), 
                      option_definition.addSetting(name='Input.alb.source', setting='ALBEDO_FROM_ALBEDO_FILE', default='NOT_DEFINED_INTEGER'), 
                      option_definition.addSetting(name='Input.alb.surf_type_map', setting=False, default=False), 
                      option_definition.addSetting(name='Input.alb.albedo', setting=-1.0, default=0.0)],
            parents=['uvspec'],  
            non_parents=brdf_non_parents, 
            plot = {'plot_type': '2D',
                'optional_args': {'column_names': (
                        "wavelength",
                        "Albedo",)
                          }
                }
        )

        brdf_cam = option_definition.option(
            name='brdf_cam',
            group='surface',
            helpstr='',
            documentation=documentation['brdf_cam'], 
            gui_inputs=(GUI_definition.ListInput(name='id',valid_range=['pcl','sal','u10','uphi']),
                        GUI_definition.FloatInput(name='Input.cm.cam[id]', default='NOT_DEFINED_FLOAT', valid_range=[0, 1e6]),),
            tokens = [option_definition.addLogical(name='id',logicals=['pcl','sal','u10','uphi'],setting='BRDF_CAM_'),
                      option_definition.addToken(name='Input.cm.param[id]', datatype=float, valid_range=[0,1e6]),
                      option_definition.addSetting(name='Input.disort2_brdf', setting='BRDF_CAM', default='BRDF_NONE')], 
            parents=['uvspec'],
            non_parents=brdf_non_parents,
                        non_parent_exceptions=['brdf_cam', 'brdf_cam_solar_wind', 'cox_and_munk_u10_map', 'mc_rossli_file', 'mc_ambrals_file'],
            speaker='rte_solver',
            enable_values=("disort","mystic"),
            non_unique=True,
        )
        
        brdf_cam_solar_wind = option_definition.option(
            name='brdf_cam_solar_wind',
            group='surface',
            helpstr='Wind azimuth identical to the incoming photon azimuth.',
            documentation=documentation['brdf_cam_solar_wind'],
            tokens=[option_definition.addSetting(name='Input.disort2_brdf', setting='BRDF_CAM', default='BRDF_NONE'), 
            option_definition.addSetting(name='Input.cm.solar_wind', setting=1, default=0)],
            parents=['uvspec'],
            non_parents=brdf_non_parents,
                        non_parent_exceptions=['brdf_cam', 'brdf_cam_solar_wind', 'cox_and_munk_u10_map', 'mc_rossli_file', 'mc_ambrals_file'],
            speaker='rte_solver',
            enable_values=("mystic",)
        )
        
        brdf_hapke = option_definition.option(
            name='brdf_hapke',
            group='surface',
            helpstr='Constant parameters for Hapke BRDF.',
            documentation=documentation['brdf_hapke'], 
            gui_inputs=(GUI_definition.ListInput(name='id',valid_range=['b0','h', 'w']), GUI_definition.FloatInput(name='Input.hapke.hapke[id]'),),
            tokens= [option_definition.addLogical(name='id', logicals=['b0','h', 'w'], setting='BRDF_HAPKE_'), 
                     option_definition.addToken(name='Input.hapke.hapke[id]', datatype=float)],
            parents=['uvspec'], 
            non_parents=brdf_non_parents,
            non_unique=True,
        )

        brdf_hapke_file = option_definition.option(
            name='brdf_hapke_file',
            group='surface',
            helpstr='File containing the Hapke BDRF parameterization.', 
            documentation=documentation['brdf_hapke_file'], 
            gui_inputs=(GUI_definition.FileInput(name='Input.filename[FN_HAPKE]'),),
            tokens = [option_definition.addToken(name='Input.filename[FN_HAPKE]', datatype=io.IOBase), 
                      option_definition.addSetting(name='Input.hapke.source', setting='HAPKE_FROM_HAPKE_FILE', default='HAPKE_CONSTANT')], 
            parents=['uvspec'], 
            non_parents=brdf_non_parents,
        )

        brdf_rossli = option_definition.option(
            name='brdf_rossli',
            group='surface',
            helpstr='Three-parameter Ross-Li BRDF fit for vegetated and non-vegetated surfaces.',
            documentation=documentation['brdf_rossli'],
            gui_inputs=(GUI_definition.ListInput(name='id', valid_range=['iso','vol','geo']), GUI_definition.FloatInput(name='Input.rossli.rossli[id]'),),
            tokens=[option_definition.addLogical(name='id', logicals=['iso','vol','geo'], setting='BRDF_ROSSLI_'),
                    option_definition.addToken(name='Input.rossli.rossli[id]', datatype=float)], 
            parents=['uvspec'],
            non_parents=brdf_non_parents, non_parent_exceptions=['brdf_rossli_hotspot'],
            non_unique=True,
            showInGui = False,
            developer = True
        )

        brdf_rossli_hotspot = option_definition.option(
            name='brdf_rossli_hotspot', 
            group='surface',
            helpstr='Turn on hot spot correction factor in Ross-Li BRDF.', 
            documentation=documentation['brdf_rossli_hotspot'], 
            tokens=option_definition.addSetting(name='Input.rossli.hotspot', setting='BRDF_ROSSLI_HOTSPOT_ON', default='BRDF_ROSSLI_HOTSPOT_OFF'),
            parents=['uvspec'],
            non_parents=brdf_non_parents, non_parent_exceptions=['brdf_rossli'],
            showInGui = False,
            developer = True
            )


        brdf_rossli_file = option_definition.option(
            name='brdf_rossli_file',
            group='surface',
            helpstr='File containing the Ross-Li BDRF parameterization.', 
            documentation=documentation['brdf_rossli_file'], 
            gui_inputs=(GUI_definition.FileInput(name='Input.filename[FN_ROSSLI]'),),
            tokens = [option_definition.addToken(name='Input.filename[FN_ROSSLI]', datatype=io.IOBase), 
                      option_definition.addSetting(name='Input.rossli.source', setting='ROSSLI_FROM_ROSSLI_FILE', default='ROSSLI_CONSTANT')], 
            parents=['uvspec'], 
            non_parents=brdf_non_parents,non_parent_exceptions=['brdf_rossli_hotspot'],
            showInGui = False,
            developer = True
        )

        brdf_ambrals = option_definition.option(
            name='brdf_ambrals',
            group='surface',
            helpstr='Three-parameter Ross-Li BRDF fit for vegetated and non-vegetated surfaces.',
            documentation=documentation['brdf_ambrals'],
            gui_inputs=(GUI_definition.ListInput(name='id', valid_range=['iso','vol','geo']), GUI_definition.FloatInput(name='Input.rossli.rossli[id]'),),
            tokens=[option_definition.addSetting(name='Input.rossli.source', setting='ROSSLI_AMBRALS_CONSTANT', default='ROSSLI_CONSTANT'),
                    option_definition.addLogical(name='id', logicals=['iso','vol','geo'], setting='BRDF_ROSSLI_'),
                    option_definition.addToken(name='Input.rossli.rossli[id]', datatype=float)], 
            parents=['uvspec'],
            non_parents=brdf_non_parents, non_parent_exceptions=['brdf_ambrals_hotspot'],
            non_unique=True,
        )

        brdf_ambrals_hotspot = option_definition.option(
            name='brdf_ambrals_hotspot', 
            group='surface',
            helpstr='Turn on hot spot correction factor in Ross-Li BRDF.', 
            documentation=documentation['brdf_ambrals_hotspot'], 
            tokens=option_definition.addSetting(name='Input.rossli.hotspot', setting='BRDF_ROSSLI_HOTSPOT_ON', default='BRDF_ROSSLI_HOTSPOT_OFF'),
            parents=['uvspec'],
            non_parents=brdf_non_parents, non_parent_exceptions=['brdf_ambrals'],
            )


        brdf_ambrals_file = option_definition.option(
            name='brdf_ambrals_file',
            group='surface',
            helpstr='File containing the Ross-Li BDRF parameterization.', 
            documentation=documentation['brdf_ambrals_file'], 
            gui_inputs=(GUI_definition.FileInput(name='Input.filename[FN_ROSSLI]'),),
            tokens = [option_definition.addToken(name='Input.filename[FN_ROSSLI]', datatype=io.IOBase), 
                      option_definition.addSetting(name='Input.rossli.source', setting='ROSSLI_FROM_AMBRALS_FILE', default='ROSSLI_CONSTANT')], 
            parents=['uvspec'], 
            non_parents=brdf_non_parents,non_parent_exceptions=['brdf_ambrals_hotspot'],
        )

        fluorescence = option_definition.option(
            name='fluorescence',
            group='surface',
            helpstr='Bottom surface isotropic fluorescence source.',
            documentation=documentation['fluorescence'], 
            gui_inputs=(GUI_definition.FloatInput(name='Input.flu.fluorescence'),),
            tokens = [option_definition.addToken(name='Input.flu.fluorescence', datatype=float, default=0.0),
                      option_definition.addSetting(name='Input.flu.source', setting='FLUORESCENCE_CONSTANT', default='NOT_DEFINED_INTEGER')], 
            parents=['uvspec'], 
            non_parents=['fluorescence_file'],
        )

        fluorescence_file = option_definition.option(
            name='fluorescence_file',
            group='surface',
            helpstr='Bottom surface isotropic fluorescence source file, wavelength dependent.',
            documentation=documentation['fluorescence_file'], 
            gui_inputs=(GUI_definition.FileInput(name='Input.filename[FN_FLUORESCENCE]'),),
            tokens = [option_definition.addToken(name='Input.filename[FN_FLUORESCENCE]', datatype=io.IOBase), 
                      option_definition.addSetting(name='Input.flu.source', setting='FLUORESCENCE_FROM_FLUORESCENCE_FILE', default='NOT_DEFINED_INTEGER')],
            parents=['uvspec'],  
            non_parents=['fluorescence'], 
            plot = {'plot_type': '2D',
                'optional_args': {'column_names': (
                        "wavelength",
                        "Fluorescence",)
                          }
                }
        )

        brdf_rpv_file = option_definition.option(
            name='brdf_rpv_file',
            group='surface',
            helpstr='File containing the RPV BDRF parameterization.', 
            documentation=documentation['brdf_rpv_file'], 
            gui_inputs=(GUI_definition.FileInput(name='Input.filename[FN_RPV]'),),
            tokens = [option_definition.addToken(name='Input.filename[FN_RPV]', datatype=io.IOBase), 
                      option_definition.addSetting(name='Input.rpv.source', setting='RPV_FROM_RPV_FILE', default='RPV_CONSTANT')], 
            parents=['uvspec'], 
            non_parents=brdf_non_parents,
            non_parent_exceptions=['brdf_rpv_file', 'brdf_rpv', 'brdf_rpv_library'],
        )
        
        brdf_rpv = option_definition.option(
            name='brdf_rpv',
            group='surface',
            helpstr='Constant RPV parameters.',
            documentation=documentation['brdf_rpv'], 
            gui_inputs=(GUI_definition.ListInput(name='id',valid_range=['k','rho0', 'theta', 'sigma', 't1', 't2', 'scale']), GUI_definition.FloatInput(name='Input.rpv.rpv[id]'),),
            tokens= [option_definition.addLogical(name='id', logicals=['k','rho0', 'theta', 'sigma', 't1', 't2', 'scale'], setting='BRDF_RPV_'), 
                     option_definition.addToken(name='Input.rpv.rpv[id]', datatype=float)],
            parents=['uvspec'], 
            non_parents=brdf_non_parents,
            non_parent_exceptions=['brdf_rpv_file', 'brdf_rpv', 'brdf_rpv_library'],
            non_unique=True,
        )


        brdf_rpv_library = option_definition.option(
            name='brdf_rpv_library',
            group='surface',
            helpstr='Collections of spectral BRDFs of different surface types.', 
            documentation=documentation['brdf_rpv_library'], 
            tokens=option_definition.addToken(name='Input.rpv.library', datatype=str, valid_range=['IGBP', io.IOBase]), 
            parents=['uvspec'], 
            non_parents=brdf_non_parents,
            non_parent_exceptions=['brdf_rpv_file', 'brdf_rpv', 'brdf_rpv_library','brdf_rpv_type'],
            childs=['brdf_rpv_type']
        )
        
        brdf_rpv_type = option_definition.option(
            name='brdf_rpv_type',
            group='surface',
            helpstr='Select surface type.',
            documentation=documentation['brdf_rpv_type'], 
            gui_inputs=(GUI_definition.IntegerInput(name='Input.alb.surface', default='NOT_DEFINED_INTEGER', valid_range=[0, 99]),),
            tokens = [option_definition.addToken(name='Input.alb.surface', datatype=int, default='NOT_DEFINED_INTEGER', valid_range=[0,99]), 
                      option_definition.addSetting(name='Input.alb.surf_type_map', setting=False)], 
            non_parents=brdf_non_parents,
            non_parent_exceptions=['brdf_rpv_file', 'brdf_rpv', 'brdf_rpv_library','brdf_rpv_type'],
            parents=['albedo_library','brdf_rpv_library'], 
        )
        
        sur_temperature = option_definition.option(
            name='sur_temperature', 
            group='surface',
            helpstr='Surface temperature.',
            documentation=documentation['sur_temperature'], 
            gui_inputs=(GUI_definition.FloatInput(name='Input.surface_temperature'),),
            tokens=option_definition.addToken(name='Input.surface_temperature', datatype=float, default='NOT_DEFINED_FLOAT'), 
            parents=['uvspec'], 
        )

        sur_temperature_file = option_definition.option(
            name='sur_temperature_file',
            group='surface', 
            helpstr='Define a MYSTIC 2D temperature input file.',
            documentation=documentation['sur_temperature_file'],
            gui_inputs=(GUI_definition.FileInput(name='Input.rte.mc.filename[FN_MC_TEMPERATURE]'),),
            tokens=option_definition.addToken(name='Input.rte.mc.filename[FN_MC_TEMPERATURE]', datatype=io.IOBase),    
            parents=['uvspec'], 
            speaker='rte_solver',
            enable_values=("mystic","montecarlo"),
            threedmystic =True
        )

        mc_albedo_file = option_definition.option(
            name='mc_albedo_file',
            group='surface', 
            helpstr='Define a MYSTIC 2D albedo input file.',
            documentation=documentation['mc_albedo_file'],
            gui_inputs=(GUI_definition.FileInput(name='Input.rte.mc.filename[FN_MC_ALBEDO]'),),
            tokens=option_definition.addToken(name='Input.rte.mc.filename[FN_MC_ALBEDO]', datatype=io.IOBase),    
            parents=['uvspec'], 
            non_parents=brdf_non_parents,
            speaker='rte_solver',
            enable_values=("mystic","montecarlo"),
            threedmystic =True
        )

        mc_albedo_type = option_definition.option(
            name='mc_albedo_type',
            group='surface', 
            helpstr='Define a MYSTIC 2D albedo type file.',
            documentation=documentation['mc_albedo_type'],
            gui_inputs=(GUI_definition.FileInput(name='Input.rte.mc.filename[FN_MC_ALBEDO_SPECTRAL]'),),
            tokens=option_definition.addToken(name='Input.rte.mc.filename[FN_MC_ALBEDO_SPECTRAL]',datatype=io.IOBase),    
            parents=['uvspec'], 
            non_parents=brdf_non_parents,
            non_parent_exceptions=['mc_albedo_spectral_file'],
            speaker='rte_solver',
            enable_values=("mystic","montecarlo"),
            threedmystic =True

        )

        mc_albedo_spectral_file = option_definition.option(
            name='mc_albedo_spectral_file',
            group='surface', 
            helpstr='Define a MYSTIC 2D spectral albedo input file.',
            documentation=documentation['mc_albedo_spectral_file'],
            gui_inputs=(GUI_definition.FileInput(name='Input.rte.mc.filename[FN_MC_ALBEDO_TYPE]'),),
            tokens=option_definition.addToken(name='Input.rte.mc.filename[FN_MC_ALBEDO_TYPE]', datatype=io.IOBase),
            parents=['uvspec'], 
            non_parents=brdf_non_parents,
            non_parent_exceptions=['mc_albedo_type'],
            speaker='rte_solver',
            enable_values=("mystic","montecarlo"),
            threedmystic =True

        )

        mc_rossli_file = option_definition.option(
            name='mc_rossli_file',
            group='surface',
            helpstr='Define a MYSTIC 2D Ross-Li BRDF input file.',
            documentation=documentation['mc_rossli_file'], 
            gui_inputs=(GUI_definition.FileInput(name='Input.rte.mc.filename[FN_MC_ROSSLI]'),),
            tokens=option_definition.addToken(name='Input.rte.mc.filename[FN_MC_ROSSLI]', datatype=io.IOBase),    
            parents=['uvspec'], 
            non_parents=brdf_non_parents,
            non_parent_exceptions=['brdf_cam', 'brdf_cam_solar_wind', 'cox_and_munk_u10_map'],
            speaker='rte_solver',
            enable_values=("mystic","montecarlo"),
            threedmystic =True,
            showInGui = False,
            developer = True
        )

        mc_ambrals_file = option_definition.option(
            name='mc_ambrals_file',
            group='surface',
            helpstr='Define a MYSTIC 2D AMBRALS (Ross-Li) BRDF input file.',
            documentation=documentation['mc_ambrals_file'], 
            gui_inputs=(GUI_definition.FileInput(name='Input.rte.mc.filename[FN_MC_AMBRALS]'),),
            tokens=option_definition.addToken(name='Input.rte.mc.filename[FN_MC_AMBRALS]', datatype=io.IOBase),    
            parents=['uvspec'], 
            non_parents=brdf_non_parents,
            non_parent_exceptions=['brdf_cam', 'brdf_cam_solar_wind', 'cox_and_munk_u10_map'],
            speaker='rte_solver',
            enable_values=("mystic","montecarlo"),
            threedmystic =True
        )

        mc_ambrals_type = option_definition.option(
            name='mc_ambrals_type',
            group='surface', 
            helpstr='Define a MYSTIC 2D ambrals spectral BRDF type file.',
            documentation=documentation['mc_ambrals_type'],
            gui_inputs=(GUI_definition.FileInput(name='Input.rte.mc.filename[FN_MC_AMBRALS_SPECTRAL]'),),
            tokens=option_definition.addToken(name='Input.rte.mc.filename[FN_MC_AMBRALS_SPECTRAL]',datatype=io.IOBase),    
            parents=['uvspec'], 
            non_parents=brdf_non_parents,
            non_parent_exceptions=['mc_ambrals_spectral_file'],
            speaker='rte_solver',
            enable_values=("mystic","montecarlo"),
            threedmystic =True

        )

        mc_ambrals_spectral_file = option_definition.option(
            name='mc_ambrals_spectral_file',
            group='surface', 
            helpstr='Define a MYSTIC 2D ambrals spectral BRDF input file.',
            documentation=documentation['mc_ambrals_spectral_file'],
            gui_inputs=(GUI_definition.FileInput(name='Input.rte.mc.filename[FN_MC_AMBRALS_TYPE]'),),
            tokens=option_definition.addToken(name='Input.rte.mc.filename[FN_MC_AMBRALS_TYPE]', datatype=io.IOBase),
            parents=['uvspec'], 
            non_parents=brdf_non_parents,
            non_parent_exceptions=['mc_ambrals_type'],
            speaker='rte_solver',
            enable_values=("mystic","montecarlo"),
            threedmystic =True

        )

        mc_rpv_file = option_definition.option(
            name='mc_rpv_file',
            group='surface',
            helpstr='Define a MYSTIC 2D RPV BRDF input file.',
            documentation=documentation['mc_rpv_file'], 
            gui_inputs=(GUI_definition.FileInput(name='Input.rte.mc.filename[FN_MC_RPV]'),),
            tokens=option_definition.addToken(name='Input.rte.mc.filename[FN_MC_RPV]', datatype=io.IOBase),    
            parents=['uvspec'], 
            non_parents=brdf_non_parents,
            non_parent_exceptions=['brdf_cam', 'brdf_cam_solar_wind', 'cox_and_munk_u10_map'],
            speaker='rte_solver',
            enable_values=("mystic","montecarlo"),
            threedmystic =True
        )

        mc_rpv_type = option_definition.option(
            name='mc_rpv_type',
            group='surface', 
            helpstr='Define a MYSTIC 2D rpv spectral BRDF type file.',
            documentation=documentation['mc_rpv_type'],
            gui_inputs=(GUI_definition.FileInput(name='Input.rte.mc.filename[FN_MC_RPV_SPECTRAL]'),),
            tokens=option_definition.addToken(name='Input.rte.mc.filename[FN_MC_RPV_SPECTRAL]',datatype=io.IOBase),    
            parents=['uvspec'], 
            non_parents=brdf_non_parents,
            non_parent_exceptions=['mc_rpv_spectral_file'],
            speaker='rte_solver',
            enable_values=("mystic","montecarlo"),
            threedmystic =True

        )

        mc_rpv_spectral_file = option_definition.option(
            name='mc_rpv_spectral_file',
            group='surface', 
            helpstr='Define a MYSTIC 2D rpv spectral BRDF input file.',
            documentation=documentation['mc_rpv_spectral_file'],
            gui_inputs=(GUI_definition.FileInput(name='Input.rte.mc.filename[FN_MC_RPV2D_TYPE]'),),
            tokens=option_definition.addToken(name='Input.rte.mc.filename[FN_MC_RPV2D_TYPE]', datatype=io.IOBase),
            parents=['uvspec'], 
            non_parents=brdf_non_parents,
            non_parent_exceptions=['mc_rpv_type'],
            speaker='rte_solver',
            enable_values=("mystic","montecarlo"),
            threedmystic =True

        )

        mc_triangular_surface_file = option_definition.option(
            name='mc_triangular_surface_file',
            group='surface',
            helpstr='Define a triangular surface for MYSTIC.',
            documentation=documentation['mc_triangular_surface_file'],
            tokens=option_definition.addToken(name='Input.rte.mc.filename[FN_MC_TRIANGULAR_SURFACE]', datatype=io.IOBase),
            parents=['uvspec'],
            non_parents=brdf_non_parents,
            speaker='rte_solver',
            enable_values=("mystic","montecarlo"),
            threedmystic =True
        )

        bpdf_tsang_u10 = option_definition.option(
            name='bpdf_tsang_u10',
            group='surface', 
            documentation=documentation['bpdf_tsang_u10'],
            gui_inputs=(GUI_definition.FloatInput(name='Input.bpdf.u10'),),
            tokens= [option_definition.addToken(name='Input.bpdf.u10', datatype=float),
                     option_definition.addSetting(name='Input.bpdf.type', setting='BPDF_TSANG')],
            parents=['uvspec'], 
            non_parents=brdf_non_parents,
            mystic=True,
        )

        self.options = [altitude,
                mc_elevation_file,
                albedo, albedo_file, mc_albedo_file,
                mc_albedo_spectral_file, mc_albedo_type, 
                brdf_rossli,brdf_rossli_hotspot,
                brdf_rossli_file,
                mc_rossli_file, 
                brdf_ambrals,brdf_ambrals_hotspot,
                brdf_ambrals_file,
                mc_ambrals_file, 
                brdf_cam, 
                brdf_cam_solar_wind, 
                bpdf_tsang_u10,
                brdf_rpv_library, 
                brdf_rpv, brdf_rpv_type, mc_rpv_type, mc_ambrals_type, 
                brdf_hapke, brdf_hapke_file,
                brdf_rpv_file, mc_rpv_file, mc_rpv_spectral_file,
                mc_ambrals_spectral_file, mc_triangular_surface_file,
                fluorescence, fluorescence_file,
                sur_temperature, sur_temperature_file, 
                ]

    def __iter__(self):
        return iter(self.options)

def get_documentation():
    return get_surface_documentation()

def get_surface_documentation():
    return {
        'brdf_hapke' : r'''

    Wavelength-independent Hapke BRDF values.
        Hapke is a parameterization of extraterrestial solid bodies such as
        the moon, asteroids or the inner planets by \cite{hapke1993}.
    \fcode{
    brdf\_hapke variable value
    }
    \code{variable} can be one of the following:
    \begin{description}
    \item[w] Single scattering albedo of surface.
    \item[B0] Hot spot factor.
    \item[h] Surface roughness parameter.
    \end{description}
        ''',
        'brdf_hapke_file' : r'''
            4 column file, containing the Hapke BDRF
            parameterization, see \code{brdf\_hapke}.
    \fcode{
    brdf\_hapke\_file file
    }
    This option is only supported with
    solvers: \code{disort}\ifmystic{ and \code{rte\_solver
    mystic}}. The columns of the input file are wavelength
    [nm], w, B0, and h.  The parameters are interpolated linearly to
    the internal wavelength grid.  To make sure that the results are
    reasonable, specify the Hapke data on a wavelength grid similar or equal
    to that used internally for the radiative transfer calculation!
        ''',
        'brdf_ambrals' : r'''
    AMBRALS (Ross-Li) BRDF, a three-parameter BRDF fit for vegetated and 
    non-vegetated surfaces \citep{roujean1992,wanner97,lucht2000,
        schaaf2002}, based on the RossThickLiSparseReciprocal
        model. The implementation is according to
        \cite{lucht2000}. See also http://www-modis.bu.edu/brdf/models.html.
    \fcode{
    brdf\_ambrals variable value 
    }
    \code{variable} can be one of the following:
    \begin{description}
        \item[iso] f\_iso, Eq. (37) \cite{lucht2000}.
    \item[vol] f\_vol, Eq. (37) \cite{lucht2000}.
    \item[geo] f\_geo, Eq. (37) \cite{lucht2000}.
    \end{description}
        The hot spot correction factor introduced by \cite{maignan2004}
        can be turned on by the option \code{brdf\_ambrals\_hotspot}.
        May be combined with \ifmystic{\code{mystic}},
    \code{disort}, and \code{fdisort2}.
        ''',
        'brdf_ambrals_hotspot' : r'''
    Turn on hot spot correction factor by \cite{maignan2004} when using the AMBRALS (Ross-Li) BRDF.
        ''',
        'brdf_ambrals_file' : r'''
            4 column file, containing the AMBRALS (Ross-Li) BDRF
    parameterization, see \code{brdf\_ambrals}.
    \fcode{
    brdf\_ambrals\_file file
    }
    This option is only supported with
    solvers: \code{disort}\ifmystic{ and \code{rte\_solver
    mystic}}. The columns of the input file are wavelength
    [nm], iso, vol, and geo.  The parameters are interpolated linearly to
    the internal wavelength grid.  To make sure that the results are
    reasonable, specify the Ambrals data on a wavelength grid similar or equal
    to that used internally for the radiative transfer calculation!
        ''',
        'brdf_rossli' : r'''
        This is the \cite{Lin2015} implementation of the Ross-Li (AMBRALS) BRDF, see also 
        \code{brdf\_ambrals} BRDF. For unknown reasons, these formulas
        differ in their overall factors of the three terms from the factors in
        \cite{lucht2000}. Identical results can be obtained by using: f\_iso=$\pi$ k\_iso, f\_geo=$\pi$ k\_geo, f\_vol = 4 k\_vol / 3.
    \fcode{
    brdf\_rossli variable value 
    }
    \code{variable} can be one of the following:
    \begin{description}
    \item[iso] k\_iso, Eq. (31) \cite{Lin2015}.
    \item[vol] k\_vol, Eq. (31) \cite{Lin2015}.
    \item[geo] k\_geo, Eq. (31) \cite{Lin2015}.
    \end{description}
        The hot spot correction factor introduced by \cite{maignan2004}
        can be turned on by the option \code{brdf\_rossli\_hotspot}.
        May be combined with \ifmystic{\code{mystic}},
    \code{disort}, and \code{fdisort2}.
        ''',
        'brdf_rossli_hotspot' : r'''
    Turn on hot spot correction factor by \cite{maignan2004} when using the Ross-Li (AMBRALS) BRDF.
        ''',
        'brdf_rossli_file' : r'''
        4 column file, containing the Ross-Li (AMBRALS) BDRF
        parameterization, see \code{brdf\_rossli} and \code{brdf\_ambrals\_file}.
        ''',
        'brdf_rpv' : r'''

    Wavelength independent RPV values, see \code{brdf\_rpv\_file}.
    \fcode{
    brdf\_rpv variable value
    }
    \code{variable} can be one of the following:
    \begin{description}
    \item[k] Overwrite the wavelength-dependent
    k value defined in \code{brdf\_rpv\_file}.
    \item[rho0] Overwrite the wavelength-dependent
    rpv0 value defined in \code{brdf\_rpv\_file}.
    \item[theta] Overwrites the wavelength-dependent
    theta value defined in \code{brdf\_rpv\_file}.
    \item[sigma] Constant RPV sigma, to be used for snow \citep{Degunther2000b}.
    \item[t1] Constant RPV t1, to be used for snow \citep{Degunther2000b}.
    \item[t2] Constant RPV t2, to be used for snow \citep{Degunther2000b}.
    \item[scale] Apply a constant scaling factor for the RPV BRDF. 
    Required e.g. if the the albedo should be set to a certain value. This
    factor is only used by \code{rte\_solver disort}, \code{rte\_solver
    fdisort2}\ifmystic{ and \code{rte\_solver mystic}}.
    \end{description}
        ''',

        'bpdf_tsang_u10' : r'''
    Wind speed for ocean BPDF (in m/s) at present only available with 
    \code{rte\_solver mystic}.
    \fcode{
    bpdf\_tsang\_u10 value
    }
    The BPDF model has been developed by \citet{Tsang1985}. 
    The wind speed is the most important parameter affecting the ocean reflectance
    matrix. The BPDF model also takes into account shadowing by surface waves. 
    The model has been implemented in a FORTRAN routine by Mishchenko 
    (http://www.giss.nasa.gov/staff/mmishchenko/brf/)
    which has been included into {\sl libRadtran}.
        ''',


        'mc_ambrals_spectral_file' : r'''
    Define a MYSTIC 2D AMBRALS spectral BRDF input file.
    \fcode{
    mc\_ambrals\_spectral\_file file
    }
        ''',

        'mc_ambrals_type' : r'''
    File containing wavelength-dependent AMBRALSs to be used in combination with \code{mc\_ambrals\_spectral\_file}.
    \fcode{
    mc\_ambrals\_type file
    }
        ''',


        
        'mc_rpv_file' : r'''
    Define a MYSTIC 2D RPV BRDF input file.
    \fcode{
    mc\_rpv\_file file
    }
        ''',

        'mc_rpv_spectral_file' : r'''
    Define a MYSTIC 2D RPV spectral BRDF input file.
    \fcode{
    mc\_rpv\_spectral\_file file
    }
        ''',

        'mc_rpv_type' : r'''
    File containing wavelength-dependent RPVs to be used in combination with \code{mc\_rpv\_spectral\_file}.
    \fcode{
    mc\_rpv\_type file
    }
        ''',

        'mc_rossli_file' : r'''
        Same as \code{mc\_ambrals\_file}, but using the \cite{Lin2015}
        implementation of the Ross-Li (AMBRALS) BRDF
        parameterization. See also \code{brdf\_rossli} and
        \code{brdf\_ambrals}.
        ''',

        'mc_ambrals_file' : r'''
        Define a MYSTIC 2D AMBRALS (Ross-Li) BRDF input file, see \code{brdf\_ambrals}.
    \fcode{
    mc\_ambrals\_file file
    }
        The format of the BRDF file is:
    \fcode{
    Nx  Ny  dx  dy
    ix  iy  iso vol geo
    ...
    }
    where Nx and Ny are the number of grid boxes in x- and y-direction,
    dx and dy are the size of the pixels in km.
    In the second and the following lines the indices in x- and y-direction
    (ix=1...Nx; iy=1...Ny) and the three parameters iso, vol, and geo. 

    Optionally, it is now possibly to mix Ross-Li with Cox and Munk
    (CaM). To this end, the format of the albedo file is:
    \fcode{
    Nx  Ny  dx  dy
    ix  iy  iso vol geo isCaM
    ...
    }
    where isCaM is 1 if CaM shall be used in the pixel, and 0 if Ross-Li
    shall be used. In case of CaM, the values of iso, vol, and geo are
    ignored. Note that this will only work if the CaM parameters are set
    (that is as minimal configuration u10).
        ''',

        'mc_albedo_file' : r'''
    Define a MYSTIC 2D albedo input file.
    \fcode{
    mc\_albedo\_file file
    }
    The format of the albedo file is:
    \fcode{
    Nx  Ny  dx  dy
    ix  iy  albedo
    ...
    }
    where Nx and Ny are the number of grid boxes in x- and y-direction,
    dx and dy are the size of the pixels in km.
    In the second and the following lines the indices in x- and y-direction
    (ix=1...Nx; iy=1...Ny) and the albedo of the pixel are specified. 
        ''',

        'mc_albedo_spectral_file' : r'''
    Define a MYSTIC 2D spectral albedo input file.
    \fcode{
    mc\_albedo\_spectral\_file file
    }
        ''',

        'mc_albedo_type' : r'''
    File containing wavelength-dependent albedos to be used in combination with \code{mc\_albedo\_spectral\_file}.
    \fcode{
    mc\_albedo\_type file
    }
        ''',

        'mc_elevation_file' : r'''
    Define a MYSTIC 2D elevation input file. 
    \fcode{
    mc\_elevation\_file file
    }
    The 
    expected format of the elevation file is:
    \fcode{
    Nx  Ny  dx  dy
    ix  iy  elevation
    ...
    }
    where Nx and Ny are the number of grid boxes in x- and y-direction,
    dx and dy are the size of the grid boxes in km.
    In the second and the following lines the indices in x- and y-direction
    (ix=1...Nx, iy=1...Ny) and the elevation in km of each point are specified.

    Attention: While the other files refer to grid boxes,
    the elevation is defined at grid points. It has to be this 
    way because each "elevation pixel" contains a surface 
    which is defined by the four cornes of the pixel. If the
    grid covers an area of 200x200 km$^2$ and the pixel sizes 
    are dx = 1km and dy = 1km, the elevation has to be defined 
    at 201x201 points (ix = 1 .. Nx, iy = 1 .. Ny).

    Also, make sure that the elevation doesn't cross or 
    even touch the 1D surface. Since the default altitude is 0
    that means that the minimum elevation must be larger than 0!
        ''',

        'altitude' : r'''
    Set the bottom level in the model atmosphere provided in 
    \code{atmosphere\_file} to be at the given altitude above sea level (km).
    \fcode{
    altitude 0.73   \# Altitude of IFU, Garmisch-Partenkirchen 
    }
    The profiles of pressure, temperature, molecular absorbers, 
    ice and water clouds are cut at the specified altitude. 
    The aerosol profile is not affected by \code{altitude} but starts
    right from the model surface. This is a convenient way for the user to calculate the 
    radiation at other altitudes than sealevel. Note that \code{altitude} is very different
    from \code{zout} where the radiation is calculated at an altitude of zout
    above the surface. E.g. to calculate the radiation field 1 km above the surface 
    at a location at 0.73 km above sealevel, one would specify '\code{altitude} 0.73' 
    and '\code{zout} 1.0'.
    If an altitude is specified which is below the lowest level in the 
    \code{atmosphere\_file}, the atmospheric profiles are extrapolated assuming a constant 
    gradient for temperature and mixing ratios.
    A second optional argument may be given to \code{altitude} as e.g.
    \fcode{
    altitude 0.73 0.5
    }
    Here the bottom level will be at 0.73 km and the vertical resolution
    of the model atmosphere will be redistributed to have a spacing
    between levels specified by the second number, here 0.5 km, starting
    however from 0km. (Levels 0.73, 1., 1.5 ... will be added to the
    original atmosphere grid and optical properties are devided into the
    new layers. In order to use interpolated properties use
    \code{zout\_interpolate}.  See verbose output for details.)  Be aware
    that specifying a fine vertical spacing will produce many layers thus
    increasing the computing time. Also the radiative transfer equation
    solvers implemented in Fortran 77 might need to have some array sizes
    increased (see \file{src\_f/DISORT.MXD}).
        ''',

        'mc_triangular_surface_file' : r'''
    Define a triangular surface for MYSTIC. 
    \fcode{
        mc\_triangular\_surface\_file file
    }
    The file includes the description of a trianglar surfaces.
    It is expected in  {\sl netcdf}-Format and
    should include the following fields: 
    \fcode{vertices (n_vert,3)} Location of vertices (x,y,z)
    \fcode{triangles (n_tris,3)} Indices of vertices for each triangle (three
    corners)
    \fcode{albedo (n_tris)} Surface albedo for each triangle. 
        ''',

        'albedo' : r'''
    The Lambertian surface albedo
    \fcode{
        albedo value
    }
    where \code{value} is a number between 0.0 and 1.0,
    constant for all wavelengths. For wavelength dependent surface albedo 
    use \code{albedo\_file}. The default albedo is 0.0.
        ''',

        'albedo_file' : r'''
    Location of surface albedo file for wavelength dependent surface albedo. 
    \fcode{
    albedo\_file file
    }
    The file must have two columns.
    Column 1 is the wavelength in nm, and column 2 the corresponding 
    Lambertian surface albedo. An arbitrary wavelength grid may be chosen as the 
    albedo will be interpolated linearely to the wavelength grid used for the 
    radiation calculation. Comments start with \code{\#}. Empty lines are ignored.
    A large collection of spectral albedos are available e.g. at http://speclib.jpl.nasa.gov/ 
    \citep{Baldridge2009}.
        ''',

        'fluorescence' : r'''
    Specifies the magnitude of a bottom surface isotropic fluorescence source.
    \fcode{
        fluorescence value
    }
    where \code{value} is a number greater or equal to 0.0, 
    constant for all wavelengths.
    Must be used together with  \code{source solar file}. The units of the fluorescence should 
    obviously be the same as for the solar source in \code{source solar file}.
    For wavelength dependent fluorescence use \code{fluorescence\_file}. The default fluorescence is 0.0. 
    Currently only works with the \code{disort} solver.
        ''',

        'fluorescence_file' : r'''
    Location of fluorescence file for wavelength dependent fluorescence emission 
    from the bottom surface. 
    \fcode{
        fluorescence\_file file
    }
    The file must have two columns.
    Column 1 is the wavelength in nm, and column 2 the corresponding 
    fluorescence. An arbitrary wavelength grid may be chosen as the 
    fluorescence will be interpolated linearely to the wavelength grid used for the 
    radiation calculation. Comments start with \code{\#}. Empty lines are ignored.
    Currently only works with the \code{disort} solver. Furthermore. if \code{raman} is 
    not set, \code{wavelength\_grid\_file} 
    must be specified with the same resolution as the \code{source solar file}, and the first value
    must be the value specified by \code{wavelength}. The units of the fluorescence should 
    obviously be the same as for the \code{source solar file}.
        ''',

        'brdf_rpv_type' : r'''
    With this option the (RPV) BRDF surface type is selected.
    This option can be used with \code{albedo\_library} in order to select a spectral albedo 
    or with \code{brdf\_rpv\_library} in order to select a BRDF function.
    \fcode{
    brdf\_rpv\_type surface\_type\_number
    }
    where \code{surface\_type\_number} is an integer starting from 0,
    where 0 refers to a black surface and the following numbers to 
    the entries in the specified \code{library}.
        ''',

        'brdf_rpv_file' : r'''
    4 to 7 column file, containing the Rahman, Pinty, and Verstraete (RPV) BDRF
    parameterization \citep{Rahman1993} and the snow extension by \citep{Degunther2000b}.
    \fcode{
    brdf\_rpv\_file file
    }
    Bidirectional reflectance distribution functions for a variety of
    surfaces are given in the paper. This option is only supported with
    solvers: \code{disort}, \code{fdisort2}\ifmystic{ and \code{rte\_solver
    mystic}}. The columns of the input file are wavelength
    [nm], rho0, k, and theta.  The parameters are interpolated linearly to
    the internal wavelength grid.  To make sure that the results are
    reasonable, specify the RPV data on a wavelength grid similar or equal
    to that used internally for the radiative transfer calculation!
    Optionally, a fifth column with a constant scaling factor may be
    defined.  If it has seven columns, the fifth to seventh are sigma, t1,
    t2, and if it has eight, the eighth is scale again.
        ''',

        'brdf_rpv_library' : r'''
    The rpv libraries are collections of spectral BRDFs of different surface types, 
    This option must be used either with \code{brdf\_rpv\_type} 
    or \code{surface\_type\_map}, in order to select the specific surface type.
                
    For using a \code{brdf\_rpv\_library} write 
    \fcode{
    brdf\_rpv\_library library\_path 
    }
    where \code{library\_path} is the path of the directory, where the BRDF data is stored. The files are 
    expected to have the names \code{IGBP.01.rpv, IGBP.02.rpv, ...} If \code{brdf\_rpv\_type 1} is specified
    the BRDF from \code{IGBP.01.rpv} will be used, and so on.
    Each file must have the structure like an \code{brdf\_rpv\_file}.
    (This option is quite the same as \code{brdf\_rpv\_file}, except that it offers you an
    easy way to use the option \code{surface\_type\_map} in combinition with your \code{brdf\_rpv\_files}.)

    %ak \emph{For our group there is also a built-in library for IGBP surface types.}
    \fcode{
    brdf\_rpv\_library IGBP
    }
    The built-in library contains the first 17 surface types see
    \code{albedo\_library}.  The data is given for the wavelengths 443nm,
    565nm, 670nm, and 865nm. Stay near this wavelength in order to get
    reasonable results. In future this the rpv-library will be NDVI
    dependent, but until now the most common NDVI class is selected
    automatically.
        ''',

        'brdf_cam' : r'''
    Set \citet{cox54a,cox54b} ocean BRDF properties.
    \fcode{
    brdf\_cam variable value
    }
    \code{variable} can be one of the following:
    \begin{description}
    \parameter{pcl} 
    Pigment concentration for \citet{cox54a,cox54b} ocean BRDF (in mg/m$^{-3}$).
    The default value is 0.01 mg/m$^{-3}$. 
    \parameter{sal}
    Salinity for \citet{cox54a,cox54b} ocean BRDF (in "per mille", 0.1\%; this 
        unit is equivalent to the other common units for salinity, 
        ppt - parts per thousand, psu - practical salinity unit).
    The default value is 34.3. 
    \parameter{u10}
    Wind speed for \citet{cox54a,cox54b} ocean BRDF (in m/s).
    The wind speed is the most important parameter affecting ocean BRDF.  The
        minimum allowed wind speed is 1 m/s because otherwise the strong
        specular reflection causes numerical problems. If a lower value is
        specified, the wind speed is automatically set to 1m/s.
    \parameter{uphi}
    Wind direction for \citet{cox54a,cox54b} ocean BRDF.
    Default value is 0 degrees,
        which is wind from the South. 90 degrees corresponds to wind from the West, etc. (Honestly, this was
        never truly validated. It could possibly be that 0 is wind from the North, 90 is wind from the East, etc.)
    The option is only implemented for the \code{mystic} solver. 
        \end{description}

    At present only available with \code{rte\_solver disort}, 
        \ifmystic{\code{rte\_solver mystic}}
        and \code{rte\_solver fdisort2}.  The
        number of streams (\code{number\_of\_streams}) is automatically increased to 16 if
        \code{cox\_and\_munk} BRDF is switched on, to avoid numerical
        problems.
    To switch on Cox and Munk BRDF,
        specify any of the \code{brdf\_cam} options and define at least
        \code{brdf\_cam u10}.
        ''',

        'brdf_cam_solar_wind' : r'''
    Use old definition of wind direction for Monte Carlo simulations. If this switch is set,
    the wind azimuth is identical to the incoming photon azimuth. Else, the wind azimuth is set by
    \code{brdf\_cam uphi} or is 0 by default.
        ''',

        'sur_temperature' : r'''
    Surface temperature, used for thermal infrared calculations. 
    \fcode{
    sur\_temperature value
    }
    If not specified, the
    temperature of the lowest atmospheric level is used as surface temperature.
        ''',

        'sur_temperature_file' : r'''
    Define a MYSTIC 2D temperature input file.
    \fcode{
    sur\_temperature\_file file
    }
    The expected format of the temperature file is:
        \fcode{
        Nx  Ny  dx  dy
        ix  iy  temperature
        }
    where Nx and Ny are the number of grid boxes in x- and y-direction,
    dx and dy are the size of the grid boxes in km.
    In the second and the following lines the indices in x- and y-direction
    and the temperature of the pixel are specified. 
        ''',
    }
