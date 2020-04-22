__author__ = 'CJ Willers'

"""Simple testing

These tests are very light and easy testing of basic set_option operation.
The functionality behind the option settings is NOT TESTED.
The tests simply invoke a set_option call for all the options found.
Some set_options do not work and these are simply commented out

Please consult with the original authors of the option_definition code
for more information and how to use it.

"""



import glob
import io

import libraddask.rad.librad as librad
import libraddask.libradtran.option_definition as option_definition
import libraddask.libradtran.aerosol_options as aerosol_options
import libraddask.libradtran.cloud_options as cloud_options
import libraddask.libradtran.general_atmosphere_options as general_atmosphere_options
import libraddask.libradtran.geometry_options as geometry_options
import libraddask.libradtran.mc_options as mc_options
import libraddask.libradtran.molecular_options as molecular_options
import libraddask.libradtran.output_options as output_options
import libraddask.libradtran.profile_options as profile_options
import libraddask.libradtran.solver_options as solver_options
import libraddask.libradtran.special_options as special_options
import libraddask.libradtran.spectral_options as spectral_options
import libraddask.libradtran.surface_options as surface_options


doAllTests = True

#==================================================
if doAllTests:
    # Find and read each uvspec example case in the examples folder
    # this test reads an input file and the existing output file
    # it does not run libRadtran, only read existing files.
    uvINPfiles = glob.glob('examples/UVSPEC_*.INP')
    for uvINPfile in uvINPfiles:
        print('Processing ' + uvINPfile)
        uvTestCase = librad.Case(filename=uvINPfile)  # Read the uvspec input file
        print(f'Outfile: {uvTestCase.outfile}')
        try:
            # Read the corresponding (existing) uvspec output file
            uvTestCase.readout()  
            # the data is fluxline format
            for key in uvTestCase.fluxline:
                print(key)
                print(getattr(uvTestCase, key).T)
        except:
            print('Reading of output failed for ' + uvINPfile)


#==================================================
if doAllTests:
    case = librad.Case(casename='TestCase')
    print(50*'='+'\n'+50*'=')
    print('aerosol_options')
    print(f'\nOptions documented:')
    for key in aerosol_options.get_documentation().keys():
        print(key)
    print('\nOptions:')
    for fn in aerosol_options.setup_aerosol_group():
        print(fn.name)
    print('\nTesting set_option:')
    case.set_option('aerosol_default',1)
    case.set_option('aerosol_file','gg','aerosol_file',1)
    case.set_option('aerosol_profile_modtran',True)
    case.set_option('aerosol_angstrom',10.,10,1)
    case.set_option('aerosol_king_byrne',1.,1.,1.,1)
    case.set_option('aerosol_modify','gg','scale',1,1)
    case.set_option('aerosol_haze',1)
    case.set_option('aerosol_set_tau_at_wvl',1.,1.,'Input.aer.tau_wvl_tau',1)
    case.set_option('aerosol_season',10,1)
    case.set_option('aerosol_species_file','continental_clean',1,1)
    case.set_option('aerosol_species_library','OPAC')
    case.set_option('aerosol_visibility',10.)
    case.set_option('aerosol_sizedist_file','aerosol_sizedist_file')
    case.set_option('aerosol_refrac_index',3.,5.,1)
    print(f'Case contents:\n{case}')

#==================================================
if doAllTests:
    case = librad.Case(casename='TestCase')
    print(50*'='+'\n'+50*'=')
    print('cloud_options')
    print(f'\nOptions documented:')
    for key in cloud_options.get_documentation().keys():
        print(key)
    print('\nOptions:')
    for fn in cloud_options.setup_cloud_group():
        print(fn.name)
    print('\nTesting set_option:')
    case.set_option('wc_file',option_definition.CaothType('wc'),'CAOTH_FROM_','str')
    case.set_option('wc_properties',option_definition.CaothType('wc'),'echam4','interpolate')
    case.set_option('wc_modify',option_definition.CaothType('wc'),'tau','scale',45.)
    case.set_option('ic_file',option_definition.CaothType('ic'),'CAOTH_FROM_',io.IOBase)
    case.set_option('ic_properties',option_definition.CaothType('ic'),'yang','interpolate')
    case.set_option('ic_modify',option_definition.CaothType('ic'),'tau550','scale',5.)
    case.set_option('ic_habit',option_definition.CaothType('ic'),'hollow-column')
    case.set_option('ic_habit_yang2013',option_definition.CaothType('ic'),'column_8elements','severe')
    case.set_option('ic_fu',option_definition.CaothType('ic'),'reff_def','off','id2')
    case.set_option('ic_raytracing_file',option_definition.CaothType('ic'),io.IOBase)
    case.set_option('cloud_fraction_file',io.IOBase)
    case.set_option('cloud_overlap','maxrand')
    case.set_option('cloudcover',option_definition.CaothType('ic'),4.,1,1)
    print(f'Case contents:\n{case}')


#==================================================
if doAllTests:
    case = librad.Case(casename='TestCase')
    print(50*'='+'\n'+50*'=')
    print('general_atmosphere_options')
    print(f'\nOptions documented:')
    for key in general_atmosphere_options.get_documentation().keys():
        print(key)
    print('\nOptions:')
    for fn in general_atmosphere_options.setup_general_atm_group():
        print(fn.name)
    print('\nTesting set_option:')
    # case.set_option('atmos_region',5,6,7,8)
    case.set_option('no_absorption','?????',option_definition.CaothoffType,1)
    case.set_option('no_scattering','?????',option_definition.CaothoffType,1)
    case.set_option('interpret_as_level',option_definition.CaothType,'FALSE')
    case.set_option('zout_interpolate','ZOUT_INTERPOLATE')
    case.set_option('z_interpolate')
    case.set_option('atm_z_grid',-5.,'ntokens')
    case.set_option('reverse_atmosphere',0)
    print(f'Case contents:\n{case}')


#==================================================
if doAllTests:
    case = librad.Case(casename='TestCase')
    print(50*'='+'\n'+50*'=')
    print('geometry_options')
    print(f'\nOptions documented:')
    for key in geometry_options.get_documentation().keys():
        print(key)
    print('\nOptions:')
    for fn in geometry_options.setup_geometry_group():
        print(fn.name)
    print('\nTesting set_option:')
    case.set_option('sza',4.,'SZA_DIRECT_INPUT')
    case.set_option('sza_file',io.IOBase,'SZA_BY_TIME_AND_LOCATION')
    case.set_option('phi0',4.)
    case.set_option('phi',-5.,6.,7.)
    case.set_option('umu',-5.,6.,7.)
    # case.set_option('mc_bw_umu_file',io.IOBase,True)
    case.set_option('earth_radius',7.)
    case.set_option('latitude','LATITUDE_SIGNUM_',5.,6.,7.)
    case.set_option('longitude','LONGITUDE_SIGNUM_',5.,6.,7.)
    case.set_option('day_of_year',7)
    case.set_option('mc_bcond','periodic')
    # case.set_option('mc_sample_grid',4,5,6.,7.)
    case.set_option('mc_spherical','1D',1,'1d')
    # case.set_option('mc_spherical3D_scene',1,5.,6.,.7,.8)
    case.set_option('mc_satellite_view','pointer', 1.,2.,3.,4.,'PAN_MODE_SATELLITE',True)
    # case.set_option('mc_satellite_position','SZA_DIRECT_INPUT',4.,5.)
    case.set_option('mc_sensordirection',1.,2.,3.,4)
    # case.set_option('mc_sensorposition',1.,2.,3.,'MC_SENSORPOSITION_CARTESIAN')
    case.set_option('mc_reference_to_nn',1)
    case.set_option('mc_panorama_forward','pointer',0,360,0,180,1,True)
    # case.set_option('mc_panorama_view','pointer',0.,360.,0.,180.,'PAN_MODE_CAMERA',True)
    case.set_option('mc_panorama_alignment','mu')
    # case.set_option('mc_panorama','no_pixel',1)
    # case.set_option('mc_blitz_position','pointer',0.,1.,2.,3.,4.,5.,'SRC_BLITZ',True)
    print(f'Case contents:\n{case}')


#==================================================
if doAllTests:
    case = librad.Case(casename='TestCase')
    print(50*'='+'\n'+50*'=')
    print('mc_options')
    print(f'\nOptions documented:')
    for key in mc_options.get_documentation().keys():
        print(key)
    print('\nOptions:')
    for fn in mc_options.setup_mc_group():
        print(fn.name)
    print('\nTesting set_option:')
    case.set_option('mc_escape','str')
    case.set_option('mc_core',1,'???','???','CAOTH_FROM_3D','filename','PROP_HU','str')
    case.set_option('mc_core_threshold',6.)
    case.set_option('mc_core_scale',-1)
    case.set_option('mc_core_savefile','str')
    case.set_option('mc_vroom',1,'on')
    # case.set_option('mc_visualize',1,'hiddenline')
    case.set_option('mc_truncate',1,4.)
    case.set_option('mc_backward',1,2,3,4,5)
    # case.set_option('mc_sample_cldprp',1,1)
    case.set_option('mc_surface_reflectalways',1)
    # case.set_option('mc_DoLE',1)
    case.set_option('mc_spectral_is',True,3.)
    # case.set_option('mc_aerosol_is',io.IOBase,True)
    case.set_option('mc_boxairmass','1D',1)
    case.set_option('mc_albedo_spectral',1)
    case.set_option('mc_azimuth_old','MCAZIMUTH_OLD')
    # case.set_option('mc_backward_heat','EMABS')
    case.set_option('mc_backward_sunshape_file',io.IOBase)
    case.set_option('mc_backward_writeback',1)
    case.set_option('mc_coherent_backscatter',1,2,3)
    # case.set_option('mc_delta_scaling',0.99,0,4,5)
    case.set_option('mc_maxscatters',3)
    case.set_option('mc_minscatters',1)
    case.set_option('mc_minphotons',4)
    case.set_option('mc_photons',4)
    case.set_option('mc_photons_file',io.IOBase)
    case.set_option('mc_polarisation',1,1,3)
    case.set_option('mc_rad_alpha',5.)
    # case.set_option('mc_radial_pathlength',5,6)
    case.set_option('mc_radial_pathlength_dt',.1)
    # case.set_option('mc_readrandomstatus',1,1)
    case.set_option('mc_randomseed',5)
    case.set_option('mc_refraction',1)
    # case.set_option('mc_ris','MC_RIS_',5.)
    print(f'Case contents:\n{case}')


#==================================================
if doAllTests:
    case = librad.Case(casename='TestCase')
    print(50*'='+'\n'+50*'=')
    print('molecular_options')
    print(f'\nOptions documented:')
    for key in molecular_options.get_documentation().keys():
        print(key)
    print('\nOptions:')
    for fn in molecular_options.setup_molecular_group():
        print(fn.name)
    print('\nTesting set_option:')
    case.set_option('atmosphere_file','str')
    # case.set_option('atmosphere_file_3D','str',True,'???','CAOTH_FROM_3D','"profile_mol3d_dummy.dat"','PROP_HU')
    case.set_option('radiosonde')
    case.set_option('radiosonde_levels_only','FALSE')
    case.set_option('mol_file','CO2',io.IOBase,'cm_3')
    case.set_option('pressure',1000)
    case.set_option('refractive_index_file','str')
    case.set_option('crs_model','o3','Bodhaine')
    case.set_option('crs_file','O3','O3')
    case.set_option('rayleigh_depol',.4)
    case.set_option('mol_abs_param','kato','str')
    case.set_option('reptran_file',io.IOBase)
    case.set_option('ck_lowtran_absorption','O4','on')
    case.set_option('ck_fu_h2o_continuum','on')
    case.set_option('mol_tau_file','FN_MOL_TAU_',io.IOBase)
    case.set_option('mol_modify','O2',5.,'CM_2')
    case.set_option('mixing_ratio','CO2',5.)
    print(f'Case contents:\n{case}')


#==================================================
if doAllTests:
    case = librad.Case(casename='TestCase')
    print(50*'='+'\n'+50*'=')
    print('output_options')
    print(f'\nOptions documented:')
    for key in output_options.get_documentation().keys():
        print(key)
    print('\nOptions:')
    for fn in output_options.setup_output_group():
        print(fn.name)
    print('\nTesting set_option:')
    case.set_option('quiet',1)
    case.set_option('verbose',1)
    # case.set_option('write_ext_to_file',1)
    case.set_option('test_optical_properties',1)
    case.set_option('print_disort_info',1)
    case.set_option('data_files_path',io.IOBase)
    case.set_option('output_process','sum')
    case.set_option('output_file','str')
    case.set_option('output_format','netCDF')
    case.set_option('output_user','str')
    case.set_option('output_quantity','brightness')
    case.set_option('heating_rate','HEAT_LAYER_CD','local')
    # case.set_option('write_output_as_netcdf',)
    case.set_option('slit_function_file',io.IOBase,1)
    case.set_option('spline',1.,1.,1.,1)
    case.set_option('spline_file',io.IOBase,1)
    case.set_option('filter_function_file',io.IOBase,'normalize')
    case.set_option('pressure_out')
    case.set_option('zout','str')
    case.set_option('zout_sea','str')
    case.set_option('mc_backward_output','edir','W_per_m3')
    case.set_option('mc_forward_output','absorption','W_per_m3')
    case.set_option('mc_basename',io.IOBase)
    # case.set_option('mc_backward_writeallpixels',1)
    case.set_option('mc_std',1)
    # case.set_option('mc_surfaceparallel',1)
    case.set_option('mc_jacobian',1)
    # case.set_option('mc_jacobian_std',1)
    print(f'Case contents:\n{case}')

#==================================================
if doAllTests:
    case = librad.Case(casename='TestCase')
    print(50*'='+'\n'+50*'=')
    print('profile_options')
    print(f'\nOptions documented:')
    for key in profile_options.get_documentation().keys():
        print(key)
    print('\nOptions:')
    for fn in profile_options.setup_cloud_group():
        print(fn.name)
    print('\nTesting set_option:')    
    case.set_option('profile_file',option_definition.CaothType,'???',io.IOBase)
    case.set_option('profile_properties',option_definition.CaothType,'hu','interpolate')
    case.set_option('profile_modify',option_definition.CaothType,'gg','set',5.)
    print(f'Case contents:\n{case}')


#==================================================
if doAllTests:
    case = librad.Case(casename='TestCase')
    print(50*'='+'\n'+50*'=')
    print('solver_options')
    print(f'\nOptions documented:')
    for key in solver_options.get_documentation().keys():
        print(key)
    print('\nOptions:')
    for fn in solver_options.setup_solver_group():
        print(fn.name)
    print('\nTesting set_option:')    
    case.set_option('rte_solver','disort')
    case.set_option('pseudospherical',1)
    case.set_option('disort_spherical_albedo',1,0)
    case.set_option('disort_intcor','phase')
    case.set_option('isotropic_source_toa',1)
    case.set_option('raman','RAMAN_CALC','CK_RAMAN','SOLVER_DISORT','PROCESS_RAMAN','original')
    case.set_option('number_of_streams',4)
    case.set_option('polradtran','aziorder',1)
    case.set_option('polradtran_quad_type','G')
    case.set_option('polradtran_max_delta_tau',5.)
    case.set_option('sdisort','ichapman',4)
    case.set_option('sos_nscat',5)
    case.set_option('deltam','on')
    # case.set_option('ipa_3d',1)
    # case.set_option('mc_tenstream',1,'str')
    case.set_option('mc_ipa',1)
    case.set_option('mc_tipa',1,'dirdiff')
    case.set_option('tipa','dirdiff')
    case.set_option('sslidar','area',5.)
    case.set_option('sslidar_nranges',5)
    case.set_option('sslidar_polarisation',1)
    case.set_option('tzs_cloud_top_height','str')
    case.set_option('mc_nca',1,'str')
    print(f'Case contents:\n{case}')

#==================================================
if doAllTests:
    case = librad.Case(casename='TestCase')
    print(50*'='+'\n'+50*'=')
    print('special_options')
    print(f'\nOptions documented:')
    for key in special_options.get_documentation().keys():
        print(key)
    print('\nOptions:')
    for fn in special_options.setup_special_group():
        print(fn.name)
    print('\nTesting set_option:')    
    case.set_option('albedo_library','str')
    case.set_option('cloud_fraction_map',io.IOBase)
    case.set_option('albedo_map','ALBEDO_FROM_ALBEDO_MAP',io.IOBase,io.IOBase,'str')
    case.set_option('altitude_map','ALT_FROM_MAP',io.IOBase,io.IOBase,5.)
    case.set_option('cox_and_munk_u10_map',io.IOBase,'BRDF_NONE')
    case.set_option('cox_and_munk_pcl_map','BRDF_CAM',io.IOBase,io.IOBase,5.)
    case.set_option('cox_and_munk_sal_map','BRDF_CAM',io.IOBase,io.IOBase,5.)
    case.set_option('emissivity_map','ALBEDO_FROM_EMISSIVITY_MAP',io.IOBase,io.IOBase,5.)
    case.set_option('surface_temperature_map',io.IOBase,io.IOBase)
    case.set_option('surface_type_map','TRUE',io.IOBase,io.IOBase)
    case.set_option('ECHAM_sza',io.IOBase,'NOT_DEFINED_FLOAT','SZA_ECHAM')
    case.set_option('ECMWF_ozone_climatology','TRUE')
    case.set_option('ECMWF_wc_file',option_definition.CaothType('wc'),'CAOTH_FROM_ECMWF',io.IOBase)
    case.set_option('ECMWF_ic_file',option_definition.CaothType('ic'),'CAOTH_FROM_ECMWF',io.IOBase)
    case.set_option('ECMWF_levels_only','FALSE')
    case.set_option('ECMWF_wind_file',io.IOBase)
    case.set_option('satellite_geometry',io.IOBase)
    case.set_option('satellite_pixel',5,5)
    print(f'Case contents:\n{case}')


#==================================================
if doAllTests:
    case = librad.Case(casename='TestCase')
    print(50*'='+'\n'+50*'=')
    print('spectral_options')
    print(f'\nOptions documented:')
    for key in spectral_options.get_documentation().keys():
        print(key)
    print('\nOptions:')
    for fn in spectral_options.setup_spectral_group():
        print(fn.name)
    print('\nTesting set_option:')
    case.set_option('wavelength',500.,600.)
    # case.set_option('wavelength_step')
    case.set_option('wavelength_index',5,6)
    case.set_option('wavelength_grid_file','wavelength_grid_file')
    case.set_option('thermal_bands_file','thermal_bands_file')
    case.set_option('thermal_bandwidth',500,'nm')
    case.set_option('source','thermal','source','per_cm-1')
    case.set_option('mc_sun_angular_size',3.14)
    # case.set_option('mc_lidar','SRC_LIDAR',0,'MCLIDAR_SPACE','polarize')
    case.set_option('mc_lidar_file','mc_lidar_file')
    # case.set_option('mc_radar','SRC_LIDAR',0,'MCRADAR')
    print(f'Case contents:\n{case}')


#==================================================
if doAllTests:
    case = librad.Case(casename='TestCase')
    print(50*'='+'\n'+50*'=')
    print('surface_options')
    print(f'\nOptions documented:')
    for key in surface_options.get_documentation().keys():
        print(key)
    print('\nOptions:')
    for fn in surface_options.setup_surface_group():
        print(fn.name)
    print('\nTesting set_option:')    
    # case.set_option('altitude','ALT_NOT_DEFINED',3.,4.)
    case.set_option('altitude',3.,4.)
    # case.set_option('mc_elevation_file',io.IOBase)
    case.set_option('albedo',.5,'NOT_DEFINED_INTEGER')
    case.set_option('albedo_file',io.IOBase,False,1)
    case.set_option('brdf_cam','pcl',6.,'BRDF_CAM')
    case.set_option('brdf_cam_solar_wind',1)
    case.set_option('brdf_hapke','b0',5.)
    case.set_option('brdf_hapke_file',io.IOBase,'HAPKE_FROM_HAPKE_FILE')
    case.set_option('brdf_rossli','iso',4.)
    case.set_option('brdf_rossli_hotspot','BRDF_ROSSLI_HOTSPOT_ON')
    case.set_option('brdf_rossli_file',io.IOBase,'ROSSLI_FROM_ROSSLI_FILE')
    case.set_option('brdf_ambrals','ROSSLI_AMBRALS_CONSTANT','vol',5.)
    case.set_option('brdf_ambrals_hotspot','BRDF_ROSSLI_HOTSPOT_ON')
    case.set_option('brdf_ambrals_file',io.IOBase,'ROSSLI_FROM_AMBRALS_FILE')
    case.set_option('fluorescence',5.,'FLUORESCENCE_CONSTANT')
    case.set_option('fluorescence_file',io.IOBase,'FLUORESCENCE_FROM_FLUORESCENCE_FILE')
    case.set_option('brdf_rpv_file',io.IOBase,'RPV_FROM_RPV_FILE')
    case.set_option('brdf_rpv','sigma',4.)
    case.set_option('brdf_rpv_library','str')
    case.set_option('brdf_rpv_type',5,False)
    case.set_option('sur_temperature',5900.)
    # case.set_option('sur_temperature_file',io.IOBase)
    # case.set_option('mc_albedo_file',io.IOBase)
    # case.set_option('mc_albedo_type',io.IOBase)
    case.set_option('mc_albedo_spectral_file',io.IOBase)
    # case.set_option('mc_rossli_file',io.IOBase)
    # case.set_option('mc_ambrals_file',io.IOBase)
    case.set_option('mc_ambrals_type',io.IOBase)
    # case.set_option('mc_ambrals_spectral_file',io.IOBase)
    # case.set_option('mc_rpv_file',io.IOBase)
    case.set_option('mc_rpv_spectral_file',io.IOBase)
    case.set_option('mc_triangular_surface_file',io.IOBase)
    case.set_option('bpdf_tsang_u10',5.,'BPDF_TSANG')
    print(f'Case contents:\n{case}')

