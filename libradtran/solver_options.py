"""--------------------------------------------------------------------
 * $Id: solver_options.py 3166 2015-08-20 13:51:47Z Claudia.Emde $
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


from . import option_definition 
from . import GUI_definition

class setup_solver_group():
	
	group_name = 'Solver'

	def __init__(self):
		documentation = get_solver_documentation()

		rte_solver = option_definition.option(
			name='rte_solver',
			group='solver',
			helpstr='Radiative transfer equation solver.',
			documentation=documentation['rte_solver'],
			gui_inputs=(GUI_definition.ListInput(name='solver', valid_range=['disort', 'twostr', 'mystic', 'rodents', 'sslidar', 'null'], optional=False),),
			tokens=option_definition.addLogical( name='Input.rte.solver', 
						logicals = [ 'disort', 'twostr', 'mystic', 'rodents', 'sslidar', 'null',
							'sdisort', 'spsdisort', 'fdisort1', 'fdisort2', 'sos', 'polradtran', 'ftwostr', 'montecarlo', 'tzs', 'sssi', 'sss', 'twostrebe' ],
						setting='SOLVER_' ),
			parents=['uvspec'],
			childs=['polradtran_quad_type','polradtran','polradtran_max_delta_tau',
			 'pseudospherical', 'disort_spherical_albedo',
			 'disort_intcor',
                         'isotropic_source_toa','raman',
			 'number_of_streams','deltam','nrefrac','brdf_cam','brdf_cam_solar_wind',
			 # Mystic options are not enabled unless solver is mystic of monte carlo, see mc_options.py
			 'mc_escape','mc_vroom','mc_backward','mc_spectral_is',
			 'mc_spectral_is_wvl','mc_forward_output',
			 'mc_albedo_file','mc_albedo_spectral_file','mc_albedo_type','mc_rossli_file',
			 'mc_backward_output', 'mc_backward_writeallpixels', 
			 'mc_basename',
			 'mc_delta_scaling', 'mc_elevation_file', 
                         'mc_tenstream',
			 'mc_ipa','mc_tipa','ipa_3d','tipa',
			 'mc_maxscatters','mc_minscatters','mc_minphotons','mc_photons', 
			 'mc_photons_file','mc_polarisation',
			 'mc_wcloud_file','mc_icloud_file',
			 'mc_randomseed', 'mc_rpv_file','mc_rpv_type','mc_spherical','mc_std',
			 'sur_temperature_file','mc_sample_grid','mc_sensordirection','mc_sensorposition',
			 'mc_reference_to_nn','mc_bw_umu_file','mc_truncate','mc_bcond','mc_coherent_backscatter',
			 'mc_surface_reflectalways','mc_refraction','mc_visualize','mc_hiddenline', 
			 'mc_surfaceparallel', 'mc_jacobian', 'mc_jacobian_std',
                         'sslidar','sslidar_nranges','sslidar_polarisation',
                         'tzs_cloud_top_height'       
                        ],
			continious_update=True,
		)

		pseudospherical = option_definition.option(
			name='pseudospherical',
			group='solver',
			helpstr='Invokes pseudospherical geometry for disort/twostr solver.',
			documentation=documentation['pseudospherical'],
			tokens=option_definition.addSetting(name='Input.rte.pseudospherical', setting=1),
			parents=['rte_solver'],
			speaker='rte_solver',
			enable_values=("disort","twostr",)
		)

		disort_spherical_albedo = option_definition.option(
			name='disort_spherical_albedo',
			group='solver',
			helpstr='Calculate spherical albedo.',
			documentation=documentation['disort_spherical_albedo'],
			tokens=[option_definition.addSetting(name='Input.rte.ibcnd', setting=1),
				  option_definition.addSetting(name='Input.rte.nphi', setting=0)],
			parents=['rte_solver'],
			speaker='rte_solver',
			enable_values=("disort",),
                        developer=True
		)

		disort_intcor = option_definition.option(
			name='disort_intcor',
			group='solver',
			helpstr='Intensity correction method', 
			documentation=documentation['disort_intcor'],
			gui_inputs=(GUI_definition.ListInput(name='disort_intcor', valid_range = ['phase', 'moments', 'off']),),
			tokens=option_definition.addLogical(name="Input.rte.disort_icm", logicals=['phase', 'moments', 'off'], setting='DISORT_ICM_'),
			parents=['rte_solver'],
			speaker='rte_solver',
			enable_values=("disort",)
		)

		isotropic_source_toa = option_definition.option(
			name='isotropic_source_toa',
			group='solver',
			helpstr='Specifies that isotropic '
			'illumination is included at top-boundary.',
			documentation=documentation['isotropic_source_toa'],
			tokens=option_definition.addSetting(name='Input.rte.fisot', setting=1),
			gui_inputs=(GUI_definition.FloatInput(name='isotropic_source_toa'),),
			parents=['rte_solver'],
			speaker='rte_solver',
			enable_values=("disort",
				      "twostr",)
		)

		raman = option_definition.option(
			name='raman',
			group='solver',
			helpstr='Include first order rotational Raman scattering.',
			documentation=documentation['raman'],
			tokens=[option_definition.addSetting(name='Input.raman', setting='RAMAN_CALC'),
				option_definition.addSetting(name='Input.ck_scheme', setting='CK_RAMAN'),
				option_definition.addSetting(name='Input.rte.solver', setting='SOLVER_DISORT'),
				option_definition.addSetting(name='Input.processing', setting='PROCESS_RAMAN'),
				option_definition.addLogical(name='Input.raman_original',logicals=['original'],optional=True) ],
			gui_inputs=(GUI_definition.ListInput(name='raman_original',valid_range=['original'],optional=True),),
			parents=['rte_solver'],
			speaker='rte_solver',
			enable_values=("disort",),
		)

		number_of_streams = option_definition.option(
			name='number_of_streams',
			group='solver',
			helpstr='Number of streams.',
			documentation=documentation['number_of_streams'],
			gui_inputs=(GUI_definition.IntegerInput(name='number_of_streams'),),
			tokens=option_definition.addToken(name='Input.rte.nstr', datatype=int, default=6),
			parents=['rte_solver'],
			speaker='rte_solver',
			enable_values=("disort","polradtran",),
		)

		polradtran = option_definition.option(
			name='polradtran',
			group='solver',
			helpstr='Set values for rte_solver polradtran.',
			documentation=documentation['polradtran'],
			tokens=[option_definition.addLogical(name='id', logicals=['aziorder','nstokes','src_code'], setting='POLRADTRAN_'),
				option_definition.addToken(name='Input.rte.polradtran[id]', datatype=int) ],
			speaker='rte_solver',
			enable_values=('polradtran',),
			showInGui=False,
			non_unique=True,
		)

		polradtran_quad_type = option_definition.option(
			name='polradtran_quad_type',
			group='solver',
			helpstr='Type of quadrature used',
			documentation=documentation['polradtran_quad_type'],
			tokens=option_definition.addToken(name='Input.rte.pol_quad_type', valid_range=['G', 'D', 'L', 'E'], datatype=str),
			parents=['rte_solver'],
			speaker="rte_solver",
			enable_values=("polradtran",),
			showInGui=False,
		)

		polradtran_max_delta_tau = option_definition.option(
			name='polradtran_max_delta_tau',
			group='solver',
			helpstr='Initial layer thickness for doubling.',
			documentation=documentation['polradtran_max_delta_tau'],
			tokens=option_definition.addToken(name='Input.rte.pol_max_delta_tau', datatype=float, default=1e-5 ),
			parents=['rte_solver'],
			speaker="rte_solver",
			enable_values=("polradtran",),
			showInGui=False,
		)


#no_sdisort		nrefrac = conditional_enable_option(
#no_sdisort			name='nrefrac',
#no_sdisort			group='solver',
#no_sdisort			helpstr='Refraction for rte_solver sdisort.',
#no_sdisort			documentation=documentation['nrefrac'],
#no_sdisort			gui_inputs=(IntegerGUI_definition.ListInput(name='Input.rte.nrefrac', default=0, valid_range=[0, 1, 2], optional=False),),
#no_sdisort			tokens=option_definition.addToken('Input.rte.nrefrac', int, default=0, valid_range=[0,1,2]),
#no_sdisort			parents=['rte_solver'],
#no_sdisort			speaker="rte_solver",
#no_sdisort			enable_values=("sdisort",)
#no_sdisort		)

#no_sdisort		nscat = option_definition.option(
#no_sdisort			name='nscat',
#no_sdisort			group='solver',
#no_sdisort			helpstr='Order of scattering.',	
#no_sdisort			documentation=documentation['nscat'],
#no_sdisort			gui_inputs=(GUI_definition.IntegerInput(name='Input.rte.nscat'),),
#no_sdisort			tokens=option_definition.addToken('Input.rte.nscat', int, default=20),
#no_sdisort			parents=['rte_solver'],
#no_sdisort			speaker="rte_solver",
#no_sdisort			enable_values=("sdisort",)
#no_sdisort		)

		sdisort = option_definition.option(
			name='sdisort',
			group='solver',
			helpstr='Set values for rte_solver sdisort.',
			documentation=documentation['sdisort'],
			tokens=[option_definition.addLogical(name='id', logicals=['ichapman','nscat','nrefrac'],setting='SDISORT_'),
				option_definition.addToken(name='Input.rte.sdisort[id]', datatype=int) ],
			speaker='rte_solver',
			enable_values=('sdisort',),
			showInGui=False,
			non_unique=True,
		)

		sos_nscat = option_definition.option(
			name='sos_nscat',
			group='solver',
			helpstr='Set order of scattering moments for rte_solver sos.',
			documentation=documentation['sos_nscat'],
			tokens=	option_definition.addToken(name='Input.rte.sos_nscat', datatype=int) ,
			speaker='rte_solver',
			enable_values=('sos',),
			showInGui=False,
                        developer=True
		)

		deltam = option_definition.option(
			name='deltam',
			group='solver',
			helpstr='Turn delta-M scaling on/off.',
			documentation=documentation['deltam'],
			gui_inputs=(GUI_definition.ListInput(name='Input.rte.deltam', valid_range=['on', 'off'], optional=False, default='on'),),
			tokens=option_definition.addLogical(name='Input.rte.deltam', logicals=['on','off'], setting='SWITCH_', default='SWITCH_ON'),
			parents=['rte_solver'],
			speaker="rte_solver",
			enable_values=("disort", "twostr",)
		)
		
		ipa_3d = option_definition.option(
			name='ipa_3d',
			group='solver',
			helpstr='',
			documentation=documentation['ipa_3d'],
			tokens=option_definition.addSetting(name='Input.ipa3d',setting=1),
			threedmystic =True,
			parents=['uvspec'],
			non_parents=['tipa',],
			speaker='rte_solver',
			enable_values=("disort","twostr","rodents",),
		)

		mc_tenstream = option_definition.option(
			name='mc_tenstream',
			group='solver', 
			helpstr='Run Tenstream approximation solver instead of Photon Tracing.', 
			documentation=documentation['mc_tenstream'],
			tokens=option_definition.addSetting(name='Input.rte.mc.tenstream', setting=1,default=0),
			parents=['uvspec'],
			speaker='rte_solver',
			enable_values=("mystic",),
			threedmystic =True
		)

		mc_ipa = option_definition.option(
			name='mc_ipa',
			group='solver', 
			helpstr='Run MYSTIC in independent pixel mode.', 
			documentation=documentation['mc_ipa'],
			tokens=option_definition.addSetting(name='Input.rte.mc.ipa', setting=1,default=0),
			parents=['uvspec'],
			non_parents=['mc_tipa',],
			speaker='rte_solver',
			enable_values=("mystic",),
			threedmystic =True
		)

		mc_tipa = option_definition.option(
			name='mc_tipa',
			group='solver',
			helpstr='Run MYSTIC in tilted independent pixel mode.',
			documentation=documentation['mc_tipa'],
			tokens=[option_definition.addSetting(name='Input.rte.mc.ipa', setting=1),
				option_definition.addLogical(name='Input.rte.mc.tipa', logicals=['dirdiff','dir3d'], setting='TIPA_') ],
			parents=['uvspec'],
			non_parents=['mc_ipa',],
			speaker='rte_solver',
			enable_values=("mystic",),
			developer=True,
		)

		tipa = option_definition.option(
			name='tipa',
			group='solver',
			helpstr='Use the tilted independent pixel approximation.',
			documentation=documentation['tipa'],
			tokens=option_definition.addLogical(name='Input.tipa',logicals=['dirdiff','dir'],setting='TIPA_'),
			parents=['uvspec'],
			non_parents=['ipa_3d',],
			speaker='rte_solver',
			enable_values=("rodents",'disort','twostr',),
			developer=True,
		)

		sslidar = option_definition.option(
			name='sslidar',
			group='solver',
			helpstr='Set single scattering lidar parameters.',
			documentation=documentation['sslidar'],
			tokens= [ option_definition.addLogical(name='id',logicals=['area','E0','eff','position','range'],setting='SSLIDAR_'),
				option_definition.addToken(name='Input.sslidar[id]', datatype=float) ],
			gui_inputs= (GUI_definition.ListInput(name='id', valid_range=['area','e0','eff','position','range']),GUI_definition.FloatInput(name='Input.sslidar[id]'),),
			parents=['uvspec'],
			speaker='rte_solver',
			enable_values=("sslidar",),
			non_unique=True,
		)

		sslidar_nranges = option_definition.option(
			name='sslidar_nranges',
			group='solver',
			helpstr='Set number fo range bins for single scattering lidar.',
			documentation=documentation['sslidar_nranges'],
			tokens= [ option_definition.addToken(name='Input.sslidar_nranges', datatype=int) ],
			gui_inputs= (GUI_definition.IntegerInput(name='Input.sslidar_nranges'),),
			parents=['uvspec'],
			speaker='rte_solver',
			enable_values=("sslidar",),
		)

		sslidar_polarisation = option_definition.option(
			name='sslidar_polarisation',
			group='solver',
			helpstr='Turn on polarisation measurement for lidar.',
			documentation=documentation['sslidar_polarisation'],
			tokens= option_definition.addSetting(name='Input.sslidar_polarisation', setting=1),
			parents=['uvspec'],
			speaker='rte_solver',
			enable_values=("sslidar",),
		)

		tzs_cloud_top_height = option_definition.not_yet_lex2py_option(
                        name='tzs_cloud_top_height',
			group='solver',
			helpstr='Cloud top height[s] of blackbody clouds',
			documentation=documentation['tzs_cloud_top_height'],
			tokens=option_definition.addToken(name="", datatype=str),
			parents=['uvspec'],
			speaker='rte_solver',
			enable_values=('tzs',),
		)

		self.options = [rte_solver, pseudospherical, disort_spherical_albedo,
				deltam, 
				isotropic_source_toa, raman, 
                                mc_tenstream,
				mc_ipa, ipa_3d, mc_tipa, tipa,
				number_of_streams, 				
				disort_intcor, 
				polradtran, 
				polradtran_quad_type, polradtran_max_delta_tau, 
				sslidar, sslidar_polarisation, sslidar_nranges, 
				sdisort, sos_nscat, tzs_cloud_top_height ]

	def __iter__(self):
		return iter(self.options)


def get_solver_documentation():
	return {
		'sslidar' : r'''
	Set single scattering lidar parameters (\code{rte_solver sslidar}).
	\fcode{
	sslidar variable value
	}
	\code{variable} can be one of the following:
	\begin{description}
	\item[area] Set area of single scattering lidar in units of square meters (default: 1.0)
	\item[E0] Set Laser pulse energy for single scattering lidar in units of (default: 0.1)
	Joule (You can also use a \code{source solar file}
	instead... not yet implemented.)
	\item[eff] Set lidar efficiency for single scattering lidar (default: 0.5)
	\item[position] Set lidar position for single scattering lidar in units of km (default: 0.0)
	\item[range] Set lidar range bin width for single scattering lidar in units of km (default: 0.1)
	\end{description}
		''',
		'sslidar_nranges' : r'''
	Set number of range bins for single scattering lidar (solver
	\code{sslidar}). Default is 100.
		''',
		'sslidar_polarisation' : r'''
	Turn on polarisation measurement for lidar (solver \code{sslidar}). Default is without polarisation.
		''',

		'ipa_3d' : r'''
	This option can be used to carry out a 3D calculation with a 1D-solver, where a 
	3D-cloud file as used for MYSTIC (flag 3) serves as an input file. Further, a sample grid needs to be specified
	via \code{mc_sample_grid}, where the grid and the resolution has to be the same as in the 3D-cloud file
	(Nx, Ny, dx, dy). To calculate heating rates, \code{mc_backward_output heat K_per_day} must stand in the input
	file. 
		''',

		'mc_tenstream' : r'''
        The tenstream solver is a mpi-parallelized finite volume method.

        First build libtenstream:\\
        Either run \\
        \code{ ./configure --with-download-tenstream } \\
        or do the checkout yourself: \\
          \code{svn co https://svn.physik.uni-muenchen.de/repos/hdcp2-rt/tenstream}\\
        Then follow the instructions to compile the tenstream lib in \\
        \code{tenstream/README} \\
        If it was successful, there should be the file \\ 
        \code{tenstream/build/lib/libtenstream.a}.

        Now configure libRadtran and it should hopefully find the tenstream-lib:\\
          \code{ ./configure  }\\
        It should say \\
          \code{ TENSTREAM solver:     yes }\\

        Input files are the same as for MYSTIC calculations with the additional 
        option mc\_tenstream.
        You can call the MPI version with:\\
        \code{mpirun -np 2 bin/uvspec_mpi <inputfile>  <PETSc options> }.

        If you are serious about solving bigger systems we highly recommend that you make yourself familiar with PETSc's options.\\
        As a start you could try:

          Solve the direct radiation with direct LU method\,(Only very small problems, 1 process)\\
              \code{-dir\_ksp\_type richardson -dir\_pc\_type lu}\\
          Solve diffuse radiation with parallel LU\\
              \code{-diff\_ksp\_type richardson -diff\_pc\_type lu -diff\_pc\_factor\_mat\_solver\_package mumps}\\
          Monitor convergence\,(recommended)\\
              \code{-dir\_ksp\_monitor -dir\_ksp\_view -dir\_ksp\_converged\_reason} \\
              \code{-diff\_ksp\_monitor -diff\_ksp\_view -diff\_ksp\_converged\_reason}\\
          Solve diffuse radiation with bi-conjugate-gradient-stabilized with incomplete LU with 2 levels of fill for each node.\,(for middle sized problems)\\
              \code{-diff\_ksp\_type bcgs -diff\_pc\_type bjacobi -diff\_sub\_pc\_factor\_levels 2} \\
          Algebraic Geometric Multigrid Preconditioner\,(for big problems)\\
              \code{-diff\_pc\_type ml} \\
        Options may also be written to \\~/.petscrc.\\

                ''',

                'mc_ipa' : r'''
	Run MYSTIC in independent pixel mode. Only meaningful with \code{rte_solver montecarlo}.
		''',
		'polradtran_max_delta_tau' : r'''
	Initial layer thickness for doubling; governs accuracy, 10E-5 should be
	adequate. Do not go beyond half the real precision, i.e. 10e-8 for REAL*8.
	Default 1.e-05.
	\fcode{
	polradtran_max_delta_tau value
	}
	This option is only relevant for \code{rte_solver polradtran}.
		''',
		'mc_tipa' : r'''
	Run MYSTIC in tilted independent pixel mode. Only meaningful with \code{rte_solver
	montecarlo}. Choose either \code{mc_tipa dir3d} or \code{mc_tipa
	dirdiff}. The latter is the analogon to \code{tipa dirdiff}, where the
	3d-cloud-matrix is tilted according to the solar zenith angle and is
	then used as input for all the following calculations. 
	In \code{mc_tipa dir3d}, photons are propagated in 3d until the first
        scattering occurs, after that the photon is captured within the current
        column.
		''',

		'tipa' : r'''
	With the option \code{tipa dirdiff}, the tilted independent
	pixel (column) approximation (TIPA, the same as TICA) is carried out,
	where a 3D-cloud file serves as input and the columns are tilted
	towards the sun (see, e.g., Wapler and Mayer, 2008).  \code{tipa
	dirdiff} works with all 1D-rte-solvers.  When choosing
	tipa (\code{tipa dir}), from the "tilted clouds" only the optical
	depth for the calculation of the direct radiation is derived and used
	later to calculate the diffuse radiation. To calculate the diffuse
	radiation, the original cloud-matrix is still necessary /
	maintained. \code{tipa dir} works only with Roberts (Buras)
	delta--Eddington two--stream model \code{rodents}.  For both tipa
	options, a \code{mc_sample_grid} has to be specified as Nx Ny dx dy. A
	3d-cloud-file (with flag 3) serves as input (see
	\code{wc_file 3D}). The output is written into mc.flx(.spc),
	however, a basename via \code{mc_basename} can be given.
		''',

		'rte_solver' : r'''
	Set the radiative transfer equation solver to be used. 
	\fcode{
	rte_solver type
	}
	If not specified the default \code{rte_solver} is \code{disort}.
	Choices for \code{type} are
	\begin{description}

	\parameter{disort}
	C-version of the disort algorithm, translated from Fortran by Tim Dowling. 
	This is the recommended discrete ordinate code in {\sl libRadtran}.
	For documentation see \file{src\_f/DISORT2.doc} as well as the papers and
	the DISORT report at
	ftp://climate1.gsfc.nasa.gov/wiscombe/Multiple\_Scatt/.  The intensity
	correction can be performed according to \citet{nakajima1988} using
	\code{disort_intcor moments} (like in the original code), or with the 
	improvements described in (Buras, Dowling, Emde, in preparation; default). Can be run in  
	plane-parallel geometry (default) or in pseudo-spherical geometry 
	(using \code{pseudospherical}).

	\parameter{twostr}
	C-version of the two--stream radiative transfer solver described by
	\citet{Kylling1995}. Can be run in plane-parallel geometry (default)
	or in pseudo-spherical geometry (using \code{pseudospherical}
	).

	\parameter{fdisort1}
	The standard plane--parallel disort algorithm by \citet{Stamnes1988c},
	version 1.3 -- provided for compatibility reasons. Use only if you have troubles with the default			         
        \code{disort} or for historical reasons. For documentation see \file{src\_f/DISORT.doc} as well 
	as the papers and the DISORT report at
	ftp://climate1.gsfc.nasa.gov/wiscombe/Multiple\_Scatt/. To optimize
	for computational time and memory, please adjust the parameters in
	\file{src\_f/DISORT.MXD} for your application and re-compile. For your
	application please use \code{rte_solver fdisort2} which is the advanced
	version, unless you e.g. want to explore how a specific feature of
	\code{fdisort2} (e.g. the \citet{nakajima1988} intensity correction)
	improves the \code{fdisort1} result.

	\parameter{fdisort2}
	Version 2 of the Fortran algorithm disort -- provided for compatibility reasons.
	Use only if you have troubles with the default \code{disort} or for historical reasons. For documentation see
	\file{src\_f/DISORT2.doc} as well as the papers and the DISORT report
	at ftp://climate1.gsfc.nasa.gov/wiscombe/Multiple\_Scatt/
	\code{fdisort2} has several improvements compared to its 'ancestor'
	\code{fdisort1} (version 1.3). To optimize for computational time and
	memory, please adjust the parameters in \code{src_f/DISORT.MXD} for
	your application and re-compile.
	Note! \code{fdisort2} is a new version of the original disort code
	which was implemented in summer 2009. It uses phase functions to
	calculate the intensity corrections by \citet{nakajima1988} instead of
	Legendre moments. Hence it needs cloud properties files which contain
	the phase functions. It is still possible to use the old version of
	disort, you need to specify \code{disort_intcor moments}.

	\parameter{sdisort}
	Pseudospherical disort as described by \citet{Dahlback1991}. Double
	precision version. To optimize for computational time and memory,
	please adjust the parameters in \file{src\_f/DISORT.MXD} for your
	application and re-compile.  

	\parameter{spsdisort} Pseudospherical
	disort as described by \citet{Dahlback1991}, single precision
	version. \strong{Warning:} it is not recommended to use
	\code{spsdisort} for really large solar zenith angles nor for cloudy
	conditions. For large optical thickness it is numerically unstable and
	may produce wrong results. To optimize for computational time and
	memory, please adjust the parameters in \file{src_f/DISORT.MXD} for
	your application and re-compile.  

	\parameter{polradtran} The plane-parallel radiative transfer solver of \citet{Evans1991}.
	Includes polarization. The full implementation of the polRadtran
	solver in \code{uvspec} is quite new (version 1.4). If you find
	unusual behaviour, please contact the {\sl libRadtran} authors.

	\parameter{ftwostr}
	Original Fortran-version of the two--stream radiative transfer solver described by
	\citet{Kylling1995}, in pseudo-spherical geometry.

	\parameter{rodents}
	Delta-Eddington two--stream code (RObert's Delta-EddingtoN Two-Stream), plane-parallel.

	\undocumented{
	 \parameter{twostrebe}
	 Delta-Eddington two--stream solver by Bernhard Mayer. \code{zout_interpolate} is activated automatically.
	}

	\parameter{sslidar}
	A simple single scattering lidar simulator by Robert Buras.

	\parameter{sos}
	A scalar pseudospherical succesive orders of scattering code. Works
	for solar zenith angles smaller than 90 degrees. Can calculate
	azimuthally averaged radiances. Set \code{sos_nscat} to specify the order
	of scattering.

	\ifmystic{
	\parameter{montecarlo}
	The MYSTIC Monte Carlo code. Monte Carlo is
	the method of choice (1) for horizontally
	inhomogeneous problems; (2) whenever
	polarization is involved; (3) for applications
	where spherical geometry plays a role; and (4)
	whenever sharp features of the scattering
	phase function play a role, like for the 
	calculation of the backscatter glory or the
	aureole.  
	\parameter{mystic}
	Same as \code{montecarlo}.
	}

	\parameter{tzs} TZS stands for "thermal, zero scattering" and is a
	very fast analytical solution for the special case of thermal emission
	in a non-scattering atmosphere. Please note that TZS does only
	radiance calculations at top of the atmosphere. "Blackbody clouds" may be included 
        using the option \code{tzs_cloud_top_height}. 

	\parameter{sss}
	SSS stands for "solar, single scattering" and is an analytical single
	scattering approximation which might be reasonable for an optically
	thin atmosphere. Please note that SSS does only radiance calculations
	at top of the atmosphere. This is an experimental solver - be careful!

	\undocumented{
	\parameter{sssi} Missing documentation!!!
	}

	\parameter{null}
	The NULL solver does not solve the radiative transfer
	equation. However, it sets up the optical properties, and does the
	post-processing; useful if you are either interested in the overhead
	time required by a particular model input or if you are simply
	interested in the optical properties, as output by \code{verbose}.
	\end{description}
	Default: \code{disort}
		''',
		'raman' : r'''
	\fcode{
	raman [original]
	}
	The \code{raman} option includes single order rotational Raman scattering in 
	the calculation. The solution treats Raman as a perturbation similar to the 
	approaches of \citet{Vountas1998} and \citet{Spurr2008}.
	
	The \code{raman} option may only be used for spectral calculation.

	The disort radiative transfer solver with a general source is needed to solve the radiative transfer
	equation including Raman scattering. This solver is automagically invoked when specifying 
	the \code{raman} option. It is thus not neccessary to set the \code{rte_solver}.

	The \code{raman} option is optimized with respect to speed. The optimized implementation
	should be just as accurate as the original version. To use the original version
	invode \code{raman original}. With the optional argument \code{original} each wavelength 
	is treated individually and is thus accurate, but computationally very expensive. 

	Please note that while the \code{raman} option has been extensively tested and verified, it 
	is nevertheless a new option, hence, use it with care.
		''',
		'pseudospherical' : r'''
	Invokes pseudo-spherical geometry in disort/twostr. Default is plane-parallel.
		''',
		'disort_spherical_albedo' : r'''
	Calculate spherical albedo using \code{disort}. When this option is enabled,
	only the spherical albedo is calculated. The output is enabled by 
	\fcode{output_user spher_alb}
		''',
		'disort_intcor' : r'''
	Intensity correction method for \code{rte_solver disort} or
	\code{rte_solver fdisort2}. Valid options are \code{phase}, i.e.~the
	phase function is used for the Nakajima intensity correction, and
	\code{moments}, i.e.~the Legendre moments are used for the
	correction. Optionally, the option \code{off} turns off the intensity
	correction. Default is \code{phase}.
		''',
		'isotropic_source_toa' : r'''
	Specifies that isotropic illumination is used at top-boundary instead
	of beam source. Useful for those who want to calculate the reflectance 
	for a homogeneous or inhomogeneous atmosphere. The intensity is still set by
	\code{source solar file}. Only works with \code{disort} and \code{twostr}.
		''',
		'number_of_streams' : r'''
	Number of streams used to solve the radiative transfer equation.
	\fcode{
	number_of_streams value
	}
	Default is 6 for fluxes and 16 for radiances.
	(For \code{rte_solver fdisort1}, \code{fdisort2} and \code{disort} only
	even number of streams are possible.)
		''',
		'polradtran_quad_type' : r'''
	Type of quadrature used:
	\fcode{
	polradtran_quad_type type
	}
	where \code{type} is one of
	\begin{description}
	\item[G] (default) gaussian
	\item[D] double gaussian, 
	\item[L] Lobatto
	\item[E] extra-angle(s), this must be used of \code{polradtran} is used in 
	combination with \code{umu}. Will internally use Gaussian scheme (\code{G}). 
	See also radtran documentation (\file{libsrc_f/README.polRadtran}).
	\end{description}
	Default G.                                  
	This option is only relevant for \code{rte_solver polradtran}.
		''',

		'polradtran' : r'''
	Specify polradtran values. This option is only relevant for \code{rte_solver polradtran}.

	\fcode{
	polradtran aziorder value
	}
	Order of Fourier azimuth series.
	The value \code{0} (default for irradiance) is the azimuthally symmetric case.
	For radiance computation a higher order is required, thus the default for radiances is 4. 

	\fcode{
	polradtran nstokes value
	}
	Number of Stokes parameters
	where \code{value} is one of
	\begin{description}
	\item[1] for I (no polarization, default)
	\item[2] for I,Q,U (Since V is very small in the atmosphere, it makes sense
	           to compute only I,Q,U. This saves computation time and
		   memory).
	\item[3] for I,Q,U,V
	\end{description}
	Default is 1.

	\fcode{
	polrdatran src_code value
	}
	Radiation sources included
	which may be
	\begin{description}
	\item[0] none
	\item[1] solar
	\item[2] thermal
	\item[3] both
	\end{description}
	Default 1.
		''',
		'sos_nscat' : r'''
	Set the order of scattering for \code{rte_solver sos}.
	\fcode{
	sos_nscat value
	}
	Default is 20.
		''',
		'sdisort' : r'''
	Specify sdisort values. This option is only relevant for \code{rte_solver sdisort}.

	%CE: Commented for ESA version
	\undocumented{
	\fcode{
	sdisort ichapman value
	}
	Undocumented option for sdisort. Default is 1.
	}

	\fcode{
	sdisort nscat value
	}
	Set the order of scattering. 
	If \code{value} is set to 1 \code{sdisort} will run in single scattering mode 
	while if set to 2, \code{sdisort} runs in full multiple scattering mode.
	Default is 2 for \code{rte_solver sdisort}. 

	\fcode{
	sdisort nrefrac value
	}
	Include refraction where \code{value} has the meaning
	\begin{description}
	\item[0] No refraction, default.
	\item[1] Refraction included using fast, but harsh method.
	\item[2] Refraction included using slow, but accurate method.
	\end{description}
	If refraction is included also set parameter \code{refraction_file}.
		''',
		'deltam' : r'''
	Turn delta-M scaling on/off. Set to either \code{on} or
	\code{off}. Note that for the \code{rte_solver disort} and
	\code{rte_solver fdisort2} delta-M scaling is hardcoded to be always
	on.
		''',
                'tzs_cloud_top_height' : r'''
	Define cloud top height of "blackbody clouds" for radiative transfer solver \code{tzs}.
		''',
	}
