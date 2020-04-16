"""--------------------------------------------------------------------
 * $Id: profile_options.py 3023 2014-04-29 13:58:13Z Claudia.Emde $
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
import io

class setup_cloud_group():

	group_name = 'Profile'

	def __init__(self):
		documentation = get_profile_documentation()

		profile_file = option_definition.option( 	#TODO: option 3D only ifthreedmystic
			name='profile_file',
			group='profile',
			helpstr='Location of file defining profile properties.',
			documentation=documentation['profile_file'],
			gui_inputs=(GUI_definition.TextInput(name='name'),GUI_definition.ListInput(name='source', valid_range=['1D', '3D', 'moments', 'ipa', 'ipa_files'], optional=False), GUI_definition.FileInput(name='Input.caoth[isp].filename'),),
			tokens = [option_definition.addToken(name='isp', datatype=option_definition.CaothType), 
				option_definition.addLogical( name='Input.caoth[isp].source', logicals=option_definition.ProfileType().get_valid_range(), setting='CAOTH_FROM_'), 
				option_definition.addToken( name='Input.caoth[isp].filename', datatype=io.IOBase ) ],

			parents=['uvspec'],
			childs=['profile_properties','profile_modify','interpret_as_level'],
			plot = {'plot_type': '2D',
				'optional_args': {'column_names': (
						"altitude",
						"liquid water content",
						"effective radius",)
						  }
				},
			non_unique=True,
			)

		profile_properties = option_definition.option(
			name='profile_properties',
			group='profile',
			helpstr='Define profile optical properties. ',
			documentation=documentation['profile_properties'],
#			gui_inputs=(GUI_definition.ListInput(name='Input.caoth[isp].properties', valid_range=['hu', 'echam4', 'mie', 'fu', 'yang', 'key', 'baum_hufit', 'baum', 'hey', 'ic_mie'], optional=False), GUI_definition.ListInput(name='Input.caoth[isp].interpolate', valid_range=['','interpolate'], optional=True),),
			tokens = [ option_definition.addToken( name='isp', datatype=option_definition.CaothType ),
				option_definition.addLogical( name='Input.caoth[isp].properties', logicals=['hu', 'echam4', 'mie', 'fu', 'yang', 'key', 'baum_v36', 'baum', 'hey', io.IOBase], setting='PROP_' ),
				option_definition.addLogical( name='Input.caoth[isp].interpolate' , logicals=['interpolate'], optional=True) ],	#TODO: BUTTON in GUI, not str!

				# option_definition.addLogical( name='Input.caoth[isp].properties', logicals=['hu', 'echam4', 'mie', 'fu', 'yang', 'key', 'baum_v36', 'baum', 'hey', file], setting='PROP_' ),

			parents=['profile_file'],
			non_unique=True,
			)

		profile_modify = option_definition.option(	#TODO: valid_ranges for GUI!!
			name='profile_modify',
			group='profile',
			helpstr='Modify profile properties.',
			documentation=documentation['profile_modify'],
#			gui_inputs=(GUI_definition.ListInput(name='id1', valid_range=['gg', 'ssa', 'tau', 'tau550'], optional=False), GUI_definition.ListInput(name='id2', valid_range=['set', 'scale'], optional=False), GUI_definition.FloatInput(name='Input.caoth[isp].modify[id1][id2]'),),
			tokens = [ option_definition.addToken( name='isp', datatype=option_definition.CaothType ),
				option_definition.addLogical( name='id1', logicals=['gg','ssa', 'tau', 'tau550'], setting='MODIFY_VAR_' ),
				option_definition.addLogical( name='id2', logicals=['set', 'scale' ], setting='MODIFY_TYPE_' ),
				option_definition.addToken( name='Input.caoth[isp].modify[id1][id2]', datatype=float ) ], 
			parents=['profile_file'],
			non_unique=True,
		)

		self.options = [profile_file, profile_properties, profile_modify]
		
	def __iter__(self):
		return iter(self.options)

def get_profile_documentation():
	return {
		'profile_file' : r'''
	Define file containing properties of clouds, aerosols, hydrometeors, etc. 
	This option is a generalization of the options \code{wc\_file} and \code{ic\_file}.
	
	Usage:
	\fcode{
	profile\_file typename type file
	}
	
	\code{typename} describes the name of the profile; typically this 
	describes what kind of particles are dealt with here. Examples
	are \code{wc} (water clouds), \code{ic} (ice clouds), \code{aer} (aerosols), drizzle. 
	The \code{typename} is needed to refer to this profile when using other options, 
	such as \code{profile\_properties}. Note that \code{typename} "wc" and "ic"
	have special effects (i.e.~default properties, and "ic" properties are
	not allowed with "wc" files, and vice versa).
	
	\code{type} defines the file type, which is identical to \code{wc\_file type}.
	Please refer to \code{wc\_file} for choises of \code{type} and a detailed description 
	on the file structures.	
		''',

		'profile_properties' : r'''
	Define how liquid/ice water content/mass concentration and effective
	particle radius are translated to optical properties for profile
	\code{typename}. This option is a generalization of the options
	\code{wc\_properties} and \code{ic\_properties}.
	
	Usage:
	\fcode{
	profile\_properties typename property [interpolate]
	}
	
	\code{typename} describes the name of the profile; it must be identical
	to the one defined in \code{profile\_file}.
	
	Please refer to \code{wc\_properties} and \code{ic\_properties} for possible choices for \code{property}.
		''',
		'profile_modify' : r'''
	Modify profile optical properties.
	\fcode{
	profile typename variable scale/set value
	}
	\code{typename} describes the name of the profile; it must be identical
	to the one defined in \code{profile\_file}.

	This option is identical to \code{wc\_modify}.
	Please refer to \code{wc\_modify} for a detailed description of \code{variable}.
		''',	
	}

