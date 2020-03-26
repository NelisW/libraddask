
"""--------------------------------------------------------------------
 * $Id: writeLex.py 2919 2013-05-17 09:24:06Z robert.buras $
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

found_lidar     = "no"
found_mystic    = "yes"
found_mystic3d  = "no"

def writeLex():
	import os

	pwd = os.getcwd()
	
	options = loadOptions()   #load python dictionary with input options
	lexParser = loadLexStarter(os.path.join(pwd,'lex_starter.l'))	#load lexStarter.l
	lexParser = insertInputOption(lexParser,options)	#add input options from python dictionary to lex parser
	lexParser = insertDocu(lexParser,options)		#add documentation to python dictionary
	lexParser = insertInitialisations(lexParser,options)	#add default initialisations

	writeUvspecLex(lexParser,os.path.join(pwd,'uvspec_lex.l'))

def loadOptionsHelper():
	"""
	Returns a list of all the groups.
	Used in the GUI.
	"""
	from . import aerosol_options, cloud_options, profile_options, geometry_options, special_options
	from . import mc_options, spectral_options, surface_options, molecular_options, solver_options, output_options, general_atmosphere_options

	profile = profile_options.setup_cloud_group()
	cloud = cloud_options.setup_cloud_group()
	aer = aerosol_options.setup_aerosol_group()
	special = special_options.setup_special_group()

	mc = mc_options.setup_mc_group()
	spectral = spectral_options.setup_spectral_group()
	gen_atm = general_atmosphere_options.setup_general_atm_group()
	sur = surface_options.setup_surface_group()
	mol = molecular_options.setup_molecular_group()
	sol = solver_options.setup_solver_group()
	out = output_options.setup_output_group()
	geom = geometry_options.setup_geometry_group()

	list_of_groups=[spectral, gen_atm, mol, aer, profile, cloud, sur, sol, mc, geom, out, special]
	for group in list_of_groups:
		for option in group:
			if found_lidar == 'no' and option['islidar']:		group.options.remove(option)
			elif found_mystic == 'no' and option['mystic']:		group.options.remove(option)
			elif found_mystic3d == 'no' and option['threedmystic']:	group.options.remove(option)

	return list_of_groups 


def loadOptions():
	options={}
	for group in loadOptionsHelper():
		for option in group:
			options[option.get('name')]=option
	return options	

def loadLexStarter(filename):
	import os

	#read lex file
	pwd = os.getcwd()
	lexstarter = open(filename,'r')
	text = lexstarter.read()
	lexstarter.close()
	return text

def writeUvspecLex(text,filename):
	newlex = open(filename,'w')
	newlex.write(text)
	newlex.close()

def getStr(string):
	if str(string) == 'True': return 'TRUE'
	elif str(string) == 'False': return 'FALSE'
	else: return str(string)
		

def getTypeScan(dtype, valid_range):
	if dtype == str or dtype == io.IOBase:
		return '{SPACE_WORD}'
	elif dtype == option_definition.Dimension:
		return '{SPACE}+{DIM}'
	elif dtype == option_definition.ProfileType:
		return '{SPACE_option_definition.ProfileType}'
	elif dtype == option_definition.SignedFloats :
		return '{SIGNED_FLOATS}'
	elif dtype == option_definition.Integers :
		return '{INTEGERS}'
	elif dtype == option_definition.CaothType or dtype == option_definition.CaothoffType :
		return '{SPACE_WORD}'
	elif dtype == int:
		return '{SPACE_INTEGER}'
	# elif dtype == long:
	# 	return '{SPACE_FLOAT}'
	elif dtype == float or dtype == option_definition.Double:
		if valid_range and valid_range[0] >= 0:
			return '{SPACE_FLOAT}'
		else:
			return '{SPACE_SIGNED_FLOAT}'
	else: print('ERROR: datatype not known: ', dtype)

def yytext2type(dtype):
	if dtype == str or dtype == io.IOBase or dtype == option_definition.ProfileType:
		return 'yytext2string'
	elif dtype == float:
		return 'yytext2float'
	elif dtype == option_definition.Double:
		return 'yytext2double'
	elif dtype == int:
		return 'yytext2int'
	# elif dtype == long:
	# 	return '(long int) yytext2double'

def getLexOptions(group):

	inputOption = ''
	
	for name, option in sorted(group.items()):
		# Do something completly different
		if hasattr(option, "writeLex"):
			inputOption += option.writeLex()
		# Do the usual
		else:
			tokens = option.get('tokens')

			#write input option name
			inputOption += name

			nTokens, rules, optional = lexRules(tokens)

			#write setting of input arguments
			variableSetting = '' 
			itok = 1
			for tok in tokens: 
				if tok.get('optional'):		variableSetting += '\tif ( ntokens >= %d )\t {\n\t' %( itok )

				if isinstance(tok, option_definition.addSetting):
					variableSetting += '\t%s = %s;\n' %(tok.get('name'),
                                                                getStr(tok.get('setting')))

				elif isinstance(tok,option_definition.addLogical):
					itok += 1
					variableSetting += getLogicals(tok,itok,nTokens)
					
				else:
					try:
						itok += 1	
						variableSetting += getTokens(tok,itok, str(nTokens) )
					except Exception as e:
						print('ERROR found in input option %s :\n\t Datatype %s of variable %s not known!' % (
							option.get('name'), tok.get('datatype'), tok.get('name') ))
						print(e)
				if  tok.get('optional'):         variableSetting += '\t}\n'

			inputOption += '%s {\n%s%s}\n' %( rules, optional, variableSetting )

	return inputOption

def lexRules(tokens):
	nTokens=1; 
	tokenStr=''; rules=''; optional=''
	for tok in tokens:
		if isinstance(tok,option_definition.addSetting):
			continue
		else:
			rules += getTypeScan(tok.get('datatype'),	   tok.get('valid_range'))
			nTokens +=1
			tokenStr = str(nTokens)
			if tok.get('optional'):
				tokenStr = 'ntokens+1'
				optional = '\tntokens = yytext2ntokens(yytext);\n'
				rules += '?'
	return [ tokenStr, rules, optional ]

def getTokens(tok,itok,nTokens):

	dtype = tok.get('datatype')
	if dtype == io.IOBase:
		inputStr = '\tyytextCstring ( %s, yytext, %d, %s);\n' % ( tok.get('name'), itok, str(nTokens) ) 
	elif dtype == option_definition.SignedFloats:
		inputStr = '\tntokens = yytext2ntokens(yytext);\n\t%s =  yytext2floats(yytext, ntokens);\n' % ( tok.get('name') )
	elif dtype == option_definition.Integers:
		inputStr = '\tntokens = yytext2ntokens(yytext);\n\t%s =  yytext2integers(yytext, ntokens);\n' % ( tok.get('name') )
	elif dtype == option_definition.CaothType:
		inputStr = '\t%s = get_caoth_index(&Input.caoth,&Input.n_caoth, yytext2string(yytext, %d, %s),0);\n' % (  tok.get('name'), itok, str(nTokens) )
	elif dtype == option_definition.CaothoffType:
		inputStr = '\t%s = get_caothoff_index(&Input.caothoff,&Input.n_caothoff, yytext2string(yytext, %d, %s));\n' % (  tok.get('name'), itok, str(nTokens) )
	else:
		inputStr = '\t%s = %s(yytext, %d, %s);\n' % ( tok.get('name'),
			  yytext2type(tok.get('datatype')),
			  itok, str(nTokens) )
	return inputStr



def getLogicals(tok,itok,nTokens):
	import re

	#dictionary to subsitute charakters like + - . to underscore _ for matching definition in C code
	rdict={'-': '_',  '.': '_',  '+':'_' }
	regex = re.compile("(%s)" % "|".join(map(re.escape, rdict.keys())))

	#uvspec_lex.l : read token for string comparison
	inputStr = '\ts = yytext2string(yytext, {0}, {1});\n\t'.format(itok,nTokens)

	#string comparison to set logical values
	logicals = tok.get('logicals')
	for ilog,log in enumerate(sorted(logicals)):
		if ilog == 0:	inputStr += 'if'
		else:		inputStr += 'else if'
		if len(logicals) == 1:  # one logical -> destination = 1
			inputStr += ' ( strcasecmp("{0}", s)==0 )\t{{ {1} = 1; }}\n\t'.format(log,tok.get('destination'))
		else: #more logicals -> destination is logical in uppercase, this variable should be defined in uvspec.h!
			inputStr += ' ( strcasecmp("{0}", s)==0 )\t{{ {1} = {2}{3}; }}\n\t'.format(log,tok.get('destination'), tok.get('setting'), ( regex.sub(lambda mo: rdict[mo.string[mo.start():mo.end()]], str(log)) ).upper() ) 
	if tok.get('logical_file'):  #save token to destination_filename
		inputStr += 'else\t{{ {0} = {1}FILE;\n\t\t strcpy ( {0}_filename, s ); }}\n\t'.format(tok.get('destination'), tok.get('setting'))
	else:	#raise error
		inputStr += r'else\t{{fprintf(stderr,"Option %s  has an invalid argument {0} on line %d: %s\\n", yytext, line_number+1, s); ierror++;}} \n\t'.format(itok-1)
	inputStr += 'free(s); \n'
	return inputStr


def insertInputOption(text,options):
	import re

	#insert input options to lexfile
	groupInput = getLexOptions(options)
	text = re.sub('\n.*PYTHON_OPTIONS.*\n' , '\n%s\n' %(groupInput), text, count=1)
	return text

def insertDocu(text,options):
	import re

	#get new documentation
	docu={}
	for key,val in sorted(options.items()):                  
		docu[key]=val['documentation']
	
	#get whole documentation
	doclist = list(docu.keys()) + re.findall('(?<=option{)\w+',text)
	doclist.sort( key = lambda member: member.lower() )
	
	#insert documentation to lexfile
	for iopt, option in enumerate(reversed(doclist)):
		i=len(doclist)-1-iopt
		if option in docu.keys() and doclist.count(option) == 1:
			if i<len(doclist)-1:
				if   text.find('\\ifthreedmystic{\n\\option{%s}' %(doclist[i+1])) > 0: n = text.find('\\ifthreedmystic{\n\\option{%s}' %(doclist[i+1])) 
				elif text.find('\\ifmystic{\n\\option{%s}'	 %(doclist[i+1])) > 0: n = text.find('\\ifmystic{\n\\option{%s}' %(doclist[i+1])) 
				elif text.find('\\iflidar{\n\\option{%s}' 	 %(doclist[i+1])) > 0: n = text.find('\\iflidar{\n\\option{%s}' %(doclist[i+1]))
				elif text.find('\\undocumented{\n\\option{%s}' 	 %(doclist[i+1])) > 0: n = text.find('\\undocumented{\n\\option{%s}' %(doclist[i+1])) 
				elif text.find('\\option{%s}' 			 %(doclist[i+1])) > 0: n = text.find('\\option{%s}' %(doclist[i+1]))
				else: n = text.find('%DOCUMENTATION END')
			else: n = text.find('%DOCUMENTATION END')

			if options[option].get('developer'):
				text = '%s\\undocumented{\n\\option{%s}%s}\n%s' %(text[:n], option, docu[option], text[n:])
			elif options[option].get('islidar'):
				text = '%s\\iflidar{\n\\option{%s}%s}\n%s' %(text[:n], option, docu[option], text[n:])
			elif options[option].get('threedmystic'):
				text = '%s\\ifthreedmystic{\n\\option{%s}%s}\n%s' %(text[:n], option, docu[option], text[n:])
			elif options[option].get('mystic'):
				text = '%s\\ifmystic{\n\\option{%s}%s}\n%s' %(text[:n], option, docu[option], text[n:])
			else:
				text = '%s\\option{%s}%s\n%s' %(text[:n], option, docu[option], text[n:])

	return text

def insertInitialisations(text,options):
	import re

	#get sorted initialisations
	inits = {}
	for val in iter(options.values()):
		#get default initialisation of tokens
		for initialize in ['tokens', 'settings']:
			if val.get(initialize):
				for tok in val.get(initialize):
					if tok.get('default') != None: 
						if tok.get('name') in inits.keys() and inits[tok.get('name')] != tok.get('default'): 

							print('ERROR: \nvariable %s has different initialisations: %s and %s\n Please check your option definitions again.\n' %(
								tok.get('name'), getStr(inits[tok.get('name')]), getStr(tok.get('default')) ))
						else: inits[tok.get('name')] = tok.get('default')
	
	initialisation = ''
	for key in sorted(inits.keys()):
		initialisation += '  %s = %s;\n' %(key,getStr(inits.get(key)) )
	text = re.sub('\n.*PYTHON INITIALISATIONS.*\n' , '\n%s\n' %(initialisation), text, count=1)
	return text

if __name__ == "__main__":
	writeLex()
