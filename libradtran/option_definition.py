
"""--------------------------------------------------------------------
 * $Id: option_definition.py 3112 2015-05-19 13:17:35Z robert.buras $
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


# import GUI_definition
from libraddask.libradtran import GUI_definition
import io

class option():
    def __init__(self, name=None, group=None, 
             helpstr="", documentation="", tokens=[],
             gui_inputs=(),
             parents=[], non_parents=[], non_parent_exceptions=[], childs=[],
             speaker=None, enable_values=None,
             mystic=False, threedmystic=False, islidar=False, developer=False,
             extra_dependencies=[], 
             plot=None, showInGui = True,
             continious_update=False,mandatory=False,non_unique=False):

        assert not name is None
        assert not group is None

        assert hasattr(gui_inputs, "__iter__"), "Oops! gui_inputs" \
            " for option {0} must be a tuple/list".format(name)
        tokens = _checkForDictionary(tokens)
        
        if tokens == -1:
            print('''ERROR: found in definition of %s:\n
            Please check your tokens and settings again''' %(name))

        self.parents = parents
        self.non_parents = non_parents
        self.non_parent_exceptions = non_parent_exceptions
        self.childs = childs
        self.name = name
        self.help = helpstr
        self.plot = plot
        self.tokens = tokens
        self.showInGui = showInGui
        self.continious_update = continious_update
        self._mandatory=mandatory
        self._extra_dependencies=extra_dependencies
        self.non_unique=non_unique

        # Please use empty lists
        if self.childs == '': self.childs = []
        if self.non_parents == '': self.non_parents = []
        if self.non_parent_exceptions == '': self.non_parent_exceptions = []

        self.dependencies = list(self.childs + self.non_parents + self.non_parent_exceptions + extra_dependencies)
        if self._mandatory and not name in self.dependencies: self.dependencies.append(self.name)
            
        for inp in gui_inputs:
            assert inp.__class__.__name__ in GUI_definition.__all__
        if not gui_inputs:    gui_inputs = self._set_gui_input(gui_inputs,tokens)
        self.gui_inputs = self._addBooleanToGuiInput(gui_inputs)

        self.dict = {
            'name'         : name,          #as defined in uvspec_lex.l
            'group'          : group,             # e.g. WC
            'help'           : helpstr,           # Help string (short), could appear as pop-up in GUI
            'documentation'  : documentation,     # Full documentation 
            'tokens'     : tokens,          # Variable in uvspec inp_struct to change
            'parents'        : parents,           # (specifies which options must also be defined together with this option)
            # One of the options must also be defined with this options
            'non_parents'    : non_parents,       # specifies which options must not be defined together with this option 
            'non_parent_exceptions'    : non_parent_exceptions,       # specifies which options inside non_parents should be ignored
            'childs'         : childs,          # (specifies which options can be defined with this option)
            # Options which will be unlocked when defining this option
            'mystic'     : mystic,          # mystic option
            'threedmystic'     : threedmystic,      # 3D mystic option
            'islidar'     : islidar,           # lidar option
            'developer'     : developer,         # developer option, undocumented for esasLight
            'plot'           : plot,              # Setup plotting for options which should be plotted
            }

        self.canEnable = self._canEnable

        if  speaker and enable_values:
            self._speaker = speaker
            assert not isinstance(enable_values, str), "Missing comma in one-item-tuple?"
            self._enable_values = enable_values
            self.canEnable = self._canEnableContinousOption
        if extra_dependencies:
            self.isMandatory = self._isMandatoryMixedOption
        

    # Pretend to be a dictionary to avoid breaking old code
    def __getitem__(self, *args, **kwargs):
        return self.dict.__getitem__(*args, **kwargs)
    def __setitem__(self, *args, **kwargs):
        return self.dict.__setitem__(*args, **kwargs)
    def __contains__(self, *args, **kwargs):
        return self.dict.__contains__(*args, **kwargs)
    def get(self, *args, **kwargs):
        return self.dict.get(*args, **kwargs)

    def _canEnableContinousOption(self, is_set, get_value):
        """
        Remember that self.speaker must be a subclass of
        continious_option.
        """
        r = self._canEnable(is_set, get_value)
        if r and is_set(self._speaker) and  \
                get_value(self._speaker)[0] in self._enable_values:
            return r
        else:
            return False

    def _canEnable(self, is_set, get_value):
        """
        Tells the GUI wether the option should be enabled or disabled.
        Returns True if the option should be enabled and False if it
        should be disabled.

        is_set is a function that returns True if an option is enabled
        and has been edited by the user, else False. It takes one
        argument, the name of an option as a string.

        get_value returns the current value of an option

        This is used to implement the logic in the GUI. If more
        complex logic than the parent, non-parent, children logic is
        needed this function should be overloaded.

        Remember to update the dependency tuple, a tuple of options
        which should be enabled or disabled depending on if this
        option is set.
        """
        parents = any([is_set(parent) for parent in self.parents]) \
            or not self.parents
        non_parents = all([(not is_set(non_parent) or self.non_parent_exceptions.count(non_parent) or non_parent==self.name) \
                       for non_parent in self.non_parents]) \
                       or not self.non_parents
        return parents and non_parents

    def isMandatory(self, is_set, get_value):
        """
        Returns True for mandatory options. Similar to canEnable.
        """
        if self._mandatory and not is_set(self.name):    return True
        return False

    def _isMandatoryMixedOption(self, is_set, get_value):
        cond = [is_set(opt) for opt in self._extra_dependencies]
        if all(cond):
            return False
        elif any(cond):
            return True
        else:
            return False
        

    def _set_gui_input(self,gui_inputs,tokens):

        if not self.showInGui:
            return gui_inputs
        for inp in tokens:
            try:
                name = inp.gui_name
            except KeyError:
                pass
            if not name:    name = inp.get('name')
            try:
                vr = inp.get('valid_range')
            except KeyError:
                vr = None

            if isinstance(inp, addSetting):
                continue
            elif isinstance(inp, addLogical):
                gui_inp = (GUI_definition.ListInput(name=name,valid_range=inp.get('valid_range'),optional=inp.get('optional'),default=inp.get('default'),logical_file=inp.get('logical_file')),)
            elif isinstance(inp, addToken):
                dtype = inp.get('datatype')
                if dtype == float or dtype==Double:
                    if not vr: vr = (-1e99, 1e99)
                    gui_inp = (GUI_definition.FloatInput(name=name,optional=inp.get('optional'),valid_range=vr,default=inp.get('default')),)
                elif dtype == int:
                    if not vr: vr = (-1e99, 1e99)
                    gui_inp = (GUI_definition.IntegerInput(name=name,valid_range=vr,optional=inp.get('optional'),default=inp.get('default')),)
                elif vr:
                    gui_inp = (GUI_definition.ListInput(name=name,valid_range=inp.get('valid_range'),optional=inp.get('optional'),default=inp.get('default'),logical_file=inp.get('logical_file')),)
                # elif dtype == file:
                elif type(inp) is type(io.IOBase):
                    gui_inp = ( GUI_definition.FileInput(name=name,optional=inp.get('optional')) ,)
                else:    gui_inp = (GUI_definition.TextInput(name=name,optional=inp.get('optional')),)
            gui_inputs = gui_inputs.__add__(gui_inp)
        return gui_inputs

    def _addBooleanToGuiInput(self,gui_inputs):
        if not self.showInGui:
            return ()
        for inp in gui_inputs:
            if not inp.optional or inp.__class__ == GUI_definition.BooleanInput:
                return gui_inputs
        return ( GUI_definition.BooleanInput(name=''), ).__add__(gui_inputs)


class Dimension():
    """
    Options which can take dimensions (number+word) as argument ( 1D, 3D )
    """
    def __init__(self):
        self.valid_range = ["1d","3d"]
    def get_valid_range(self):
        return self.valid_range

class ProfileType():
    """
    Options which can take several profile files as argument (e.g. 1D, 3D, moments, ipa_files)
    """
    def __init__(self):
        self.valid_range = ["1d","3d","ipa_files","moments"]
    def get_valid_range(self):
        return self.valid_range

class CaothType():
    """
    Options which can take several profile as argument (e.g. wc, ic or any other profile)
    """
    def __init__(self,caoth=None):
        self.caoth = caoth
    def get_caoth(self):
        return self.caoth
class CaothoffType():
    """
    Quick fix for new option names no_scattering and no_absorption
    """

class Double(float):
    """Double for c allocation double"""

class SignedFloats():
    """Signed floats for c allocation multiple floats"""

class Integers():
    """Integers for c allocation multiple integers"""

class VariableNumberOfLines():
    pass

# valid_datatypes = ( # i.e. datatypes the GUI support (soon)
#     ProfileType,
#     CaothType,
#     CaothoffType,
#     Double,
#     SignedFloats,
#     VariableNumberOfLines,
    
#     # In my opinion these should be replaced by DatatFile, Integer etc., since as
#     # far as I know, the only property used is that file, int etc. are distinct
#     # and globally avilable variables. It is completely irrelevant, that they
#     # happen to be internal datatypes of python.
#     file, int, float, long, str
# )

class not_yet_lex2py_option(option):
    """
    Quick fix for options which are currently hard to implement in
    the python-structures.
    """
    def writeLex(self):
        return ""

class addInput():
    def __init__(self, name='', default=None, gui_name=None):
        self.name = name
        self.default = default
        self.gui_name = gui_name

        self.dict = {
            'name'        : self.name,        # variable name to set as defined in uvspec_lex.l
            'default'    : self.default,    # default value for GUI AND uvspec_lex.l
        }

    def __getitem__(self, *args, **kwargs):
        return self.dict.__getitem__(*args, **kwargs)

    def get(self, *args, **kwargs):
        return self.dict.get(*args, **kwargs)


class addToken(addInput):
    def __init__(self, datatype=None, valid_range=None, optional=False, **kwargs):
        addInput.__init__(self,**kwargs)

        logical_file = False
        if valid_range:
            # if isinstance(valid_range[-1], io.IOBase):
            # if valid_range[-1] == file:
            if type(valid_range[-1]) == type(io.IOBase):
                valid_range.pop(-1)
                logical_file = True
        self.datatype = datatype
        self.valid_range = valid_range
        self.optional = optional

        self.dict['datatype']    = self.datatype    # datatype of input token, possible types are float, int, str, file (for GUI)
        self.dict['valid_range'] = self.valid_range    # important for GUI and unsigned floats
        self.dict['optional']     = self.optional    # optional argument
        self.dict['logical_file']= logical_file        #if not logical argument, then input is filename
        
class addLogical(addToken):
    def __init__(self, logicals=[], setting='', logical_file=False, datatype=str, **kwargs):
        
        if type(logicals[-1]) is type(io.IOBase):
            logicals.pop(-1)
            logical_file = True
        elif isinstance( logicals[-1], CaothType ):
            logicals.pop(-1)
            caoth = True

        if   sorted([ dim.lower() for dim in logicals ]) == ProfileType().get_valid_range():
            datatype = ProfileType
        elif sorted([ dim.lower() for dim in logicals ]) ==   Dimension().get_valid_range():
            datatype = Dimension

        addToken.__init__(self,datatype=datatype,valid_range=logicals,**kwargs)
            
        self.setting = setting
        self.logical_file = logical_file

        if isinstance(logicals, ProfileType): self.datatype=ProfileType

        self.dict['logicals']        = self.valid_range    # variable name to set as defined in uvspec_lex.l
        self.dict['setting']        = self.setting         # value variable should be set to as in uvspec_lex.l
        self.dict['destination']    = self.name         # default value for uvspec_lex.l
        self.dict['logical_file']    = self.logical_file    #if not logical argument, then input is filename
        

class addSetting(addInput):
    def __init__(self, setting='', **kwargs):
        addInput.__init__(self,**kwargs)

        if isinstance( setting, CaothType ):
            if setting.get_caoth():
                self.setting = 'get_caoth_index(&Input.caoth,&Input.n_caoth,"%s",1)' %( setting.get_caoth() )
            else:
                print("Error: Must define caoth wc/ic")
        else:
            self.setting = setting

        self.dict['setting']    = self.setting     # value variable should be set to as in uvspec_lex.l
        

def addPlot(plot_type, optional_args={}):
    return {
        'plot_type'     : plot_type,     # One of '2D' or 'map'
        'optional_args' : optional_args # A dictionary of other information required by the specific plot_types. Currently it should be empty for map and contain a list of column_names for 2D.
        }

def _checkForDictionary(argument):
    
    if not argument: return []
    elif isinstance(argument,list): 
        for arg in argument:
            if   isinstance(arg,addToken):    continue
            elif isinstance(arg,addSetting):    continue
            elif isinstance(arg,addLogical):    continue
            else:    return -1 
            
        return argument
                
    elif isinstance(argument, addToken):    return [argument]
    elif isinstance(argument, addSetting):    return [argument]
    elif isinstance(argument, addLogical):    return [argument]
    else: return -1    
