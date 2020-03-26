__all__ = ["FileInput", "FloatInput", "TextInput", "IntegerInput", "ListInput", "IntegerListInput", "BooleanInput", "VariableNumberOfLinesInput"]


class Input():
    def __init__(self, name=None, optional=False):
        """
        name = displayed above input in GUI for options with multiple inputs
        """
        assert not name is None
        self.name = name
        self.optional = optional

class NumberInput(Input):
    def __init__(self, default=None, valid_range=(-1e99, 1e99),**kwargs):
        Input.__init__(self, **kwargs)

        # This should be removed when/if the option files are cleaned up
        if default in ("NOT_DEFINED_INTEGER", "NOT_DEFINED_FLOAT"):
            default = None

        self.default = default
        self.valid_range = valid_range

class FileInput(Input):
    def __init__(self, **kwargs):
        Input.__init__(self, **kwargs)

class FloatInput(NumberInput):
    pass

class TextInput(Input):
    def __init__(self, default=None, **kwargs):
        Input.__init__(self, **kwargs)
        self.default = default

class IntegerInput(NumberInput):
    def __init__(self, default=None, **kwargs):
        # This should be removed when/if the option files are cleaned up
        if default in ("NOT_DEFINED_INTEGER", "NOT_DEFINED_FLOAT"):
            default = None

        if not default is None:
            assert type(default) == int, \
                "Default of integer input must be an integer!"
        NumberInput.__init__(self, default=default,  **kwargs)

class ListInput(Input):
    """ Valid inputs are one among a list of strings """

    def __init__(self, default=None, valid_range=None, optional=False,logical_file=False,**kwargs):
        Input.__init__(self, optional=optional, **kwargs)
        
        assert not valid_range is None, "You must provide a range of choices!"

        self.valid_range = []

        for val in valid_range:
            if isinstance(val,str):
                self.valid_range.append( val.lower() )
        else:
            self.valid_range.append( val )

        if optional:
            if self.valid_range.count(""):
                self.valid_range.remove("")
            self.valid_range.insert(0,"")

        if isinstance(default,str): 
            default=default.lower()
        if default is None: 
            default = self.valid_range[0]
        assert default in self.valid_range, "Default not among valid options!"
        self.default = default
        self.logical_file=logical_file

class IntegerListInput(ListInput):
    def __init__(self, **kwargs):
        ListInput.__init__(self, **kwargs)

        self.default = str(self.default)
        self.valid_range = tuple([str(i) for i in self.valid_range])
    
class BooleanInput(Input):
    pass

class VariableNumberOfLinesInput(Input):
    def __init__(self, valid_range=None, **kwargs):
        Input.__init__(self, **kwargs)
        assert not valid_range is None
        self.valid_range = valid_range

