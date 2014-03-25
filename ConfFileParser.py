class ConfFileParser:
    """This can be used to parse file containing a list of variables and its associated values
    The separator beetween vaiable and value must be the same for all field and can be customized by users (by default =)
    the class store the list of variable in a dictionnary named conf_dict"""


########################################################################################################################
#   FONDAMENTAL METHODS
########################################################################################################################

    # Object constructor initialized with the path of the configuration file to parse
    # In the configuration files parameter names must start the line and be follow by the value. Both filed needs
    # to be separated by a blank space (space or tab)
    # Empty lines and lines starting by # will be filtered out

    def __init__ (self, filename):
        self.d = self._make_conf_dict (filename)


    def __repr__ (self):
        """Long description string used by interpreter and repr"""
        key_list = self.d.keys()
        key_list.sort()

        result = "<Instance of ConfFileParser\n"
        for key in key_list:
            result += "{0} :\t{1}\n".format(key,self.d[key])
        result += ">"
        return result


    def __str__(self):
        """Short representation"""
        return  "<Instance of ConfFileParser\n"


########################################################################################################################
#   GETERS
########################################################################################################################

    # Grant acces to the complete dictionary
    def getDict (self):
        return self.d

    # Give acces to individual values in conf_dict by using its key name.
    def get (self, key ):
        return self.d[key]


########################################################################################################################
#   PRIVATE METHODS
########################################################################################################################


    # Return a dictionnary containing all parameters assosiated with its value
    def _make_conf_dict (self, filename):
        return {name : value for name, value in self._make_conf_list(filename)}

    # Return a 2 element list for each value containing lines
    def _make_conf_list (self, filename):
        return [line[0:2] for line in self._open_file(filename)]

    # Open the file containing configurations and return the list of splitted value containing lines
    def _open_file (self, filename):
        try:
            with open (filename) as filename:
                return [line.split() for line in filename if line[0] != '#' and line[0] != '\n']

        except IOError:
            print ("Error : " + filename + " is not readable. Can't parse the configuration file")
            return {}

