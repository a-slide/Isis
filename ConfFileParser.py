####################################################################################################

class ConfFileParser:
    """This class can be used to parse a configuration, file containing
    a list of variables names and its associated values separated by a
    blank space. the class store the list of variable in a dictionnary
    named conf_dict"""

###    FONDAMENTAL METHODS    ###

    def __init__ (self, filename):
        """ Parse a configuration file. Each parameter names must start
        the line and be follow by the value separated by a blank space.
        Empty lines and lines starting by # are filtered out.  
        """
        self.d = self._make_conf_dict (filename)

    def __repr__ (self):
        """Long description string used by interpreter and repr
        """
        # A key list is created to output a sorted list of dict entries
        key_list = self.d.keys()
        key_list.sort()

        result = self.__str__()
        for key in key_list:
            result += "{} :\t{}\n".format(key,self.d[key])
        result += ">"
        return result

    def __str__(self):
        """Short representation
        """
        return "<Instance of " + self.__module__ + ">\n"

###    GETERS    ###

    def getDict (self):
        return self.d

    def get (self, key ):
        return self.d[key]

###    PRIVATE METHODS    ###

    def _make_conf_dict (self, filename):
        """ Return a dictionnary containing all parameters associated
        with its value
        """
        return {name : value for name, value in self._make_conf_list(filename)}

    def _make_conf_list (self, filename):
        """ Return a 2 element list for each value containing lines
        """
        return [line[0:2] for line in self._open_file(filename)]

    def _open_file (self, filename):
        """ Open the file containing configurations and return the list
        of splitted value containing lines.
        """
        try:
            with open (filename) as filename:
                return [line.split() for line in filename if line[0] != '#' and line[0] != '\n']

        except IOError:
            print ("Error : " + filename + " is not readable. Can't parse the configuration file")
            return {}
            
