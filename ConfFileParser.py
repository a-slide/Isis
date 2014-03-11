class ConfFileParser:
    """This can be used to parse file containing a list of variables and its associated values
    The separator beetween vaiable and value must be the same for all field and can be customized by users (by default =)
    the class store the list of variable in a dictionnary named conf_dict"""
    
# FONDAMENTAL METHODS
    
    # Object constructor initialized the path of the field to parse and a custom separator
    # Create dictionnary of parameters
    def __init__ (self, filename):
        self.conf_dict = self._makeConfDict (filename)
    
    # Short description string returned by print and str 
    def __str__(self):
        return "Instance of ConfFileParser\n"
    
    # Long description string used by interpreter and repr
    def __repr__ (self):
        result = "Instance of ConfFileParser\n\n"
        for key,value in self.conf_dict.items():
            result = result + key + "\t" + value + "\n"
            
        return result
        
# GETERS
    
    # Grant acces to the complete dictionary
    def getDict (self):
        return self.conf_dict
    
    # Give acces to individual values in conf_dict by using its key name.
    def getVar (self, varkey ):
        return self.conf_dict[varkey]
        
# PRIVATE SUPPORT METHODS
    
    # Return a dictionnary containing all parameters assosiated with its value
    def _makeConfDict (self, filename):
        return {name : value for name, value in self._makeConfList(filename)}
    
    # Return a 2 element list for each value containing lines
    def _makeConfList (self, filename):
        return [line[0:2] for line in self._openFile(filename)]
    
    # Open the file containing configurations and return the list of splitted value containing lines
    # TODO = define a alternative method of initialisation in case of failure to read the file
    def _openFile (self, filename):
        try: # bloc try pour gerer l'erreur d'ouverture de fichier
            with open (filename) as filename:
                return [line.split() for line in filename if line[0] != '#' and line[0] != '\n']
            
        except IOError:
            print '\n', filename, 'is not readable. The file will be ignored\n'
            return {}
        
