from ConfFileParser import ConfFileParser
from sys import maxsize
import optparse

class IsisConf:
    """Class description"""

########################################################################################################################
#   FONDAMENTAL METHODS
########################################################################################################################

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
    # TODO Define possible remaining exceptions?
    def __init__(self):

        # dictionnary to store validated parameters from command line arguments and conf_file parsing
        self.d = {}

        #### Parse command line argument ####
        # add hg_filename, vg_filename, conf_filename and output_prefix entries in dictionnary (dict concatenation operator)
        self._optparser("hg_filename", "vg_filename", "conf_filename", "output_prefix")

        #### Extract configuration parameters from the conf file by creating a instance of ConfFileParser ####
        self.conf = ConfFileParser(self.d["conf_filename"])

        # Overall quantity of reads
        self.d["read_num"] = self._import_integer("read_num",1,None,100)
        # Read length in pb
        self.d["read_len"] = self._import_integer("read_len",1,None,150)

        freq_hg, freq_vg, freq_tj, freq_fj =  self._import_freq(0.45, 0.45, 0.09, 0.01)
        # host genome frequency
        self.d["nread_hg"] = int (freq_hg * self.d["read_num"])
        # Virus genome frequency
        self.d["nread_vg"] = int (freq_vg * self.d["read_num"])
        # True junctions frequency
        self.d["nread_tj"] = int (freq_tj * self.d["read_num"])
        # False junctions frequency
        self.d["nread_fj"] = int (freq_fj * self.d["read_num"])

        # Mean number of sampling in true junctions and false junction
        samp_tj = self._import_integer("samp_tj",0,self.d["nread_tj"], 100) # Should be more than 10
        samp_fj = self._import_integer("samp_fj",0,self.d["nread_fj"], 1) # Should be 1 or 0

        # Number of uniq true junction in the pool of junction to be generated
        self.d["uniq_tj"] = int (float(self.d["nread_tj"])/samp_tj)
        # Number of uniq true junction in the pool of junction to be generated
        self.d["uniq_fj"] = int (float(self.d["nread_fj"])/samp_fj)

        ## Allow or forbid sampling in repeat regions (lowercase characters), ambigous DNA (Ambigous IUPAC code) and duplicated reads.
        self.d["repeats"] = self._import_boolean("repeats", False)
        self.d["ambigous"] = self._import_boolean("ambigous", False)
        self.d["duplicate"] = self._import_boolean("duplicate", False)

        # Sequencing mode. If true = pair end alse = single end
        self.d["pair_end"] = self._import_boolean("pair_end", True)

        ## Options of fragment sonication distribution to set up if paired end mode
        # Maximal Minimal and mode of sonication size max > mode > min
        self.d["sonic_min"] = self._import_integer ("sonic_min", self.d["read_len"], None, 200)
        self.d["sonic_mode"] = self._import_integer ("sonic_mode", self.d["sonic_min"], None, 400)
        self.d["sonic_max"] = self._import_integer ("sonic_max", self.d["sonic_mode"], None, 1000)

        # Certainty of the sonication smire: 5 (wide peak) to 50 (thin peak). Recommanded = 10
        self.d["sonic_certainty"] = self._import_integer("sonic_certainty", 5, 50, 10)

        # Quality score scale of fastq output file
        self.d["qual_scale"] = self._import_string("qual_scale", ["illumina", "solexa", "sanger"])

        # Mean quality score of fastq output file at start, middle and end of reads between 0 and 40
        self.d["qual_start"] = self._import_integer ("qual_start",0, 40, 30)
        self.d["qual_mid"] = self._import_integer ("qual_mid",0, 40, 35)
        self.d["qual_end"] = self._import_integer ("qual_end",0, 40, 30)

        # Range of quality variation at start, middle and end of reads between 0 and 20
        self.d["qvar_start"] = self._import_integer ("qvar_start", 0, 20, 5)
        self.d["qvar_mid"] = self._import_integer ("qvar_mid", 0, 20, 2)
        self.d["qvar_end"] = self._import_integer ("qvar_end", 0, 20, 10)

        # Number of available threads for distributed computing
        self.d["nb_thread"] = self._import_integer ("nb_thread", 1, None, 1)


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
    def __repr__ (self):
        """Long description string used by interpreter and repr"""
        # A key list is created to output a sorted list of dict entries
        key_list = self.d.keys()
        key_list.sort()

        result = "<Instance of IsisConf>\n"
        for key in key_list:
            result += "{0} :\t{1}\n".format(key,self.d[key])
        result += ">"
        return result


    def __str__(self):
        """Short representation"""
        return  "<Instance of IsisConf>\n"

########################################################################################################################
#   GETERS
########################################################################################################################

    def getDict (self):
        """Grant acces to the complete dictionary"""
        return self.d


    def get (self, key ):
        """Give acces to individual values in conf_dict by using its key name"""
        return self.d[key]


########################################################################################################################
#   PRIVATE METHODS
########################################################################################################################

    def _optparser(self, hg, vg, conf, output):
        """Parse command line arguments prompt and verify the filename"""
        usage_string = "%prog -H Host_genome.fasta[.gz] -V Viral_genome.fasta[.gz] -C Conf_file.txt [-o Output_prefix]"
        optparser = optparse.OptionParser(usage = usage_string)
        optparser.add_option( '-H', '--Host_genome', dest = hg,
        help = "Path of the fasta file containing the host genome sequence. The source can be gziped")
        optparser.add_option( '-V', '--Virus_genome', dest = vg,
        help = "Path of the fasta file containing the viruq genome sequence. The source can be gziped")
        optparser.add_option( '-C', '--Conf_file', dest = conf,
        help = "Path of the configuration text file. Option values can be changed but the template should remain unmodified")
        optparser.add_option( '-o', '-O', '--output', default ="out", dest = output,
        help = "Facultative option to indicate the name of the output prefix. By default : out")

        # Parse arg and return a dictionnary_like object of options
        options, args = optparser.parse_args()
        # Verify each path and enter them in the object conf dictionnary
        self.d[hg] = self._check_file (options.hg_filename, hg)
        self.d[vg] = self._check_file (options.vg_filename, vg)
        self.d[conf] = self._check_file (options.conf_filename, conf)
        # No need to check this options
        self.d[output] = options.output_prefix


    def _check_file (self, path, key):
        """Try to path from the opt dictionnary"""
        if not path:
            path = raw_input ("{0} is a mandatory parameter. Please enter a valid path : ".format(key))

        while True:
            try:
                with open (path) as filename:
                    print (path + " is valid\n")
                    return path
            except IOError:
                print ("Error : " + path + " can't be read")
                path = raw_input ("Please enter a valid path :  ")


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
    def _import_integer(self, name, min = None, max = None, default = None):
        """Import an integer from conf dict and validate its value"""
        try:
            val = int(self.conf.get(name))
            if self._valid_int(val, min, max):
                return val

        except KeyError:
            print ("Error : {0} was not found in {1})".format(name, self.d["conf_filename"]))
            return self._define_int (min, max, default, name)
        except ValueError:
            print ("Error : {0} value is not a invalid integer or is outside of boundaries (min = {1}, max = {2})".format(name, min, max))
            return self._define_int (min, max, default, name)


    def _define_int (self, min, max, default, name):
        while (True):
            if raw_input ("Would you like to use the default value < {0} > (Y/N)?\n".format(default)) in 'Yy':
                return default
            try:
                val = int(raw_input ("Please enter a new value :  "))
                if self._valid_int(val, min, max):
                    return val
            # if the value is still not suitable print a meassage and try again
            except ValueError:
                print ("Error : {0} value is not a invalid integer or is outside of boundaries (min = {1}, max = {2})".format(name, min, max))


    def _valid_int(self, val, min, max):
        """Predicate raising a ValueError exception if a value is not in the given range"""
        max = maxsize if max == None else max
        min = -maxsize if min == None else min

        if val >= min and val <= max:
            return True
        else:
            raise ValueError


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
    def _import_freq(self,default_hg, default_vg, default_tj, default_fj):
        """Import frequencies from conf dict and validate their values"""
        try:
            hg = float(self.conf.get("freq_hg"))
            vg = float(self.conf.get("freq_vg"))
            tj = float(self.conf.get("freq_tj"))
            fj = float(self.conf.get("freq_fj"))
            if self._valid_freq([hg, vg, tj, fj], 1):
                return (hg, vg, tj, fj)

        except KeyError:
            print ("Error : Frequency values were not found in {0})".format(self.d["conf_filename"]))
            return self._define_freq (default_hg, default_vg, default_tj, default_fj)
        except ValueError:
            print ("Error : Frequencies are not numeric values or the sum of frequencies is not egual to 1")
            return self._define_freq (default_hg, default_vg, default_tj, default_fj)


    def _define_freq (self, default_hg, default_vg, default_tj, default_fj):
        while (True):
            if raw_input ("Would you like to use the default values < {0}, {1}, {2}, {3} > (Y/N)?\n".format(default_hg, default_vg, default_tj, default_fj)) in 'Yy':
                return (default_hg, default_vg, default_tj, default_fj)
            try:
                hg = float(raw_input ("Please enter host genome frequency :  "))
                vg = float(raw_input ("Please enter viral genome frequency :  "))
                tj = float(raw_input ("Please enter true junction frequency :  "))
                fj = float(raw_input ("Please enter false junction frequency :  "))

                if self._valid_freq([hg, vg, tj, fj], 1):
                    return (hg, vg, tj, fj)
            # if the value is still not suitable print a meassage and try again
            except ValueError:
                print ("Error : Frequencies are not numeric values or the sum of frequencies is not egual to 1")


    def _valid_freq (self, val_list, target_sum):
        """Predicate raising a ValueError exception if the sum of values in val_list is not equal to target_sum"""
        val_sum = 0
        for val in val_list:
           val_sum += val

        if val_sum == target_sum:
            return True
        else:
            raise ValueError


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
    def _import_string(self,name,allowed_entries):
        """Import a string from conf dict and validate its value"""
        try:
            val = self.conf.get(name).lower()
            if self._valid_string(val, allowed_entries):
                return val
        except KeyError:
            print ("Error : {0} was not found in {1})".format(name, self.d["conf_filename"]))
            return self._define_string (allowed_entries, name)
        except ValueError:
            print ("Error : {0} value is not in the list of allowed values:\n {1}".format(name, allowed_entries))
            return self._define_string (allowed_entries, name)


    def _define_string (self, allowed_entries, name):
        while (True):
            if raw_input ("Would you like to use the default value < {0} > (Y/N)?\n".format(allowed_entries[0])) in 'Yy':
                return allowed_entries[0]
            try:
                val = raw_input ("Please enter a new value :  ")
                if self._valid_string(val, allowed_entries):
                    return val
            # if the value is still not suitable print a meassage and try again
            except ValueError:
                print ("Error : {0} value is not in the list of allowed values:\n {1}".format(name, allowed_entries))


    def _valid_string(self, val, allowed_entries):
        """Predicate raising a ValueError exception if the string is not in the list of allowed entries"""
        if val in allowed_entries:
            return True
        else:
            raise ValueError


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
    def _import_boolean (self, name, default):
        """Import a boolean from conf dict and validate its value"""
        true_list = ['True', 'true', 1, 't', 'T']
        false_list = ['False', 'false', 0, 'f', 'F']
        try:
            val = self.conf.get(name)
            if self._valid_string(val,true_list+false_list):  # Use the same predicate as for strings
                return True if val in true_list else False
        except KeyError:
            print ("Error : {0} was not found in {1})".format(name, self.d["conf_filename"]))
            return self._define_bool (true_list, false_list, default, name)
        except ValueError:
            print ("Error : {0} value is not in the list of allowed values:\n {1}".format(name, true_list+false_list))
            return self._define_bool (true_list, false_list, default, name)


    def _define_bool (self, true_list, false_list, default, name):
        while (True):
            if raw_input ("Would you like to use the default value < {0} > (Y/N)?\n".format(default)) in 'Yy':
                return default
            try:
                val = raw_input ("Please define manualy your own value :  ")
                if self._valid_string(val,true_list+false_list): # Use the same
                    return True if val in true_list else False
            # if the value is still not suitable print a meassage and try again
            except ValueError:
                print ("Error : {0} value is not in the list of allowed values:\n {1}".format(name, true_list+false_list))


