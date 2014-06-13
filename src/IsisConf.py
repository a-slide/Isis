"""
@package    IsisConf
@brief      **Configuration file and command line argument parser/verifier**
Import command line arguments with optparse and parse configuration file with
ConfigParser. Values are imported in a dictionnary after verification. If a
value is invalid, a message describing the error will be printed before program
exit.
@copyright  [GNU General Public License v2](http://www.gnu.org/licenses/gpl-2.0.html)
@author     Adrien Leger - 2014
* <adrien.leger@gmail.com>
* <adrien.leger@inserm.fr>
* <adrien.leger@univ-nantes.fr>
* [Github](https://github.com/a-slide)
* [Atlantic Gene Therapies - INSERM 1089] (http://www.atlantic-gene-therapies.fr/)
"""

#~~~~~~~PACKAGE IMPORTS~~~~~~~#

# Standard library packages
from sys import exit as sys_exit
from imp import find_module
import ConfigParser
import optparse

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class IsisConfException(Exception):
    """
    @class IsisConfException
    @brief Custom Exception Class to handle error while parsing Isis options
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    def __init__(self, msg): # Object constructor initialized with a custom user message
        """
        Catch an IsisConfException, print a description of the error and exit without python
        error printing traceback
        @param msg spefific error message
        """
        err_msg = "IsisConfException : An error occured while parsing Isis options!\n"
        err_msg += "\t{}.\n\tPlease ajust your settings\n".format(msg)
        print (err_msg)
        sys_exit ()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class IsisConf(object):
    """
    @class IsisConf
    @brief Import command line arguments with optparse and parse configuration file with
    ConfigParser. Values are imported in a dictionnary after verification. If a value is invalid,
    a message describing the error will be printed before program exit.
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~FONDAMENTAL METHODS~~~~~~~#

    def __init__(self, program_name, program_version):
        """
        Import command line arguments with optparse and parse configuration file
        with ConfigParser. Values are imported in a dictionnary after verification.
        If a value is invalid, a message describing the error will be printed
        before program exit.
        @param program_name Name of the program
        @param program_version Version of the program
        """
        print "Set-up of the program parameters..."
        ## Dictionnary to store valid parameters from command line arguments and configuration file
        self.d = {}

        #~~~Parse command line argument~~~#
        print "\tParsing command-line arguments"
        # add hg_filename, vg_filename, conf_filename and output_prefix entries to self.d
        self.d.update (self._optparser(program_name, program_version))

        #~~~Extract configuration parameters from the conf file~~~#
        print "\tParsing configuration file"
        ## Instance of RawConfigParser
        self.config = ConfigParser.RawConfigParser(allow_no_value = False)

        self.config.read(self.d["conf_file"])

        #~GENERAL SECTION~#
        self.d.update (self._get_int("General", "read_num", 1, None))
        self.d.update (self._get_int("General", "read_len", 1, None))
        self.d.update (self._get_float("General", "mut_freq", 0, 1))
        self.d.update (self._get_bool("General", "repeats"))
        self.d.update (self._get_bool("General", "ambiguous"))
        self.d.update (self._get_bool("General", "graph"))
        self.d.update (self._get_bool("General", "report"))

        #~FREQUENCY SECTION~#
        freq_host = self._get_float ("Frequency", "freq_host", 0, 1).values()[0]
        freq_virus = self._get_float ("Frequency", "freq_virus", 0, 1).values()[0]
        freq_tjun = self._get_float ("Frequency", "freq_tjun", 0, 1).values()[0]
        freq_fjun = self._get_float ("Frequency", "freq_fjun", 0, 1).values()[0]
        if (freq_host+freq_virus+freq_tjun+freq_fjun) != 1:
            raise IsisConfException ("The sum of read number frequencies is not equal to 1")

        # Calculate the number of reads from validated frequencies
        self.d["nread_host"] = int (freq_host * self.d["read_num"])
        self.d["nread_virus"] = int (freq_virus * self.d["read_num"])
        self.d["nread_tjun"] = int (freq_tjun * self.d["read_num"])
        self.d["nread_fjun"] = int (freq_fjun * self.d["read_num"])

        #~JUNCTION SECTION~#
        if self.d["pair"]:
            self.d.update (self._get_int ("Junction", "min_chimeric", 0, self.d["read_len"]))
        else:
            self.d.update (self._get_int ("Junction", "min_chimeric", 0, self.d["read_len"]/2))

        samp_tjun = self._get_float ("Junction", "samp_tjun", 0, self.d["nread_tjun"]).values()[0]
        samp_fjun = self._get_float ("Junction", "samp_fjun", 0, self.d["nread_fjun"]).values()[0]

        # Calculate the number of uniq junctions in ReferenceJunctions to be generated
        self.d["uniq_tjun"] = int (float(self.d["nread_tjun"])/samp_tjun)
        self.d["uniq_fjun"] = int (float(self.d["nread_fjun"])/samp_fjun)

        #~SONICATION SECTION~#
        min_size = self.d["read_len"] + self.d["min_chimeric"]
        self.d.update (self._get_int ("Sonication", "sonic_min", min_size, None))
        self.d.update (self._get_int ("Sonication","sonic_mode", self.d["sonic_min"], None))
        self.d.update (self._get_int ("Sonication","sonic_max", self.d["sonic_mode"], None))
        self.d.update (self._get_int ("Sonication","sonic_certainty", 5, 50))

         #~QUALITY SECTION~#
        self.d.update (self._get_str ( "Quality", "qual_scale", ["fastq-sanger", "fastq-solexa", "fastq-illumina"]))
        self.d.update (self._get_str ( "Quality", "qual_range", ["very-good", "good", "medium", "bad", "very-bad"]))
        
        #~~~Check third party dependencies~~~#
        print "\tChecking third party dependencies"
        
        try:
            find_module('Bio')
            print("\t\tBiopython package is available")
        except ImportError:
            raise IsisConfException ("Biopython package is required")
        
        if self.d["graph"]:
            try:
                find_module('matplotlib')
                print("\t\tmatplotlib package is available")
            except ImportError:
                print ("\t\tmatplotlib package is required for graphical output")
                print ("\t\t\tSwitching to non graphical mode")
                self.d["graph"] = False
                
                
    def __repr__ (self):
        # A key list is created to output a sorted list of dict entries
        key_list = self.d.keys()
        key_list.sort()

        result = self.__str__()
        for key in key_list:
            result += "{} :\t{}\n".format(key,self.d[key])
        return result

    def __str__(self):
        return "<Instance of {} from {} >".format(self.__class__.__name__, self.__module__)

    #~~~~~~~ACCESS METHODS~~~~~~~#

    def getDict (self):
        return self.d

    def get (self, key):
        return self.d[key]

    #~~~~~~~PRIVATE METHODS~~~~~~~#

    def _optparser(self, program_name, program_version):
        """
        Parse command line arguments with optparse and verify the file validity
        @param program_name Name of the program
        @param program_version Version of the program
        @return A dictionnary containing:
        * Host and virus fasta file paths
        * Configuration file path
        * Sequencing mode : pair or single
        * Basename for output files
        """
        # Usage and version strings
        usage_string = "%prog -H Host_genome.fa[.gz] -V Viral_genome.fa[.gz] -C Conf_file.txt [-o Output_prefix] [-p |-s]"
        version_string = program_name + program_version
        optparser = optparse.OptionParser(usage = usage_string, version = version_string)

        # Define optparser options
        hstr = "Path of the fasta file containing the host genome sequence (can be gziped)"
        optparser.add_option( '-H', '--host_genome', dest="hg", help=hstr)
        hstr = "Path of the fasta file containing the viral genome sequence (can be gziped)"
        optparser.add_option( '-V', '--virus_genome', dest="vg", help=hstr)
        hstr = "Path of the configuration text file"
        optparser.add_option( '-C', '--conf_file', dest="conf", help=hstr)
        hstr = "Facultative option to indicate the name of the output prefix (default = out)"
        optparser.add_option( '-o', '--output', default="out", dest="output", help=hstr)
        htr = "Single end mode overwriten if -p option is also indicated"
        optparser.add_option( '-s', '--single', dest="single", action='store_true', help=hstr)
        hstr = "Pair end mode incompatible with -s option (default mode)"
        optparser.add_option( '-p', '--pair', dest="pair", action='store_true', help=hstr)

        # Parse arg and return a dictionnary_like object of options
        options, args = optparser.parse_args()

        # Validate option and generate a dictionnary
        arg_dict = {'host_genome' : self._check_file (options.hg, "host_genome"),
                    'virus_genome' : self._check_file (options.vg, "virus_genome "),
                    'conf_file' : self._check_file (options.conf, "conf_file"),
                    'basename' : options.output,
                    'pair' : self._check_mode (options.single, options.pair)}
                    
        return arg_dict

    def _check_file (self, path, descr):
        """
        Try to open the file at a given path. In case of impossibility an IsisConfException is raised
        @param path     Path of the file to verify
        @param descr    name of the option associated with the file
        @return         Valid path
        """
        if not path:
            raise IsisConfException ("{} is a mandatory command line argument".format(descr))

        try:
            handle = open(path, "r")
            handle.close()
            print ("\t\tValid file for {}".format (descr))
            return path

        except IOError:
            raise IsisConfException ("Error : " + path + " can't be read\n\
            Please enter a valid path")

    def _check_mode (self, single, pair):
        """
        Verify the compatibility of mode option. If both options had been selected an
        IsisConfException is raised
        @param single Boolean value of the single option
        @param pair Boolean value of the pair option
        @return A boolean indicating if pair end mode should be used
        """
        if pair and single:
            raise IsisConfException ("-p and -s are incompatible options")
        return False if single else True


    def _get_int(self, section, name, min=None, max=None):
        """
        Import an integer from self.config dict and verify its value. If an error occur, an
        IsisConfException is raised.
        @param section Name of the section in the configuration file where the option is (string)
        @param name Name of the option (string)
        @param min Minimal value (integer)
        @param max Maximal value (integer)
        @return An integer within the requested range
        """
        try:
            # Parse file with RawConfigParser
            val = self.config.getint (section, name)
            # Verify if the value is valid
            if min and val < min:
                raise IsisConfException ("{} value must be greater than {}.".format(name, min))
            elif max and val > max:
                raise IsisConfException ("{} value must be lower than {}.".format(name, max))
            else:
                return {name : val}

        # Handle ConfigFileParser errors
        except ConfigParser.NoOptionError:
            raise IsisConfException ("Option {} was not found in conf file.".format(name))
        except ConfigParser.NoSectionError:
            raise IsisConfException ("Section {} was not found in conf file.".format(section))
        except ValueError:
            raise IsisConfException ("{} value is not a valid integer".format(name))


    def _get_float(self, section, name, min=None, max=None):
        """
        Import an float from self.config dict and verify its value. If an error occur, an
        IsisConfException is raised.
        @param section Name of the section in the configuration file where the option is (string)
        @param name Name of the option (string)
        @param min Minimal value (float)
        @param max Maximal value (float)
        @return A float within the requested range
        """
        try:
            # Parse file with RawConfigParser
            val = self.config.getfloat (section, name)
            # Verify if the value is valid
            if min and val < min:
                raise IsisConfException ("{} value must be greater than {}.".format(name, min))
            elif max and val > max:
                raise IsisConfException ("{} value must be lower than {}.".format(name, max))
            else:
                return {name : val}

        # Handle ConfigFileParser errors
        except ConfigParser.NoOptionError:
            raise IsisConfException ("Option {} was not found in conf file.".format(name))
        except ConfigParser.NoSectionError:
            raise IsisConfException ("Section {} was not found in conf file.".format(section))
        except ValueError:
            raise IsisConfException ("{} value is not a valid float".format(name))


    def _get_bool(self, section, name):
        """
        Import a boolean from self.config dict and verify its value. If an error occur, an
        IsisConfException is raised
        @param section Name of the section in the configuration file where the option is (string)
        @param name Name of the option (string)
        @return A valid boolean value
        """
        try:
            # Parse file with RawConfigParser
            val = self.config.getboolean (section, name)
            return {name : val}

        # Handle ConfigFileParser errors
        except ConfigParser.NoOptionError:
            raise IsisConfException ("Option {} was not found in conf file.".format(name))
        except ConfigParser.NoSectionError:
            raise IsisConfException ("Section {} was not found in conf file.".format(section))
        except ValueError:
            raise IsisConfException ("{} value is not a valid boolean".format(name))


    def _get_str (self, section, name, allowed):
        """
        Import an string from self.config dict and verify that its value is in the list of
        autorized entries. If an error occur, an IsisConfException is raised
        @param section Name of the section in the configuration file where the option is (string)
        @param name Name of the option (string)
        @param allowed List of allowed entries (list of string)
        @return A valid string
        """
        try:
            # Parse file with RawConfigParser
            val = self.config.get (section, name)
            # Verify if the value is valid
            if val not in allowed:
                raise IsisConfException ("{} must be in this list: {}.".format(name, allowed))
            else:
                return {name : val}

        # Handle ConfigFileParser errors
        except ConfigParser.NoOptionError:
            raise IsisConfException ("Option {} was not found in conf file.".format(name))
        except ConfigParser.NoSectionError:
            raise IsisConfException ("Section {} was not found in conf file.".format(section))
        except ValueError:
            raise IsisConfException ("{} value is not a valid float".format(name))
