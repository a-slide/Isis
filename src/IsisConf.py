"""
@package    IsisConf
@brief
@copyright  [GNU General Public License v2](http://www.gnu.org/licenses/gpl-2.0.html)
@author     Adrien Leger <adrien.leger@gmail.com>
"""

#~~~~~~~PACKAGE IMPORTS~~~~~~~#

# Standard library packages
from sys import exit as sys_exit
import ConfigParser
import optparse

####################################################################################################

class IsisConfException(Exception):
    """Custom Exception Class to handle error while parsing Isis options
    """
####################################################################################################

    def __init__(self, msg): # Object constructor initialized with a custom user message
        err_msg = "IsisConfException : An error occured while parsing Isis options!\n"
        err_msg += "\t{}.\n\tPlease ajust your settings\n".format(msg)
        print (err_msg)
        sys_exit ()

####################################################################################################

class IsisConf(object):
    """Import configuration parameters from command lines arguments and
    verify parameters value range and validity. Users are invited to
    corect their values if a parameters is invalid.
    """
####################################################################################################

#####    FONDAMENTAL METHODS    #####

    def __init__(self, program_name, program_version):
        """ Import parameters and verify their values.
        """

        print "Parsing and verification of config file with IsisConf"
        # dictionnary to store validated parameters from command line
        # arguments and conf_file parsing
        self.d = {}

        ## Parse command line argument ####
        # add hg_filename, vg_filename, conf_filename and output_prefix
        # entries to self.d
        self.d.update (self._optparser(program_name, program_version))

        ## Extract configuration parameters from the conf file
        # creating a instance of RawConfigParser ####
        self.config = ConfigParser.RawConfigParser(allow_no_value = False)
        self.config.read(self.d["conf_file"])

        ## GENERAL SECTION ##
        self.d.update (self._get_int("General", "read_num", 1, None))
        self.d.update (self._get_int("General", "read_len", 1, None))
        self.d.update (self._get_float("General", "mut_freq", 0, 1))
        self.d.update (self._get_bool("General", "repeats"))
        self.d.update (self._get_bool("General", "ambiguous"))
        self.d.update (self._get_bool("General", "graph"))
        self.d.update (self._get_bool("General", "report"))

        ## FREQUENCY SECTION ##
        freq_host = self._get_float ("Frequency", "freq_host", 0, 1).values()[0]
        freq_virus = self._get_float ("Frequency", "freq_virus", 0, 1).values()[0]
        freq_tjun = self._get_float ("Frequency", "freq_tjun", 0, 1).values()[0]
        freq_fjun = self._get_float ("Frequency", "freq_fjun", 0, 1).values()[0]
        self._verify_freq_sum ([freq_host, freq_virus, freq_tjun, freq_fjun])

        # Calculate the number of reads from validated frequencies
        self.d["nread_host"] = int (freq_host * self.d["read_num"])
        self.d["nread_virus"] = int (freq_virus * self.d["read_num"])
        self.d["nread_tjun"] = int (freq_tjun * self.d["read_num"])
        self.d["nread_fjun"] = int (freq_fjun * self.d["read_num"])

        ## JUNCTION SECTION ##
        if self.d["pair"]:
            self.d.update (self._get_int ("Junction", "min_chimeric", 0, self.d["read_len"]))
        else:
            self.d.update (self._get_int ("Junction", "min_chimeric", 0, self.d["read_len"]/2))

        samp_tjun = self._get_float ("Junction", "samp_tjun", 0, self.d["nread_tjun"]).values()[0]
        samp_fjun = self._get_float ("Junction", "samp_fjun", 0, self.d["nread_fjun"]).values()[0]

        # Calculate the number of uniq source junctions to be generated
        self.d["uniq_tjun"] = int (float(self.d["nread_tjun"])/samp_tjun)
        self.d["uniq_fjun"] = int (float(self.d["nread_fjun"])/samp_fjun)

        ## SONICATION SECTION ##
        min_size = self.d["read_len"] + self.d["min_chimeric"]
        self.d.update (self._get_int ("Sonication", "sonic_min", min_size, None))
        self.d.update (self._get_int ("Sonication","sonic_mode", self.d["sonic_min"], None))
        self.d.update (self._get_int ("Sonication","sonic_max", self.d["sonic_mode"], None))
        self.d.update (self._get_int ("Sonication","sonic_certainty", 5, 50))

        ## QUALITY SECTION ##
        self.d.update (self._get_str ( "Quality", "qual_scale",
            ["fastq-sanger", "fastq-solexa", "fastq-illumina"]))
        self.d.update (self._get_str ( "Quality", "qual_range",
            ["very-good", "good", "medium", "bad", "very-bad"]))

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

#####    GETERS    #####

    def getDict (self):
        return self.d

    def get (self, key):
        return self.d[key]

#####    PRIVATE METHODS    #####

    def _optparser(self, program_name, program_version):
        """Parse command line arguments prompt and verify the filename
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
        arg_dict = {'host_genome' : self._check_file (options.hg, "--host_genome | -H"),
                    'virus_genome' : self._check_file (options.vg, "--virus_genome | -V "),
                    'conf_file' : self._check_file (options.conf, "conf_file | -C"),
                    'basename' : options.output,
                    'pair' : self._check_mode (options.single, options.pair)}

        return arg_dict

    def _check_file (self, path, descr):
        """Try to path from the opt dictionnary"""

        if not path:
            raise IsisConfException ("{} is a mandatory command line argument".format(descr))

        try:
            handle = open(path, "r")
            handle.close()
            return path

        except IOError:
            raise IsisConfException ("Error : " + path + " can't be read\n\
            Please enter a valid path")

    def _check_mode (self, single, pair):
        """Try to path from the opt dictionnary"""

        if pair and single:
            raise IsisConfException ("-p and -s are incompatible options")

        return False if single else True


    def _get_int(self, section, name, min=None, max=None):
        """Import an integer from self.config dict and verify its value
        If an error occur, an IsisConfException is raised
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
        """Import an float from self.config dict and verify its value
        If an error occur, an IsisConfException is raised
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
        """Import a boolean from self.config dict and verify its value
        If an error occur, an IsisConfException is raised
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
        """Import an float from self.config dict and verify its value
        If an error occur, an IsisConfException is raised
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


    def _verify_freq_sum(self, val_list):
        """ Verify if the sum of frequencies is equal to one
        """
        if sum(val_list) != 1:
            raise IsisConfException ("The sum of read number frequencies is not equal to 1")
