from ConfFileParser import ConfFileParser

class IsisConf:
    """Class description"""

########################################################################################################################
#   FONDAMENTAL METHODS
########################################################################################################################

    def __init__(self):
        try:
            #### Parse command line argument ####

            self.host_genome_file, self.virus_genome_file, self.conf_file, self.output_prefix = _self.parse_arg()
            #Except if not openable

            #### Extract configuration parameters from the conf file by creating a instance of ConfFileParser ####
            self.conf = ConfFileParser(self.conf_file)
            # Except if non readable

            # Overall quantity of reads
            self.read_num       = self._import_numeric("read_num", 1, None, 100)
            # Except if value not suitable
            # Except if not found in dict
            
            # Read length
            self.read_len       = self._import_numeric("read_len", 1, None, 150)

            # Host genome frequency
            self.freq_hg        = self._import_numeric("freq_hg",0,1,0.45)
            # Virus genome frequency
            self.freq_vg        = self._import_numeric("freq_vg",0,1,0.45)
            # True junctions frequency
            self.freq_tj        = self._import_numeric("freq_tj",0,1,0.09)
            # False junctions frequency
            self.freq_fj        = self._import_numeric("freq_fj",0,1,0.01)
            
            self._verify_freq_sum()
            # Except if != 1

            # Mean number of sampling in true junctions and false junction (should be 1 for false junctions)
            self.samp_tj        = self._import_numeric("samp_tj", 0, self.freq_tj*self.read_num, None, 100) # Should be more than 1
            self.samp_fj        = self._import_numeric("samp_fj", 0, self.freq_fj*self.read_num, None, 100) # Should be 1 or 0 

            ## Allow or forbid sampling in repeat regions (lowercase characters), ambigous DNA (Ambigous IUPAC code) and duplicated reads.
            self.repeats        = self._import_boolean("repeats", False)
            self.ambigous       = self._import_boolean("ambigous", False)
            self.duplicate      = self._import_boolean("duplicate", False)
            # Except if != Boolean

            # Sequencing mode. If true = pair end alse = single end
            self.pair_end       = self._import_boolean("pair_end", True)

            ## Options of fragment sonication distribution to set up if paired end mode
            # Maximal Minimal and mode of sonication size max > mode > min
            self.sonic_min      = self._import_numeric (sonic_min, self.read_len, None, 200)
            self.sonic_mode     = self._import_numeric (sonic_mode, self.sonic_min, None, 400)
            self.sonic_max      = self._import_numeric (sonic_max, self.sonic_mode, None, 1000)

            # Certainty of the sonication smire: 5 (wide peak) to 50 (thin peak). Recommanded = 10
            sonic_certainty    = self._import_numeric(sonic_certainty, 5, 50, 10)

            # Quality score scale of fastq output file = solexa (CASAVA < 1.3), illumina (CASAVA 1.3 to 1.7), sanger (which is CASAVA >= 1.8)
            self.qual_scale     = self._import_string(qual_scale, ["solexa", "illumina", "sanger"])

            # Mean quality score of fastq output file at start, middle and end of reads between 0 and 40
            self.qual_start     = self._import_numeric (qual_start,0, 40, 30)
            self.qual_mid       = self._import_numeric (qual_mid,0, 40, 35)
            self.qual_end       = self._import_numeric (qual_end,0, 40, 30)

            # Range of quality variation at start, middle and end of reads between 0 and 20
            self.qvar_start     = self._import_numeric (qvar_start, 0, 20, 5)
            self.qvar_mid       = self._import_numeric (qvar_mid, 0, 20, 2)
            self.qvar_end       = self._import_numeric (qvar_end, 0, 20, 10)

            # Number de available threads for distributed computing
            self.nb_thread      = self._import_numeric (nb_thread, 1, None, 1)
        
        except

    def __repr__(self):
        """Long representation"""
        return

    def __str__(self):
        """Short representation"""
        return  "<Instance of IsisConf>\n"

########################################################################################################################
#   GETERS
########################################################################################################################

    def get_variable (self):
        """return variable"""
        return self.variable
    
    def get_qual_scale(self)
        # return a string of possible value according to 


########################################################################################################################
#   PRIVATE METHODS
########################################################################################################################

    def _do_something (self):
        """Do something inside the object"""
        pass
        
    def _self.parse_arg(self): 
        #return a tuple containing the arg
        pass
    
    def _import_numeric(name, min, max, default):
        
        try:
            var = self.conf.getVar(name)
        except NameError:
            print ("Error : " + name + " was not found in " + self.conf_file)
            print ("Would you like to use the default value (y/Y/n/N) ?\t<" + default + ">")
        
        
