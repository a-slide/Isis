
class FastqGenerator(object):
    """Class description"""

#####    FONDAMENTAL METHODS    #####
    
    def __repr__(self):
        """Long representation"""
        return self.__str__()

    def __str__(self):
        """Short representation"""
        return "<Instance of " + self.__module__ + ">\n"
    
#####    PUBLIC METHODS    #####

    def fastq_single (self, slice_picker, qual_generator, source):
        """"""
        # Ask a seqReccord slice to a reference sequence source
        read = s.pick_slice(source)

        # Add a phred quality string the seqReccord
        read.letter_annotations["phred_quality"] = qual_generator.random_qual_string()

        # Write a proper id to the read
            read.id = "{}|{}|loc_{}_{}_{}|".format(
                read.annotations["source"],
                read.annotations["refseq"],
                read.annotations["location"][0],
                read.annotations["location"][1],
                read.annotations["orientation"])
            
        return read

    def fastq_pair (self, slice_picker, qual_generator, source, read_len, sonic_min, sonic_max,
                    sonic_mode, sonic_certainty, repeats, ambiguous, mut_freq = 0):
        """"""
        # Ask a pair of seqReccord slice to a reference sequence source
        read = s.pick_pair(source, read_len, sonic_min, sonic_max, sonic_mode, sonic_certainty,
                           repeats, ambiguous, mut_freq)
        
        
#####    PRIVATE METHODS    #####

    def _do_something (self):
        """Do something inside the object"""
        pass
