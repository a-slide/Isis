
####################################################################################################

class FastqGenerator(object):
    """Base class for FastqGeneratorSingle and FastqGeneratorPair
    """
####################################################################################################

#####    FONDAMENTAL METHODS    #####

    def __init__(self, slicer, qualgen):
        """ Init object with link to SlicePicker and QualGenerator and
        determine the number of available thread for multiprocessing
        """
        # Store link to SlicePicker and QualGenerator object
        self.slicer = slicer
        self.qualgen = qualgen

    def __repr__(self):
        return "{}\n QualGenerator :\n{}\nSlicePicker :\n{}\n".format(
            self.__str__(),
            self.qualgen.repr(),
            self.slicer.repr())

    def __str__(self):
        return "<Instance of {} from {} >".format(self.__class__.__name__, self.__module__)

#####    GETTERS    #####

    def get_slicer(self):
        return self.slicer

    def get_qualgen(self):
        return self.qualgen

####################################################################################################

class FastqGeneratorSingle(FastqGenerator):
    """ Generate single end fastq reads
    """
####################################################################################################

#####    PUBLIC METHODS    #####

    def generate_fastq (self,source):
        """ Ask a read to slicer and add a quality string generated by
        qualgen
        """
        # Ask a seqReccord slice to a reference sequence source
        try:
            read = self.slicer.pick_slice(source)
        except Exception as e:
            print e
            exit (0)

        # Add a phred quality string the seqReccord
        read.letter_annotations["phred_quality"] = self.qualgen.qual_score()

        return read

####################################################################################################

class FastqGeneratorPair(FastqGenerator):
    """ Generate paired end fastq reads
    """
####################################################################################################

#####    PUBLIC METHODS    #####

    def generate_fastq (self,source):
        """ Ask a pair of reads to slicer and add a quality string
        generated by qualgen to both of them.
        """
        # Ask a seqReccord slice to a reference sequence source
        try:
            read1, read2 = self.slicer.pick_slice(source)
        except Exception as e:
            print e
            exit (0)

        # Add a phred quality string the seqReccord
        read1.letter_annotations["phred_quality"] = self.qualgen.qual_score()
        read2.letter_annotations["phred_quality"] = self.qualgen.qual_score()

        return read1, read2

