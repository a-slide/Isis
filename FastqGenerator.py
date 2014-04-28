# Standard library packages
import multiprocessing as mp

####################################################################################################

class FastqGenerator(object):
    """Base class for FastqGeneratorSingle and FastqGeneratorPair
    """
####################################################################################################

#####    FONDAMENTAL METHODS    #####
    
    def __init__(self, slicer, qualgen, qual_scale):
        """ Init object with link to SlicePicker and QualGenerator and
        determine the number of available thread for multiprocessing
        """
        
        # Store link to SlicePicker and QualGenerator object
        self.slicer = slicer
        self.qualgen = qualgen
        
        # Store the quality scale that will be used to generate fastq
        self.qual_scale = qual_scale
        
        # Automatically determine the number of available thread
        try:
            self.nb_thread = mp.cpu_count()
        except NotImplementedError:
            print "cpu_count method is not available on your system" 
            self.nb_thread = 1
            
    def __repr__(self):
        """Long representation
        """
        return "{}\n QualGenerator :\n{}\nSlicePicker :\n{}\n".format(
            self.__str__(),
            self.qualgen.repr(),
            self.slicer.repr())
        
    def __str__(self):
        """Short representation
        """
        return "<Instance of " + self.__module__ + ">"
        
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
    
    def write_fastq (self, source, nread, basename):
        """ 
        """
        filename = basename + ".fastq"
        
        try:
            with (open(filename, "w")) as f: # RENAME READS TO BE UNIQLY IDENTIFIED HERE ####################################
                for i in range (nread):
                    f.write(self.generate_fastq(source).format(self.qual_scale))
        
        except IOError:
            print ("Error : " + filename + " can't be write")
    
    def generate_fastq (self,source):
        """ Ask a read to slicer and add a quality string generated by
        qualgen 
        """
        # Ask a seqReccord slice to a reference sequence source
        read = self.slicer.pick_slice(source)

        # Add a phred quality string the seqReccord
        read.letter_annotations["phred_quality"] = self.qualgen.qual_score()

        return read
        
####################################################################################################

class FastqGeneratorPair(FastqGenerator):
    """ Generate paired end fastq reads 
    """
####################################################################################################
    
#####    PUBLIC METHODS    #####

    def write_fastq (self, source, nread, filename):
        """ 
        """
        filename1 = basename + "_R1.fastq"
        filename2 = basename + "_R2.fastq"
        
        try:
            with (open(filename1, "w")) as f1:
                with (open(filename2, "w")) as f2:
                    for i in range (nread):
                        r1, r2 = self.generate_fastq(source)
                        f1.write(r1.format(self.qual_scale))
                        f2.write(r2.format(self.qual_scale))
                        
        except IOError:
            print ("Error : " + filename1 + " can't be write")

    def generate_fastq (self,source):
        """ Ask a pair of reads to slicer and add a quality string 
        generated by qualgen to both of them. 
        """
        # Ask a seqReccord slice to a reference sequence source
        read1, read2 = self.slicer.pick_slice(source)

        # Add a phred quality string the seqReccord
        read1.letter_annotations["phred_quality"] = self.qualgen.qual_score()
        read2.letter_annotations["phred_quality"] = self.qualgen.qual_score()

        return read1, read2
