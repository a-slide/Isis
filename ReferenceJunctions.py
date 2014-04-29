# Standard library packages
from random import sample, randint
import csv

# Local packages
from SlicePicker import SlicePickerSingle

####################################################################################################

class ReferenceJunctions(object):
    """ Create de dictionnary of sequence by merging from 2 references 
    and by creating a chimeric junction between both for n junctions
    """

###    FONDAMENTAL METHODS    ###

    def __init__(self, name, min_chimeric, size, njun, ref1, ref2, repeats, ambiguous):
        """Create a dictionnary of junctions by merging 2 SeqRecords
        """
        print "Initialisation of {}...".format(name)
         
        # Store object variables
        self.name = name
        self.min_chimeric = min_chimeric
        self.lenjun = 2 * size
        
        # Initialise junction dictionnary 
        self.d = self._create_junctions_dict(size, njun, ref1, ref2, repeats, ambiguous)
        
        # Initialize a counter for each reference that will be
        # incremented each time _random_slice choose this reference
        self.reset_samp_counter()
        
    def __repr__(self):
        """Long description string used by interpreter and repr
        """
        result = "{}\n".format(self.__str__())
        for entry in self.d.values():
            result += "{}\nLenght:{}\n".format(entry, len(entry))
        return result

    def __str__(self):
        """Short representation
        """
        return "{} : Instance of {}".format(self.name, self.__module__)

###    GETERS    ###

    def get(self, varkey):
        return self.d[varkey]

    def getDict(self):
        return self.d
        
    def getName(self):
        return self.name

###    PUBLIC METHODS    ####

    def get_slice (self, size):
        """Generate a slice overlapping a junction and return it
        """
        # Guard conditions
        if  size < 2 * self.min_chimeric:
            raise Exception ("The size of the slice is too short.\n\
            Cannot define a fragment with the require minimal number of bases\
            overlapping each reference sequence")
        if size > self.lenjun:
            raise Exception ("The size of the slice is longer than the reference")
        
        # Pick a random reference junction and return a slice of it
        refseq = sample(self.d, 1)[0]
        return self._random_slice (refseq, size)

    def reset_samp_counter(self):
        """ Initialize or reset counter for each reference that will be
        incremented each time _random_slice choose this reference
        """
        for name, record in self.d.items():
            record.annotations ["nb_samp"] = 0

    def write_samp_report (self):
        """ Create a simple csv report 
        """
        # Open a file for writting with python csv module
        with open(self.name+"_samp_report.csv", 'w') as csvfile:
            writer = csv.writer(csvfile, delimiter='\t', quoting=csv.QUOTE_NONE)
            writer.writerow([
                "junction_id",
                "ref1_source",
                "ref1_chr",
                "ref1_loc",
                "ref2_source",
                "ref2_chr",
                "ref2_loc",
                "nb_samp"])
            
            # Create a sorted list of refseq
            ref_list = self.d.keys()
            ref_list.sort()
            
            # Export each refseq characteristics in a file
            for ref in ref_list:
                writer.writerow([
                ref,
                self.d[ref].annotations["ref1"].getName(),
                self.d[ref].annotations["ref1_refseq"],
                self.d[ref].annotations["ref1_location"],
                self.d[ref].annotations["ref2"].getName(),
                self.d[ref].annotations["ref2_refseq"],
                self.d[ref].annotations["ref2_location"],
                self.d[ref].annotations["nb_samp"]])

###    PRIVATE METHODS    ###

    def _create_junctions_dict (self, size, njunctions, ref1, ref2, repeats, ambiguous):
        """Fill a dictionary with junctions between sequences from ref1
        and ref2
        """
        
        print("\tCreating junctions between {} and {}".format(ref1.getName(), ref2.getName()))
        
        slicer = SlicePickerSingle(size, repeats, ambiguous)
        junctions_dict = {}

        for i in range(njunctions):
            
            # Pick 1 slice in each reference
            s1 = slicer.pick_slice(ref1)
            s2 = slicer.pick_slice(ref2)
            
            # Create a new junction from the 2 slice
            junction = s1 + s2
            
            # Define id, name and description to the SeqReccord object
            junction.name = junction.description = ""
            junction.id = "Junction_{:06}".format(i)
            
            # Add informations to the annotations dictionnary
            junction.annotations = {
                "ref1" : ref1,
                "ref2" : ref2,
                "ref1_refseq" : s1.annotations["refseq"],
                "ref2_refseq" : s2.annotations["refseq"],
                "ref1_location" : s1.annotations["location"],
                "ref2_location" : s2.annotations["location"]}
            
            junctions_dict[junction.id] = junction

        return junctions_dict

    def _random_slice (self, refseq, size):
        """Return a slice overlapping a junction from a biopython Seqrecord in d
        """
        # Randomly choose the slice start position in the autorized area
        start = randint((self.lenjun/2 + self.min_chimeric - size), (self.lenjun/2 - self.min_chimeric))
        end = start + size

        # Randomly choose an orientation reverse or forward and sample 
        # a slice
        if randint(0,1):
            s = self.d[refseq][start:end]
            s.annotations["location"] = [start, end]
        else:
            s = self.d[refseq][start:end].reverse_complement()
            s.annotations["location"] = [end, start]
        
        # Add informations to the annotations dictionnary
        s.annotations["refseq"] = refseq
        
        # Undefine description and define id and name  
        s.name = s.description = ""
        s.id = "{}|{}:{}-{}".format(
            self.name,
            s.annotations["refseq"],
            s.annotations["location"][0],
            s.annotations["location"][1])
        
        # Increment the sampling counter of the refseq
        self.d[refseq].annotations["nb_samp"] += 1

        return s

