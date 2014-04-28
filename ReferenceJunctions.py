# Standard library packages
from random import sample, randint

# Local packages
from SlicePicker import SlicePickerSingle

####################################################################################################

class ReferenceJunctions(object):
    """ Create de dictionnary of sequence by merging from 2 references 
    and by creating a chimeric junction between both for n junctions
    """

###    FONDAMENTAL METHODS    ###

    def __init__(self, min_chimeric, size, njun, ref1, ref2, repeats, ambiguous):
        """Create a dictionnary of junctions by merging 2 SeqRecords
        """
        # Store global variable
        self.min_chimeric = min_chimeric
        # Store the len of junctions
        self.lenjun = 2 * size
        
        # Initialise junction dictionnary 
        self.d = self._create_junctions_dict(size, njun, ref1, ref2, repeats, ambiguous)

    def __repr__(self):
        """Long description string used by interpreter and repr
        """
        result = "{}\n".format(self.__str__())
        for entry in self.d.values():
            result += "<{}\nLenght:{}>\n".format(entry, len(entry))
        return result

    def __str__(self):
        """Short representation
        """
        return "<Instance of " + self.__module__ + ">"

###    GETERS    ###

    def get(self, varkey):
        return self.d[varkey]

    def getDict(self):
        return self.d

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
    

###    PRIVATE METHODS    ###

    def _create_junctions_dict (self, size, njunctions, ref1, ref2, repeats, ambiguous):
        """Fill a dictionary with junctions between sequences from ref1
        and ref2
        """
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
            junction.id = "#{:010}".format(i)
            
            # Add informations to the annotations dictionnary
            junction.annotations = {
                "ref1" : ref1,
                "ref2" : ref2,
                "ref1_refseq" : s1.annotations["refseq"],
                "ref2_refseq" : s2.annotations["refseq"],
                "ref1_location" : s1.annotations["location"],
                "ref2_location" : s2.annotations["location"],
                "nb_samp" : 0,}
                # nb_samp will be incremented each time a junction is
                # sampled in the reference
            
            junctions_dict[junction.id] = junction

        return junctions_dict

    def _random_slice (self, refseq, size):
        """Return a slice overlapping a junction from a biopython Seqrecord in d
        """
        # Example of the strategy with size = 20 and min chimeric = 6
        # Junction      -----------------------------|-----------------------------
        # Chimeric bases                       /////////////
        # Start area    ///////////////ooooooo////////////////////////////////////
        # End area      ////////////////////////////////////ooooooo////////////////

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
        s.id = "{}:{}-{}".format(
            s.annotations["refseq"],
            s.annotations["location"][0],
            s.annotations["location"][1])
        
        # Increment the sampling counter of the refseq
        self.d[refseq].annotations["nb_samp"] += 1

        return s

