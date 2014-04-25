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

    def __init__(self, source, min_chimeric, size, njun, ref1, ref2, repeats, ambiguous):
        """Create a dictionnary of junctions by merging 2 SeqRecords
        """
        # Store global variable
        self.source = source
        self.min_chimeric = min_chimeric
        # initialise junction dictionnary 
        self.d = self._create_junctions_dict(size, njun, ref1, ref2, repeats, ambiguous)

    def __repr__(self):
        """Long description string used by interpreter and repr
        """
        result = "{}<Source = {}>\n\n".format(self.__str__(), self.source)
        for entry in self.d.values():
            result += "<{}\nLenght:{}>\n".format(entry, len(entry))
        return result

    def __str__(self):
        """Short representation
        """
        return "<Instance of " + self.__module__ + ">\n"

###    GETERS    ###

    def getSource(self):
        """Give acces to the name of the source
        """
        return self.source

    def get(self, varkey):
        """ Give acces to individual values in d by using its key name
        """
        return self.d[varkey]

    def getDict(self):
        """ Give acces to the complete dictionary
        """
        return self.d

###    PUBLIC METHODS    ####

    def get_slice (self, size):
        """Generate a slice overlapping a junction and return it
        """
        # Guard condition
        for count in range (100):
            # Pick a random reference junction
            refseq = sample(self.d, 1)[0]
            
            # If the size is valid return a slice
            if 2*self.min_chimeric <= size <= len(self.d[refseq]):
                return self._random_slice (refseq, size, lenjun)
        
        # if no valid size was found
        raise Exception ("No valid slice was found")

###    PRIVATE METHODS    ###

    def _create_junctions_dict (self, size, njunctions, ref1, ref2, repeats, ambiguous):
        """Fill a dictionary with junctions between sequences from ref1
        and ref2
        """
        slicer = SlicePickerSingle(size, repeats, ambiguous)
        junctions_dict = {}

        for i in range(njunctions):
            
            # Pick 1 slice in each reference
            seq1 = slicer.pick_single(ref1)
            seq2 = slicer.pick_single(ref2)
            
            # Create a new junction from the 2 slice
            junction = seq1 + seq2
            
            # Define id, name and description to the SeqReccord object
            junction.id = "#{:010}".format(i)
            junction.name = junction.description = None
            
            # Add informations to the annotations dictionnary
            junction.annotations = {
            "seq1_source"       : seq1.annotations["source"],
            "seq2_source"       : seq2.annotations["source"],
            "seq1_refseq"       : seq1.annotations["refseq"],
            "seq2_refseq"       : seq2.annotations["refseq"],
            "seq1_location"     : seq1.annotations["location"],
            "seq2_location"     : seq2.annotations["location"],
            "seq1_orientation"  : seq1.annotations["orientation"],
            "seq2_orientation"  : seq2.annotations["orientation"]}

            junctions_dict[junction.name] = junction

        return junctions_dict

    def _random_slice (self, refseq, size, lenjun):
        """Return a slice overlapping a junction from a biopython Seqrecord in d
        """
        # Example of the strategy with size = 20 and min chimeric = 6
        # Junction      -----------------------------|-----------------------------
        # Chimeric bases                       /////////////
        # Start area    ///////////////ooooooo////////////////////////////////////
        # End area      ////////////////////////////////////ooooooo////////////////

        # Randomly choose the slice start position in the autorized area
        start = randint((lenjun/2 + self.min_chimeric - size), (lenjun/2 - self.min_chimeric))
        end = start + size

        # Randomly choose an orientation reverse or forward for the fragment
        if randint(0,1):
            slice = self.d[refseq][start:end]
            slice.annotations["location"] = [start, end]
        else:
            slice = self.d[refseq][start:end].reverse_complement()
            slice.annotations["location"] = [end, start]
        
        # Add informations to the annotations dictionnary
        slice.annotations["source"] = self.source
        slice.annotations["refseq"] = refseq
        slice.annotations["size"] = size
        
        # Undefine Name and description and define id  
        slice.name = slice.description = None
        slice.id = "{}|{}:{}-{}|size:{}".format(
            self.d["source"],
            self.d["refseq"],
            self.d["location"][0],
            self.d["location"][1],
            self.d["size"])

        return slice

