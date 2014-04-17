from SlicePicker import SlicePicker
from random import sample
from random import randint

class ReferenceJunctions:
    """"""

########################################################################################################################
#   FONDAMENTAL METHODS
########################################################################################################################

    def __init__(self, source, max_chimeric, size, njunctions, ref1, ref2, repeats = False, ambiguous = False):
        """Create a dictionnary of junctions by merging 2 Biopython records"""
        self.source = source
        self.max_chimeric = max_chimeric
        self.d = self._create_junctions_dict (size, njunctions, ref1, ref2, repeats, ambiguous)

    def __repr__(self):
        """Long description string used by interpreter and repr"""
        result = self.__str__()
        for entry in self.d.values():
            result += "<{0}\nLenght:{1}>\n\n".format(entry, len(entry))
        return result

    def __str__(self):
        """Short representation"""
        return "<Instance of " + self.__module__ + ">\n"

########################################################################################################################
#   GETERS
########################################################################################################################

    # Give acces to the name of the source
    def getSource (self):
        return self.source

    # Give acces to individual values in d by using its key name.
    def get (self, varkey ):
        return self.d[varkey]

    # Give acces to the complete dictionary
    def getDict (self):
        return self.d

########################################################################################################################
#   PUBLIC METHODS
########################################################################################################################

    def get_slice (self, size):
        """Generate a slice overlapping a junction and return it"""

        # No need to use a try block. Junction are mandatory longer than size of region to be sample
        # Pick a random reference junction
        refseq = sample(self.d, 1)[0]

        # Determine the size of the junction
        lenjun = len(self.d[refseq])

        # Calculate the size of left and right boundaries
        offset = lenjun/2 - size + self.max_chimeric

        # Randomly choose the slice start position
        start = randint(offset,lenjun - offset - size)
        end = start+size

        # Randomly choose an orientation reverse or forward for the fragment
        forward = randint(0,1)

        slice = self.d[refseq][start:end] if forward else self.d[refseq][start:end].reverse_complement()

        # Adding informations to the biopython record
        slice.id = slice.name = slice.description = ""
        slice.annotations = {
        "orientation"   : "+" if forward else "-",
        "source"        : self.source,
        "refseq"        : refseq,
        "location"      : [start, end]}

        return slice

########################################################################################################################
#   PRIVATE METHODS
########################################################################################################################

    def _create_junctions_dict (self, size, njunctions, ref1, ref2, repeats, ambiguous):
        """Fill a dictionary with junctions between sequences from ref1 and ref2"""
        slicer = SlicePicker()
        junctions_dict = {}

        for i in range(njunctions):
            seq1 = slicer.peak_single(ref1, size, repeats, ambiguous)
            seq2 = slicer.peak_single(ref2, size, repeats, ambiguous)

            junction = seq1 + seq2
            junction.name = junction.id = "#{:010}".format(i)

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




