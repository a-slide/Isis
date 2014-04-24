# Standard library packages
from random import random, sample, betavariate

# Third party packages
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA as Ambiguous
from Bio.Alphabet.IUPAC import IUPACUnambiguousDNA as Unambiguous

####################################################################################################
####################################################################################################

class SlicePicker(object):
    """Base class for SlicePickerSingle and SlicePickerPair 
    """

########################################################################
#   FONDAMENTAL METHODS
########################################################################
    
    def __init__(self, repeats, ambiguous, mut_freq):
        """Initialize the class with generic variables
        """
        # Definition of a generic mutation frequency
        self.mut_freq = mut_freq
        # Definition of a generic alphabet of allowed DNA bases
        self.alphabet = self._IUPAC_alphabet(repeats, ambiguous)

    def __repr__(self):
        """Long representation"""
        return "{}\n Alphabet : {}\n Mutation frequency : {}\n". format(
            self.__str__(),
            self.alphabet,
            self.mut_freq)

    def __str__(self):
        """Short representation"""
        return "<Instance of " + self.__module__ + ">\n"

########################################################################
#   PRIVATE METHODS
########################################################################

    def _IUPAC_alphabet(self, repeats, ambiguous):
        """ Define an DNA letters alphabet according to the user 
        specifications (allow repeats, allow ambiguous)
        """
        alphabet = Ambiguous.letters if ambiguous else Unambiguous.letters
        alphabet = alphabet + alphabet.lower() if repeats else alphabet
        return alphabet

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

    def _valid_sequence (self, seq):
        """ Verify if the users sequence 
        """
        return all ([(base in self.alphabet) for base in seq])

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

    def _mutate_sequence(self, read):
        """Introduce mutations in a biopython record"""
        
        read.annotations ["Mutations"] = []
        
        if self.mut_freq != 0:
            
            # Change seqreccord to mutable object
            seq = read.seq.tomutable()

            for i in range(len(seq)):
                if random() < self.mut_freq:
                    original_base = seq[i]
                    seq[i] = self._mutate_base(seq[i])
                    read.annotations ["Mutations"].append ("Pos {} {} -> {}".format(i+1, original_base, seq[i]))
            
            # Change back to non-mutable object
            read.seq = seq.toseq()
        
        # if no mutation was introduced
        if not read.annotations ["Mutations"]:
            read.annotations ["Mutations"] = "No mutation"

        return read

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

    def _mutate_base(self, base):
        """return a DNA base different from the given base"""
        return sample([mut_base for mut_base in Unambiguous.letters if mut_base not in base.upper()],1)[0]

####################################################################################################
####################################################################################################

class SlicePickerSingle(SlicePicker):
    """For single read sampling"""
    
########################################################################
#   FONDAMENTAL METHODS
########################################################################

    def __init__(self, read_len, repeats, ambiguous, mut_freq):
        """Initialize the class with super class init method and the 
        read len.  
        """
        # Use the super class init method
        super(self.__class__, self).__init__(repeats, ambiguous, mut_freq)
        
        # Store read length
        self.read_len = read_len

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

    def __repr__(self):
        """Long representation
        """
        return "{} Read Lenght {}\n". format(
            super(self.__class__, self).__repr__(),
            self.read_len)

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

    def __str__(self):
        """Short representation
        """
        return "<Instance of " + self.__module__ + ">\n"

########################################################################
#   ACTION METHODS
########################################################################

    def pick_slice (self, source):
        """Generate a candidate read. The presence or absence of repeats 
        and ambiguous DNA bases can be checked and a given frequency of 
        bases can be randomly mutated
        """
        # Guard condition if not possible to find a valid pair after 100 tries
        for count in range (100):

            # Ask a random sequence to the source
            try:
                read = source.get_slice(self.read_len)
            except Exception:
                continue

            # Verify the validity of the candidate sequence in terms of repeats and ambiguity
            if self._valid_sequence(str(read.seq)):
                return self._mutate_sequence(read)

        # If not candidate sequence if found after all tries a generic Exception is raised
        raise Exception ("No valid slice was found")


####################################################################################################
####################################################################################################


class SlicePickerPair(SlicePicker):
    """For paired read sampling"""
    
########################################################################
#   FONDAMENTAL METHODS
########################################################################

    def __init__(self, read_len, sonic_min, sonic_max, sonic_mode, sonic_certainty, repeats,
                 ambiguous, mut_freq):
        """Initialize the class with super class init method, the 
        read len and sonication parameters
        """
        # Use the super class init method
        super(self.__class__, self).__init__(repeats,  ambiguous, mut_freq)
        
        if not sonic_min <= sonic_mode <= sonic_max:
            raise Exception ("Wrong sonication parameters")
        
        # Store read length and sonication parameters
        self.read_len = read_len
        self.sonic_min = sonic_min
        self.sonic_max = sonic_max
        self.sonic_mode = sonic_mode
        self.sonic_certainty = sonic_certainty
        self.alpha, self.beta = self._beta_shape()
        
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

    def __repr__(self):
        """Long representation
        """
        
        return "{} Read Lenght {}\n Sonication minimum {}\n Sonication mode {}\n Sonication max {}\n Sonication certainty {}\n alpha {}\n beta {}\n". format(
                super(self.__class__, self).__repr__(), self.read_len, self.sonic_min, self.sonic_max, self.sonic_mode, self.sonic_certainty, self.alpha, self.beta)

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

    def __str__(self):
        """Short representation
        """
        return "<Instance of " + self.__module__ + ">\n"

########################################################################
#   ACTION METHODS
########################################################################

    def pick_slice (self, source):
        """Generate a candidate read pair from a fragment of a length
        following a distribution mimicking sonication smire. The
        presence or absence of repeats and ambiguous DNA bases can be
        checked and a given frequency of bases can be randomly mutated
         """

        # Guard condition if not possible to find a valid pair after 100 tries
        for count in range (100):
            # Generate a random size following a beta distribution (mimick sonication smire)
            frag_len = self._beta_distrib()

            # Ask a random sequence to the source
            try:
                fragment = source.get_slice(frag_len)
            except Exception:
                continue

            # Extract pair reads from slice
            read1, read2 = self._extract_pair (fragment, self.read_len, frag_len)

            # Verify the validity of the candidate sequence in terms of repeats and ambiguity
            if self._valid_sequence(str(read1.seq)) and self._valid_sequence(str(read2.seq)):
                return (self._mutate_sequence(read1), self._mutate_sequence(read2))

        # If not candidate sequence if found after all tries a generic Exception is raised
        raise Exception ("No valid slice was found")

########################################################################
#   PRIVATE METHODS
########################################################################

    def _beta_shape(self):
        """Calculate shape parameters alpha and beta to fit experimental indication from user"""
        
        alpha = float(self.sonic_mode-self.sonic_min) / (self.sonic_max-self.sonic_min) * (self.sonic_certainty-2) + 1
        beta = self.sonic_certainty - alpha
        return (alpha, beta)

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

    def _beta_distrib(self):
        """Define a pseudorandom size according to a beta distribution giving alpha and beta"""
        return int (betavariate(self.alpha, self.beta) * (self.sonic_max - self.sonic_min) + self.sonic_min)

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

    def _extract_pair(self, fragment, read_len, frag_len):
        """Extract reads forward and reverse from a fragment sequence and return a list with name, and both records"""

        forward = fragment[:read_len]
        reverse = fragment[-read_len:].reverse_complement()

        forward.annotations = {
        "orientation"   : fragment.annotations["orientation"],
        "source"        : fragment.annotations["source"],
        "refseq"        : fragment.annotations["refseq" ],
        "location"      : fragment.annotations["location"],
        "frag_length"   : frag_len,
        "pair_overlap"  : self._pair_overlap(read_len, frag_len),
        "read"          : "R1"}
        
        reverse.annotations = {
        "orientation"   : fragment.annotations["orientation"],
        "source"        : fragment.annotations["source"],
        "refseq"        : fragment.annotations["refseq" ],
        "location"      : fragment.annotations["location"],
        "frag_length"   : frag_len,
        "pair_overlap"  : self._pair_overlap(read_len, frag_len),
        "read"          : "R2"}
        
        forward.id = forward.name = forward.description = reverse.id = reverse.name = reverse.description=  ""

        return (forward, reverse)

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

    def _pair_overlap(self,read_len, frag_len):
        return "No_overlap" if frag_len > 2*read_len else "{0}_bp_overlap".format(2*read_len-frag_len)
