# Standard library packages
from random import random, sample, betavariate

# Third party packages
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA as Ambiguous
from Bio.Alphabet.IUPAC import IUPACUnambiguousDNA as Unambiguous

####################################################################################################

class SlicePicker(object):
    """Base class for SlicePickerSingle and SlicePickerPair 
    """
    
#####    FONDAMENTAL METHODS    #####
    
    def __init__(self, read_len, repeats, ambiguous, mut_freq=0):
        """Initialize the class with generic variables
        """
        # Definition of a generic mutation frequency
        self.mut_freq = mut_freq
        # Definition of a generic alphabet of allowed DNA bases
        self.alphabet = self._IUPAC_alphabet(repeats, ambiguous)
        # Store read length
        self.read_len = read_len

    def __repr__(self):
        """Long representation"""
        
        descr = self.__str__(),
        descr += "Read Lenght {}\n".format(self.read_len)
        descr += "Alphabet {}\n".format(self.alphabet)
        descr += "Mutation frequency {}\n".format(self.mut_freq)
        return descr

    def __str__(self):
        """Short representation"""
        return "<Instance of " + self.__module__ + ">\n"

#####    PRIVATE METHODS    #####

    def _IUPAC_alphabet(self, repeats, ambiguous):
        """ Define an DNA letters alphabet according to the user 
        specifications (allow repeats, allow ambiguous)
        """
        alphabet = Ambiguous.letters if ambiguous else Unambiguous.letters
        alphabet = alphabet + alphabet.lower() if repeats else alphabet
        return alphabet

    def _valid_sequence(self, seq):
        """ Verify if the users sequence 
        """
        return all([(base in self.alphabet) for base in seq])

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
                    read.annotations ["Mutations"].append("Pos {} {} -> {}".format(i+1, original_base, seq[i]))
            
            # Change back to non-mutable object
            read.seq = seq.toseq()
        
        # if no mutation was introduced
        if not read.annotations ["Mutations"]:
            read.annotations ["Mutations"] = "No mutation"

        return read

    def _mutate_base(self, base):
        """return a DNA base different from the given base"""
        return sample([mut_base for mut_base in Unambiguous.letters if mut_base not in base.upper()],1)[0]

####################################################################################################

class SlicePickerSingle(SlicePicker):
    """For single read sampling"""
    
#####    FONDAMENTAL METHODS    #####

    # User base class __init__, __repr__ and __str__ methods

#####    ACTION METHODS    #####

    def pick_slice(self, source):
        """Generate a candidate read. The presence or absence of repeats 
        and ambiguous DNA bases can be checked and a given frequency of 
        bases can be randomly mutated
        """
        # Guard condition
        for count in range(100):
            # Ask a random sequence to the source
            try:
                read = source.get_slice(self.read_len)
            except Exception:
                continue

            # Verify the validity of the candidate sequence
            if self._valid_sequence(str(read.seq)):
                return self._mutate_sequence(read)

        # If no candidate sequence was found an Exception is raised
        raise Exception("No valid slice was found")

####################################################################################################

class SlicePickerPair(SlicePicker):
    """For paired read sampling"""
    
#####    FONDAMENTAL METHODS    #####

    def __init__(self, read_len, sonic_min, sonic_max, sonic_mode, sonic_certainty, repeats,
                 ambiguous, mut_freq):
        """Surdefine the super class init by adding sonication
        parameters
        """
        # Use the super class init method
        super(self.__class__, self).__init__(read_len, repeats, ambiguous, mut_freq)
        
        if not sonic_min <= sonic_mode <= sonic_max:
            raise Exception("Wrong sonication parameters")
        
        # Store read length and sonication parameters
        self.sonic_min = sonic_min
        self.sonic_max = sonic_max
        self.sonic_mode = sonic_mode
        self.sonic_certainty = sonic_certainty
        self.alpha, self.beta = self._beta_shape()

    def __repr__(self):
        """Long representation
        """
        descr = super(self.__class__, self).__repr__()
        descr += "Sonication minimum {}\n".format(self.sonic_min)
        descr += "Sonication mode {}\n".format(self.sonic_mode)
        descr += "Sonication max {}\n".format(self.sonic_max)
        descr += "Sonication certainty {}\n".format(self.sonic_certainty)
        descr += "alpha {}\n".format(self.alpha)
        descr += "beta {}\n".format(self.beta)
        return descr

    def __str__(self):
        """Short representation
        """
        return "<Instance of " + self.__module__ + ">\n"

#####    ACTION METHODS    #####

    def pick_slice(self, source):
        """Generate a candidate read pair from a fragment of a length
        following a distribution mimicking sonication smire. The
        presence or absence of repeats and ambiguous DNA bases can be
        checked and a given frequency of bases can be randomly mutated
         """

        # Guard condition if not possible to find a valid pair after 100 tries
        for count in range(100):
            # Generate a random size following a beta distribution
            frag_len = self._beta_distrib()

            # Ask a random sequence to the source
            try:
                fragment = source.get_slice(frag_len)
            except Exception:
                continue

            # Extract pair reads from slice
            read1, read2 = self._extract_pair(fragment, self.read_len, frag_len)

            # Verify the validity of the candidate sequence
            if self._valid_sequence(str(read1.seq)) and self._valid_sequence(str(read2.seq)):
                return(self._mutate_sequence(read1), self._mutate_sequence(read2))

        # If no candidate sequence was found an Exception is raised
        raise Exception("No valid slice was found")

#####    PRIVATE METHODS    #####

    def _beta_shape(self):
        """Calculate shape parameters alpha and beta to fit experimental
        indication from user
        """
        alpha = float(self.sonic_mode-self.sonic_min) / (self.sonic_max-self.sonic_min) * (self.sonic_certainty-2) + 1
        beta = self.sonic_certainty - alpha
        return(alpha, beta)

    def _beta_distrib(self):
        """Define a pseudorandom size according to a beta distribution
        giving alpha and beta
        """
        return int(betavariate(self.alpha, self.beta) * (self.sonic_max - self.sonic_min) + self.sonic_min)

    def _extract_pair(self, fragment, read_len, frag_len):
        """Extract reads forward and reverse from a fragment sequence
        and return a list with name, and both records
        """
        # Extract forward read and add information to SeqReccord
        forward = fragment[:read_len]
        forward.annotations = fragment.annotations
        forward.annotations ["read"] = "R1"
        forward.annotations ["pair_overlap"] = self._pair_overlap(read_len, frag_len)
        forward.id = "{}|{}|{}".format(
            fragment.id
            forward.annotations ["pair_overlap"],
            forward.annotations ["read"])
        
        # Extract Reverse read and add information to SeqReccord
        reverse = fragment[-read_len:].reverse_complement()
        reverse.annotations = fragment.annotations
        reverse.annotations ["read"] = "R2"
        reverse.annotations ["pair_overlap"] = self._pair_overlap(read_len, frag_len)
        reverse.id = "{}|{}|{}".format(
            fragment.id
            reverse.annotations ["pair_overlap"],
            reverse.annotations ["read"])

        return(forward, reverse)

    def _pair_overlap(self,read_len, frag_len):
        """Simple fonction returning a string the overlap between the 2
        mates of a read
        """
        return "No_overlap" if frag_len > 2*read_len else "{}_bp_overlap".format(2*read_len - frag_len)
