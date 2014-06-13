"""
@package    SlicePicker
@brief
@copyright  [GNU General Public License v2](http://www.gnu.org/licenses/gpl-2.0.html)
@author     Adrien Leger - 2014
* <adrien.leger@gmail.com>
* <adrien.leger@inserm.fr>
* <adrien.leger@univ-nantes.fr>
* [Github](https://github.com/a-slide)
* [Atlantic Gene Therapies - INSERM 1089] (http://www.atlantic-gene-therapies.fr/)
"""

#~~~~~~~PACKAGE IMPORTS~~~~~~~#

# Standard library packages
from random import random, sample, betavariate

# Third party packages
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA as Ambiguous
from Bio.Alphabet.IUPAC import IUPACUnambiguousDNA as Unambiguous

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

class SlicePicker(object):
    """Base class for SlicePickerSingle and SlicePickerPair
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~FONDAMENTAL METHODS~~~~~~~#

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
        descr = "{}\n".format(self.__str__())
        descr += "Read Lenght {}\n".format(self.read_len)
        descr += "Alphabet {}\n".format(self.alphabet)
        descr += "Mutation frequency {}\n".format(self.mut_freq)
        return descr

    def __str__(self):
        return "<Instance of {} from {} >".format(self.__class__.__name__, self.__module__)

    #~~~~~~~ACCESS METHODS~~~~~~~#

    def get_mut_freq (self):
        return self.mut_freq

    def get_alphabet (self):
        return self.alphabet

    def get_read_len (self):
        return self.read_len

    #~~~~~~~PRIVATE METHODS~~~~~~~#

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

            # Change Seq to mutable object
            seq = read.seq.tomutable()

            for i in range(len(seq)):
                if random() < self.mut_freq:
                    original_base = seq[i]
                    seq[i] = self._mutate_base(seq[i])
                    read.annotations ["Mutations"].append("Pos {} {} -> {}".format(i+1, original_base, seq[i]))

            # Change back to non-mutable object
            read.seq = seq.toseq()

        return read

    def _mutate_base(self, base):
        """return a DNA base different from the given base"""
        return sample([mut_base for mut_base in Unambiguous.letters if mut_base not in base.upper()],1)[0]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

class SlicePickerSingle(SlicePicker):
    """For single read sampling
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~PUBLIC METHODS~~~~~~~#

    def pick_slice(self, source):
        """Generate a candidate read. The presence or absence of repeats
        and ambiguous DNA bases can be checked and a given frequency of
        bases can be randomly mutated
        """
        # Guard condition if not possible to find a valid pair after 100 tries
        for count in range(100):
            # Ask a random sequence to the source
            try:
                read = source.slicer(self.read_len)
            except Exception as e:
                print e
                exit (0)

            # Verify the validity of the candidate sequence
            if self._valid_sequence(str(read.seq)):
                read.annotations["frag_len"] = len(read.seq)
                return self._mutate_sequence(read)

        # If no candidate sequence was found an Exception is raised
        raise Exception("ERROR. Unable to find a valid slice after 100 tries.\n\
        Please review repetition and ambiguity parameters and verify your reference sequences")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

class SlicePickerPair(SlicePicker):
    """For paired read sampling
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~FONDAMENTAL METHODS~~~~~~~#

    def __init__(self, read_len, sonic_min, sonic_mode, sonic_max, sonic_certainty, repeats,
                 ambiguous, mut_freq):
        """Surdefine the super class init by adding sonication
        parameters
        """
        # Use the super class init method
        super(self.__class__, self).__init__(read_len, repeats, ambiguous, mut_freq)

        if not read_len <= sonic_min <= sonic_mode <= sonic_max:
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

    #~~~~~~~ACCESS METHODS~~~~~~~#

    def get_sonic_min (self):
        return self.sonic_min

    def get_sonic_mode (self):
        return self.sonic_mode

    def get_sonic_max (self):
        return self.sonic_max

    def get_sonic_certainty (self):
        return self.sonic_certainty

    def get_alpha (self):
        return self.alpha

    def get_beta (self):
        return self.beta

    #~~~~~~~PUBLIC METHODS~~~~~~~#

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
                fragment = source.slicer(frag_len)
            except Exception as e:
                print e
                exit (0)

            # Extract pair reads from slice
            read1, read2 = self._extract_pair(fragment, source)

            # Verify the validity of the candidate sequence
            if self._valid_sequence(str(read1.seq)) and self._valid_sequence(str(read2.seq)):
                return(self._mutate_sequence(read1), self._mutate_sequence(read2))

        # If no candidate sequence was found an Exception is raised
        raise Exception("ERROR. Unable to find a valid slice after 100 tries.\n\
        Please review repetition and ambiguity parameters and verify your reference sequences")

    #~~~~~~~PRIVATE METHODS~~~~~~~#

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

    def _extract_pair(self, fragment, source):
        """Extract reads forward and reverse from a fragment sequence
        and return a list with name, and both records
        """
        # Extract forward and reverse reads
        forward = fragment[:self.read_len]
        reverse = fragment[-self.read_len:].reverse_complement()

        # Add information to SeqReccord
        forward.annotations = self._annotate_read (fragment, "R1")
        reverse.annotations = self._annotate_read (fragment, "R2")

        # Postfix names of reads with R1 or R2 + reset name and dscr
        forward.id = reverse.id = fragment.id
        forward.name = reverse.name = forward.description = reverse.description = ""

        return(forward, reverse)

    def _annotate_read (self, fragment, mate):
        """Extract annotations for each of the reads from a pair
        """
        # Define shortcuts
        f_start = fragment.annotations["location"][0]
        f_end = fragment.annotations["location"][1]
        f_len = len(fragment.seq)
        r_len = self.read_len

        # Empty dict
        d = {}

        # Specific filling according the mate and fragment orientation
        if mate == "R1":
            d["mate"] = "R1"
            if fragment.annotations["orientation"] == "+":
                d["location"] = [f_start, f_start + r_len]
                d["orientation"] = "+"
            else:
                d["location"] = [f_end - r_len, f_end]
                d["orientation"] = "-"

        elif mate == "R2":
            d["mate"] = "R2"
            if fragment.annotations["orientation"] == "+":
                d["location"] = [f_end - r_len, f_end]
                d["orientation"] = "-"
            else:
                d["location"] = [f_start, f_start + r_len]
                d["orientation"] = "+"

        # Generic filling
        d["refseq"] = fragment.annotations["refseq"]
        d["frag_len"] = f_len
        d["source"] = fragment.annotations["source"]
        d["pair_overlap"] = 0 if f_len > 2*r_len else 2*r_len - f_len

        return d

