from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA as Ambiguous
from Bio.Alphabet.IUPAC import IUPACUnambiguousDNA as Unambiguous
from random import random
from random import sample
from random import betavariate

class SlicePicker:
    """"""

########################################################################################################################
#   FONDAMENTAL METHODS
########################################################################################################################

    def __repr__(self):
        """Long representation"""
        return self.__str__()

    def __str__(self):
        """Short representation"""
        return "<Instance of " + self.__module__ + ">\n"

########################################################################################################################
#   ACTION METHODS
########################################################################################################################

    def pick_single (self, source, read_len, repeats, ambiguous, mut_freq = 0):

        # Guard condition if not possible to find a valid pair after 100 tries
        for count in range (100):

            # Ask a random sequence to the source
            try:
                read = source.get_slice(read_len)
            except Exception:
                continue

            # Verify the validity of the candidate sequence in terms of repeats and ambiguity
            if self._valid_sequence(str(read.seq), repeats, ambiguous):
                return self._mutate_sequence(read, mut_freq)

        # If not candidate sequence if found after all tries a generic Exception is raised
        raise Exception ("No valid slice was found")

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

    def pick_pair (self, source, read_len, sonic_min, sonic_max, sonic_mode, sonic_certainty, repeats, ambiguous, mut_freq = 0):
        """Generate a candidate read or read pair of a given lenght with or without repeats and ambiguous DNA bases"""

        # Calculate the parameters of shape for the beta distribution to mimick DNA shearing distribution by sonication
        alpha, beta = self._beta_shape(sonic_min, sonic_max, sonic_mode, sonic_certainty)
        ######################## TEST SONICATION VALUES

        # Guard condition if not possible to find a valid pair after 100 tries
        for count in range (100):
            # Generate a random size following a beta distribution (mimick sonication smire)
            frag_len = self._beta_distrib(alpha, beta, sonic_min, sonic_max)

            # Ask a random sequence to the source
            try:
                fragment = source.get_slice(frag_len)
            except Exception:
                continue

            # Extract pair reads from slice
            read1, read2 = self._extract_pair (fragment, read_len, frag_len)

            # Verify the validity of the candidate sequence in terms of repeats and ambiguity
            if self._valid_sequence(str(read1.seq), repeats, ambiguous) and self._valid_sequence(str(read2.seq), repeats, ambiguous):
                return (self._mutate_sequence(read1, mut_freq), self._mutate_sequence(read2, mut_freq))

        # If not candidate sequence if found after all tries a generic Exception is raised
        raise Exception ("No valid slice was found")

########################################################################################################################
#   PRIVATE METHODS
########################################################################################################################

    def _valid_sequence (self, seq, repeats, ambiguous):
        """Define if the candidate region is valid according to user specifications (allow repeats, allow ambiguous)"""
        valid_base = self._IUPAC(repeats, ambiguous)
        return all ([(base in valid_base) for base in seq])

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

    def _IUPAC(self, repeats, ambiguous):
        """Return the IUPAC DNA alphabet corresponding to users requirements"""
        return self._repeats (Ambiguous.letters, repeats) if ambiguous else self._repeats (Unambiguous.letters, repeats)

    def _repeats(self, alphabet, repeats):
        return alphabet + alphabet.lower() if repeats else alphabet

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

    def _beta_shape(self, min, max, mode, certainty):
        """Calculate shape parameters alpha and beta to fit experimental indication from user"""

        alpha = float((mode - min)) / (max - min) * (certainty-2) + 1
        beta = certainty - alpha
        return (alpha, beta)

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

    def _beta_distrib(self, alpha, beta, min, max):
        """Define a pseudorandom size according to a beta distribution giving alpha and beta"""
        return int (betavariate(alpha, beta) * (max - min) + min)

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

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

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

    def _pair_overlap(self,read_len, frag_len):
        return "No_overlap" if frag_len > 2*read_len else "{0}_bp_overlap".format(2*read_len-frag_len)

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

    def _mutate_sequence(self, read, mut_freq):
        """Introduce mutations in a biopython record"""
        
        read.annotations ["Mutations"] = []
        
        if mut_freq != 0:
            
            # Change seqreccord to mutable object
            seq = read.seq.tomutable()

            for i in range(len(seq)):
                if random() < mut_freq:
                    original_base = seq[i]
                    seq[i] = self._mutate_base(seq[i])
                    read.annotations ["Mutations"].append ("Pos {} {} -> {}".format(i+1, original_base, seq[i]))
            
            # Change back to non-mutable object
            read.seq = seq.toseq()
        
        # if no mutation was introduced
        if not read.annotations ["Mutations"]:
            read.annotations ["Mutations"] = "No mutation"

        return read

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
    def _mutate_base(self, base):
        """return a DNA base different from the given base"""
        return sample([mut_base for mut_base in Unambiguous.letters if mut_base not in base.upper()], 1)[0]

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
