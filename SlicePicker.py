from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA as Ambiguous
from Bio.Alphabet.IUPAC import IUPACUnambiguousDNA as Unambiguous

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

    def peak_single (self, source, size, repeats, ambiguous, mut_freq = 0):
        
        # Guard condition if not possible to find a valid pair after 100 tries
        for count in range (100):
            try:
                # Pick a random sequence in dictionnary 
                slice = source.get_slice(size)
                # Verify the validity of the candidate sequence in terms of repeats and ambiguity
                if self._valid_sequence(str(slice.seq), repeats, ambiguous):
                    return slice
            
            except Exception:
                continue
        
        # If not candidate sequence if found after all tries a NoCandidateFound error is raised
        raise Exception ("No valid slice was found")
            
            
    def peak_pair (self, read_len, repeats, ambiguous, mut_freq, sonic_min, sonic_max, sonic_mode, sonic_certainty):
        pass

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

