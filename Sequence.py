from random import random
from random import randint
from random import betavariate
from random import sample

import matplotlib.pyplot as plt

from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA as Ambiguous
from Bio.Alphabet.IUPAC import IUPACUnambiguousDNA as Unambiguous
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

class Sequence:
    """"""

########################################################################################################################
#   FONDAMENTAL METHODS
########################################################################################################################

    def __init__ (self, nseq, min_len, max_len,repeats = False, ambiguous = False):
        """Create a simple dictionnary containing randomly generated sequences"""
        self.d = self._create_simple_dict(nseq, min_len, max_len, repeats, ambiguous)

    def __repr__(self):
        """Long description string used by interpreter and repr"""
        result = self.__str__()
        for entry in self.d.values():
            result += "<{0}\nLenght:{1}>\n\n".format(entry, len(entry))
        return result

    def __str__(self):
        """Short representation"""
        return "<Instance of " + self.__module__ + ">\n"

    # TODO Write this documentation
    def __doc__ (self):
        pass

########################################################################################################################
#   GETERS
########################################################################################################################

    # Give acces to the complete dictionary
    def getDict (self):
        return self.d

    # Give acces to individual values in conf_dict by using its key name.
    def getVar (self, varkey ):
        return self.d[varkey]


########################################################################################################################
#   ACTION METHODS
########################################################################################################################

    # TODO PARRALELIZE READ Picking
    def generate_read_dict (self, nread, read_len, repeats, ambiguous, duplicate, mut_freq, pair, min = None, max = None, mean = None, certainty = None):
        """Generate a list of sequence containg the name (with localisation in mother sequence) size and DNA sequence"""

        read_dict = {}
        alpha = beta = ndup = nfail = i = 0 # Counter allowing to iterate in case

        # Calculate the parameters of shape for the beta distribution to mimick DNA shearing distribution by sonication
        if pair:
            alpha, beta = self._beta_shape(min, max, mean, certainty)

        while i < nread:
            # Try to generate a valid read. Return a list with name, length [+frag size if pair], sequence [both if pair]
            read = self._generate_read (read_len, repeats, ambiguous, mut_freq, pair, alpha, beta, min, max)

            # In case it is impossible to generate a valid read in d
            if not read:
                nfail +=1
                # Try to resample until the maximal number of failure to sample a valid read is reached
                if nfail < nread:
                    continue
                break

            # If the dictionnary already contains this entry (same sequence name)
            if read[0] in read_dict:
                ndup +=1

                # If duplicates are not allowed
                if not duplicate:
                # Try to resample until the maximal number of failure to find a non duplicated read is reached
                    if ndup < nread:
                        continue
                    break

                # If duplicates are allowed
                else:
                    # Adding a number to the name (starting at 1) and checking if the new name is non in dict
                    for j in range (ndup):
                        name = '{0}_{1}'.format(read[0], j+1)
                        if name not in read_dict:
                            read[0] = name
                            break

            # Finally if the read was already in dict or if it was renamed add a new entry in dict
            read_dict[read[0]] = read[1:]
            i+=1

        len_dict = len(read_dict)
        self._out_message(duplicate, ndup, nread, len_dict)

        if pair and len_dict != 0:
            self._draw_distribution(read_dict, (max-min)/5 )

        return read_dict


########################################################################################################################
#   PRIVATE METHODS
########################################################################################################################

    def _create_simple_dict(self, nseq, max_len, min_len, repeats = False, ambiguous = False):
        """Simple method that can be used to generate a test dictionnary of Seq record with randomly generated sequences"""
        d={}
        for i in range (nseq):
            name_seq = "Seq#"+str(i)
            len_seq = randint(max_len, min_len)
            seq = Seq(self._random_seq(len_seq,repeats, ambiguous))
            d[name_seq] = SeqRecord(seq, id=name_seq, name=name_seq )

        return d

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

    def _random_seq (self, len_seq, repeats, ambiguous):
        """Return a sequence of lenght len_seq including a 1% chance to add repeats or ambiguous DNA if selected"""

        alphabet = self._IUPAC(repeats, ambiguous)
        seq =''
        for i in range (len_seq):
            if random() <= 0.01:
                seq += alphabet[randint(0,len(alphabet)-1)]
            else:
                seq += 'ATCG'[randint(0,3)]
        return seq

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

    def _IUPAC(self, repeats, ambiguous):
        """Return the IUPAC DNA alphabet corresponding to users requirements"""
        return self._repeats (Ambiguous.letters, repeats) if ambiguous else self._repeats (Unambiguous.letters, repeats)

    def _repeats(self, alphabet, repeats):
        return alphabet + alphabet.lower() if repeats else alphabet

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

    # TODO create an Exception in case no read was found
    def _generate_read (self, read_len, repeats, ambiguous, mut_freq, pair, alpha = None, beta = None, min = None, max = None):
        """Generate a candidate read or read pair of a given lenght with or without repeats and ambiguous DNA bases"""

        # Guard condition if not possible to find a valid pair after 10 tries
        for count in range (10):

            # Pick a random sequence in dictionnary
            ref_seq = self._random_ref()

            # Generate a pseudo-random size resulting from sonication for paired-end reads
            if pair:
                frag_len = self._beta_distrib(alpha, beta, min, max)
            # For single read just use the size of one read
            else:
                frag_len = read_len

            # Start again if the choosen reference sequence is shorter than the lenght of the fragment to be sampled
            if len(self.d[ref_seq]) < frag_len :
                continue

            # Define a random position in the reference sequence and sample a candidate region
            candidate = self._random_candidate (ref_seq, frag_len)

            # Verify the validity of the candidate sequence in terms of repeats and ambiguity
            if self._valid_sequence(str(candidate.seq), repeats, ambiguous):

                # Both extremities of candidate are sampled for pair end
                if pair:
                    return self._extract_pair(candidate, read_len, frag_len, mut_freq)
                # For single end the whole fragment is returned
                else:
                    return [candidate.name, frag_len, self._mutate(candidate, mut_freq)]

        #print ("\nNo valid read found")
        return None

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

    def _random_ref (self):
        """Simple random reference sequence sampler in d"""
        return sample(self.d, 1)[0]

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

    def _random_candidate(self, refseq, size):
        """Return a splice from a biopython Seqrecord in d"""
        startpos = randint(0, len(self.d[refseq])-size)

        # Randomly choose an orientation reverse or forward for the fragment
        if randint(0,1):
            rec = self.d[refseq][startpos:startpos+size]
            rec.id = rec.name = '{0}_{1}-{2}'.format(refseq, startpos, startpos+size)
            rec.description= "Forward"
        else:
            rec = self.d[refseq][startpos:startpos+size].reverse_complement()
            rec.id = rec.name = '{0}_{1}-{2}'.format(refseq, startpos+size, startpos)
            rec.description= "Reverse"

        return rec

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

    def _valid_sequence (self, seq, repeats, ambiguous):
        """Define if the candidate region is valid according to user specifications (allow repeats, allow ambiguous)"""
        valid_base = self._IUPAC(repeats, ambiguous)
        return all ([(base in valid_base) for base in seq])


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

    def _beta_shape(self,min, max, mean, k):
        """Calculate shape parameters alpha and beta to fit experimental indication from user"""
        mode = float((mean - min)) / (max - min)
        a = mode *(k-2) + 1
        b = k - a
        return (a, b)

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

    def _beta_distrib(self, a, b, min, max):
        """Define a pseudorandom size according to a beta distribution giving alpha and beta"""
        return int(betavariate(a,b)*(max-min) + min)

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

    def _extract_pair(self, candidate, read_len, frag_len):
        """Extract reads forward and reverse from a candidate sequence and return a list with name, and both records"""

        forward = candidate[:read_len]
        reverse = candidate[-read_len:].reverse_complement()
        forward.id = forward.name = reverse.id = reverse.name = candidate.name

        forward.description = self._pair_overlap(read_len, frag_len) + "__R1"
        reverse.description = self._pair_overlap(read_len, frag_len) + "__R2"

        return [candidate.name, frag_len, self._mutate(forward, mut_freq), self._mutate(reverse, mut_freq)]

    def _mutate(self, candidate, mut_freq):
        if mut_freq == 0:
             candidate.description = self._pair_overlap(read_len, frag_len) + "__R1"





## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

    def _pair_overlap(self,read_len, frag_len):
        return "No_overlap" if frag_len > 2*read_len else "{0}_bp_overlap".format(2*read_len-frag_len)

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

    def _out_message (self, duplicate, ndup, nread, len_dict):
        """ Print an output massage according to the success or failure to generate the required sequences"""

        if len_dict == 0:
            print ('\nThe dictionnary is empty. Parameters or reference sequences are not suitable')
        elif len_dict < nread:
            print ('\nImpossible to generate the required number of reads ({0}). The dictionnary contains only {1} entries'.format(nread, len_dict))
        elif duplicate and ndup >=1 :
            print ('\nThe dictionnary contains the required number of sequences ({0}) but includes {1} duplicate(s)'.format(len_dict, ndup))
        else:
            print ('\nThe dictionnary contains the required number of sequences ({0}) and no duplicate'.format(len_dict))

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

    def _draw_distribution(self, read_dict, bin_range):
        """Use Pyplot histogram function to draw an histogram of fragment size distribution"""
        list_len = []
        # Extract frag_len from read_dict and append them in the liste
        for val in read_dict.values():
            list_len.append(val[0])

        # Represent data with an histogramm using pyplot
        h = plt.hist(list_len, bins = bin_range, normed = True)
        plt.show()

