from random import random
from random import randint
import gzip

from Bio import SeqIO

class ReferenceGenome:

########################################################################################################################
#   FONDAMENTAL METHODS
########################################################################################################################

    def __init__(self, source, filename):
        """Object constructor importing reference sequences from fasta file"""
        # Dictionnary of bioPython record created from fasta file
        self.source = source
        self.d = self._import_fasta (filename)
        # List cummulative probabilities of each sequence to be picked calculated from to their respective size.
        self.proba_list = self._calculate_proba()

    def __repr__(self):
        """Long description string used by interpreter and repr"""
        result = "{}<Source = {}>\n\n".format(self.__str__(), self.source)
        for entry in self.d.values():
            result += "<{}\nLenght:{}>\n".format(entry, len(entry))
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
    
    # Give acces to the name of the source
    def getSource (self):
        return self.source
    
    # Give acces to individual values in d by using its key name.
    def get (self, varkey ):
        return self.d[varkey]
    
    # Give acces to the complete dictionary
    def getDict (self):
        return self.d
    
    # Acces to the cumulative frequency list
    def getProba (self):
        return self.proba_list

########################################################################################################################
#   PUBLIC METHODS
########################################################################################################################

    def get_slice (self, size):
        """Generate a candidate slice and return it"""
        
        # Guard condition if not possible to find a valid slice after 100 tries
        for count in range (100):
            # Pick a random sequence in dictionnary
            refseq = self._random_refseq()

            # Start again if the choosen reference sequence is shorter than the lenght of the fragment to be sampled
            if len(self.d[refseq]) >= size :
                return self._random_slice (refseq, size)
        
        # if no valid size was found
        raise Exception ("The size of the slice seems to be too long to be sampled in reference sequences")

########################################################################################################################
#   PRIVATE METHODS
########################################################################################################################

    def _import_fasta (self, filename):
        """Import fasta files in a dictionary of biopython SeqRecord"""

        try: # try to open the file
            if filename.rpartition(".")[-1] == "gz":
                print ("Uncompressing and extracting data")
                handle = gzip.open(filename, "r")
            else:
                print ("Extracting data")
                handle = open(filename, "r")

            d = SeqIO.to_dict(SeqIO.parse( handle, "fasta"))
            handle.close()

            return d

        except IOError:
               print ('CRITICAL ERROR. The fasta file ' + filename + ' is not readable. Exit')
               exit

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

    def _calculate_proba(self):
        """Return a 2 entries list / 1 = name of the sequence / 2 = cumulative frequency of the sequence"""

        l = []
        cumulative_len = 0.0
        total_len = sum([len(record) for record in self.d.values()])

        for record in self.d.values():
            cumulative_len += len(record)
            l.append([record.name, cumulative_len/total_len])

        return l

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

    def _random_refseq (self):
        """Return a random sequence from d according the respective size of references"""

        # Define a pseudo-random decimal frequency
        rand_freq = random()

        # Attibute this frequency to a sequence from d based on proba_list
        for name, freq  in self.proba_list:
            if freq > rand_freq :
                return name

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

    def _random_slice(self, refseq, size):
        """Return a slice from a biopython Seqrecord in d"""
        
        # Randomly choose the slice start position
        start = randint(0, len(self.d[refseq])-size)
        end = start+size
        
        # Randomly choose an orientation reverse or forward for the fragment
        forward = randint(0,1)
        
        slice = self.d[refseq][start:end] if forward else self.d[refseq][start:end].reverse_complement()
        
        slice.id = slice.name = slice.description = ""
        slice.annotations = {
        "orientation"   : "+" if forward else "-",
        "source"        : self.source,
        "refseq"        : refseq,
        "location"      : [start, end]}
        
        return slice




