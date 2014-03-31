import gzip
from Sequence import Sequence
from Bio import SeqIO

class Reference (Sequence):

########################################################################################################################
#   FONDAMENTAL METHODS
########################################################################################################################

    def __init__(self, filename):
        """Object constructor importing reference sequences from fasta file"""
        # Dictionnary of bioPython record created from fasta file
        self.d = self._import_fasta (filename)
        # List cummulative probabilities of each sequence to be picked calculated from to their respective size.
        self.proba_list = self._calculate_proba()


    # TODO Write this documentation
    def __doc__ (self):
        pass


########################################################################################################################
#   GETERS
########################################################################################################################

    # Acces to the cumulative frequency list
    def getProba (self):
        return self.proba_list


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


    def _calculate_proba(self):
        """Return a 2 entries list / 1 = name of the sequence / 2 = cumulative frequency of the sequence"""

        l = []
        cumulative_len = 0.0
        total_len = sum([len(record) for record in self.d.values()])

        for record in self.d.values():
            cumulative_len += len(record)
            l.append([record.name, cumulative_len/total_len])

        return l


    def _random_seq (self):
        """Return a random sequence from d according the respective size of references"""

        # Define a pseudo-random decimal frequency
        rand_freq = random()

        # Attibute this frequency to a sequence from d based on proba_list
        for name, freq  in self.proba_list:
            if freq > rand_freq :
                return name

