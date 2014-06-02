# Standard library packages
from random import random, randint
import gzip
import csv

# Third party packages
from Bio import SeqIO # Require Biopython

####################################################################################################

class ReferenceGenome(object):
    """Import reference sequences from a fasta files and store each of
    them in a dictionary. The class provide a method get a random
    biopython seqrecord slice from one of the reference sequences,
    proportionally to the size of all reference sequences.
    """

###    FONDAMENTAL METHODS    ###

    def __init__(self, name, filename):
        """Import reference sequences from fasta file and create a list
        of probability to sample in each ref sequence
        """

        print "Initialisation of {}...".format(name)

        # Store object variables
        self.name = name

        # Dictionnary of bioPython record created from fasta file
        self.d = self._import_fasta(filename)

        # Initialize a counter for each reference that will be
        # incremented each time _random_slice choose this reference
        self.reset_samp_counter()

        # List cummulative probabilities of each sequence to be picked
        # calculated from to their respective size.
        self.proba_list = self._calculate_proba()

    def __repr__(self):
        result = "{}\n".format(self.__str__())
        for entry in self.d.values():
            result += "{}\nLenght:{}\n".format(entry, len(entry))
        return result

    def __str__(self):
        return "<Instance of {} from {} >".format(self.__class__.__name__, self.__module__)

###    GETERS    ###

    def get(self, varkey):
        return self.d[varkey]

    def getDict(self):
        return self.d

    def getProba(self):
        return self.proba_list

    def getName(self):
        return self.name

###    PUBLIC METHODS    ####

    def get_slice(self, size):
        """Generate a candidate slice and return it
        """
        # Guard condition
        for count in range(100):
            # Pick a random sequence in dictionnary
            refseq = self._random_refseq()

            # If the size is valid return a slice
            if size <= len(self.d[refseq]):
                return self._random_slice(refseq, size)

        # if no valid size was found
        raise Exception("ERROR. Unable to find a slice of a valid length after 100 tries.\n\
        Please review the size of your references or the size of the reads.")

    def reset_samp_counter(self):
        """ Initialize or reset counter for each reference that will be
        incremented each time _random_slice choose this reference
        """
        for name, record in self.d.items():
            record.annotations ["nb_samp"] = 0

    def write_samp_report (self):
        """ Create a simple csv report
        """
        # Open a file for writting with python csv module
        with open(self.name+"_report.csv", 'w') as csvfile:
            writer = csv.writer(csvfile, delimiter='\t', quoting=csv.QUOTE_NONE)
            writer.writerow(["chr", "lenght", "nb_samp"])

            # Create a sorted list of refseq
            ref_list = self.d.keys()
            ref_list.sort()

            # Export each refseq characteristics in a file
            for ref in ref_list:
                writer.writerow([ref, len(self.d[ref]), self.d[ref].annotations["nb_samp"]])

###    PRIVATE METHODS    ###

    def _import_fasta(self, filename):
        """Import fasta files in a dictionary of biopython SeqRecord
        """
        # Try to open the file fist gz compressed and uncompressed
        try:
            if filename.rpartition(".")[-1] == "gz":
                print("\tUncompressing and extracting data")
                handle = gzip.open(filename, "r")
            else:
                print("\tExtracting data")
                handle = open(filename, "r")

            d = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
            handle.close()

            return d

        # Try to open the file fist gz compressed and uncompressed
        except IOError:
               print('CRITICAL ERROR. The fasta file ' + filename + ' is not readable. Exit')
               exit

    def _calculate_proba(self):
        """Return a 2 entries list with the name of the sequence and a
        cumulative frequency of the sequence
        """
        proba_list = []
        cumulative_len = 0.0
        total_len = sum([len(record) for record in self.d.values()])

        for record in self.d.values():
            cumulative_len += len(record)
            proba_list.append([record.name, cumulative_len/total_len])

        return proba_list

    def _random_refseq(self):
        """Return a random sequence from d according the respective size
         of references
        """
        # Define a pseudo-random decimal frequency
        rand_freq = random()

        # Attibute this frequency to a sequence from d based on proba_list
        for name, freq  in self.proba_list:
            if freq > rand_freq :
                return name

    def _random_slice(self, refseq, size):
        """Return a slice from a biopython Seqrecord in d
        """
        # Randomly choose the slice start position
        start = randint(0, len(self.d[refseq])-size)
        end = start+size

        # Randomly choose an orientation reverse or forward and sample
        # a slice
        if randint(0,1):
            s = self.d[refseq][start:end]
            s.annotations["orientation"] = "+"
        else:
            s = self.d[refseq][start:end].reverse_complement()
            s.annotations["orientation"] = "-"

        # Add informations to the annotations dictionnary
        s.annotations["refseq"] = refseq
        s.annotations["location"] = [start, end]
        s.annotations["source"] = self

        # Undefine description, id and name
        s.name = s.description = s.id = ""

        # Increment the sampling counter of the refseq
        self.d[refseq].annotations["nb_samp"] += 1

        return s
