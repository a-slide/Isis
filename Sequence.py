from random import gauss
from random import random

class Sequence:
    """This can be used to parse file containing a list of variables and its associated values
    The separator beetween vaiable and value must be the same for all field and can be customized by users (by default =)
    the class store the list of variable in a dictionnary named conf_dict"""

# DECLARATION OF VARIABLES
    # Dictionnary to store sequences dictionary
    seq_dict = {}
    # Store a probability list for each sequence in seq_dict which is calculated according to the respective size
    # of each sequences in the seq_dict. First field =
    proba_list = []

# FONDAMENTAL METHODS

    # Simple object constructor to test the mother class initialized with a name and a sequence
    def __init__ (self):
        self.seq_dict = {
            "chr1" : "ATAGCGTCGCTAGCTCGTCAGCTGCATCAGTCAGCTCA",
            "chr2" : "TAGAGCGCGCGGGGGCTATTAGATATTGGA",
            "chr3" : "TATAATGCTCGTCGCGCGGTTTTTCGTAGCTCAGACTGCATCAGTAGCTCGACTCGATCGCTCTTTTTTTTTTTTT",
            "chr4" : "ATTATTAAAAAAAAAAA"}

        self.proba_list = self._calculate_proba()

# GETERS

    # Grant acces to the complete dictionary
    def getDict (self):
        return self.seq_dict

    # Give acces to individual values in conf_dict by using its key name.
    def getVar (self, varkey ):
        return self.seq_dict[varkey]

    # Give acces to individual values in conf_dict by using its key name.
    def getProba (self):
        return self.proba_list

# ACTION METHODS

    # Create a generator
    def pair_generator (self, read_len = 150 , min_sonic = 150, max_sonic = 1000, mean_sonic = 500, avoid_repeats = True):


          while (guard condition):

            return ... if ...

        refseq = random_seq (self)
        #size = valid_region_size (self, min_sonic, max_sonic, mean_sonic)
        #readseq = random_region (self, size)
        #UNTIL valid_region(self, avoid_repeats,candidate)
        #yield read_pair from readseq with end and start positions
        pass

    def read_generator (self, read_len = 150 , avoid_repeats = True):
        #DO
            #seq = random_seq (self)
            #readseq = random_region (self, read_len)
        #UNTIL valid_region(self, avoid_repeats, read)
        #yield read with end and start positions
        pass

# PRIVATE METHODS

    # Return a 3 entries list / 1 = name of the sequence / 2 = size of the sequence / 3 = cumulative frequency of the sequence
    def _calculate_proba(self):
        list_len = []
        cumulative_len = 0

        # Fill the list with name, length and cumulative length for each seq in seq_dict
        for name, sequence in self.seq_dict.items():
            seq_len = len(sequence)
            cumulative_len += seq_len
            list_len.append([name, seq_len, cumulative_len])

        # Convert cummulative length in cum
        for item in list_len:
            item[2] =  float(item[2])/cumulative_len

        return list_len

    # Return a random sequence from seq_dict according to the probability vector
    def _random_seq (self):

        print(self.getProba())
        rand_freq = random()
        print rand_freq

        for name, length, freq  in self.proba_list:
            if freq > rand_freq :
                return name


    def _valid_region_size (self, min_sonic, max_sonic, mean_sonic):
        #Determine a random size of sonication based on a gausian distribution
        #gauss(mu, sigma) # mu is the mean, and sigma is the standard deviation
        #remove extreme size = bellow min or above max
        pass

    def _random_region (self, size): # choose random region of size length
        pass

    def _valid_region (self, avoid_repeats): # no N and or repeats in reads
        pass



#class Reference (Sequence):
    #pass


#class Junction (Sequence):
    #pass
