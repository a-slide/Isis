from random import random
from random import randint
from random import triangular

class Sequence:
    """This can be used to parse file containing a list of variables and its associated values
    The separator beetween vaiable and value must be the same for all field and can be customized by users (by default =)
    the class store the list of variable in a dictionnary named conf_dict"""

# DECLARATION OF VARIABLES
    # Dictionnary to store sequences dictionary
    seq_dict = {}
    # Store a probability list for each sequence in seq_dict which is calculated according to the respective size           
    # of each sequences in the seq_dict.
    proba_list = []

# FONDAMENTAL METHODS

    # Simple object constructor to test the mother class initialized with a name and a sequence
    def __init__ (self):
        self.seq_dict = {
            "chr1" : [96, "ATAGCGTCGCTAGCTCGTCAGCTGCATACGTAGCTATCgcgcgatcgattacgCGGGGGGGCATGCTAGTAGCATATACATGCGACAGTCAGCTCA"],
            "chr2" : [90, "TCGCGGGGGCTATTAGATATTCGTCCTNNNNNNNNNNGCGCGCGATTACTATTTTTTAGCGCGGACTCGCGGCGGCGCGGCGGCGCTGGA"],
            "chr3" : [92, "TATCTCGTCGCGCGGTTTTTCGTAGCTAGagtactATCAGTAGCTCGACTCGATCGCGCGCATGTCGTCGCGGGGGTTCTTTTTTTTTTTTT"],
            "chr4" : [80, "ATTGTCGCTGCTCGATCGCATGCGCGGCGCGGCGCGCGGCGATNNTATATTTTTTAATATAAAATAGCGGCGCTTCGCCC"]}

        self.proba_list = self._calculate_proba()
        
        self._DNA_bases_strict = "ATCG"
        self._DNA_bases_repeats = "ATCGatcg"
        self._DNA_bases_ambigous = "ATCGKMRYSWBVHDN"
        self._DNA_bases_ambigous_repeats = "ATCGKMRYSWBVHDNatcgkmryswbvhdn"

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

    def pair_generator (self, read_len = 150 , min_sonic = 150, max_sonic = 1000, mean_sonic = 500, allow_repeats = False, allow_ambigous = False):
        pass
        
        #while (guard condition):
            #refseq = random_seq (self)
            #size = valid_region_size (self, min_sonic, max_sonic, mean_sonic)
            #readseq = random_region (self, size)
            #return ... if valid_region(self, avoid_repeats,candidate)


    #TODO raise no valid pair exception.......
    def read_generator (self, read_len = 150 , allow_repeats = False, allow_ambigous = False):
        """Generate a candidate read choosen proportionally in seq_dict and satisfying user requirements"""
        
        # Guard condition if not possible to find a valid pair after 10000 tries
        for count in range (10000):
            
            # Pick a random sequence in dictionnary proportionally to its length
            ref_seq = self._random_seq ()
            
            # Start again if the reference choosen reference sequence is shorter than the lenght of reads
            if self.seq_dict[ref_seq][0] < read_len :
                print "The reference is too short to sample a valid read. New try "
                continue
                
            # Define a random position in the reference sequence and sample a candidate region
            candidate_seq = self._random_candidate_sequence (ref_seq, read_len)
            
            # Return the candidate if it mets user request in terms of repeats and ambiguity
            if self._valid(candidate_seq[2], allow_repeats, allow_ambigous):
                return candidate_seq
            
            print "Fail to sample a valid read. New try"
            
        # else if impossible to find a valid pair = raise no valid pair exception.......
        print ("No valid pair found")
        return None
       

# PRIVATE METHODS


    # Return a 2 entries list / 1 = name of the sequence / 2 = cumulative frequency of the sequence
    def _calculate_proba(self):
        proba_list = []
        cumulative_len = 0
        # Fill the list with name and cumulative length for each seq in seq_dict
        for name, info in self.seq_dict.items():
            cumulative_len += info[0]
            proba_list.append([name, cumulative_len])
        # Convert cumulative lengths in cumulative frquencies
        for item in proba_list:
            item[1] =  float(item[1])/cumulative_len
        return proba_list
        
        
    # Return a random sequence from seq_dict according to the probability vector
    def _random_seq (self):
        rand_freq = random()
        for name, freq  in self.proba_list:
            if freq > rand_freq :
                return name
                
                
    # Return a list with [refseq name+startpos+endpos+orientation, size, sequence]
    # TODO = Random orientation of reads
    def _random_candidate_sequence(self, refseq, size):
        startpos = randint (0, self.seq_dict[refseq][0]-size)
        sequence = self.seq_dict[refseq][1][startpos:startpos+size]
        name = '{0}_{1}-{2}_forward'.format(refseq, str(startpos), str(startpos+size))
        return [name, size, sequence]
        
        
    def _valid (self, candidate_seq, allow_repeats, allow_ambigous):
        """Define if the candidate region is valid according to user specifications (allow repeats, allow ambigous)"""
        if not allow_repeats:
            if not allow_ambigous:
                # no repeats and no ambigous DNA base
                return self._valid_region(candidate_seq, self._DNA_bases_strict)
            else:
                # no repeats but ambigous DNA base
                return self._valid_region(candidate_seq, self._DNA_bases_ambigous)
        else:
            if not allow_ambigous:
                # Repeats but no ambigous DNA base
                return self._valid_region(candidate_seq, self._DNA_bases_repeats)
            else:
                # Repeats and ambigous DNA base
                return self._valid_region(candidate_seq, self._DNA_bases_ambigous_repeats)
                
                
    def _valid_region(self, candidate_seq, alphabet):
        """Define if the sequence contains only the letters in alphabet"""
        for i in candidate_seq:
            if i not in alphabet:
                print ("Invalid Sequence")
                return False
            else:
                print ("Valid Sequence")
                return True
                
                
    def _region_size (self, min_sonic, max_sonic, mean_sonic):
        #Determine a random size of sonication based on a gausian distribution
        #gauss(mu, sigma) # mu is the mean, and sigma is the standard deviation
        #remove extreme size = bellow min or above max
        #int(random.triangular(low, high, mode))
        pass

#class Reference (Sequence):
    #pass


#class Junction (Sequence):
    #pass
