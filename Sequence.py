from random import random
from random import randint
from random import betavariate
from string import maketrans

#To visualyse beta distribution
#import matplotlib.pyplot as plt
#liste = []
#for i in range (10000):
#   liste.append(int(betavariate(2.125, 2.875)* 800 + 200))
#h = plt.hist(liste, bins = 100, normed = True)
#plt.show()

class Sequence:
    """"""


########################################################################################################################
#   FONDAMENTAL METHODS
########################################################################################################################

    # Simple object constructor to test the mother class
    def __init__ (self):
        
        # Dictionnary to store name, size and sequences
        self.seq_dict = {
            "chr1" : [168, "ATAGCGTCGCTAGCTCGTCAGCTGCATACGTAGCTATCgcgcgatcgattacgCGGGGGGGCATGCTAGTAGGGGCATGCTAGTAGCATATACATGCGACAGTCAGCTGGGCATGCTAGTAGCATATACATGCGACAGTCAGCTCATATACATGCGACAGTCAGCTCA"],
            "chr2" : [162, "TCGCGGGGGCTATTGGGCATGCTAGTAGCATATACATGCGACAGTCAGCTGGGCATGCTAGTAGCATATACATGCGACAGTCAGCTAGATATTCGTCCTNNNNNNNNNNGCGCGCGATTACTATTTTTTAGCGCGGACTCGCGGCGGCGCGGCGGCGCTGGA"],
            "chr3" : [164, "TATCTCGTCGCGCGGTTTTTCGTAGCTAGagtactATCAGTAGCTCGACTCGGGGCATGCTAGTAGCATATACATGCGACAGTCAGCTGGGCATGCTAGTAGCATATACATGCGACAGTCAGCTATCGCGCGCATGTCGTCGCGGGGGTTCTTTTTTTTTTTTT"],
            "chr4" : [116, "ATTGTCGCTGCTCGATCGCATGCGCGGCGCGGCGCGCGGCGATNNGGGCATGCTAGTAGCATATACATGCGACAGTCAGCTTATATTTTTTAATATAAAATAGCGGCGCTTCGCCC"]}
            
        # List cummulative probabilities of each sequence to be picked calculated from to their respective size.
        self.proba_list = self._calculate_proba()
        
        # Strings containing allowed DNA bases
        self._DNA_strict = "ATCG"
        self._DNA_repeats = "ATCGatcg"
        self._DNA_ambigous = "ATCGKMRYSWBVHDN"
        self._DNA_ambigous_repeats = "ATCGKMRYSWBVHDNatcgkmryswbvhdn"
        
        # Complement table to generate DNA complementary strand
        self._complement_table = maketrans("AGCTYRWSKMDVHBNagctyrwskmdvhbn", "TCGARYWSMKHBDVNtcgarywsmkhbdvn")


########################################################################################################################
#   GETERS
########################################################################################################################

    # Grant acces to the complete dictionary
    def getDict (self):
        return self.seq_dict

    # Give acces to individual values in conf_dict by using its key name.
    def getVar (self, varkey ):
        return self.seq_dict[varkey]

    # Give acces to individual values in conf_dict by using its key name.
    def getProba (self):
        return self.proba_list


########################################################################################################################
#   ACTION METHODS
########################################################################################################################

    # TODO Catch exception in case no read was found

    def generate_read_dict (self, nread, read_len, repeats = False, ambigous = False, duplicate = False,
                            pair = False, min = None, max = None, mean = None):
        """Generate a dictionnay of sequence containg the name (with localisation in mother sequence) size and DNA sequence"""
        
        read_dict = {}
        ndup = i = 0
        alpha = beta = None
        k = 10 # To be ajusted. Lower value = larger with, Higher value = thiner peak. Should not be bellow 2
        
        # Calculate the parameters of shape for the beta distribution to mimick DNA shearing distribution by sonication
        if pair:
            alpha, beta = self._beta_shape(min, max, mean, k)
            print ("APLHA = ", alpha, "BETA = ", beta)
        
        while i < nread:
            read = self._generate_read (read_len, repeats, ambigous, pair, alpha, beta, min, max)
            
            # In case it is impossible to generate a valid read in seq_dict
            if not read:
                return None
            
            # If the dictionnary already contains this entry (same sequence name)
            if read[0] in read_dict:
                ndup +=1
                
                # If duplicates are not allowed
                if not duplicate:
                    # If the maximal number of tries to generate non duplicated reads was not yet exceded = cancel count and try to resample
                    if ndup < nread:
                        #print (read[0]+" already exists in the dictionnary. Resampling")
                        continue
                    # If the maximal number of tries to generate non duplicated reads was exceded = return the 
                    else:
                        print ('\nImpossible to generate the required number of reads. The dictionnary contains only {0} entries'.format(len(read_dict)))
                        return read_dict
                        
                # If duplicates are allowed
                else:
                    #print (read[0]+" already exists in the dictionnary. Renaming and adding")
                    # Adding a number to the name (starting at 2) and checking if the new name is non in dict
                    for j in range (ndup):
                        name = '{0}_{1}'.format(read[0], j+1)
                        if name not in read_dict:
                            read[0] = name
                            break
                            
            # Finaly if the read was already in dict or if it was renamed add a new entry in dict
            read_dict[read[0]] = [read[1], read[2]]
            i+=1
        
        if duplicate:
            print ('\nThe dictionnary contains the required number of sequences but includes {0} duplicate(s)'.format(ndup))
        else:
            print ('\nThe dictionnary contains the required number of sequences and no duplicate')
        
        return read_dict
        

########################################################################################################################
#   PRIVATE METHODS
########################################################################################################################

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
        
        
    # TODO create an Exception in case no read was found
    def _generate_read (self, read_len, repeats, ambigous, pair, alpha = None, beta = None, min = None, max = None):
        """Generate a candidate read or read pair of a given lenght with or without repeats and ambigous DNA bases"""
        
        # Guard condition if not possible to find a valid pair after 100 tries
        for count in range (100):
            
            # Pick a random sequence in dictionnary proportionally to its length
            ref_seq = self._random_seq ()
            
            # Generate a pseudo-random size resulting from sonication for paired-end reads
            if pair:
                frag_len = self._beta_distrib(alpha, beta, min, max)
            # For single read just use the size of one read
            else:
                frag_len = read_len
            
            # Start again if the choosen reference sequence is shorter than the lenght of the fragment to be sampled
            if self.seq_dict[ref_seq][0] < frag_len :
                continue
                
            # Define a random position in the reference sequence and sample a candidate region
            candidate_seq = self._random_candidate_sequence (ref_seq, frag_len)
            
            # Verify the valibility of the candidate sequence in terms of repeats and ambiguity
            if self._valid(candidate_seq[2], repeats, ambigous):
                
                # Both extremities of candidate are sampled for pair end, but for single end candidate
                if pair:
                    return self._extract_pair(candidate_seq, read_len, frag_len)
                else:
                    return candidate_seq
                    
        print ("\nNo valid read found")
        return None
    
    
    def _random_seq (self):
        """Return a random sequence from seq_dict according the respective size of references"""

        # Define a pseudo-random decimal frequency
        rand_freq = random()
        
        # Attibute this frequency to a sequence from seq_dict based on proba_list
        for name, freq  in self.proba_list:
            if freq > rand_freq :
                return name
                
                
    def _random_candidate_sequence(self, refseq, size):
        """Return a list with [refseq name+startpos+endpos, size, sequence]"""
        startpos = randint (0, self.seq_dict[refseq][0]-size)
        sequence = self.seq_dict[refseq][1][startpos:startpos+size]
        
        # Randomly choose an orientation reverse or forward for the fragment
        if randint(0,1):
            name = '{0}_{1}-{2}'.format(refseq, startpos, startpos+size)
        else:
            name = '{0}_{1}-{2}'.format(refseq, startpos+size, startpos)
            sequence = self._rc(sequence)
            
        return [name, size, sequence]
        
        
    def _valid (self, seq, repeats, ambigous):
        """Define if the candidate region is valid according to user specifications (allow repeats, allow ambigous)"""
        if not repeats:
            if not ambigous:
                # no repeats and no ambigous DNA base
                return self._valid_region(seq, self._DNA_strict)
            else:
                # no repeats but ambigous DNA base
                return self._valid_region(seq, self._DNA_ambigous)
        else:
            if not ambigous:
                # Repeats but no ambigous DNA base
                return self._valid_region(seq, self._DNA_repeats)
            else:
                # Repeats and ambigous DNA base
                return self._valid_region(seq, self._DNA_ambigous_repeats)
                
                
    def _valid_region(self, seq, alphabet):
        # determine if the sequence contains only the letters in alphabet
        
        # if any of seq characters are not in the alphabet = return false
        for i in seq:
            if i not in alphabet:
                return False
        
        # If all characters of seq are in the alphabet = return true
        return True
        
        
    def _rc (self, sequence):
        """Return the reverse complementary of any DNA sequence"""
        return sequence.translate(self._complement_table)[::-1]
        
        
    def _beta_shape(self,min, max, mean, k):
        """Calculate shape parameters alpha and beta to fit experimental indication from user"""
        mode = (mean - min) / (max - min)
        a = mode *(k-2) + 1
        b = k - a
        return (a, b)
        
        
    def _beta_distrib(self, a, b, max, min):
        """Define a pseudorandom size according to a beta distribution giving alpha and beta"""
        return int(betavariate(a, b)*(max-min) + min)
        
        
    def _extract_pair(self, candidate_seq, read_len, frag_len):
        """Extract reads forward and reverse from a candidate sequence and return a list with name,read size, frag length and both sequences"""
        forward_read = candidate_seq[2][0:read_len]
        reverse_read = self._rc(candidate_seq[2][-read_len:])
        return [candidate_seq[0],read_len, frag_len, forward_read, reverse_read]

###################################################
#
# TO DO =   beta dis seems inverted
#           forward_read and reverse_read are not added to the list
#
###################################################




#class Reference (Sequence):
    #pass


#class Junction (Sequence):
    #pass
