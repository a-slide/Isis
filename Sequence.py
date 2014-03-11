from random import gauss

class Sequence:
    """This can be used to parse file containing a list of variables and its associated values
    The separator beetween vaiable and value must be the same for all field and can be customized by users (by default =)
    the class store the list of variable in a dictionnary named conf_dict"""
    
    seq_dict{}
    probability vector[]

# GETERS
    
    # Grant acces to the complete dictionary
    def getDict (self):
        return self.conf_dict
    
    # Give acces to individual values in conf_dict by using its key name.
    def getVar (self, varkey ):
        return self.seq_dict[varkey]
        
# ACTION METHODS
    
    # Create a generator 
    def pair_generator (self, read_len = 150 , min_sonic = 150, max_sonic = 1000, mean_sonic = 500, avoid_repeats = True):
        DO
            refseq = random_seq (self)
            size = valid_region_size (self, min_sonic, max_sonic, mean_sonic)
            readseq = random_region (self, size)
        UNTIL valid_region(self, avoid_repeats,,candidate)
        
        yield read_pair from readseq with end and start positions
    
    def read_generator (self, read_len = 150 , avoid_repeats = True):
        DO
            seq = random_seq (self)
            readseq = random_region (self, read_len)
        UNTIL valid_region(self, avoid_repeats, read)
        
        yield read with end and start positions
        
# PRIVATE METHODS
    
    def _random_seq (self):
        choose random seq in seq_dict based probability_vector
    
    def _valid_region_size (self, min_sonic, max_sonic, mean_sonic)
        Determine a random size of sonication based on a gausian distribution
        gauss(mu, sigma) # mu is the mean, and sigma is the standard deviation
        remove extreme size = bellow min or above max
    
    def _random_region (self, size): # choose random region of size length
    

    def _valid_region (self, avoid_repeats) # no N and or repeats in reads
    
    
class Reference (Sequence):
    pass
    
    
class Junction (Sequence):
    pass
