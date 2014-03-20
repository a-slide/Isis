from random import random
from random import randint
from random import betavariate
from string import maketrans
from pprint import pprint as pp
import matplotlib.pyplot as plt

class Sequence:
    """"""

########################################################################################################################
#   FONDAMENTAL METHODS
########################################################################################################################

    def __init__ (self):
        
        # Empty dictionnary to store sequences and description
        self.seq_dict = {}
        
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


########################################################################################################################
#   ACTION METHODS
########################################################################################################################

    # TODO Catch exception in case no read was found
    # PARRALELIZE READ Picking
    def generate_read_dict (self, nread, read_len, repeats = False, ambigous = False, duplicate = False,
                            pair = False, min = None, max = None, mean = None, certainty = None):
        """Generate a dictionnay of sequence containg the name (with localisation in mother sequence) size and DNA sequence"""
        
        read_dict = {}
        ndup = i = 0
        alpha = beta = None
        
        # Calculate the parameters of shape for the beta distribution to mimick DNA shearing distribution by sonication
        if pair:
            alpha, beta = self._beta_shape(min, max, mean, certainty)
        
        while i < nread:
            read = self._generate_read (read_len, repeats, ambigous, pair, alpha, beta, min, max)
            
            # In case it is impossible to generate a valid read in seq_dict
            if not read:
                break
            
            # If the dictionnary already contains this entry (same sequence name)
            if read[0] in read_dict:
                ndup +=1
                
                # If duplicates are not allowed
                if not duplicate:
                    # If the maximal number of tries to generate non duplicated reads was not yet exceded = cancel count and try to resample
                    if ndup < nread:
                        continue
                    # If the maximal number of tries to generate non duplicated reads was exceded break out the loop  
                    break
                        
                # If duplicates are allowed
                else:
                    # Adding a number to the name (starting at 2) and checking if the new name is non in dict
                    for j in range (ndup):
                        name = '{0}_{1}'.format(read[0], j+1)
                        if name not in read_dict:
                            read[0] = name
                            break
                            
            # Finaly if the read was already in dict or if it was renamed add a new entry in dict
            read_dict[read[0]] = [read[1], read[2]]
            i+=1
        
        len_dict = len(read_dict)
        self._out_message(duplicate, ndup, nread, len_dict)
        
        if pair and len_dict != 0:
            self._draw_distribution(read_dict)
        
        return read_dict
        

########################################################################################################################
#   PRIVATE METHODS
########################################################################################################################      
    
    
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
        """Pick a random sequence in seq_dict. Usage shall be different for subclasses"""
        pass
                
                
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
        mode = float((mean - min)) / (max - min)
        a = mode *(k-2) + 1
        b = k - a
        return (a, b)
        
        
    def _beta_distrib(self, a, b, min, max):
        """Define a pseudorandom size according to a beta distribution giving alpha and beta"""
        return int(betavariate(a,b)*(max-min) + min)
        
        
    def _extract_pair(self, candidate_seq, read_len, frag_len):
        """Extract reads forward and reverse from a candidate sequence and return a list with name,read size, frag length and both sequences"""
        forward_read = candidate_seq[2][0:read_len]
        reverse_read = self._rc(candidate_seq[2][-read_len:])
        return [candidate_seq[0],[read_len, frag_len], [forward_read, reverse_read]]
        
                
    def _out_message (self, duplicate, ndup, nread, len_dict):
        """ Print an output massage according to the success or failure to generate the required sequences"""
        
        if len_dict == 0:
            print ('\nThe dictionnary is empty. Paramaters or reference sequences are not suitable')
        elif len_dict < nread:
            print ('\nImpossible to generate the required number of reads ({0}). The dictionnary contains only {1} entries'.format(nread, len_dict))
        elif duplicate and ndup >=1 :
            print ('\nThe dictionnary contains the required number of sequences ({0}) but includes {1} duplicate(s)'.format(len_dict, ndup))
        else:
            print ('\nThe dictionnary contains the required number of sequences ({0}) and no duplicate'.format(len_dict))
        
        
    def _draw_distribution(self, read_dict):
        """Use Pyplot histogram function to draw an histogram of fragment size distribution"""
        list_len = []
        # Extract frag_len from read_dict and append them in the liste
        for val in read_dict.values():
            list_len.append(val[0][1])
        
        # Represent data with an histogramm using pyplot
        h = plt.hist(list_len, bins = 100, normed = True)
        plt.show()


########################################################################################################################
########################################################################################################################
########################################################################################################################


class Reference (Sequence):

########################################################################################################################
#   FONDAMENTAL METHODS
########################################################################################################################

    def __init__(self):
        """Object constructor importing reference sequences from fasta file"""
        super().__init__()
        
        # Dictionnary of sequences storing name, size and sequence string
        self.seq_dict = self._create_seq_dict(filename)

        # List cummulative probabilities of each sequence to be picked calculated from to their respective size.
        self.proba_list = self._calculate_proba()
    
    def __repr__(self):
        result = "<Instance of Reference(Sequence)>\n"
        for key, value in self.seq_dict.items():
            result += "<{0}\tSize : {1}\tSequence : {2}...{3}>\n".format(key, value[0], value[1][:10], value[1][-10:])
        
        return result

    def __str__(self):
        return "<Instance of Reference(Sequence)>\n"
    
    
########################################################################################################################
#   GETERS
########################################################################################################################

    # Give acces to the cummulative frequency list
    def getProba (self):
        return self.proba_list


########################################################################################################################
#   ACTION METHODS
########################################################################################################################


########################################################################################################################
#   PRIVATE METHODS
########################################################################################################################
    
    def _import_seq_dict(self, filename):
        """Ouvre un fichier contenant des séquence au format fasta et les ajoute au dictionaire de séquence en entrée/sortie"""
        seqDict = {}
        header = None



### readlines()
### convert to list
### convert to dict


        try: # bloc try pour gerer l'erreur d'ouverture de fichier
            with open(filename) as fasta:
                for line in fasta: # itère sur tt les lignes du fasta
                    if line[0] == '>':  # si début de ligne
                        if header != None: # si ce n'est pas la première itération
                            seqDict [header] = sequence # remplissage du dictionnaire
                        
                        header = line[1:].replace('\n', '').replace('\r', '') # stockage du header
                        sequence = '' # reinitialisation de la séquence
                    else:
                        sequence += line[0:].replace('\n', '').replace('\r', '') # Obtention itérative de la séquence
                
                seqDict [header] = sequence     # ajout de la dernière séquence trouvée
                return seqDict

        except IOError:
            print '\n', filename, 'is not readable. The file will be ignored\n'
            return None


    def _calculate_proba(self):
        """Return a 2 entries list / 1 = name of the sequence / 2 = cumulative frequency of the sequence"""
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


    def _random_seq (self):
        """Return a random sequence from seq_dict according the respective size of references"""

        # Define a pseudo-random decimal frequency
        rand_freq = random()
        
        # Attibute this frequency to a sequence from seq_dict based on proba_list
        for name, freq  in self.proba_list:
            if freq > rand_freq :
                return name



class Junction (Sequence):
    pass
