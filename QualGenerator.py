from random import randint
from random import gauss

class QualGenerator:
    """Accessory class generating quality strings following a given pattern"""

########################################################################################################################
#   FONDAMENTAL METHODS
########################################################################################################################

    def __init__(self, length, quality):
        """Create the class by asigning a lenght and a pattern of quality mean, sd variation along the read"""
    
        self.length = length
        self.qual_pattern = self._quality_pattern (length, quality)
        
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

    def __repr__(self):
        """Long representation"""
        return self.__str__()
        
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

    def __str__(self):
        """Short representation"""
        return "<Instance of " + self.__module__ + ">\n"

########################################################################################################################
#   GETERS
########################################################################################################################

    def get_variable (self):
        """return variable"""
        return self.variable


########################################################################################################################
#   PUBLIC METHODS
########################################################################################################################

    def random_qual_string (self):

        # Create an empty list to store scores 
        qual_string = []

        for i in range (self.length):
            
            # Calculation of the mean qual score evolution based on the previous score + an attraction factor toward
            # the mean value of quality in qual_pattern (a quarter of the interval seems realistic) 
            # In the case of the first loop round, the score is the mean value for the pattern at position 0
            if i == 0:
                qual_mean = self.qual_pattern[i][0]
            else: 
                qual_mean = qual_string[i-1] + int((self.qual_pattern[i][0] - qual_string[i-1])/2)
            
             # Adding a new valid score to the qual string list based on the score mean a stochastic finite variation
            qual_string.append (self._valid_score (qual_mean, self.qual_pattern[i][1]))

        return qual_string

########################################################################################################################
#   PRIVATE METHODS
########################################################################################################################

    def _valid_score(self, qual_mean, sigma): # sigma >= 1
        """Calculate a new valid score to based on the mean score + a stochastic gaussian fluctuation following sigma"""
        
        # Generation of a score including a variation following a gaussian distribution
        # except if the value if lower than 0 or higher than 40
        score = qual_mean + int (gauss (0, sigma))
        
        return score if 0 <= score <= 40 else qual_mean
        
        
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

    def _quality_pattern (self, length, quality):
        
        qual_param = self._quality_parameters (length, quality)
        
        # Create a smoothed pattern of quality mean and sd for each positions in all areas of the read.
        qual_pattern = []
        
        # For all 4 areas defined in the 
        for i in range (4):
            start_mean  = qual_param[0][i]
            end_mean    = qual_param[0][i+1]
            start_sd    = qual_param[1][i]
            end_sd      = qual_param[1][i+1]
            len_area    = qual_param[2][i+1]
            
            # For each j position of the i area
            for j in range (len_area):
                # Calculate the smoothed mean and sd at the current position
                mean =  int (start_mean + (float (end_mean - start_mean)/ len_area * j))
                sd =  int (start_sd + (float (end_sd - start_sd)/ len_area * j))
                # Add mean and sd to the list
                qual_pattern.append([mean,sd])
        
        return qual_pattern


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

    def _quality_parameters (self, length, quality):
        
        # Calculate the lenght of each quality area in the read 
        a1 = int(0.1*length)
        a2 = int(0.2*length)
        a3 = int(0.3*length)
        a4 = int(0.4*length)
        
        # Verify partitioning of the size and if needed correct it
        sum_len = a1 + a2 + a3 + a4
        if sum_len != length:
            print "Invalid integer partitioning. Correction of the area lenght" 
            a3 += length - sum_area
        
        # Dictionnary storing mean and sd at start, 3 intermediate points and end of a model read
        qual_param_dict = {
        "very_good" : [[33,38,39,38,36], [2,1,1,1,5],  [0,a1,a2,a3,a4]],
        "good"      : [[32,37,38,37,32], [2,2,1,3,10],  [0,a1,a2,a3,a4]],
        "medium"    : [[30,36,37,32,25], [4,3,2,10,10], [0,a1,a2,a3,a4]],
        "bad"       : [[25,30,32,28,20], [4,3,2,10,10], [0,a1,a2,a3,a4]],
        "very_bad"  : [[15,20,20,15,5],  [4,3,2,10,10], [0,a1,a2,a3,a4]]}
        
        return qual_param_dict[quality]
