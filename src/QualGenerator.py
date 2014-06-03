# Standard library packages
from random import gauss

####################################################################################################

class QualGenerator(object):
    """Accessory class generating quality strings following a given
    pattern"""

###    FONDAMENTAL METHODS    ###

    def __init__(self, length, quality):
        """Create the class by asigning a lenght and a pattern of
        quality mean, sd variation for all positions along a given len
        """
        self.length = length
        self.qual_pattern = self._quality_pattern(length, quality)

    def __repr__(self):
        return "{}\n Mean qual pattern :\n{}\nStandard dev pattern :\n{}".format(
            self.__str__(),
            "".join([str(i[0])+" " for i in self.qual_pattern]),
            "".join([str(i[1])+" " for i in self.qual_pattern]))

    def __str__(self):
        return "<Instance of {} from {} >".format(self.__class__.__name__, self.__module__)

###    GETERS    ###

    def get_qual_pattern(self):
        return self.qual_pattern

    def get_mean_pattern(self):
        return [i[0] for i in self.qual_pattern]

    def get_sd_pattern(self):
        return [i[1] for i in self.qual_pattern]

    def get_length(self):
        return self.length

###    PUBLIC METHODS    ####

    def qual_score(self):
        """Generate a quality string mimicking a real fastq qual track
        Based on the pattern defined an class instanciation.
        """
        # Create an empty list to store scores
        qual_string = []

        for i in range(self.length):

            # Calculation of the mean qual score evolution based on the
            # previous score and an attraction factor toward the mean
            # value of quality in qual_pattern. In the case of the first
            # round, the score is the first mean value for the pattern
            if i == 0:
                qual_mean = self.qual_pattern[i][0]
            else:
                # Factor attracting values toward the pattern mean
                # Lower divider = higher attraction
                mean_attraction_factor = int((self.qual_pattern[i][0] - qual_string[i-1])/2)
                qual_mean = qual_string[i-1] + mean_attraction_factor

             # Adding a new valid score to the qual string list based on
             # the score mean a stochastic finite variation
            qual_string.append(self._valid_score(qual_mean, self.qual_pattern[i][1]))

        return qual_string

###    PRIVATE METHODS    ###

    def _valid_score(self, qual_mean, sigma): # sigma >= 1
        """Calculate a new valid score to based on the mean score and a
        stochastic gaussian fluctuation following standard deviation.
        """
        score = qual_mean + int(gauss(0, sigma))
        # if the score is outside of Phred quality borders (0 to 40).
        # the part coresponding to the gaussian fluctuation will not be
        # took into acount
        return score if 0 <= score <= 40 else qual_mean

    def _quality_pattern(self, length, quality):
        """Define a quality pattern for every positions in a given
        lenght based on quality parameters depending of the user defined
        quality range.
        """
        qual_param = self._quality_parameters(length, quality)

        # Create a smoothed pattern of quality mean and sd for each
        # positions in all areas of the read.
        qual_pattern = []

        # For all 4 areas defined
        for i in range(4):
            start_mean = qual_param[0][i]
            end_mean = qual_param[0][i+1]
            start_sd = qual_param[1][i]
            end_sd = qual_param[1][i+1]
            len_area = qual_param[2][i+1]

            # For each j position of the i area
            for j in range(len_area):
                # Calculate the smoothed mean and sd at the cur position
                mean =  int(start_mean +(float(end_mean - start_mean)/ len_area * j))
                sd =  int(start_sd +(float(end_sd - start_sd)/ len_area * j))
                # Add mean and sd to the list
                qual_pattern.append([mean, sd])

        return qual_pattern

    def _quality_parameters(self, length, quality):
        """Return quality parameters at borders of 4 areas within the
        lenght (mean, sd and length).
        """
        # Calculate the lenght of each quality area in the read
        a1 = int(0.1*length)
        a2 = int(0.2*length)
        a3 = int(0.3*length)
        a4 = int(0.4*length)

        # Verify partitioning of the size and if needed correct it
        sum_len = a1 + a2 + a3 + a4
        if sum_len != length:
            print "Invalid integer partitioning: {} instead of {}".format(sum_len, length)
            print "a1 = {}, a2 = {}, a3 = {}, a4 = {}\n".format(a1,a2,a3,a4)
            a4 += (length - sum_len)
            print "Correction of areas lenghts"
            print "a1 = {}, a2 = {}, a3 = {}, a4 = {}\n".format(a1,a2,a3,a4)

        # Dictionnary storing mean sd and lenght of each area boundaries
        # (start, 3 intermediate points and end of a model read for
        # different quality categories
        qual_param_dict = {
        "very-good": [[33, 38, 39, 38, 36], [2, 1, 1, 1, 5], [0, a1, a2, a3, a4]],
        "good": [[32, 37, 38, 37, 32], [2, 2, 1, 3, 8], [0, a1, a2, a3, a4]],
        "medium": [[30, 36, 37, 32, 25], [4, 3, 2, 6, 8], [0, a1, a2, a3, a4]],
        "bad": [[20, 25, 27, 23, 15], [4, 3, 2, 7, 8], [0, a1, a2, a3, a4]],
        "very-bad": [[10, 15, 15, 10, 5], [4, 3, 2, 8, 8], [0, a1, a2, a3, a4]]}

        return qual_param_dict[quality]
