from random import randint

class QualGenerator:
    """Accessory class generating quality strings following a given pattern"""

########################################################################################################################
#   FONDAMENTAL METHODS
########################################################################################################################

    def __init__(self, len, Q0m, Q1m, Q2m, Q3m, Q4m, Q0r, Q1r, Q2r, Q3r, Q4r, Q0f, Q1f, Q2f, Q3f, Q4f):
        """Create a quality template for mean, min and max value"""

        l   = [[Q0m, Q0m - Q0r/2, Q0m + Q0r/2, int(Q0f * len)],
               [Q1m, Q1m - Q1r/2, Q1m + Q1r/2, int(Q1f * len)],
               [Q2m, Q2m - Q2r/2, Q2m + Q2r/2, int(Q2f * len)],
               [Q3m, Q3m - Q3r/2, Q3m + Q3r/2, int(Q3f * len)],
               [Q4m, Q4m - Q4r/2, Q4m + Q4r/2, int(Q4f * len)]]

        self.qual_pattern = []

        for area in range (4):
            start_mean  = l [area] [0]
            delta_mean  = int(l [area + 1] [0] - l [area] [0])
            start_min   = l [area] [1]
            delta_min   = int(l [area + 1] [1] - l [area] [1])
            start_max   = l [area] [2]
            delta_max   = int(l [area + 1] [2] - l [area] [2])
            len_area    = l [area + 1] [3]

            for pos in range(len_area):
                mean_qual =  int( start_mean + (float (delta_mean)/ (len_area)*pos))
                min_qual  =  int( start_min  + (float (delta_min) / (len_area)*pos))
                max_qual  =  int( start_max  + (float (delta_max) / (len_area)*pos))
                self.qual_pattern.append([mean_qual, min_qual, max_qual])

    def __repr__(self):
        """Long representation"""
        return self.__str__()

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

    def random_qual_string (self, length):

        # Chose a random first score within ranges
        qual_string = [randint (self.qual_pattern[0][1], self.qual_pattern[0][2])]


        for i in range (1, length):

            # Evolution of the mean quality between the previous position and the current
            qual_evol = self.qual_pattern[i][0] - self.qual_pattern[i-1][0]
            # Half of the range of variation
            qual_var = int((self.qual_pattern[i][2] - self.qual_pattern[i][1]) / 2)
            # Score for the current position based on the previous position + the evolution of mean + a fluctuation in the range
            score = qual_string [i-1] + qual_evol + randint (- qual_var, qual_var)
            ############################################# Introduce a lot of variation = gaussian or beta law instead ?

            # If the score is under the minimun append the minimum
            if score < self.qual_pattern[i][1]:
                qual_string.append(self.qual_pattern[i][1])
                continue

            # If the score is above the maximun append the maximum
            if score > self.qual_pattern[i][2]:
                qual_string.append(self.qual_pattern[i][2])
                continue

            # Else add the calculated score
            else:
                qual_string.append(score)

        return qual_string

########################################################################################################################
#   PRIVATE METHODS
########################################################################################################################
