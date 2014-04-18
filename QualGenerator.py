from random import gauss

class QualGenerator:
    """Accessory class generating quality strings following a given pattern"""

########################################################################################################################
#   FONDAMENTAL METHODS
########################################################################################################################
    
    def __init__(self, len, Q0m, Q1m, Q2m, Q3m, Q4m, Q0r, Q1r, Q2r, Q3r, Q4r, Q0f, Q1f, Q2f, Q3f, Q4f): 
        """Create a quality template for mean value and  """
        
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
        
# TODO POUR ATTRIBUER string faire variation qui suit le pas du pattern +1  0 -1 en suivant un loi normale (ou beta)


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
        
    def get_variable (self):
        """return variable"""
        return self.variable

########################################################################################################################
#   PUBLIC METHODS
########################################################################################################################

    def random_qual_string (self, mean_qual, sd_qual, length):
        return [valid_qual(mean_qual[i], sd_qual[i]) for i in range (length)]

########################################################################################################################
#   PRIVATE METHODS
########################################################################################################################

    def _qual_template (self, qual_start, qual_mid, qual_end, qvar_start, qvar_mid, qvar_end, length):
        """Generate a mean quality string and sd quality string that will be used as as a template qual strings"""
        # Define boundaries dividing lenght in  parts start, middle and end
        # These 3 parts ressemble to illumina quality strings starting by a short area of quality increase followed
        # by a intermediate area of good and stable quality and finally an area of quality decrease
        # transition beetween areas are progressive.
        start_border = int(0.2*length)  # First 20% of read lenght
        mid_border = int(0.66*length)   # From 20 to 66 percent
        end_border = length             # Until the end
        
        # Empty list to store quality means and quality 
        self.mean = []
        self.sd =[]
        
        for i in range (start_border) :
            mean_qual.append    (int(m_start + (float (m_mid - m_start))/start_border*i))
            sd_qual.append      (int(sd_start + (float (sd_mid - sd_start))/start_border*i))

        for i in range (mid_border - start_border) :
            mean_qual.append(m_mid)
            sd_qual.append(sd_mid)

        for i in range(end_border - mid_border):
            mean_qual.append    (int(m_mid + (float (m_end - m_mid)/ (end_border - mid_border)*i)))
            sd_qual.append      (int(sd_mid + (float (sd_end - sd_mid)/ (end_border - mid_border)*i)))

        return (mean_qual, sd_qual)

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

    def _valid_qual(self):
        while True:
            qual = int(gauss(mean, sd))
            if 0 <= qual <= 40:
                return qual

