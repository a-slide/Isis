"""
@package    ReferenceJunctions
@brief
@copyright  [GNU General Public License v2](http://www.gnu.org/licenses/gpl-2.0.html)
@author     Adrien Leger - 2014
* <adrien.leger@gmail.com>
* <adrien.leger@inserm.fr>
* <adrien.leger@univ-nantes.fr>
* [Github](https://github.com/a-slide)
* [Atlantic Gene Therapies - INSERM 1089] (http://www.atlantic-gene-therapies.fr/)
"""

#~~~~~~~PACKAGE IMPORTS~~~~~~~#

# Standard library packages
from random import sample, randint

# Local packages
from SlicePicker import SlicePickerSingle

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

class ReferenceJunctions(object):
    """ Create de dictionnary of sequence by merging from 2 references
    and by creating a chimeric junction between both for n junctions
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~FONDAMENTAL METHODS~~~~~~~#

    def __init__(self, name, min_chimeric, size, njun, ref1, ref2, repeats, ambiguous):
        """Create a dictionnary of junctions by merging 2 SeqRecords
        """
        print "Initialisation of {}...".format(name)

        # Store object variables
        self.name = name
        self.min_chimeric = min_chimeric
        self.half_len = size

        # Initialise junction dictionnary
        self.d = self._create_junctions_dict(size, njun, ref1, ref2, repeats, ambiguous)

        # Initialize a counter for each reference that will be
        # incremented each time _random_slice choose this reference
        self.reset_samp_counter()

    def __repr__(self):
        result = "{}\n".format(self.__str__())
        for entry in self.d.values():
            result += "{}\nLenght:{}\n".format(entry, len(entry))
        return result

    def __str__(self):
        return "<Instance of {} from {} >".format(self.__class__.__name__, self.__module__)

    #~~~~~~~ACCESS METHODS~~~~~~~#

    def get(self, varkey):
        return self.d[varkey]

    def getDict(self):
        return self.d

    def getName(self):
        return self.name

    def getLenDict(self):
        return len(self.d)

###    PUBLIC METHODS    ####

    def get_slice (self, size):
        """Generate a slice overlapping a junction and return it
        """
        # Guard conditions
        if  size < 2 * self.min_chimeric:
            raise Exception ("ERROR. The size of the slice is too short to define a fragment\n\
            with the require minimal number of chimeric bases in each reference sequence.\n\
            Review your size of min chimeric base or read length")
        if size > self.half_len:
            raise Exception ("ERROR. The size of the slice is longer than the reference.\n\
            Review your maximal sonication size or the size of junctions")

        # Pick a random reference junction and return a slice of it
        refseq = sample(self.d, 1)[0]
        return self._random_slice (refseq, size)

    def reset_samp_counter(self):
        """ Initialize or reset counter for each reference that will be
        incremented each time _random_slice choose this reference
        """
        for name, record in self.d.items():
            record.annotations ["nb_samp"] = 0

    def samp_report (self):
        """ Create a simple list report
        """
        # Add column header
        samp_list = [["junction_id", "ref1_source", "ref1_chr", "ref1_loc",
                "ref1_orientation", "ref2_source", "ref2_chr", "ref2_loc",
                "ref2_orientation", "nb_samp"]]

        # Create a sorted list of refseq
        ref_list = self.d.keys()
        ref_list.sort()

        # Add values for each reference
        for ref in ref_list:
            samp_list.append([
                self.d[ref].id,
                self.d[ref].annotations["ref1_source"].getName(),
                self.d[ref].annotations["ref1_refseq"].id,
                self.d[ref].annotations["ref1_location"],
                self.d[ref].annotations["ref1_orientation"],
                self.d[ref].annotations["ref2_source"].getName(),
                self.d[ref].annotations["ref2_refseq"].id,
                self.d[ref].annotations["ref2_location"],
                self.d[ref].annotations["ref2_orientation"],
                self.d[ref].annotations["nb_samp"]])

        return samp_list

    def origin_coord (self, refseq, start, end):
        """ Return a string describing the the origin of a sequence from a
        junction from the dictionnary
        """
        # If the give coord overlap only the left reference of a junction
        if end < self.half_len:
            return ("{}-{}={}".format(
            1, end-start,
            self._coord_to_str(refseq, "ref1", start, end+1)))

        # If the give coord overlap only the right reference of a junction
        elif start >= self.half_len:
            return ("{}-{}={}".format(
            1, end-start,
            self._coord_to_str(refseq, "ref2", start-self.half_len, end-self.half_len+1)))

        # If the give coord overlap both references of a junction
        else :
            return ("{}-{}={}|{}-{}={}".format(
            1,self.half_len-start,
            self._coord_to_str(refseq, "ref1", start, self.half_len),
            self.half_len-start+1, end-start,
            self._coord_to_str(refseq, "ref2", 0, end-self.half_len+1)))

    #~~~~~~~PRIVATE METHODS~~~~~~~#

    def _create_junctions_dict (self, size, njunctions, ref1, ref2, repeats, ambiguous):
        """Fill a dictionary with junctions between sequences from ref1
        and ref2
        """

        print("\tCreating junctions between {} and {}".format(ref1.getName(), ref2.getName()))

        slicer = SlicePickerSingle(size, repeats, ambiguous)
        junctions_dict = {}

        # Calculate the number of digits in njunctions for to generate junctions ids
        max_len = len(str(njunctions))

        for i in range(njunctions):

            # Pick 1 slice in each reference
            s1 = slicer.pick_slice(ref1)
            s2 = slicer.pick_slice(ref2)
            # Reset the sampling counter of references
            ref1.reset_samp_counter()
            ref2.reset_samp_counter()

            # Create a new junction from the 2 slice
            junction = s1 + s2

            # Define id, name and description to the SeqReccord object
            junction.name = junction.description = ""
            junction.id = "Junction_{0:0{1}}".format(i, max_len)

            # Add informations to the annotations dictionnary
            junction.annotations = {
                "ref1_source" : ref1,
                "ref2_source" : ref2,
                "ref1_refseq" : s1.annotations["refseq"],
                "ref2_refseq" : s2.annotations["refseq"],
                "ref1_location" : s1.annotations["location"],
                "ref2_location" : s2.annotations["location"],
                "ref1_orientation" : s1.annotations["orientation"],
                "ref2_orientation" : s2.annotations["orientation"]}

            junctions_dict[junction.id] = junction

        return junctions_dict

    def _random_slice (self, refseq, size):
        """Return a slice overlapping a junction from a biopython Seqrecord in d
        """
        # Randomly choose the slice start position in the autorized area
        start = randint((self.half_len + self.min_chimeric - size), (self.half_len - self.min_chimeric))
        end = start + size

        # Randomly choose an orientation reverse or forward and sample
        # a slice
        if randint(0,1):
            s = self.d[refseq][start:end]
            s.annotations["orientation"] = "+"
        else:
            s = self.d[refseq][start:end].reverse_complement()
            s.annotations["orientation"] = "-"

        # Add informations to the annotations dictionnary
        s.annotations["refseq"] = self.d[refseq]
        s.annotations["location"] = [start, end]
        s.annotations["source"] = self

        # Undefine description, id and name
        s.name = s.description = s.id = ""

        # Increment the sampling counter of the refseq
        self.d[refseq].annotations["nb_samp"] += 1

        return s


    def _coord_to_str (self, refseq, refsource, start, end):
        """ Return the coordinate on ref1 or ref2 from a junction
        @param refseq SeqReccord object cointaining a junction from the class
        junction dict
        @param refsource Indicates if the source reference sequence is the left
        part ("ref1") or the right part ("ref2") of the junction
        @param start Start coordinate relative to the start of the reference
        @param start End coordinate relative to the start of the reference
        @return A string composed of the id of the source reference followed by
        start and end coordinates on the source reference (ex : chr12:124-434 )
        """

        ref_id = refseq.annotations[refsource+"_refseq"].id

        if refseq.annotations[refsource+"_orientation"] == '+':
            ref_start = refseq.annotations[refsource+"_location"][0] + start
            ref_end = refseq.annotations[refsource+"_location"][0] + end

        else:
            ref_start = refseq.annotations[refsource+"_location"][1] - end
            ref_end = refseq.annotations[refsource+"_location"][1] - start

        return ("{}:{}-{}".format(ref_id, ref_start, ref_end))
