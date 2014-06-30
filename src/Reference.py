"""
@package    Reference
@brief      Allow to create and manipulate dictionaries of Biopython SeqRecord objects
representing either a Reference Genome (for instance mm10) or a Pool of junctions between 2 other
references intended to mimick integration of exogenous DNA sequences in a another reference DNA.
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
from random import random, randint, sample
import gzip

# Local packages
from Utilities import import_seq
from SlicePicker import SlicePickerSingle

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class Reference(object):
    """
    @class Reference
    @brief Super class containing basic functions to reduce code redondancy in ReferenceGenome and
    ReferenceJunctions. NOT intended to be instanciated.
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~FONDAMENTAL METHODS~~~~~~~#

    def __init__(self, name):
        """
        Print a status message and initialize the name of the Reference object
        @param name Name of the reference
        """
        print "Initialisation of {}...".format(name)
        # Store object variables
        self.name = name

    def __repr__(self):
        result = "{}\n".format(self.__str__())
        for entry in self.d.values():
            result += "{}\nLenght:{}\n".format(entry, len(entry))
        return result

    def __str__(self):
        return "<Instance of {} from {} >".format(self.__class__.__name__, self.__module__)

    #~~~~~~~ACCESS METHODS~~~~~~~#

    def get(self, varkey):
        """
        @param varkey Name of a sequence contained in the SeqRecord dictionnary (string)
        @return A reference to the seqRecord coresponding to the Name required (*SeqRecord object)
        """
        return self.d[varkey]

    def getDict(self):
        return self.d

    def getName(self):
        return self.name

    def getLenDict(self):
        return len(self.d)

    def getSlice(self, refseq, start, end):
        """
        Return a slice from start to end in a given SeqRecord in a random orientation.
        The SeqRecord annotations dict is updated to contain relevant informations.
        @param refseq Name of the sequence (string)
        @param start Start position of the slice (int)
        @param end End position of the slice (int)
        @return A slice of the refseq at the given positions (SeqRecord object)
        """

        # Randomly choose an orientation reverse or forward and sample a slice
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

    #~~~~~~~PUBLIC METHODS~~~~~~~#

    def reset_samp_counter(self):
        """
        Initialize or reset counter for each Seqrecord that will be incremented each time
        _random_slice will sample a slice in this Seqrecord
        """
        for name, record in self.d.items():
            record.annotations ["nb_samp"] = 0

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class ReferenceGenome(Reference):
    """
    @class ReferenceGenome
    @brief **Import reference sequences from a fasta files and store each of them in a dict**
    The class provide a method to get a random seqrecord slice from one of the reference sequences,
    proportionally to the size of all reference sequences. This class require the Third party
    package Biopython via Utilities library.
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~FONDAMENTAL METHODS~~~~~~~#

    def __init__(self, name, fasta_path):
        """
        Import reference sequences from fasta file in a dictionnary and create a list
        of probability to allow a random sampling in SeqRecords proportionally to their size
        @param name Name of the reference (string)
        @param fasta_path Path of the fasta or fasta.gz file containing sequences (string)
        """
        # Use the super class init method
        super(self.__class__, self).__init__(name)

        # Dictionnary of bioPython record created from fasta file
        self.d = import_seq(fasta_path, "dict", "fasta")

        # Initialize a counter for each reference that will be
        # incremented each time _random_slice choose this reference
        self.reset_samp_counter()

        # List cummulative probabilities of each sequence to be picked
        # calculated from to their respective size.
        self.proba_list = self._calculate_proba()


    #~~~~~~~ACCESS METHODS~~~~~~~#

    def getProba(self):
        return self.proba_list

    #~~~~~~~PUBLIC METHODS~~~~~~~#

    def slicer (self, size):
        """
        Generate a random candidate slice in a random Seqrecord from the dictionnary. If not
        SeqRecord of a sufficient size was found after 100 tries, a generic Exception will be
        raised
        @param size Size of the slice to sample (int)
        @return Random slice of a valid size (SeqRecord object)
        """
        # Guard condition
        for count in range(100):
            # Pick a random sequence in dictionnary
            refseq = self._random_refseq()

            # If the size is valid return a slice
            if size <= len(self.d[refseq]):
                return self._random_slice(refseq, size)

        # if no valid size was found
        raise Exception("ERROR. Unable to find a slice of a valid length after 100 tries.\n\
        Please review the size of your references or the size of the reads.")

    def samp_report (self):
        """
        Create a simple list report containing the name, size and number of time a slice was sample
        into each SeqRecord in the dictionnary.
        @return A sampling report for each SeqRecord as a list (list)
        """
        # Add column header
        samp_list = [["chr", "lenght", "nb_samp"]]

        # Create a sorted list of refseq
        ref_list = self.d.keys()
        ref_list.sort()

        # Add values for each reference
        for ref in ref_list:
            samp_list.append([self.d[ref].id, len(self.d[ref]), self.d[ref].annotations["nb_samp"]])

        return samp_list

    #~~~~~~~PRIVATE METHODS~~~~~~~#

    def _calculate_proba(self):
        """
        Create a 2 entries list with the name of the sequence and the cumulative lenght of
        the sequence normalize to 1 (same as a frequency).
        @return A list of cumulative frequencies for all SeqRecord in the dictionnary (list)
        """
        proba_list = []
        cumulative_len = 0.0
        total_len = sum([len(record) for record in self.d.values()])

        for record in self.d.values():
            cumulative_len += len(record)
            proba_list.append([record.name, cumulative_len/total_len])

        return proba_list

    def _random_refseq(self):
        """
        Select a random sequence from the RefSeq dictionary according the cumulative probability
        list. A short sequence have a lower chance to be selected that a long one.
        @return ID of a SeqRecord (string)

        """
        # Define a pseudo-random decimal frequency
        rand_freq = random()

        # Attibute this frequency to a sequence from d based on proba_list
        for name, freq  in self.proba_list:
            if freq > rand_freq :
                return name

    def _random_slice(self, refseq, size):
        """
        Return a random slice in a given SeqRecord reference sequence in a random orientation.
        The SeqRecord annotations dict is updated to contain relevant informations.
        @param refseq Name of the sequence (string)
        @param size Size of the slice to sample (int)
        @return A random slice in refseq (SeqRecord object)
        """
        # Randomly choose the slice start position
        start = randint(0, len(self.d[refseq])-size)
        end = start+size

        return self.getSlice(refseq, start, end)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

class ReferenceJunctions(Reference):
    """
    @class ReferenceJunctions
    @brief Create a dictionnary of sequences by merging 2 SeqRecord slices from 2 Reference object
    for n junctions. The class provide a method to get a random seqrecord slice from one of the
    SeqRecord in the junction dictionnary. These slices overlap the middle of the chosen junction
    with a possible asymmetry given by the min_chimeric parameter.
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~FONDAMENTAL METHODS~~~~~~~#

    def __init__(self, name, min_chimeric, size, njun, ref1, ref2, repeats, ambiguous):
        """
        Create a dictionnary of junctions by merging SeqRecords from 2 Reference objects and
        initialize other class parameters.
        @param name Name of the reference (string)
        @param min_chimeric Minimal number of pb from one of the 2 references in a read or read
        pair overlaping a junction (int)
        @param size Size of the slice to sample in each reference to constitute a junction
        @param njun Number of junctions to generate
        @param ref1 Reference to the first source Reference
        @param ref2 Reference to the second source Reference
        @param repeats Allow lowercase characters (repeats) in the junctions
        @param ambiguous Allow Ambiguous DNA bases in the junctions
        """
        # Use the super class init method
        super(self.__class__, self).__init__(name)

        # Store object variables
        self.min_chimeric = min_chimeric
        self.half_len = size

        # Initialise junction dictionnary
        self.d = self._create_junctions_dict(size, njun, ref1, ref2, repeats, ambiguous)
        # Reset the sampling counter of references
        ref1.reset_samp_counter()
        ref2.reset_samp_counter()

        # Initialize a counter for each reference that will be
        # incremented each time _random_slice choose this reference
        self.reset_samp_counter()


###    PUBLIC METHODS    ####

    def slicer (self, size):
        """
        Generate a random candidate slice ovelaping a junction in a random Seqrecord from the
        dictionnary. If the require size is too short or too long to be sample, a generic
        Exception will be raised
        @param size Size of the slice to sample (int)
        @return Random slice of a valid size (SeqRecord object)

        Generate a slice overlapping a junction and return it
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

    def samp_report (self):
        """
        Create a simple list report containing the name, origin and number of time a slice was
        sample into each SeqRecord in the dictionnary.
        @return A sampling report for each SeqRecord as a list (list)
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
        """
        Return a string describing the the origin of a sequence correponding to a slice of a
        SeqRecord junction in the dictionnary.
        @param refseq Name of the sequence (string)
        @param start Start position of the slice (int)
        @param end End position of the slice (int)
        @return A string detailing the origin of a sequence (string)
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

    def _create_junctions_dict (self, size, njun, ref1, ref2, repeats, ambiguous):
        """
        Fill a dictionary with junctions between sequences from ref1 and ref2. An annotation
        dictionnary associated with each junction allows to retrieve easily its origin.
        @param size Size of the slice to sample in each reference to constitute a junction
        @param njun Number of junctions to generate
        @param ref1 Reference to the first source Reference
        @param ref2 Reference to the second source Reference
        @param repeats Allow lowercase characters (repeats) in the junctions
        @param ambiguous Allow Ambiguous DNA bases in the junctions
        @return A dictionnary containing junctions as SeqRecord objects
        """
        print("\tCreating junctions between {} and {}".format(ref1.getName(), ref2.getName()))

        slicer = SlicePickerSingle(size, repeats, ambiguous)
        junctions_dict = {}

        # Calculate the number of digits in njunctions for to generate junctions ids
        max_len = len(str(njun))

        for i in range(njun):

            # Pick 1 slice in each reference
            s1 = slicer.pick_slice(ref1)
            s2 = slicer.pick_slice(ref2)

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
        """
        Return a slice overlapping a junction from a given SeqRecord reference sequence in a
        random orientation. The SeqRecord annotations dict is updated to contain relevant
        informations.
        @param refseq Name of the sequence (string)
        @param size Size of the slice to sample (int)
        @return A random slice in refseq (SeqRecord object)
        """
        # Randomly choose the slice start position in the autorized area
        start = randint((self.half_len + self.min_chimeric - size), (self.half_len - self.min_chimeric))
        end = start + size

        return self.getSlice(refseq, start, end)


    def _coord_to_str (self, refseq, refsource, start, end):
        """
        Return a string containg the coordinate on ref1 or ref2 from a junction
        @param refseq SeqReccord object containing a junction from the dict
        @param refsource Indicates if the source reference sequence is the left
        part ("ref1") or the right part ("ref2") of the junction
        @param start Start coordinate relative to the start of the reference
        @param end End coordinate relative to the start of the reference
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
