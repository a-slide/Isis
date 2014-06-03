#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
@package    IsisMain
@brief      Main package of Isis program. Store configuration parameters,
instanciate reference genome and junctions, ask fastq sequences and write them
in compressed file(s). No function parameters are required but several command
line arguments are needed as explained in the documentation or with -h option
* -V virus/vector genome fasta file path
* -H host genome fasta file path
* -C Configuration file path
* (-o outupt name)
@copyright  [GNU General Public License v2](http://www.gnu.org/licenses/gpl-2.0.html)
@author     Adrien Leger <adrien.leger@gmail.com>
"""

#~~~~~~~GLOBAL IMPORTS~~~~~~~#

# Standard library packages
import gzip

# Local packages
from IsisConf import IsisConf, IsisConfException
from ReferenceGenome import ReferenceGenome as RefGen
from ReferenceJunctions import ReferenceJunctions as RefJun
from SlicePicker import SlicePickerSingle, SlicePickerPair
from QualGenerator import QualGenerator
from FastqGenerator import FastqGeneratorSingle, FastqGeneratorPair


#~~~~~~~FUNCTIONS FOR SINGLE END MODE~~~~~~~#

def IsisSingle():
    """
    Specific for single end reads. Instantiate accessory classes and call the
    single end fastq writer. Need global parameters to be executed correctly.
    * READ_LEN  Lenght of the reads to be gerenated (integer)
    * REPEATS   Allow lowercase characters in the DNA sequence (bool)
    * AMBIGUOUS Allow IUPAC ambigous characters in the DNA sequence (bool)
    * MUT_FREQ  Frequency of single nucleotide mutation to be randomly introduced in reads (float)
    * QUAL_RANGE Range of quality for reads (string)
    """
    # Instantiate accessory classes
    slicer = SlicePickerSingle (READ_LEN, REPEATS, AMBIGUOUS, MUT_FREQ)
    qualgen = QualGenerator (READ_LEN, QUAL_RANGE)
    fastgen = FastqGeneratorSingle (slicer, qualgen)

    # Write fastq file for all source in source_list
    write_fastq_single (fastgen)

    return 1

def write_fastq_single (fastgen):
    """ Function printing the request number of single end reads per reference
    in a fastq.gz file
    """
    try:
        f = gzip.open (BASENAME + ".fastq.gz", 'w')

        for source, nread in SOURCE_LIST:
            print ("\tWritting {} read(s) in Fastq file from {}".format (nread, source.getName()))
            # Calculate the number of digits in nread for read naming
            max_len = len (str (nread))

            for i in range (nread):
                # Ask a read to the source throught fastgen
                read = fastgen.generate_fastq (source)
                # Generate a uniq identifier
                read.id = generate_id (i, max_len, read.annotations)
                # Write the fastq formated read
                f.write (read.format (QUAL_SCALE))

                # Add read coverage over junction to JUN_COV list if the source is a junction
                if GRAPH and isinstance (source, RefJun):
                    update_jun_cov (read.annotations["location"][0], read.annotations["location"][1])
        f.close()

    except IOError as E:
        print (E)
        exit (0)


#~~~~~~~FUNCTIONS FOR PAIR END MODE~~~~~~~#

def IsisPair():
    """
    """
    # Instantiate accessory classes
    slicer = SlicePickerPair (READ_LEN, SONIC_MIN, SONIC_MODE, SONIC_MAX, SONIC_CERTAINTY, REPEATS, AMBIGUOUS, MUT_FREQ)
    qualgen = QualGenerator (READ_LEN, QUAL_RANGE)
    fastgen = FastqGeneratorPair (slicer, qualgen)

    # Write paired fastq files for all source in source_list
    write_fastq_pair (fastgen)

    return 1

def write_fastq_pair (fastgen):
    """ Function printing the request number of pair end reads per reference
    in a fastq.gz file
    """
    try:
        f1 = gzip.open(BASENAME + ".R1.fastq.gz", 'w')
        f2 = gzip.open(BASENAME + ".R2.fastq.gz", 'w')

        for source, nread in SOURCE_LIST:
            print ("\tWritting {} read(s) in Fastq file from {}".format(nread, source.getName()))
            # Calculate the number of digits in nread for read naming
            max_len = len(str(nread))

            for i in range (nread):
                # Ask a read to the source throught fastgen
                read1, read2 = fastgen.generate_fastq(source)
                # Uniq read identifier
                read1.id = generate_id(i, max_len, read1.annotations)
                read2.id = generate_id(i, max_len, read2.annotations)
                # Write the fastq formated read
                f1.write(read1.format(QUAL_SCALE))
                f2.write(read2.format(QUAL_SCALE))

                # Add "sonication" fragment lenght to the FRAG_LEN list
                if GRAPH:
                    update_frag_len (read1.annotations["frag_len"])
                    # Add read coverage over junction to JUN_COV list if the source is a junction
                    if isinstance (source, RefJun):
                        update_jun_cov (read1.annotations["location"][0], read1.annotations["location"][1])
                        update_jun_cov (read2.annotations["location"][0], read2.annotations["location"][1])
        f1.close()
        f2.close()

    except IOError as E:
        print (E)
        exit (0)

#~~~~~~~HELPER FUNCTIONS~~~~~~~#

def generate_id(i,max_len, d):
    """ Generate an uniq identifier indicating from where the read was sampled
    """
    # Add source name
    id_string = "{}".format(d["source"].getName())
    # Add an uniq numeric identifier per source with zero padding
    id_string += "|{0:0{1}}".format(i, max_len)

    # If the source is a junction coordinate along original references are asked
    # to the source reference junction in which the read was sampled
    if isinstance (d["source"], RefJun):
        id_string += "|1|{}".format(d["source"].origin_coord(d["refseq"], d["location"][0], d["location"][1]))
    # Just enter strored coordinates if the reference is a ReferenceGenome obj
    else:
        id_string += "|0|All={}:{}-{}".format(d["refseq"].id, d["location"][0], d["location"][1])

    return id_string


#~~~~~~~FUNCTIONS FOR GRAPHICAL OUTPUT~~~~~~~#

def update_jun_cov (start, end):
    for i in range (start, end):
        JUN_COV [i] += 1

def update_frag_len (frag_len):
    FRAG_LEN.append(frag_len)

def frag_len_graph ():
    """ Output a graphical representation of the fragment len distribution
    """
    print ("\tCreating a graphical output of fragment length distribution")

    # Create a figure object and adding details
    fig = plt.figure(figsize=(15, 10), dpi=100)
    plt.title("Distribution of fragment length")
    plt.ylabel('Relative Count')
    plt.xlabel('Size of fragment')

    # Plot value from FRAG_LEN list in an histogram reprensentation
    plt.hist(FRAG_LEN, bins=(SONIC_MAX-SONIC_MIN)/5, range=(SONIC_MIN, SONIC_MAX),
    normed=1, facecolor='green', alpha=0.5, align='mid')

    # Tweak spacing to prevent clipping of ylabel
    plt.subplots_adjust(left=0.15)

    # Export figure to file
    fig.savefig(BASENAME+'_distribution.png')

def jun_cov_graph ():
    """ Output a graphical representation of read coverage over junction
    """
    print ("\tCreating a graphical output of read coverage over junctions")

    # Create a figure object and adding details
    fig = plt.figure(figsize=(15, 10), dpi=100)
    plt.title("Coverage of reads over junctions")
    plt.ylabel('Count')
    plt.xlabel('Position (mid = junction)')

    # List of numbers for x axis positions
    x = [i+1 for i in range (len(JUN_COV))]

    # Plot a vertical line indicating the position of the junction
    plt.axvline(len(JUN_COV)/2, color="red")
    # Plot an area representing the coverage depth
    plt.fill(x,JUN_COV, facecolor='green', alpha=0.5)

    # Tweak spacing to prevent clipping of ylabel
    plt.subplots_adjust(left=0.15)

    # Export figure to file
    fig.savefig(BASENAME+'_junction_coverage.png')

#~~~~~~~FUNCTIONS FOR REPORT WRITING~~~~~~~#

def write_samp_report ():
    """
    """
    with open(BASENAME+'_sampling_report.txt', 'w') as report:

        # Ask report list to all source
        for source, nread in SOURCE_LIST:
            # Write general source informations
            report.write("Source = {}\n".format(source.getName()))
            report.write("Total number of read = {}\n".format(nread))
            report.write("Number of sequences = {}\n".format(source.getLenDict()))

            # Write a csv like report for all items in the list returned by source
            for i in source.samp_report():
                line = ''
                for j in i:
                    line += '{}\t'.format(j)
                report.write(line+"\n")
            report.write("\n")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~MAIN FUNCTION~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

if __name__ == '__main__':

    PROGRAM_NAME = "Isis"
    PROGRAM_VERSION = "0.01"

    # Instantiate the configuration file parser/verifier
    CONFIG = IsisConf(PROGRAM_NAME, PROGRAM_VERSION)

    # Import configurations in local variables
    HOST_GENOME = CONFIG.get("host_genome")
    VIRUS_GENOME = CONFIG.get("virus_genome")
    BASENAME = CONFIG.get("basename")
    PAIR = CONFIG.get("pair")
    READ_NUM = CONFIG.get("read_num")
    READ_LEN = CONFIG.get("read_len")
    MUT_FREQ = CONFIG.get("mut_freq")
    REPEATS = CONFIG.get("repeats")
    AMBIGUOUS = CONFIG.get("ambiguous")
    GRAPH = CONFIG.get("graph")
    REPORT = CONFIG.get("report")
    N_HOST = CONFIG.get("nread_host")
    N_VIRUS = CONFIG.get("nread_virus")
    N_TJUN = CONFIG.get("nread_tjun")
    N_FJUN = CONFIG.get("nread_fjun")
    MIN_CHIMERIC = CONFIG.get("min_chimeric")
    UNIQ_TJUN = CONFIG.get("uniq_tjun")
    UNIQ_FJUN = CONFIG.get("uniq_fjun")
    SONIC_MIN = CONFIG.get("sonic_min")
    SONIC_MODE = CONFIG.get("sonic_mode")
    SONIC_MAX = CONFIG.get("sonic_max")
    SONIC_CERTAINTY = CONFIG.get("sonic_certainty")
    QUAL_SCALE = CONFIG.get("qual_scale")
    QUAL_RANGE = CONFIG.get("qual_range")
    JUN_LEN = SONIC_MAX if PAIR else READ_LEN

    # Import reference genomes
    VIRUS = RefGen ("virus", VIRUS_GENOME)
    HOST = RefGen ("host", HOST_GENOME)
    # Create reference junctions
    TJUN = RefJun("True_Junction", MIN_CHIMERIC, JUN_LEN, UNIQ_TJUN, VIRUS, HOST, REPEATS, AMBIGUOUS)
    FJUN = RefJun("False_Junction", MIN_CHIMERIC, JUN_LEN, UNIQ_FJUN, VIRUS, HOST, REPEATS, AMBIGUOUS)

    # Import tools and Initialize
    if GRAPH:
        # Third party package matplotlib imported only if needed
        from matplotlib import pyplot as plt
        # Init an empty list to store "sonication" fragment lenght
        FRAG_LEN = []
        # Init a zero filled list to store coverage over junctions
        JUN_COV = [0 for i in range (JUN_LEN*2)]
        # Define a convenient ref list to simplify subsequent operations. The 3rd
        # position bool if used for sources coverage graph will be output
        SOURCE_LIST = [[VIRUS, N_VIRUS], [HOST, N_HOST], [TJUN, N_TJUN], [FJUN, N_FJUN]]
        # Initialisation of the coverage list

    else :
        # Define a convenient ref list to simplify subsequent operations.
        SOURCE_LIST = [[VIRUS, N_VIRUS], [HOST, N_HOST], [TJUN, N_TJUN], [FJUN, N_FJUN]]

    # Lauch main fastq writing functions according to the mode selected
    if PAIR:
        # Pair end mode specific imports
        print ("Start pair end mode fastq sampling")
        IsisPair()
    else:
        # Single end mode specific imports
        print ("Start single end mode fastq sampling")
        IsisSingle()

    # Write sampling reports if requested
    if REPORT:
        write_samp_report ()

    # Draw graphs if requested
    if GRAPH:
        # Graph of read coverage over junctions
        jun_cov_graph ()
        if PAIR:
            # Graph of "sonication" fragment len distribution
            frag_len_graph ()

    print "DONE"
    exit (1)
