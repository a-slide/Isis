#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
@package    IsisMain
@brief      **Main package of Isis program**.
Store configuration parameters,instanciate reference genome and junctions, ask
for fastq sequences and write themin compressed file(s). No function parameters
are required but several command line arguments are needed as explained in the
documentation or with -h option :
* -V virus/vector genome fasta file path
* -H host genome fasta file path
* -C Configuration file path
* (-o outupt name)

@copyright  [GNU General Public License v2](http://www.gnu.org/licenses/gpl-2.0.html)
@author     Adrien Leger - 2014
* <adrien.leger@gmail.com>
* <adrien.leger@inserm.fr>
* <adrien.leger@univ-nantes.fr>
* [Github](https://github.com/a-slide)
* [Atlantic Gene Therapies - INSERM 1089] (http://www.atlantic-gene-therapies.fr/)
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
    * MUT_FREQ  Frequency of single nucleotide mutation to be randomly
    introduced in reads (float)
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
    """
    Write the requested number of single end reads per reference in a fastq.gz
    file. Need global parameters to be executed correctly.
    * BASENAME  basename for the name of the output file (string)
    * QUAL_SCALE Quality scale in which fastq quality string will be converted
    (sanger, illumina, solexa...) (string)
    * SOURCE_LIST   List containing source sublists with a reference to an
    instance of ReferenceGenome or ReferenceJunction and the requested number of
    reads to be generated (list)
    @param  fastgen Instance of FastqGenerator
    """
    try:
        f = gzip.open (BASENAME + ".fastq.gz", 'w')

        for source, nread in SOURCE_LIST:
            print ("\tWritting {} read(s) in Fastq file from {}".format (nread, source.getName()))
            # Calculate the number of digits in nread for read id
            id_len = len (str (nread))

            for i in range (nread):
                # Ask a read to the source throught fastgen
                read = fastgen.generate_fastq (source)
                # Generate a uniq identifier
                read.id = generate_id (i, id_len, read.annotations)
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
    Specific for pair end reads. Instantiate accessory classes and call the pair
    end fastq writer. Need global parameters to be executed correctly.
    * READ_LEN  Lenght of the reads to be generated (integer)
    * REPEATS   Allow lowercase characters in the DNA sequence (bool)
    * AMBIGUOUS Allow IUPAC ambigous characters in the DNA sequence (bool)
    * MUT_FREQ  Frequency of single nucleotide mutation to be randomly
    introduced in reads (float)
    * QUAL_RANGE Range of quality for reads (string)
    * SONIC_MIN Minimal size of sonication fragments (integer)
    * SONIC_MODE Modal size of sonication fragments (integer)
    * SONIC_MAX Maximal size of sonication fragments (integer)
    * SONIC_CERTAINTY Thickness of the sonication peak (integer)
    """
    # Instantiate accessory classes
    slicer = SlicePickerPair (READ_LEN, SONIC_MIN, SONIC_MODE, SONIC_MAX, SONIC_CERTAINTY, REPEATS, AMBIGUOUS, MUT_FREQ)
    qualgen = QualGenerator (READ_LEN, QUAL_RANGE)
    fastgen = FastqGeneratorPair (slicer, qualgen)

    # Write paired fastq files for all source in source_list
    write_fastq_pair (fastgen)

    return 1

def write_fastq_pair (fastgen):
    """
    Write the requested number of pair end reads per reference in a fastq.gz
    file. Need global parameters to be executed correctly.
    * BASENAME  basename for the name of the output file (string)
    * QUAL_SCALE Quality scale in which fastq quality string will be converted
    (sanger, illumina, solexa...) (string)
    * SOURCE_LIST   List containing source sublists with a reference to an
    instance of ReferenceGenome or ReferenceJunction and the requested number of
    reads to be generated (list)
    @param  fastgen Instance of FastqGenerator
    """
    try:
        f1 = gzip.open(BASENAME + ".R1.fastq.gz", 'w')
        f2 = gzip.open(BASENAME + ".R2.fastq.gz", 'w')

        for source, nread in SOURCE_LIST:
            print ("\tWritting {} read(s) in Fastq file from {}".format(nread, source.getName()))
            # Calculate the number of digits in nread for read id
            id_len = len(str(nread))

            for i in range (nread):
                # Ask a read to the source throught fastgen
                read1, read2 = fastgen.generate_fastq(source)
                # Uniq read identifier
                read1.id = generate_id(i, id_len, read1.annotations)
                read2.id = generate_id(i, id_len, read2.annotations)
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

def generate_id(i,id_len, d):
    """
    Generate an identifier indicating from where the read was sampled. The id
    contains the name of the source an uniq numeric identifier, a bolean set to
    True is the original sequence is a junction and a description of the origin
    sequence by read positions. Example: 1-75=Bact:330-405 means that bases
    position 1 to 75 from the read originates from positions 333-405 from bact
    reference sequences. Need global parameters to be executed correctly.
    * READ_LEN  Lenght of the reads to be generated (integer)
    @param  i   Num of the read (integer)
    @param  id_len Max number of digit
    @param  d   annotation dictionnary from a seq record object
    @return A descriptive and uniq string.
    """
    # Add source name
    id_string = "{}".format(d["source"].getName())
    # Add an uniq numeric identifier per source with zero padding
    id_string += "|{0:0{1}}".format(i, id_len)

    # If the source is a junction coordinate along original references are asked
    # to the source reference junction in which the read was sampled
    if isinstance (d["source"], RefJun):
        id_string += "|1|{}".format(d["source"].origin_coord(d["refseq"], d["location"][0], d["location"][1]))
    # Just enter strored coordinates if the reference is a ReferenceGenome obj
    else:
        id_string += "|0|1-{}={}:{}-{}".format(READ_LEN, d["refseq"].id, d["location"][0], d["location"][1])

    return id_string


#~~~~~~~FUNCTIONS FOR GRAPHICAL OUTPUT~~~~~~~#

def update_jun_cov (start, end):
    for i in range (start, end):
        JUN_COV [i] += 1

def update_frag_len (frag_len):
    FRAG_LEN.append(frag_len)

def frag_len_graph ():
    """
    Output a graphical representation of the fragment len distribution. Require
    the third party package pyplot from matplotlib. Need global parameters to be
    executed correctly.
    * BASENAME  basename for output files (string)
    * SONIC_MIN Minimal size of sonication fragments (integer)
    * SONIC_MAX Maximal size of sonication fragments (integer)
    * FRAG_LEN  List of size from sonication fragments (list)
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

    return 1

def jun_cov_graph ():
    """
    Output a graphical representation of read coverage over junction. Require
    the third party package pyplot from matplotlib. Need global parameters to be
    executed correctly.
    * BASENAME  basename for output files (string)
    * JUN_COV  Read coverage list for each position overlaping junctions (list)
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

    return 1

#~~~~~~~FUNCTIONS FOR REPORT WRITING~~~~~~~#

def write_samp_report ():
    """
    Output a text file containin csv like tabs containing informations from
    source references and the number of reads sampled in each of its sequences.
    Need global parameters to be executed correctly.
    * BASENAME  basename for output files (string)
    * SOURCE_LIST   List containing source sublists with a reference to an
    instance of ReferenceGenome or ReferenceJunction and the requested number of
    reads to be generated (list)
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

    return 1


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~MAIN FUNCTION~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

if __name__ == '__main__':

    ## Name of the main program
    PROGRAM_NAME = "Isis"
    ## Version of the main program
    PROGRAM_VERSION = "0.01"

    #~~~~~~~Import configurations in global variables~~~~~~~#

    ## Instance of the configuration file parser/verifier IsisConf
    CONFIG = IsisConf(PROGRAM_NAME, PROGRAM_VERSION)

    ## Host genome fasta file path (string)
    HOST_GENOME = CONFIG.get("host_genome")
    ## Virus genome fasta file path (string)
    VIRUS_GENOME = CONFIG.get("virus_genome")
    ## Basename for output files (string)
    BASENAME = CONFIG.get("basename")
    ## Pair end mode. If false single end mode (bolean)
    PAIR = CONFIG.get("pair")
    ## Total number of reads to be generated (integer)
    READ_NUM = CONFIG.get("read_num")
    ## Lenght of the reads to be generated (integer)
    READ_LEN = CONFIG.get("read_len")
    ## Frequency of single nucleotide mutation to be randomly introduced in reads (float)
    MUT_FREQ = CONFIG.get("mut_freq")
    ## Allow lowercase characters in the DNA sequence (bool)
    REPEATS = CONFIG.get("repeats")
    ## Allow IUPAC ambigous characters in the DNA sequence (bool)
    AMBIGUOUS = CONFIG.get("ambiguous")
    ## Activate graphical output (bool)
    GRAPH = CONFIG.get("graph")
    ## Write out a sampling report (bool)
    REPORT = CONFIG.get("report")
    ## Number of reads to be sampled in host genome (integer)
    N_HOST = CONFIG.get("nread_host")
    ## Number of reads to be sampled in virus genome (integer)
    N_VIRUS = CONFIG.get("nread_virus")
    ## Number of reads to be sampled in true junctions (integer)
    N_TJUN = CONFIG.get("nread_tjun")
    ## Number of reads to be sampled in false junctions (integer)
    N_FJUN = CONFIG.get("nread_fjun")
    ## Minimal number of pb from one of the 2 references in a read or read pair overlaping a junction (integer)
    MIN_CHIMERIC = CONFIG.get("min_chimeric")
    ## Number of uniq source junctions to be generated in true junction dict (integer)
    UNIQ_TJUN = CONFIG.get("uniq_tjun")
    ## Number of uniq source junctions to be generated in false junction dict (integer)
    UNIQ_FJUN = CONFIG.get("uniq_fjun")
    ## Minimal size of sonication fragments (integer)
    SONIC_MIN = CONFIG.get("sonic_min")
    ## Modal size of sonication fragments (integer)
    SONIC_MODE = CONFIG.get("sonic_mode")
    ## Maximal size of sonication fragments (integer)
    SONIC_MAX = CONFIG.get("sonic_max")
    ## Thickness of the sonication peak (integer)
    SONIC_CERTAINTY = CONFIG.get("sonic_certainty")
    ## Range of quality for reads (string)
    QUAL_SCALE = CONFIG.get("qual_scale")
    ## Quality scale in which fastq quality string will be converted (sanger, illumina, solexa...) (string)
    QUAL_RANGE = CONFIG.get("qual_range")
    ## Lenght of junctions to be generated (integer)
    JUN_LEN = SONIC_MAX if PAIR else READ_LEN

    #~~~~~~~Instanciate reference genomes and junctions~~~~~~~#

    ## Viral genome reference object (ReferenceGenome)
    VIRUS = RefGen ("virus", VIRUS_GENOME)
    ## Host genome reference object (ReferenceGenome)
    HOST = RefGen ("host", HOST_GENOME)
    ## True junctions reference object (ReferenceJunctions)
    TJUN = RefJun("True_Junction", MIN_CHIMERIC, JUN_LEN, UNIQ_TJUN, VIRUS, HOST, REPEATS, AMBIGUOUS)
    ## False junctions reference object (ReferenceJunctions)
    FJUN = RefJun("False_Junction", MIN_CHIMERIC, JUN_LEN, UNIQ_FJUN, VIRUS, HOST, REPEATS, AMBIGUOUS)

    ## Convenient list of reference sequence source associated with the number of read to generate for each of them
    SOURCE_LIST = [[VIRUS, N_VIRUS], [HOST, N_HOST], [TJUN, N_TJUN], [FJUN, N_FJUN]]

    #~~~~~~~Prepare tools and variables for graphical output~~~~~~~#

    if GRAPH:
        # Third party package matplotlib imported only if needed
        from matplotlib import pyplot as plt
        ## Lenghts of sonication fragments (list)
        FRAG_LEN = []
        ## Read depth over all positions of junctions (list)
        JUN_COV = [0 for i in range (JUN_LEN*2)]

    #~~~~~~~Lauch main fastq writing functions~~~~~~~#

    if PAIR: # Pair end mode
        print ("Start pair end mode fastq sampling")
        IsisPair()

    else: # Single end mode specific imports
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
