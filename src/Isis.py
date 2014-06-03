#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
@package    Isis
@brief      Main package of Isis program. Contains functions to call sequence generation and write fastq in file
@copyright  [GNU General Public License v2](http://www.gnu.org/licenses/gpl-2.0.html)
@author     Adrien Leger <adrien.leger@gmail.com>
"""

#~~~~~~~GLOBAL IMPORTS~~~~~~~#

# Standard library packages
import gzip

# Third party package
from matplotlib import pyplot as plt

# Local packages
from IsisConf import IsisConf, IsisConfException
from ReferenceGenome import ReferenceGenome as RefGen
from ReferenceJunctions import ReferenceJunctions as RefJun
from SlicePicker import SlicePickerSingle, SlicePickerPair
from QualGenerator import QualGenerator
from FastqGenerator import FastqGeneratorSingle, FastqGeneratorPair

#~~~~~~~GLOBAL VARIABLES~~~~~~~#

# Store lenght of fragment for pair end sonication profile graphical output
FRAG_LEN = []
# Store read coverage of over junctions
JUN_COV = []

#~~~~~~~MAIN~~~~~~~#

def main ():
    """main function"""

    # Instantiate the configuration file parser/verifier
    try:
        config = IsisConf()
    except IsisConfException as E:
        print (E)
        exit (0)

    # Import configurations in local variables
    host_genome = config.get("host_genome")
    virus_genome = config.get("virus_genome")
    basename = config.get("basename")
    pair = config.get("pair")
    read_num = config.get("read_num")
    read_len = config.get("read_len")
    mut_freq = config.get("mut_freq")
    repeats = config.get("repeats")
    ambiguous = config.get("ambiguous")
    graph = config.get("graph")
    report = config.get("report")
    n_host = config.get("nread_host")
    n_virus = config.get("nread_virus")
    n_tjun = config.get("nread_tjun")
    n_fjun = config.get("nread_fjun")
    min_chimeric = config.get("min_chimeric")
    uniq_tjun = config.get("uniq_tjun")
    uniq_fjun = config.get("uniq_fjun")
    sonic_min = config.get("sonic_min")
    sonic_mode = config.get("sonic_mode")
    sonic_max = config.get("sonic_max")
    sonic_certainty = config.get("sonic_certainty")
    qual_scale = config.get("qual_scale")
    qual_range = config.get("qual_range")
    jun_len = sonic_max if pair else read_len

    # Import reference genomes
    virus = RefGen ("virus", virus_genome)
    host = RefGen ("host", host_genome)
    # Create reference junctions
    tjun = RefJun("True_Junction", min_chimeric, jun_len, uniq_tjun, virus, host, repeats, ambiguous)
    fjun = RefJun("False_Junction", min_chimeric, jun_len, uniq_fjun, virus, host, repeats, ambiguous)

    # Define a convenient ref list to simplify subsequent operations. The 3rd
    # position bool if used for sources coverage graph will be output
    if graph :
        source_list = [[virus, n_virus, 0], [host, n_host, 0], [tjun, n_tjun, 1], [fjun, n_fjun, 1]]
        # Initialisation of the coverage list
        init_cov_list (jun_len*2)
    else :
        source_list = [[virus, n_virus, 0], [host, n_host, 0], [tjun, n_tjun, 0], [fjun, n_fjun, 0]]

    # Pair end mode
    if pair:
        print ("Start pair end mode fastq sampling...")
        # Instantiate accessory classes
        slicer = SlicePickerPair(read_len, sonic_min, sonic_mode, sonic_max,
            sonic_certainty, repeats, ambiguous, mut_freq)
        qualgen = QualGenerator (read_len, qual_range)
        fastgen = FastqGeneratorPair(slicer, qualgen)
        # Write paired fastq files for all source in source_list
        write_fastq_pair (fastgen, source_list, basename, qual_scale)
        # Output a graphical representation of the fragment len distribution
        if graph:
            frag_len_graph (sonic_min, sonic_max, basename)

    # Single end mode
    else:
        print ("Start single end mode fastq sampling")
        # Instantiate accessory classes
        slicer = SlicePickerSingle (read_len, repeats, ambiguous, mut_freq)
        qualgen = QualGenerator (read_len, qual_range)
        fastgen = FastqGeneratorSingle (slicer, qualgen)
        # Write paired fastq files for all source in source_list
        write_fastq_single (fastgen, source_list, basename, qual_scale)

    # Write sampling reports if requested
    if report:
        write_report ([virus, host, tjun, fjun])
    if graph:
        jun_cov_graph (basename)

    print "DONE"
    exit (1)

#~~~~~~~FUNCTIONS~~~~~~~#

def write_fastq_single (fastgen, source_list, basename, qual_scale):
    """ Function printing the request number of single end reads per reference
    in a fastq.gz file
    """
    try:
        f = gzip.open(basename + ".fastq.gz", 'w')

        for source, nread, graph in source_list:
            print ("\tWritting {} read(s) in Fastq file from {}".format(nread, source.getName()))
            # Calculate the number of digits in nread for read naming
            max_len = len(str(nread))

            for i in range (nread):
                # Ask a read to the source throught fastgen
                read = fastgen.generate_fastq(source)
                # Generate a uniq identifier
                read.id = generate_id(i, max_len, read.annotations)
                # Write the fastq formated read
                f.write(read.format(qual_scale))

                # if graph is True
                if graph:
                    update_jun_cov (read.annotations["location"][0], read.annotations["location"][1])
        f.close()

    except IOError as E:
        print (E)
        exit (0)

def write_fastq_pair (fastgen, source_list, basename, qual_scale):
    """ Function printing the request number of pair end reads per reference
    in a fastq.gz file
    """
    try:
        f1 = gzip.open(basename + ".R1.fastq.gz", 'w')
        f2 = gzip.open(basename + ".R2.fastq.gz", 'w')

        for source, nread, graph in source_list:
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
                f1.write(read1.format(qual_scale))
                f2.write(read2.format(qual_scale))

                # if graph is True
                if graph:
                    update_jun_cov (read1.annotations["location"][0], read1.annotations["location"][1])
                    update_jun_cov (read2.annotations["location"][0], read2.annotations["location"][1])
                    update_frag_len (read1.annotations["frag_len"])
        f1.close()
        f2.close()

    except IOError as E:
        print (E)
        exit (0)

def init_cov_list (size):
    for i in range (size):
        JUN_COV.append (0)

def update_jun_cov (start, end):
    for i in range (start, end):
        JUN_COV [i] += 1

def update_frag_len (frag_len):
    FRAG_LEN.append(frag_len)

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


def frag_len_graph (min, max, basename):
    """ Output a graphical representation of the fragment len distribution
    """
    print ("\tCreating a graphical output detailling fragment length distribution")

    # Create a figure object and adding details
    fig = plt.figure(figsize=(15, 10), dpi=100)
    plt.title("Distribution of fragment length")
    plt.ylabel('Relative Count')
    plt.xlabel('Size of fragment')

    # Plot value from FRAG_LEN list in an histogram reprensentation
    plt.hist(FRAG_LEN, bins=(max-min)/5, range=(min,max), normed=1,
        facecolor='green', alpha=0.5, align='mid')

    # Tweak spacing to prevent clipping of ylabel
    plt.subplots_adjust(left=0.15)

    # Export figure to file
    fig.savefig(basename+'_distribution.png')

def jun_cov_graph (basename):
    """ Output a graphical representation of read coverage over junction
    """
    print ("\tCreating a graphical output detailling read coverage over junctions")

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
    fig.savefig(basename+'_junction_coverage.png')

def write_report (source_list):
    for source in source_list:
        source.write_samp_report()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

if __name__ == '__main__':
    main()
