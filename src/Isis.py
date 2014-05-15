#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Standard library packages
from multiprocessing import Pool, Manager, cpu_count
import gzip

# Third party package
from matplotlib import pyplot

# Local packages
from IsisConf import IsisConf, IsisConfException
from ReferenceGenome import ReferenceGenome as RefGen
from ReferenceJunctions import ReferenceJunctions as RefJun
from SlicePicker import SlicePickerSingle, SlicePickerPair
from QualGenerator import QualGenerator
from FastqGenerator import FastqGeneratorSingle, FastqGeneratorPair

##### Global Variables #####

# Store lenght of fragment for pair end sonication profile graphical output
FRAG_LEN = []

#####    MAIN   #####

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
    nread_host = config.get("nread_host")
    nread_virus = config.get("nread_virus")
    nread_tjun = config.get("nread_tjun")
    nread_fjun = config.get("nread_fjun")
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
    tjun = RefJun("True_Junction", min_chimeric, jun_len, uniq_tjun, virus,
        host, repeats, ambiguous)
    fjun = RefJun("False_Junction", min_chimeric, jun_len, uniq_fjun, virus,
        host, repeats, ambiguous)

    # Reset sampling counters of Virus and Host reference object
    virus.reset_samp_counter()
    host.reset_samp_counter()

    # Define a convenient ref list to simplify subsequent operations
    ref_list = [[virus, nread_virus],
                [host, nread_virus],
                [tjun, nread_tjun],
                [fjun, nread_fjun]]

    # Pair end mode
    if pair:
        print ("Start pair end mode fastq sampling...")
        # Instantiate accessory classes
        slicer = SlicePickerPair(read_len, sonic_min, sonic_mode, sonic_max,
            sonic_certainty, repeats, ambiguous, mut_freq)
        qualgen = QualGenerator (read_len, qual_range)
        fastgen = FastqGeneratorPair(slicer, qualgen)
        # Reset eventual files and write paired fastq files
        reset_file (basename+"_R1.fastq.gz")
        reset_file (basename+"_R2.fastq.gz")
        for ref, nread in ref_list:
            write_fastq_pair (fastgen, ref, nread, basename, qual_scale)
        # Output a graphical representation of the fragment len distribution
        frag_len_graph (sonic_min, sonic_max, basename)

    # Single end mode
    else:
        print ("Start single end mode fastq sampling")
        # Instantiate accessory classes
        slicer = SlicePickerSingle (read_len, repeats, ambiguous, mut_freq)
        qualgen = QualGenerator (read_len, qual_range)
        fastgen = FastqGeneratorSingle (slicer, qualgen)
        # Reset eventual file and write the single end fastq file
        reset_file (basename+".fastq.gz")
        for ref, nread in ref_list:
            write_fastq_single (fastgen, ref, nread, basename, qual_scale)

    # Write sampling reports
    for ref, nread in ref_list:
        ref.write_samp_report()

    print "DONE"
    exit (1)

#####    FUNCTIONS   #####

def reset_file (filename):
    """ Simple function to initialize the file and reset its content if
    existing
    """
    try:
        handle = gzip.open(filename, 'w')
        handle.close()
        print "\t{} initialized".format(filename)

    except IOError:
        print('CRITICAL ERROR. {} cannot by open for writing'.format(filename))
        exit

def write_fastq_single (fastgen, source, nread, basename, qual_scale):
    """
    """
    if nread == 0:
        print ("\tNo read sampled in {}".format(source.getName()))

    else:
        print ("\tWritting {} reads in Fastq file from {}".format(nread, source.getName()))
        filename = basename + ".fastq.gz"

        try:
            f = gzip.open(filename, 'a')
            for i in range (nread):
                # Ask a read to the source throught fastgen
                read = fastgen.generate_fastq(source)
                # Uniq read identifier
                read.id += "|#{:010}".format(i)
                # Write the fastq formated read
                f.write(read.format(qual_scale))
            f.close()

        except IOError:
            print('CRITICAL ERROR. {} cannot by open for writing'.format(filename))
            exit

def write_fastq_pair (fastgen, source, nread, basename, qual_scale):
    """
    """
    if nread == 0:
        print ("\tNo read sampled in {}".format(source.getName()))

    else:
        print ("\tWritting {} reads in Fastq files from {}".format( nread, source.getName()))
        filename1 = basename + "_R1.fastq.gz"
        filename2 = basename + "_R2.fastq.gz"

        try:
            f1 = gzip.open(filename1, 'a')
            f2 = gzip.open(filename2, 'a')
            for i in range (nread):
                # Ask a pair of read to the source throught fastgen
                read1, read2 = fastgen.generate_fastq(source)
                # Uniq read identifier
                read1.id += "|#{:010}".format(i)
                read2.id += "|#{:010}".format(i)
                # Write the fastq formated reads
                f1.write(read1.format(qual_scale))
                f2.write(read2.format(qual_scale))
                # Append the size of sonication frag to the global list
                FRAG_LEN.append(read1.annotations["frag_len"])

            f1.close()
            f2.close()

        except IOError:
            print('CRITICAL ERROR. {} cannot by open for writing'.format(filename))
            exit

def frag_len_graph (min, max, basename):
    """ Output a graphical representation of the fragment len
    distribution
    """
    print ("\tCreating a graphical output detailling fragment length distribution")

    # Create a figure object and adding details
    fig = pyplot.figure(figsize=(15, 10), dpi=100)
    pyplot.title("Distribution of fragment length")
    pyplot.ylabel('Relative Count')
    pyplot.xlabel('Size of fragment')

    # Plot value from FRAG_LEN list in an histogram reprensentation
    pyplot.hist(FRAG_LEN, bins=(max-min)/5, range=(min,max), normed=1,
        facecolor='green', alpha=0.5, align='mid')

    # Tweak spacing to prevent clipping of ylabel
    pyplot.subplots_adjust(left=0.15)

    # Export figure to file
    fig.savefig(basename+'_distribution.png')


#def write_single_mp (fastgen, source, nread, filename, qual_scale):

    ## Define a queue manager and associate a queue
    #manager = Manager()
    #queue = manager.Queue()

    ## Automatically determine the number of available thread + 1
    #try:
        #nb_thread = cpu_count() + 1
    #except NotImplementedError:
        #print "cpu_count method is not available on your system"
        #nb_thread = 2

    ## Define a pool of workers
    #pool = Pool(nb_thread)

    ## Put writter to work first
    #writer = pool.apply_async(writer, (basename+".fastq", queue))

    ## Start all workers
    #jobs = []
    #print "Done\n"
    #for i in range(nread):
        #job = pool.apply_async(self.worker, (i, queue))
        #jobs.append(job)

#def worker(self, i,queue):
    ## Ask a fastq to generate_fastq
    ##read = self.generate_fastq(source)
    ## Add Uniq read identifier
    ##read.id += "|#{:010}".format(i)
    ## Add the read to the queue
    #queue.put(i)

#def writer(self, filename, q):
    #'''listens for messages on the q, writes to file. '''
    #f = open(filename, 'w')
    #while True:
        #read = q.get()
        #if m == 'kill':
            #break

        #f.write(read)#.format(self.qual_scale))
        #f.flush()
    #f.close()


########################################################################################################################

if __name__ == '__main__':
    main()
