#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Standard library packages
from multiprocessing import Pool, Manager, cpu_count
from sys import argv
import gzip

# Third party package
from matplotlib import pyplot

# Local packages
#from IsisConf import IsisConf
from ReferenceGenome import ReferenceGenome
from ReferenceJunctions import ReferenceJunctions
from SlicePicker import SlicePickerSingle
from SlicePicker import SlicePickerPair
from QualGenerator import QualGenerator
from FastqGenerator import FastqGeneratorSingle
from FastqGenerator import FastqGeneratorPair

##### Global Variables #####

# Store lenght of each fragment for pair end sonication profile
# graphical output

FRAG_LEN = []

#####    MAIN   #####

def main ():
    """main function"""

    #try:
        #fasta_filelist_import()
    #except IsisConfException as E:
        #print (E)
        #exit (1)

    pair = True
    basename = "test"
    path_vg = "../datasets/Bacterial_backbone.fa"
    path_hg = "../datasets/Helper_plasmid.fa"
    nread_vg = 2000
    nread_hg = 1000
    nread_tj = 15000
    nread_fj = 10000
    sonic_min = 300
    sonic_mode = 350
    sonic_max = 700
    sonic_certainty = 10
    read_len = 150
    ambiguity = 1
    repeats = 1
    mut_freq = 0.05
    quality = "medium"
    quality_scale = "fastq-sanger"
    min_chimeric = 50
    size_junction = sonic_max if pair else read_len
    n_tj = 150
    n_fj = 10000

    # Import reference genomes and create Reference Junctions
    virus = ReferenceGenome ("virus", path_vg)
    host = ReferenceGenome ("host", path_hg)
    tj = ReferenceJunctions("True_Junction",min_chimeric, size_junction, n_tj, virus, host, repeats, ambiguity,)
    fj = ReferenceJunctions("False_Junction",min_chimeric, size_junction, n_fj, virus, host, repeats, ambiguity,)

    # Reset sampling counters of Virus and Host reference object
    virus.reset_samp_counter()
    host.reset_samp_counter()

    # Paired end mode
    if pair:
        print ("Start pair end mode fastq sampling...")
        # Instantiate accessory classes
        slicer = SlicePickerPair(read_len, sonic_min, sonic_mode, sonic_max, sonic_certainty, repeats, ambiguity, mut_freq)
        qualgen = QualGenerator (read_len, quality)
        fastgen = FastqGeneratorPair(slicer, qualgen)
        # Reset eventual files and write paired fastq files
        reset_file (basename+"_R1.fastq.gz")
        reset_file (basename+"_R2.fastq.gz")
        write_fastq_pair (fastgen, virus, nread_vg, basename, quality_scale)
        write_fastq_pair (fastgen, host, nread_hg, basename, quality_scale)
        write_fastq_pair (fastgen, tj, nread_tj, basename, quality_scale)
        write_fastq_pair (fastgen, fj, nread_fj, basename, quality_scale)
        # Output a graphical representation of the fragment len distribution
        frag_len_graph (sonic_min, sonic_max, basename)

    # Single end mode
    else:
        print ("Start single end mode fastq sampling")
        # Instantiate accessory classes
        slicer = SlicePickerSingle (read_len, repeats, ambiguity, mut_freq)
        qualgen = QualGenerator (read_len, quality)
        fastgen = FastqGeneratorSingle (slicer, qualgen)
        # Reset eventual file and write the single endfastq file
        reset_file (basename+".fastq.gz")
        write_fastq_single (fastgen, virus, nread_vg, basename, quality_scale)
        write_fastq_single (fastgen, host, nread_hg, basename, quality_scale)
        write_fastq_single (fastgen, tj, nread_tj, basename, quality_scale)
        write_fastq_single (fastgen, fj, nread_fj, basename, quality_scale)

    # Write sampling reports
    virus.write_samp_report()
    host.write_samp_report()
    tj.write_samp_report()
    fj.write_samp_report()

    print "END"
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
    fig = pyplot.figure(figsize=(15, 10), dpi=100)
    pyplot.title("Distribution of fragment length")
    pyplot.ylabel('Relative Count')
    pyplot.xlabel('Size of fragment')

    pyplot.hist(FRAG_LEN, bins=(max-min)/5, range=(min,max), normed=0, facecolor='green', alpha=0.5, align='mid')

    # Tweak spacing to prevent clipping of ylabel
    pyplot.subplots_adjust(left=0.15)

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
