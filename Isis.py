#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Standard library packages
from multiprocessing import Pool, Manager, cpu_count
from sys import argv
import gzip

# Local packages
#from IsisConf import IsisConf
from ReferenceGenome import ReferenceGenome
from ReferenceJunctions import ReferenceJunctions
from SlicePicker import SlicePickerSingle
from SlicePicker import SlicePickerPair
from QualGenerator import QualGenerator
from FastqGenerator import FastqGeneratorSingle
from FastqGenerator import FastqGeneratorPair

#####    MAIN   #####

def main ():
    """main function"""     
    
    pair = True
    basename = "test"
    nread_vg = 200
    nread_hg = 100
    nread_tj = 1500
    nread_fj = 1000
    
    # Import reference genomes and create Reference Junctions
    virus = ReferenceGenome ("virus", "../datasets/Bacterial_backbone.fa")
    host = ReferenceGenome ("host","../datasets/Helper_plasmid.fa")
    tj = ReferenceJunctions("True_Junction",50, 700, 10, virus, host, 1, 0)
    fj = ReferenceJunctions("False_Junction",50, 700, 10, virus, host, 1, 1)
    
    # Reset sampling counters of Virus and Host reference object 
    virus.reset_samp_counter()
    host.reset_samp_counter()
    
    # Paired end mode
    if pair:
        print ("Start pair end mode fastq sampling...") 
        # Instantiate accessory classes
        slicer = SlicePickerPair(150, 300, 350, 700, 10, 1,1,0.01)
        qualgen = QualGenerator (150, "medium")
        fastgen = FastqGeneratorPair(slicer, qualgen, "fastq-sanger")
        # Reset eventual files and write paired fastq files
        reset_file (basename+"_R1.fastq.gz")
        reset_file (basename+"_R2.fastq.gz")
        write_fastq_pair (fastgen, virus, nread_vg, basename, "fastq-sanger")
        write_fastq_pair (fastgen, host, nread_hg, basename, "fastq-sanger")
        write_fastq_pair (fastgen, tj, nread_tj, basename, "fastq-sanger")
        write_fastq_pair (fastgen, fj, nread_fj, basename, "fastq-sanger")
    
    # Single end mode
    else:
        print ("Start single end mode fastq sampling") 
        # Instantiate accessory classes
        slicer = SlicePickerSingle (150,1,1,0.01)
        qualgen = QualGenerator (150, "medium")
        fastgen = FastqGeneratorSingle (slicer, qualgen, "fastq-sanger")
        # Reset eventual file and write the single endfastq file
        reset_file (basename+".fastq.gz")
        write_fastq_single (fastgen, virus, nread_vg, basename, "fastq-sanger")
        write_fastq_single (fastgen, host, nread_hg, basename, "fastq-sanger")
        write_fastq_single (fastgen, tj, nread_tj, basename, "fastq-sanger")
        write_fastq_single (fastgen, fj, nread_fj, basename, "fastq-sanger")

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
            f1.close()
            f2.close()
            
        except IOError:
            print('CRITICAL ERROR. {} cannot by open for writing'.format(filename))
            exit


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

        ## Start all workers
        #jobs = []
        #print "Done\n"
        #for i in range(nread):
            #job = pool.apply_async(self.worker, (i, queue))
            #jobs.append(job)

        ## Collect results from the workers through the pool result queue
        #for job in jobs: 
            #job.get()

        ## Now we are done, kill the writter
        #q.put('kill')
        #pool.close()
    
