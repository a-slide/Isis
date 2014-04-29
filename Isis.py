#!/usr/bin/env python
# -*- coding: utf-8 -*-

from multiprocessing import Pool, Manager, cpu_count
from sys import argv

#from IsisConf import IsisConf
from ReferenceGenome import ReferenceGenome
from ReferenceJunctions import ReferenceJunctions
from SlicePicker import SlicePickerSingle
from SlicePicker import SlicePickerPair
from QualGenerator import QualGenerator
from FastqGenerator import FastqGeneratorSingle
from FastqGenerator import FastqGeneratorPair

########################################################################################################################

def main ():
    """main function"""     
    
    pair = False
    basename = "test"
    
    # Import reference genomes and create Reference Junctions
    virus = ReferenceGenome ("virus", "datasets/virusl_backbone.fa")
    host = ReferenceGenome ("host","datasets/host_plasmid.fa")
    tj = ReferenceJunctions("True_Junction",50, 300, 10, virus, host, 1, 1)
    fj = ReferenceJunctions("False_Junction",50, 300, 10, virus, host, 1, 1)

    # Reset sampling counters of Virus and Host reference object 
    virus.reset_samp_counter()
    host.reset_samp_counter()

    # Different usage for single end or paired end
    if pair:
        slicer = SlicePickerPair(200, 300, 700, 350, 10, 1,1,0.1)
    else:
        slicer = SlicePickerSingle(150,1,1,0.1)
    
    # 
    qualgen = QualGenerator (150, "medium")
    fastgen = FastqGenerator(slicer, qualgen, "fastq-sanger")
    reset_file (basename)
    
    write_fastq (fastgen, virus, 100, basename, "fastq-sanger")
    write_fastq (fastgen, host, 100, basename, "fastq-sanger")
    write_fastq (fastgen, tj, 100, basename, "fastq-sanger")
    write_fastq(fastgen, fj, 100, basename, "fastq-sanger")

    else:
        slicer = SlicePickerSingle(150,1,1,0.1)
        qualgen = QualGenerator (150, "medium")
        fastgen = FastqGeneratorSingle (slicer, qualgen, "fastq-sanger")
        
        reset_single (basename)
        
        write_single (fastgen, virus, 100, basename, "fastq-sanger")
        write_single (fastgen, host, 100, basename, "fastq-sanger")
        write_single (fastgen, tj, 100, basename, "fastq-sanger")
        write_single (fastgen, fj, 100, basename, "fastq-sanger")

    # Write sampling reports
    virus.write_samp_report()
    host.write_samp_report()
    tj.write_samp_report()
    fj.write_samp_report()
    
    exit (1)
    
def reset_single (basename):
    with open (basename+".fastq", "w"):
        pass
    return 1
        
def reset_pair(basename):
    with open (basename+"R1.fastq", "w"):
        pass
    with open (basename+"R2.fastq", "w"):
        pass
    return 1

def write_single_mp (fastgen, source, nread, filename, qual_scale):

    # Define a queue manager and associate a queue 
    manager = Manager()
    queue = manager.Queue()
           
    # Automatically determine the number of available thread + 1
    try:
        nb_thread = cpu_count() + 1
    except NotImplementedError:
        print "cpu_count method is not available on your system" 
        nb_thread = 2
    
    # Define a pool of workers
    pool = Pool(nb_thread)

    # Put writter to work first
    writer = pool.apply_async(writer, (basename+".fastq", queue))
        
    # Start all workers
    jobs = []
    print "Done\n"
    for i in range(nread):
        job = pool.apply_async(self.worker, (i, queue))
        jobs.append(job)

def worker(self, i,queue):
    # Ask a fastq to generate_fastq
    #read = self.generate_fastq(source)
    # Add Uniq read identifier
    #read.id += "|#{:010}".format(i)
    # Add the read to the queue
    queue.put(i)

def writer(self, filename, q):
    '''listens for messages on the q, writes to file. '''
    f = open(filename, 'w') 
    while True:
        read = q.get()
        if m == 'kill':
            break
            
        f.write(read)#.format(self.qual_scale))
        f.flush()
    f.close()


def write_single (fastgen, source, nread, filename, qual_scale):
    
    with (open(filename, "a")) as f:
        for i in range (nread):
            # Ask a read to 
            read = fastgen.generate_fastq(source)
            # Uniq read identifier
            read.id += "|#{:010}".format(i)
            # Write the a fastq formated read
            f.write(read.format(qual_scale))
            

########################################################################################################################

if __name__ == '__main__':
    main()
    
if __name__ == '__write_fastq__':
    write_fastq()
    
if __name__ == '__reset_file__':
    reset_file()

if __name__ ==  '__QualGenerator__'
    if argv[1] == "pe":
        QualGeneratorPair()
    else if argv[1] == "se":
        QualGeneratorSingle()

if __name__ ==  '__QualGenerator__'
    if argv[1] == "pe":
        QualGeneratorPair()
    else if argv[1] == "se":
        QualGeneratorSingle()



   def write_fastq_mp(self, source, nread, basename):
        """
        """


        # Start all workers
        jobs = []
        print "Done\n"
        for i in range(nread):
            job = pool.apply_async(self.worker, (i, queue))
            jobs.append(job)

        # Collect results from the workers through the pool result queue
        for job in jobs: 
            job.get()

        # Now we are done, kill the writter
        q.put('kill')
        pool.close()
    


 #def write_fastq (self, source, nread, basename):
        #""" 
        #"""
        #("\tWritting Fastq file from {}".format(source.getName()))
        
        #filename1 = basename + "_R1.fastq"
        #filename2 = basename + "_R2.fastq"
        
        #with (open(filename1, "w")) as f1:
            #with (open(filename2, "w")) as f2:
                #for i in range (nread):
                    #read1, read2 = self.generate_fastq(source)
                    #read1.id += "|#{:010}".format(i)
                    #read2.id += "|#{:010}".format(i)
                    #f1.write(read1.format(self.qual_scale))
                    #f2.write(read2.format(self.qual_scale))
