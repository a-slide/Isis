
from Bio import SeqIO

from ReferenceGenome import ReferenceGenome
from SlicePicker import SlicePicker
from QualGenerator import QualGenerator

def FastqGenerator (source, file, quality, length, mutfreq):

    g = ReferenceGenome (source, file)
    print(repr(g))
    
    s = SlicePicker ()
    print(repr(s))
    
    q = QualGenerator (length, quality)
    print(repr(q))

    read = s.pick_single(g,length,1,1,mutfreq)
    read.letter_annotations["phred_quality"] = q.random_qual_string()
    read.id = "{}|{}|loc_{}_{}_{}".format(
        read.annotations["source"],
        read.annotations["refseq"],
        read.annotations["location"][0],
        read.annotations["location"][1],
        read.annotations["orientation"])
    
    print(read)
    print(read.format("qual"))
    print(read.format("fastq-illumina"))

    return read
    

    
