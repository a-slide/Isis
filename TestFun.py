
from Bio import SeqIO

from ReferenceGenome import ReferenceGenome
from SlicePicker import SlicePickerPair
from QualGenerator import QualGenerator

def FastqGenerator (source, file, quality, read_len, sonic_min, sonic_max, sonic_mode, sonic_certainty, repeats, ambiguous, mut_freq):

    g = ReferenceGenome (source, file)
    print(repr(g))
    
    s = SlicePickerPair(read_len, sonic_min, sonic_max, sonic_mode, sonic_certainty, repeats, ambiguous, mut_freq)
    print(repr(s))
    
    q = QualGenerator (read_len, quality)
    print(repr(q))

    read1, read2 = s.pick_slice(g)
    read1.letter_annotations["phred_quality"] = q.random_qual_string()
    read2.letter_annotations["phred_quality"] = q.random_qual_string()
    
    read1.id = "R1"
    read2.id = "R2"

    print(read1)
    print(read2)
    
    print(read1.format("fastq-illumina"))
    print(read2.format("fastq-illumina"))
    
    return (read1, read2)
    

    
