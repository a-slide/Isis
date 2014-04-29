from ReferenceGenome import ReferenceGenome
from ReferenceJunctions import ReferenceJunctions

from SlicePicker import SlicePickerSingle
from SlicePicker import SlicePickerPair

from QualGenerator import QualGenerator

from FastqGenerator import FastqGenerator
from FastqGenerator import FastqGeneratorSingle
from FastqGenerator import FastqGeneratorPair


bacteria = ReferenceGenome ("Bacteria", "datasets/Bacterial_backbone.fa")
helper = ReferenceGenome ("Helper","datasets/Helper_plasmid.fa")
junction = ReferenceJunctions("Junction",50, 300, 10, bacteria, helper, 1, 1)

slicer1 = SlicePickerSingle(200,1,1,0.1)
slicer2 = SlicePickerPair(200, 220, 300, 350, 10, 1,1,0.1)

qualgen = QualGenerator (200, "medium")

fast1 = FastqGeneratorSingle (slicer1, qualgen, "fastq-sanger")
fast2 = FastqGeneratorPair (slicer2, qualgen, "fastq-sanger")

fast1.write_fastq_mp(helper, 10000, "test")

helper.write_samp_report ()

