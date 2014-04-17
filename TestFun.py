from ReferenceJunctions import ReferenceJunctions
from ReferenceGenome import ReferenceGenome
from SlicePicker import SlicePicker

r1 = ReferenceGenome("Mouse","datasets/small_multi.fa")
r2 = ReferenceGenome("AAV","datasets/Bacterial_backbone.fa")
j = ReferenceJunctions ("True_junctions", 20 , 40, 200, r1, r2, True, False)
