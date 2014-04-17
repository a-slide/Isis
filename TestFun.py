from ReferenceJunctions import ReferenceJunctions
from ReferenceGenome import ReferenceGenome
from SlicePicker import SlicePicker
from pprint import pprint as pp

r1 = ReferenceGenome("Mouse","datasets/small_multi.fa")
r2 = ReferenceGenome("AAV","datasets/Bacterial_backbone.fa")
tj = ReferenceJunctions ("True_junctions", 10, 40, 20, r1, r2, True, False)
fj = ReferenceJunctions ("False_junctions", 10, 40, 20, r1, r2, True, False)
s = SlicePicker()

liste =[]
