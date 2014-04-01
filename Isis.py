#!/usr/bin/env python
# -*- coding: utf-8 -*-

from multiprocessing import Pool
from pprint import pprint as pp

from IsisConf import IsisConf
from ReferenceGenome import ReferenceGenome
from ReferenceJunctions import ReferenceJunctions
from SlicePicker import SlicePicker


########################################################################################################################

def main ():
    """main function"""
    conf = IsisConf()
    print(repr(conf))
    
    genome_ref = ReferenceGenome (conf.get("hg_filename"))
    print(repr(genome_ref))
    virus_ref = ReferenceGenome (conf.get("vg_filename"))
    print(repr(virus_ref))
    true_junctions = ReferenceJunctions (conf.get("sonic_max"), conf.get("uniq_tj"), genome_ref, virus_ref, conf.get("repeats"), conf.get("ambigous"))
    print(repr(true_junctions))
    false_junctions = ReferenceJunctions (conf.get("sonic_max"), conf.get("uniq_fj"), genome_ref, virus_ref, conf.get("repeats"), conf.get("ambigous"))
    print(repr(false_junctions))
    exit (0)

########################################################################################################################

if __name__ == '__main__':
    main()
