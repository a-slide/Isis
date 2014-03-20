#!/usr/bin/env python
# -*- coding: utf-8 -*-

import gexf
from sys import argv as argv
from Bio import SeqIO
from multiprocessing import Pool
from pprint import pprint
from Sequence import Sequence
from ConfFileParser import ConfFileParser


########################################################################################################################

def main ():
    """main function"""
    filename = argv[0]
    conf = ConfFileParser()

    exit (0)

########################################################################################################################


def usage():
    """Simple usage function"""
    print "Usage: ", argv[0], "<Number of available threads><Threshold> <output> <fasta file 1> [<fasta file 2, ...]"
    print "\tExample : ", argv[0], " 4  1000  my_alignment my_seq.fasta  my_other_seq.fasta"

########################################################################################################################

if __name__ == '__main__':
    if len(argv) < 1:        # if not enought arg call usage function
        usage()
    else:
        main()              # else call the main function

