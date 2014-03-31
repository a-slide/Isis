#!/usr/bin/env python
# -*- coding: utf-8 -*-

from sys import argv
from Bio import SeqIO
from multiprocessing import Pool
from pprint import pprint as pp
from IsisConf import IsisConf
#from Sequence import Sequence
#from ConfFileParser import ConfFileParser



########################################################################################################################

def main ():
    """main function"""
    filename = argv[0]
    conf = IsisConf()
    print(repr(conf))




    exit (0)

########################################################################################################################

if __name__ == '__main__':
    main()
