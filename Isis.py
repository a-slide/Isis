#!/usr/bin/env python
# -*- coding: utf-8 -*-

import gexf
from sys import argv as argv
from Bio import SeqIO
from Bio import pairwise2
from multiprocessing import Pool

########################################################################

# Global variable to simplify communication between functions
HOSTGENOME = []
VIRUSGENOME = []

########################################################################

def main ():
    """main function"""

    # Try bloc to catch errors if input data files are empty are incorectly formated
    try:
        fasta_filelist_import()
    except MissingDataException, ex:
        print (ex)
        exit (1)

    print_fasta ()

    exit (0)


def fasta_filelist_import():
    """Iterative filling of a list of sequence from the list of file given in argument"""

    print '\nStarting to import data...\n'

    # For each filepath in arguments
    for filename in argv[4:]:
        fasta_file_import(filename)

    # Raising exception if the sequence list is empty after parsing
    if not LISTSEQ:
        raise MissingDataException ("File(s) are either empty or not in fasta format. Please verify your files")


def fasta_file_import(filename):
    """Open an parse a fasta containing file with SeqIO.parse and return a list of Biopython sequence objects"""

    try: # Bloc try to manage opening errors
        # Iterative call of SeqIO.parse from Biopython and concatenation of sequences
        for seq_record in SeqIO.parse(filename, "fasta"):
            LISTSEQ.extend([seq_record])
        print '\n', filename, 'had been imported\n'

    except IOError:
        print '\n', filename, 'is not readable. The file will be ignored\n'


def print_fasta():
    """May be used to verify imported sequences"""
    for seq in LISTSEQ:
        print 'ID        ', (seq.id)
        print 'Sequence  ', (repr(seq.seq))
        print 'Lenght    ', (len(seq)), 'bp\n'

        print (len(LISTSEQ)), 'fasta sequences had been imported\n'







def usage():
    """Simple usage function"""
    print "Usage: ", argv[0], "<Number of available threads><Threshold> <output> <fasta file 1> [<fasta file 2, ...]"
    print "\tExample : ", argv[0], " 4  1000  my_alignment my_seq.fasta  my_other_seq.fasta"
    print "\tNumber of thread is used to distribute processing thanks to the multiprocessing package"
    print "\tThreshold value is used to filter low score interactions (use 0 to obtain a complete graph)"
    print "\tOutput will be used as a prefixe for the gexf output file"
    print "\One or sereval files containing fasta sequence(s) are needed as input data\n"
    print "\tA sequence in FASTA format consists of: "
    print "\t\tOne line starting with a '>' sign, followed by a sequence identification code."
    print "\t\tIt is optionally be followed by a textual description of the sequence."
    print "\t\tOne or more lines containing the sequence itself (Usually 80 characters blocks).\n"
    print "\tA file in FASTA format may comprise more than one sequence.\n"
    print "\tRequire the following python packages :"
    print "\tBiopython, multiprocessing and gefx (http://pythonhosted.org/pygexf/users.html)"



########################################################################

class ConfFileParser:








########################################################################

class MissingDataException(Exception):
    """Custom Exception Class to handle empty or incorrectly formated data files"""
    def __init__(self, msg):    # Object constructor initialized with a custom user message
        self.msg = msg

    def __str__(self):          # String returned by print
        return "MissingDataException\n" + self.get_msg()

    def get_msg (self):         # msg getter
        return self.msg

########################################################################

if __name__ == '__main__':
    if len(argv) < 4:        # if not enought arg call usage function
        usage()
    else:
        main()              # else call the main function

if __name__ == '__fasta_filelist_import__':
    fasta_filelist_import()
if __name__ == '__fasta_file_import__':
    fasta_file_import()
if __name__ == '__print_fasta__':
    print_fasta()
