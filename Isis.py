#!/usr/bin/env python
# -*- coding: utf-8 -*-

import gexf
from sys import argv as argv
from Bio import SeqIO
from multiprocessing import Pool
from pprint import pprint

########################################################################################################################



def main ():
    """main function"""
    filename = argv[0]
    conf = ConfFileParser()

    exit (0)
    


########################################################################################################################

class ConfFileParser:
"""This can be used to parse file containing a list of variables and its associated values
The separator beetween vaiable and value must be the same for all field and can be customized by users (by default =)
the class store the list of variable in a dictionnary named conf_dict"""
    
    # FONDAMENTAL METHODS
    
    def __init__ (self, filename, separator="="):
        # Object constructor initialized the path of the field to parse and a custom separator
        self.conf_dict = makemake_conf_dict (self, filename, separator)
    
    def __str__(self):
        # Short description string returned by print and str 
        print "ConfFileParser ... \n"
        
    def __repr__ (self):
        # Long description string used by interpreter and repr
        print "ConfFileParser ... \n"
        pprint (conf_dict)
        
    # GETERS
    
    def get_conf_dict (self):
        return conf_dict
    
    # Give acces to individual values in conf_dict by using its key name.
    def get_conf_var (self, varkey ):
        return conf_dict[varkey]
        
    # PRIVATE SUPPORT METHODS
    
    def make_conf_dict (self, filename, separator)
        return {name : value for name, value in conf_list()}
        
    def make_conf_list (self, filename, separator)
        with open (filename)
            return [
        
        except IOError:
            print '\n', filename, 'is not readable. The file will be ignored\n'

    
    
    
    
    
    
    
    def 


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

