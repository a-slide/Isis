from pprint import pprint as pp


def fasta_open (filename):

    try: # Try block to manage 
        with open(filename) as fasta:
            return fasta.readlines()

    except IOError:
        print (filename + 'is not readable')
        return None

def fasta_parse (filename):
    for line in fasta_open(filename):
        if line[0:1] == ">":
            
