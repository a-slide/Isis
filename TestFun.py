from pprint import pprint as pp

def fasta_open (filename):

    try: # Try block to manage
        with open(filename) as fasta:
            return fasta.readlines()

    except IOError:
        print (filename + 'is not readable')
        return None

def fasta_parse (filename):
    return[line.rstrip() for line in fasta_open (filename)]



#def fasta_parse (filename):
#    seq_dict = {}
#   for line in fasta_open (filename):
#        if line[0] == ">":
#            key = line[1:].rstrip().replace(" ", "_")
#            seq_dict[key] = ""
#        else:
#            seq_dict[key] += line.rstrip()
#    return seq_dict

