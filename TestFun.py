import gzip
from Bio import SeqIO

def fasta_open (filename):

    try: # Try block to manage
        if filename.rpartition(".")[-1] == "gz":
            print ("Uncompressing and extracting data")
            handle = gzip.open(filename, "r")
        else:
            print ("Extracting data")
            handle = open(filename, "r")

        seq_dict = SeqIO.to_dict(SeqIO.parse( handle, "fasta"))
        handle.close()
        return seq_dict
        
    except IOError:
           print (filename + 'is not readable')
           return None

def _calculate_proba(seq_dict):
    """Return a 2 entries list / 1 = name of the sequence / 2 = cumulative frequency of the sequence"""
    cumulative_len = 0

    # Calculate the cumulative length for all seq in seq_dict
    for record in seq_dict.values():
        cumulative_len += len(record.seq)

    # Calculate a cumulative frequency for all reference sequences and return the object
    return [[record.name, float(len(record.seq))/cumulative_len] for record in seq_dict.values()]
