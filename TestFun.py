from random import random
from random import randint
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pprint import pprint as pp

def _random_seq (len_seq):
    seq =''
    for i in range (len_seq):
        seq += 'ATCG'[randint(0,3)]
    return seq

def _create_simple_dict(nseq, len_seq):
    d = {}
    for i in range (nseq):
        name_seq = "Seq#"+str(i)
        seq = Seq(_random_seq(len_seq),IUPAC.IUPACAmbiguousDNA())
        d[name_seq] = SeqRecord(seq, id=str(i), name=name_seq )
    pp(d)
