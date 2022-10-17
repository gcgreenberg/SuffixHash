import os
import numpy as np
from multiprocessing import Pool
from numpy.random import randint
from sympy import nextprime

base_to_bin = {'A':'00', 'T':'01', 'G':'10', 'C':'11'}
# rc_base = {'A':'T','T':'A','G':'C','C':'G'}
# revcomp = lambda seq: ''.join([rc_base[b] for b in seq[::-1]])

class random_hash_func():
    def __init__(self,n_bits,k):
        self.a = randint(2**n_bits)
        self.b = randint(2**n_bits)
        self.max_hash = nextprime(2**n_bits)
        self.k = k
    
    def seq_to_bin(self, seq):
#         if rc: seq = revcomp(seq)
        bin_repr = ''.join([base_to_bin[b] for b in seq])
        return bin_repr 
    
    def hash_val(self, bin_repr):
        hash_val = (int(bin_repr,2)*self.a + self.b) % self.max_hash
        return hash_val
    
    def __call__(self, seq):
        n_kmers = len(seq)-self.k+1; assert n_kmers>0
        bin_repr = self.seq_to_bin(seq)
        hashes = [self.hash_val(bin_repr[2*i:2*(i+self.k)]) for i in range(n_kmers)]
        return hashes
    
def get_seq_sketches(sketcher,seq,hs): 
    return [sketcher(seq,h) for h in hs]
#     return [sketcher(seq,h) for h in hs], [sketcher(revcomp(seq),h) for h in hs]

def get_all_sketches(sketcher, seqs, k, n_hash, n_bits):
    hs = [random_hash_func(n_bits=n_bits, k=k) for _ in range(n_hash)]
    args = ((sketcher,seqs[i],hs) for i in range(len(seqs)))
    with Pool() as pool:
        sketches = pool.starmap(get_seq_sketches, args)
    return sketches


