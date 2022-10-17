import numpy as np
from sympy import nextprime
import sys
from multiprocessing import Pool, cpu_count
from time import time

sys.path.append('/home/gcgreen2/alignment/SequenceAlignmentAndSketching/utils')
import seq_utils as su

base_to_bin = {'A':'00', 'T':'01', 'G':'10', 'C':'11'}
# rc_base = {'A':'T','T':'A','G':'C','C':'G'}
# revcomp = lambda seq: ''.join([rc_base[b] for b in seq[::-1]])

class random_hash_func():
    def __init__(self,n_bits,k):
        self.a = np.random.randint(2**n_bits)
        self.b = np.random.randint(2**n_bits)
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
    
def get_seq_sketches(i, rc):
    if rc:
        return [min(h(seqs_rc[i])) for h in hash_funcs]
    else:
        return [min(h(seqs[i])) for h in hash_funcs]

def get_all_sketches(rc):
    args = ((i,rc) for i in range(len(seqs)))
    start = time()
    with Pool() as pool:
        sketches = pool.starmap(get_seq_sketches, args)
    print('sketching took {} minutes'.format(int((time()-start)/60)))
    return sketches    

###########################################################################

def est_overlap(i, rc):
    n_hash = len(sketches[i])
    ests = {}
    pairs = [j%len(seqs) for j in range(i+1,i+1+len(seqs)//2)]
    for j in pairs:
        alpha = sum([hash_i==hash_j for hash_i,hash_j in zip(sketches[i],sketches_rc[j] if rc else sketches[j])])
        if alpha >= 2: 
            theta_hat = (len(seqs[i])+len(seqs[j])) * alpha/(n_hash+alpha)
            ests[j] = round(theta_hat, 1)
    return ests

def pairwise_overlap_ests(rc):
    n = len(seqs)
    args = ((i,rc) for i in range(n))
    chunksize = int(np.ceil(n / cpu_count()))
    start = time()
    with Pool() as pool:
        ests = pool.starmap(est_overlap, args, chunksize)
    print('pairwise alignment took {} minutes'.format(int((time()-start)/60)))
    return ests


###############################################################################

def write_overlaps(aln_path, estimates, ids, modifier='w', direc='+', use_ids=False):
    print(f'# overlaps: {sum([len(estimates[i]) for i in range(len(estimates))])}')
    ids = ids if use_ids else [str(i) for i in range(len(ids))]
    with open(aln_path, modifier) as fh:
        for i in range(len(estimates)):
            for j in estimates[i]:
                line = [ids[i], ids[j], str(estimates[i][j]), direc]
                fh.write('\t'.join(line)+'\n')

def find_overlaps(fasta_file, aln_paths, **args):
    global seqs, seqs_rc, hash_funcs, sketches, sketches_rc
    seqs,ids = su.get_seqs(fasta_file, return_ids=True)
    
    for aln_path,k in zip(aln_paths, args['ks']):
        hash_funcs = [random_hash_func(n_bits=args['n_bits'], k=k) for _ in range(args['n_hash'])]
        sketches = get_all_sketches(rc=False)
        estimates = pairwise_overlap_ests(rc=False)
        write_overlaps(aln_path, estimates, ids, modifier='w', direc='+', use_ids=args['use_ids'])

        if args['rc']:
            del estimates
            seqs_rc = su.revcomp_seqs(seqs)
            sketches_rc = get_all_sketches(rc=True)
            estimates = pairwise_overlap_ests(rc=True)
            write_overlaps(aln_path, estimates, ids, modifier='a', direc='-', use_ids=args['use_ids'])    
    
#######################################################
#ARCHIVE
    
# def pairwise_overlap_ests(seq_lens, sketches):
#     pairwise_ests = []
#     n = len(seq_lens)
#     for i in tqdm(range(n), desc='Estimating pairwise overlaps', leave=True):
#         for j in range(n):
#             if i==j: continue
#             theta_hat = est_overlap(sketches[i], sketches[j], seq_lens[i], seq_lens[j])
#             if theta_hat > 0: 
#     pairwise_ests.append((i+1,j+1,theta_hat,seq_lens[i],seq_lens[j]))
#     return pairwise_ests

# def find_overlaps(fasta_file, aln_file, k, n_hash=100, n_bits=20, **args):
#     seqs = su.get_seqs(fasta_file); seq_lens = [len(s) for s in seqs]
#     minhashes = hu.get_all_sketches(minhash, seqs, k, n_hash, n_bits)
#     au.find_overlaps(aln_file, minhashes, seq_lens, est_overlap)

