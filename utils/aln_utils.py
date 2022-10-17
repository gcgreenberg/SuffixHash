import numpy as np
from tqdm import tqdm

def pairwise_overlap_ests(seq_lens, sketches, est_overlap):
    pairwise_ests = []
    n = len(seq_lens)
    for i in tqdm(range(n), desc='Estimating pairwise overlaps', leave=True):
        for j in range(n):
            if i==j: continue
            theta_hat = est_overlap(sketches[i], sketches[j], seq_lens[i], seq_lens[j])
            if theta_hat > 0: 
                pairwise_ests.append((i+1,j+1,theta_hat,seq_lens[i],seq_lens[j]))
    return pairwise_ests

def find_overlaps(aln_file, sketches, seq_lens, est_overlap):
    pairwise_ests = pairwise_overlap_ests(seq_lens, sketches, est_overlap)
    print(f'# overlaps: {len(pairwise_ests)}')
    with open(aln_file, 'w') as fh:
        for aln in pairwise_ests:
            fh.write('\t'.join([str(x) for x in aln])+'\n')