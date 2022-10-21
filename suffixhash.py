import numpy as np
from scipy.optimize import minimize
from scipy.stats import norm
from tqdm import tqdm
from os.path import exists
from time import time
import sys
from multiprocessing import Pool, cpu_count, set_start_method

# sys.path.append('/home/gcgreen2/alignment/SequenceAlignmentAndSketching/utils')
from utils import seq_utils as su, utils
# import utils

base_to_int = {'A':0, 'T':1, 'G':2, 'C':3}

# def get_stat_inflec(arr,k_thresh): return sum([max(0,v-k_thresh) for v in arr])
# def get_stat_inflec2(arr,k_thresh): return sum([v for v in arr if v>k_thresh])
def get_stat_max(arr,k_thresh): return max(arr)
# def get_stat_max2(arr,k_thresh): return round(np.dot(np.sort(np.partition(arr,-5)[-5:])[::-1], [2**(-x) for x in range(5)]),1)    


def obtain_masks(string_size,number_of_masks): 
    masks = [np.random.choice([0,1,2,3],string_size) for _ in range(number_of_masks)]
    return masks

def lexicographic_first(seq,mask):
    best,idxs = [],[]
    for i,b in enumerate(seq):
        idxs.append(i)
        best.append((seq[i]-mask[i-idxs[0]]) % 4) # append to end of best (best idx is at idxs[0])
        j = 1
        while(len(idxs) > j):
            bm = (b-mask[i-idxs[j]]) % 4
            if bm > best[i-idxs[j]]:
                del idxs[j]
            elif bm < best[i-idxs[j]]:
                best = [*best[:i-idxs[j]], bm]
                idxs = idxs[j:]
                j = 1
            else: 
                j += 1
    return idxs[0]
    
# def get_seq_sketches(seq, masks): 
def get_seq_sketches(i, rc): 
    if rc:
        seq_int_repr = [base_to_int[b] for b in seqs_rc[i]]
    else:
        seq_int_repr = [base_to_int[b] for b in seqs[i]]
    sketches = [lexicographic_first(seq_int_repr,mask) for mask in masks]
    return sketches


def get_all_sketches(rc):
    args = ((i,rc) for i in range(len(seqs)))
    chunksize = int(np.ceil(len(seqs)/cpu_count()/4))
    start = time()
    with Pool() as pool:
        sketches = pool.starmap(get_seq_sketches, args, chunksize)
    print('sketching took {} minutes'.format(int((time()-start)/60)))
    return sketches


######################################################

def get_matching_bases(i, j, seq1, seq2):
    count = 0
    while(i<len(seq1) and j<len(seq2) and seq1[i]==seq2[j]):
        count += 1
        i += 1
        j += 1
    return count

def est_overlap_top_matching(sketch1, sketch2, seq1, seq2):
    n_matching = [get_matching_bases(i,j,seq1,seq2) for i,j in zip(sketch1, sketch2)]
#     total_match = sum([m for m in n_matching if m>=thresh])
    return n_matching

#######################################################

def get_thresh(sample_ks):
    all_ks = np.concatenate(sample_ks)
    unique_vals, counts = np.unique(all_ks, return_counts=True)
    thresh = [np.log(counts[i])-np.log(counts[i-1]) for i in range(1,len(counts))]
    thresh = [thresh[i]-thresh[i-1] for i in range(1,len(thresh))]
    thresh = [unique_vals[i+2] for i in range(1,len(thresh)) if np.sign(thresh[i])!=np.sign(thresh[i-1])][0]
    print('inflection point threshold: {}'.format(thresh))
    return thresh


def get_sample_ks(sample_size=10):
    sample_ks = []
    n = len(seqs)
    for i in np.random.choice(np.arange(n), size=sample_size, replace=False):
        for j in range(n):
            if i==j: continue
            arr = est_overlap_top_matching(sketches[i], sketches[j], seqs[i], seqs[j])
            sample_ks.append(arr)
    return sample_ks

def ks_to_stat(method, gt_path=None):
    if method == 'inflec':
        return get_stat_inflec
    elif method == 'inflec2':
        return get_stat_inflec2
    elif method == 'max':
        return get_stat_max
    elif method == 'max2':
        return get_stat_max2

def get_stat_threshs(stat_funcs, k_thresh, sample_ks, quantile):
    method_stats = [[func(ks,k_thresh) for ks in sample_ks] for func in stat_funcs]
    stat_threshs = [np.quantile(stats, quantile) for stats in method_stats]
    return stat_threshs

def get_threshs_and_funcs(methods, rc):
    sample_ks = get_sample_ks() 
    k_thresh = get_thresh(sample_ks)
    stat_funcs = [ks_to_stat(method) for method in methods]
    stat_threshs = get_stat_threshs(stat_funcs, k_thresh, sample_ks, quantile=0.5 if rc else 0.25)
    return stat_threshs, stat_funcs, k_thresh

######################################################


def estimate_pair(i,rc):
#     i,direc = x
#     with open('/home/gcgreen2/alignment/tmp/pairs.txt','a') as fh: fh.write('{}\n'.format(str(i)))
    ests = [{} for _ in range(len(stat_funcs))]
    pairs = [j%len(seqs) for j in range(i+1,i+1+len(seqs)//2)]
    for j in pairs:
        if rc:
            ks = est_overlap_top_matching(sketches[i], sketches_rc[j], seqs[i], seqs_rc[j])
        else:
            ks = est_overlap_top_matching(sketches[i], sketches[j], seqs[i], seqs[j])
        for m in range(len(stat_funcs)):
            stat = stat_funcs[m](ks,k_thresh)
            if stat > stat_threshs[m]:
                ests[m][j] = stat
    return ests

def pairwise_overlap_ests(rc):
    n = len(seqs)
    args = ((i,rc) for i in range(n))
    chunksize = int(np.ceil(n/cpu_count()/4))
    start = time()
    with Pool() as pool:
        method_ests = pool.starmap(estimate_pair, args, chunksize)
    print('pairwise alignment took {} minutes, rc={}'.format(int((time()-start)/60), str(rc)))
    return method_ests

####################################################################

def write_overlaps(aln_path, estimates, ids, m=0, modifier='w', direc='+', use_ids=False):
    print(f'# overlaps: {sum([len(estimates[i][m]) for i in range(len(estimates))])}')
    ids = ids if use_ids else [str(i) for i in range(len(ids))]
    with open(aln_path, modifier) as fh:
        for i in range(len(estimates)):
            for j in estimates[i][m]:
                line = [ids[i], ids[j], str(estimates[i][m][j]), direc]
                fh.write('\t'.join(line)+'\n')

def find_overlaps(fasta_file, aln_path, n_hash, **args):
    global seqs, seqs_rc, masks, sketches, sketches_rc, stat_threshs, stat_funcs, k_thresh
    seqs,ids = su.get_seqs(fasta_file, return_ids=True)
    seqs_rc = su.revcomp_seqs(seqs) if args['rc'] else None
    
    utils.print_banner('SKETCHING READS')
    if 'sketches_path' in args and exists(args['sketches_path']):
        sketches, sketches_rc = np.load(args['sketches_path'], allow_pickle=True)
        sketches = [sketch[:n_hash] for sketch in sketches]
        if sketches_rc is not None: sketches_rc = [sketch[:args['n_hash']] for sketch in sketches_rc] 
    else:
        masks = obtain_masks(max([len(seq) for seq in seqs]), n_hash)
        sketches = get_all_sketches(rc=False)
        sketches_rc = get_all_sketches(rc=True) if args['rc'] else None
    if 'sketches_path' in args and not exists(args['sketches_path']):
        np.save(args['sketches_path'], [sketches, sketches_rc]) 
        print(f'saved sketches at {args["sketches_path"]}')
        
    stat_threshs, stat_funcs, k_thresh = get_threshs_and_funcs(['max'], rc=args['rc'])
    
    utils.print_banner('PERFORMING PAIRWISE COMPARISON')
    method_ests = pairwise_overlap_ests(rc=False)
    write_overlaps(aln_path, method_ests, ids)
        
    if args['rc']:
        del method_ests
        method_ests = pairwise_overlap_ests(rc=True)
        write_overlaps(aln_path, method_ests, ids, modifier='a', direc='-')
                
    
    