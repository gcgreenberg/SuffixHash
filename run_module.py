import sys
import os
import argparse
from os.path import join
import minhash, suffixhash
from utils import utils
  
def check_args(args):
    assert args['method'] in ['suffixhash', 'minhash'], 'method argument must be either "suffixhash" or "minhash"'
    assert args['k'] > 8 and args['k'] < 24, 'k-value must be between 8 and 24'

def add_args(args):
    args['use_ids'] = False
    args['n_bits'] = 24
    args['rc'] = True
    args['aln_path'] = join(args['out_dir'], 'aln.tsv')
    args['sketch_path'] = join(args['out_dir'], 'tmp', 'sketches.npy')
    return args

def parse_args():
    parser = argparse.ArgumentParser(description='Sequence Alignment via Variable-Length Substring Matching')
    parser.add_argument('--out', type=str, dest='out_dir', help='output directory', required=True)
    parser.add_argument('--method', type=str, help='alignment method (suffixhash/minhash)', required=True)
    parser.add_argument('--fasta', type=str, dest='fasta_file', help='path to reads fasta', required=True)
    parser.add_argument('--n_hash', type=int, default=500, help='number of hash functions to use')
    parser.add_argument('--k', type=int, help='k-value (for minhash only)', required=False, default=14)
    return vars(parser.parse_args())

if __name__ == "__main__":
    args = parse_args()
    args = add_args(args)
    check_args(args)
    utils.setup_out_dir(**args)
    utils.print_header()
    
    print(f'running {args["method"]} on {args["fasta_file"]}')
    if args['method'] == 'suffixhash':
        suffixhash.find_overlaps(**args)
    else:
        suffixhash.find_overlaps(**args)
   


    
        
