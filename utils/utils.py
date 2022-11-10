import numpy as np
import os
from os.path import join

def print_banner(x):
    print('\n=================== {} ==================='.format(str(x)))

def print_header():
    print_banner('SuffixHash Sequence Alignment Estimator')
    for _ in range(2): print('=======================================')
        
def setup_out_dir(out_dir,  **args):
    tmp_dir = join(out_dir, 'tmp')
    os.makedirs(tmp_dir, exist_ok=True)
    save_file(join(tmp_dir, 'args.npy'), args)
        
def save_file(filepath, x):
    x = np.array(x, dtype=object)
    np.save(filepath, x, allow_pickle=True)
