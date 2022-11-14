
# SuffixHash

## Sequence Alignment via Variable-Length Substring Matching

For help running suffixhash/minhash, run the following in SuffixHash directory:

`python3 run_method.py -h`

Example on NCTC1080 dataset:

`python3 run_module.py --out ../out --method suffixhash --fasta data/NCTC1080_reads.fasta --n_hash 500`

The groundtruth and output alignmnet files have the following column format

1. First read index (zero-based)
2. Second read index (zero-based)
3. Alignment size in base-pairs
4. Second read alignment orientation (+/-)
