import numpy as np
import argparse
import pyfastx
import os
#from joblib import Parallel, delayed

tab = str.maketrans("ACTG", "TGAC") 
def reverse_comp(seq):
    return seq.translate(tab)[::-1]

def count_unitigs(kmerfile, unitig_index, unitig_count, outdir):
    # check if we've already processed this file
    prefix = os.path.basename(kmerfile).split(".")[0]
    if os.path.exists(outdir + prefix + ".csv"):
        return
    
    # iterate through unitigs
    # the count for each unitig is the kmer with the maximum value
    counts = np.zeros(unitig_count+1, dtype=int)

    # load kmers
    for name, seq in pyfastx.Fasta(kmerfile, build_index=False):
        if seq in unitig_index:
            index = unitig_index[seq]
            counts[index] = max(counts[index], int(name))
        else:
            rev = reverse_comp(seq)
            if rev in unitig_index:
                index = unitig_index[reverse_comp(seq)]
                counts[index] = max(counts[index], int(name))
            # else:
            #     print("missing!")
            #     sys.exit()
        

    with open(outdir + prefix + ".csv", 'w') as outfile:
        outfile.write(prefix + "\n")
        for c in counts:
            outfile.write(str(c) + "\n")

    return


def main():
    
    parser = argparse.ArgumentParser(description=(
        'Summarise counts of unitigs by aligning kmers from dsk output. ' +
        'Assumes the header of the each kmer contains the count of that ' +
        'kmer in the original dataset.'))

    parser.add_argument("--unitig",
                        dest="unitigfile",
                        help="the location of a fasta file containing unitigs.",
                        type=str)

    parser.add_argument("--kmers",
                        dest="input_files",
                        required=True,
                        help=("location of kmer files."),
                        type=str,
                        nargs='+')

    parser.add_argument("-o",
                        "--out",
                        dest="outdir",
                        required=True,
                        help="name of output file.",
                        type=str)

    parser.add_argument("--quiet",
                        dest="quiet",
                        help="supress additional output",
                        action='store_true',
                        default=False)

    parser.add_argument("-t",
                        "--threads",
                        dest="n_cpu",
                        help="number of threads to use (default=1)",
                        type=int,
                        default=1)

    parser.add_argument("-K",
                        "--kmersize",
                        dest="K",
                        help="number of threads to use (default=1)",
                        type=int)

    args = parser.parse_args()
    args.outdir = os.path.join(args.outdir, "")

    files = []
    if len(args.input_files)==1:
        with open(args.input_files[0], 'r') as infile:
            for line in infile:
                files.append(line.strip())
        args.input_files = files

    # First read in unitigs to a list
    if not args.quiet: print("reading in unitigs...")
    unitig_count=0
    unitig_index = {}
    for name, seq in pyfastx.Fasta(args.unitigfile, build_index=False):
        for p in range(len(seq)-args.K+1):
            unitig_index[seq[p:p+args.K]]=unitig_count
        unitig_count+=1
    if not args.quiet: print("reading completed...")

    # Align each set of kmers to the contigs using exact matches
    if not args.quiet: print("counting unitigs...")

    for kfile in args.input_files:
        count_unitigs(kfile, unitig_index, unitig_count, args.outdir)

    #Parallel(n_jobs=args.n_cpu)(
    #    delayed(count_unitigs)(kfile, unitigs, args.outdir
    #        ) for kfile in args.input_files)

    return


if __name__ == '__main__':
    main()
