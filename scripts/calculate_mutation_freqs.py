import os, sys
import argparse
import pyfastx
import numpy as np
from collections import Counter, defaultdict
import ete3 as et
import itertools
import tempfile
import subprocess
import shutil

def run_iqtree(msa, outputdir, model='GTR+G', ncpu=1, quiet=False):

    prefix = os.path.splitext(os.path.basename(msa).replace('.gz', ''))[0]

    cmd = "iqtree2"
    cmd += " -s " + msa
    cmd += " --prefix " + outputdir + prefix
    cmd += " -m " + model
    cmd += " -T " + str(ncpu)
    cmd += ' --ancestral'
    cmd += ' -fast'

    if quiet:
        cmd += ' --quiet'

    if not quiet:
        print('running cmd: ', cmd)
    subprocess.run(cmd, shell=True, check=True)

    return (outputdir + prefix + '.state',
            outputdir + prefix + '.treefile')

def generate_chars(length):
    a = np.chararray(length)
    a[:] = '0'
    return(a)


def count_muts(ancs_file, treefile, msa, quiet=False):
    # initially work out how long the alignment is and set up dictionary
    for name, seq in pyfastx.Fasta(msa, build_index=False):
        seq_len = len(seq)
        break

    nodes_seqs = defaultdict(lambda: generate_chars(seq_len))

    for name, seq in pyfastx.Fasta(msa, build_index=False):
        name = name.replace('#', '_')
        name = name.replace(';', '_')
        seq = seq.upper()
        seq = seq.replace('-', 'N')
        nodes_seqs[name] = np.array(list(seq), dtype=bytes)
        # print(name)
        # print(nodes_seqs[name])
    
    # Load sequences
    with open(ancs_file, 'r') as infile:
        for i in range(9):
            next(infile)
        for line in infile:
            line = line.split()
            if (line[3]=='-nan') | (line[3]=='nan') | (line[2]=='-'):
                nodes_seqs[line[0]][int(line[1])-1] = 'N'
            else:
                nodes_seqs[line[0]][int(line[1])-1] = line[2]

    tree = et.Tree(treefile, format=1)

    # iterate through tree counting mutations
    mutations = []
    prev = tree 
    root = nodes_seqs[tree.name]
    for node in tree.traverse("preorder"):
        if node in prev.children:
            parent = prev.name
            child = node.name
            if parent not in nodes_seqs: raise ValueError('missing node!', parent)
            if child not in nodes_seqs: raise ValueError('missing node!', child)
            # print(parent, child)
            for i in np.argwhere(nodes_seqs[parent]!=nodes_seqs[child]):
                i = int(i)
                mut = (str(i), 
                    nodes_seqs[parent][i].decode('utf-8'), 
                    nodes_seqs[child][i].decode('utf-8'), 
                    nodes_seqs[parent][(i-1):(i+2)].tostring().decode('utf-8'))
                if (mut[1]!='N') and  (mut[2]!='N') and not ('N' in mut[3]):
                    mutations.append(mut)
        prev = node

    ref = Counter([root[(i-1):(i+2)].tostring().decode('utf-8') for i in range(1,len(root)-1)])

    return(mutations, ref)

def get_options():
    description = 'Runs the pairsnp and transcluster algorithms.'
    parser = argparse.ArgumentParser(description=description,
                                     prog='fasttranscluster')

    io_opts = parser.add_argument_group('Input/output')
    io_opts.add_argument(
        "--msa",
        dest="msa",
        help="Location of fasta formatted multiple sequence alignment")


    io_opts.add_argument("-o",
                         "--out_dir",
                         dest="output_dir",
                         required=True,
                         help="location of an output directory",
                         type=str)

    #Iqtree options
    iqtree_opts = parser.add_argument_group('iqtree options')
    iqtree_opts.add_argument(
        "--iqtree-model",
        dest="model",
        help="passed to the '-m' option when running iqtree (default='GTR+G')",
        type=str,
        default="GTR+G")

    # MSA filter options
    # clean = parser.add_argument_group('filter options')
    # clean.add_argument(
    #     "--min-sequences",
    #     dest="minseq",
    #     help='min sequences passing qc to proceed with analysis (default=5)',
    #     type=int,
    #     default=5)

    # clean.add_argument(
    #     "--min-occupancy",
    #     dest="minocc",
    #     help="min proportion of N's or missing sites (default=0.1)",
    #     type=float,
    #     default=0.1)

    # Other options
    parser.add_argument("-t",
                        "--threads",
                        dest="n_cpu",
                        help="number of threads to use (default=1)",
                        type=int,
                        default=1)
    parser.add_argument("--quiet",
                        dest="quiet",
                        help="suppress additional output",
                        action='store_true',
                        default=False)

    args = parser.parse_args()

    return (args)


def main():
    args = get_options()

    # set up output directory
    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir) 
    # make sure trailing forward slash is present
    args.output_dir = os.path.join(args.output_dir, "")

    # generate temporary directory
    temp_dir = os.path.join(os.path.abspath(tempfile.mkdtemp(dir=args.output_dir)), "")

    prefix = os.path.splitext(os.path.basename(args.msa).replace('.gz', ''))[0]

    ancs_file, treefile = run_iqtree(args.msa, temp_dir, model=args.model, ncpu=args.n_cpu, quiet=args.quiet)

    mutations, ref = count_muts(ancs_file, treefile, args.msa, quiet=args.quiet)

    with open(args.output_dir + prefix + "_mutations.csv", 'w') as outfile:
        for mut in mutations:
            outfile.write(prefix + ',' + ','.join(mut) + '\n')

    with open(args.output_dir + prefix + "_ref_tri_context.csv", 'w') as outfile:
        for r in itertools.product('ACGT', repeat=3):
            r = ''.join(r)
            outfile.write(prefix + ',' + r + ',' + str(ref[r]) + '\n')

    # remove temporary directory
    shutil.rmtree(temp_dir)

    return


if __name__ == '__main__':
    main()