import subprocess
import tempfile
import os
import sys
import argparse
from collections import defaultdict
import numpy as np

def get_mash_dist_medoid(input_gffs, n_cpu=1, quiet=True):

    # build mash sketch
    mash_cmd = "mash triangle"
    mash_cmd += " -p " + str(n_cpu)

    #Set up two lists of gffs to input into mash. This is a bit messy but
    # allows for an arbitary number of gffs.
    temp_output_file1 = tempfile.NamedTemporaryFile(delete=False)
    temp_output_file2 = tempfile.NamedTemporaryFile(delete=False)
    temp_output_dist = tempfile.NamedTemporaryFile(delete=False)
    temp_output_file1.close()
    temp_output_file2.close()
    temp_output_dist.close()
    
    with open(temp_output_file1.name, 'w') as outfile:
        outfile.write(input_gffs[0])
    with open(temp_output_file2.name, 'w') as outfile:
        for gff in input_gffs[1:]:
            outfile.write(gff + "\n")
    
    mash_cmd += " -l " + temp_output_file1.name
    mash_cmd += " " + temp_output_file2.name
    mash_cmd += " > " + temp_output_dist.name

    if not quiet:
        print("running cmd: " + mash_cmd)

    subprocess.run(mash_cmd, shell=True, check=True)

    # load distance matrix
    dist_mat = np.zeros((len(input_gffs), len(input_gffs)), dtype=float)
    with open(temp_output_dist.name, 'r') as infile:
        next(infile)
        for i, line in enumerate(infile):
            line = line.strip().split()
            for j, d in enumerate(line[1:]):
                dist_mat[i, j] = float(d)
                dist_mat[j, i] = dist_mat[i, j]

    # get simplified file names
    file_names = [
        os.path.splitext(os.path.basename(gff))[0] for gff in input_gffs
    ]

    # find medoid
    medoid = int(np.argmin(np.sum(dist_mat, 0)))

    # clean up
    os.remove(temp_output_file1.name)
    os.remove(temp_output_file2.name)
    os.remove(temp_output_dist.name)

    return file_names[medoid],  input_gffs[medoid]

def get_options():
    

    description = 'Finds the medoids of a set of genome clusters using Mash'
    parser = argparse.ArgumentParser(description=description,
                                     prog='find_medoids')

    io_opts = parser.add_argument_group('Input/output')
    io_opts.add_argument("-i",
                         "--input",
                         dest="input_files",
                         required=True,
                         help="input fasta files",
                         type=str,
                         nargs='+')
    io_opts.add_argument("-c",
                         "--clusters",
                         dest="clusters",
                         required=True,
                         help="tab sepearated file with 'genome cluster'",
                         type=str)
    io_opts.add_argument("-o",
                         "--out_dir",
                         dest="output_dir",
                         required=True,
                         help="location of an output directory",
                         type=str)

    # Other options
    parser.add_argument("-t",
                        "--threads",
                        dest="n_cpu",
                        help="number of threads to use (default=1)",
                        type=int,
                        default=1)

    args = parser.parse_args()
    return (args)


def main():
    args = get_options()

    # make sure trailing forward slash is present
    args.output_dir = os.path.join(args.output_dir, "")

    # collect inoput file base names
    files = {}
    for f in args.input_files:
        files[os.path.splitext(os.path.basename(f))[0]] = os.path.abspath(f)

    # collect clusters into dictionary
    clusters = defaultdict(list)
    with open(args.clusters, 'r') as infile:
        for line in infile:
            line = line.strip().split()
            clusters[line[1]].append(files[line[0]])


    # For each cluster calculate the pairwise distance matrix and find the medoid
    for cluster in clusters:
        print('cluster:', cluster)
        prefix, file = get_mash_dist_medoid(clusters[cluster], n_cpu=args.n_cpu, quiet=True)
        print('centroid prefix:', prefix)
        print('centroid file:', file)
        cmd = 'any2fasta ' + file + ' > '
        cmd += args.output_dir + cluster + '_' + prefix + '.fasta'
        subprocess.run(cmd, shell=True, check=True)

    return



if __name__ == '__main__':
    main()

