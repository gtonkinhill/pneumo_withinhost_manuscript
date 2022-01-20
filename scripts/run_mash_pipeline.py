import os
import argparse
import subprocess

MASH = "mash"
MASH_DB = "./mash/split_reference_combined_gps.msh"
SEQ_TO_GROUPS = "./mash/themisto_gpsc_groups.tab"

def get_options():
    import argparse

    description = 'Run the MASH algorithm on deep sequenced S. pneumoniae samples using the GPS project as reference.'
    parser = argparse.ArgumentParser(description=description,
                                     prog='run-mash')

    io_opts = parser.add_argument_group('Input/output')

    io_opts.add_argument(
        "--r1",
        dest="read1",
        required=True,
        type=str,
        help="forward reads (fastq)")

    io_opts.add_argument(
        "--r2",
        dest="read2",
        required=True,
        type=str,
        help="reverse reads (fastq)")

    io_opts.add_argument("-o",
                         "--out_dir",
                         dest="output_dir",
                         required=True,
                         type=str,
                         help="location of an output directory")

    parser.add_argument(
        "-t",
        "--threads",
        dest="ncpu",
        type=int,
        default=5,
        help=("number of threads to use."))

    args = parser.parse_args()

    return (args)


def main():
    args = get_options()

    # make sure trailing forward slash is present
    args.output_dir = os.path.join(os.path.abspath(args.output_dir), "")

    # get prefix and create output sub-directory
    prefix = os.path.basename(args.read1).split(".")[0]
    if prefix[-2:]=="_1": 
        prefix = prefix[:-2]
    elif prefix[-2:]=="_2":
        prefix = prefix[:-2]

    if not os.path.exists(args.output_dir + prefix):
        os.mkdir(args.output_dir + prefix)
    args.output_dir = args.output_dir + prefix + "/"

    # run mash
    print("running mash")
    
    # generate mash command and run
    cmd = MASH + " screen "
    cmd += "-p " + str(args.ncpu) + " "
    cmd += MASH_DB + " "
    cmd += args.read1 + " "
    cmd += args.read2 + " "
    cmd += "> " + args.output_dir + "mash_screen.tab"

    subprocess.run(cmd, shell=True, check=True)

    # summarise results
    print("summarising results")

    seq_to_group = {}
    with open(SEQ_TO_GROUPS, 'r') as infile:
        for line in infile:
            line=line.strip().split()
            seq_to_group[line[0]] = line[1]

    top_hits = {}
    with open(args.output_dir + "mash_screen.tab", 'r') as infile:
        for line in infile:
            line = line.strip().split()
            group = seq_to_group[line[-1]]
            if group not in top_hits:
                top_hits[group] = line + [group]
            else:
                if float(top_hits[group][0])<float(line[0]):
                    top_hits[group] = line + [group]

    top_hits = sorted([(float(top_hits[g][0]), top_hits[g]) for g in top_hits], reverse=True)

    with open(args.output_dir + "mash_screen_summary.tab", 'w') as outfile:
        for hit in top_hits:
            o = outfile.write("\t".join([prefix] + hit[1])+"\n")

    return


if __name__ == '__main__':
    main()
