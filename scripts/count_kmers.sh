
#!/bin/bash

# exit when any command fails
set -e
# keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG

# A POSIX variable
OPTIND=1         # Reset in case getopts has been used previously in the shell.

# Initialize our own variables:
output_file=""
verbose=0

function print_help { echo "Usage: count_kmers.sh -1 read1.fq -2 read2.fq -t n_threads -p out_prefix " >&2 ; }

while getopts "h?::1:2:p:t:v:" opt; do
    case "$opt" in
    h|\?)
        print_help
        exit 0 || true
        ;;
    1)  read1_file=$OPTARG
        ;;
    2)  read2_file=$OPTARG
        ;;
    p)  prefix=$OPTARG
        ;;
    t)  threads=$OPTARG
        ;;
    v)  verbose=1
        ;;
    esac    
done

shift $((OPTIND-1))

[ "${1:-}" = "--" ] && shift

echo $read1_file
echo $read2_file

#count kmers
dsk -file ${read1_file},${read2_file} -kmer-size 31 -abundance-min 3 -nb-cores $threads -out $prefix > ${prefix}_dsk.log

#convert to fasta for input to bifrost
gzip -c <(head -n -8 <(awk '{print ">"$2"\n"$1}' <(dsk2ascii -file ${prefix}.h5 -nb-cores $threads -out "/dev/stdout"))) > "${prefix}_kmers.fasta.gz"

# remove intermediate file
rm ${prefix}.h5