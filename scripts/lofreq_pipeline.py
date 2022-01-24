#Takes .gff output from prokka and outputs combined gene/protien sequences from each isolate

import subprocess
import tempfile
import os
import sys
import argparse
from collections import defaultdict
import shutil

def get_options():
    

    description = 'Call variants with lofreq pipeline'
    parser = argparse.ArgumentParser(description=description,
                                     prog='lofreq_pipe')

    io_opts = parser.add_argument_group('Input/output')

    io_opts.add_argument("--r1",
                         dest="read1",
                         required=True,
                         help="first fastq read file",
                         type=str)
    io_opts.add_argument("--r2",
                         dest="read2",
                         required=True,
                         help="second fastq read file",
                         type=str)
    io_opts.add_argument("--ref",
                         dest="ref",
                         required=True,
                         help="reference fasta file",
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

def run_bwa(fq1, fq2, ref, outputdir, ncpu, quiet=False):
    outfile = outputdir + os.path.basename(fq1).replace('_1.fastq.gz', '') + '.bam'
    cleaned_outfile = outputdir + os.path.basename(fq1).replace('_1.fastq.gz', '') + '_cleaned.bam'
    cleaned_cram = outputdir + os.path.basename(fq1).replace('_1.fastq.gz', '') + '_cleaned.cram'
    out_pileup = outputdir + os.path.basename(fq1).replace('_1.fastq.gz', '') + '_pileup.txt.gz'

    # run bwa index
    if not os.path.isfile(ref + '.bwt'):
        cmd = 'bwa'
        cmd += " index " + ref
        if not quiet:
            print("running cmd: " + cmd)
        subprocess.run(cmd, shell=True, check=True)

    # run samtools index
    if not os.path.isfile(ref + '.fai'):
        cmd = 'samtools'
        cmd += " faidx " + ref
        if not quiet:
            print("running cmd: " + cmd)
        subprocess.run(cmd, shell=True, check=True)

    # run bwa
    cmd = 'bwa mem'
    cmd += " -R '@RG\\tID:group\\tSM:sample\\tPL:Illumina'"
    # @RG\tID:foo\tSM:bar
    cmd += " -M"
    cmd += " " + ref
    cmd += " " + fq1
    cmd += " " + fq2
    cmd += " -t " + str(ncpu)
    cmd += " | samtools view -S -b "
    cmd += " | samtools sort > " 
    cmd += outfile
    if not quiet:
        print("running cmd: " + cmd)
    subprocess.run(cmd, shell=True, check=True)

    # clean bam file using picard tools
    cmd = 'picard CleanSam'
    cmd += ' INPUT=' + outfile
    cmd += ' OUTPUT=' + cleaned_outfile
    if not quiet:
        print("running cmd: " + cmd)
    subprocess.run(cmd, shell=True, check=True)

    # run pileup and allelecounter
    cmd = 'pysamstats'
    cmd += ' --type variation_strand'
    cmd += ' --min-mapq 60 --min-baseq 13'
    cmd += ' --pad'
    cmd += ' --fasta ' + ref
    cmd += ' ' + cleaned_outfile
    cmd += ' | gzip'
    cmd += ' > ' +  out_pileup
    if not quiet:
        print("running cmd: " + cmd)
    subprocess.run(cmd, shell=True, check=True)

    # convert bam to cram for keeping
    cmd = 'samtools view'
    cmd += ' -T ' + ref
    cmd += ' -C'
    cmd += ' -o ' + cleaned_cram
    cmd += ' ' + cleaned_outfile
    if not quiet:
        print("running cmd: " + cmd)
    subprocess.run(cmd, shell=True, check=True)

    return(cleaned_outfile, cleaned_cram, out_pileup)

def run_lofreq(bam, ref, outputdir, ncpu, minCov=1, indels=False, quiet=False):

    output = outputdir + os.path.basename(bam).split('.')[0] + '.vcf'

    # add indel quality scores is requested
    temp_bam = bam.split('.')[0] + '_indel.bam'
    if indels:
        cmd = 'lofreq indelqual'
        cmd += ' --dindel'
        cmd += ' --ref ' + ref
        cmd += ' -o ' + temp_bam
        cmd += ' ' + bam
        if not quiet:
            print("running cmd: " + cmd)
        subprocess.run(cmd, shell=True, check=True)
    else:
        temp_bam = bam
    
    # run samtools index
    cmd = 'samtools'
    cmd += " index " + temp_bam
    if not quiet:
        print("running cmd: " + cmd)
    subprocess.run(cmd, shell=True, check=True)

    # run lofreq
    cmd = 'lofreq call-parallel'
    cmd += ' --pp-threads ' + str(ncpu)
    # cmd = 'lofreq call'
    cmd += ' -f ' + ref
    cmd += ' -C ' + str(minCov)
    cmd += ' -o ' + output
    cmd += ' -m 60 -q 13 --min-alt-bq 13'
    
    if indels:
        cmd += ' --call-indels'
    
    cmd += ' ' + temp_bam

    if not quiet:
        print("running cmd: " + cmd)
    subprocess.run(cmd, shell=True, check=True)
    
    return(output)

def run_base_recab(vcf, bam, ref, outputdir, ncpu, quiet=False):

    # create dictionary for gatk recalibration
    if not os.path.isfile(ref.replace('.fasta', '.dict')):
        cmd = 'gatk CreateSequenceDictionary'
        cmd += " -R " + ref
        if not quiet:
            print("running cmd: " + cmd)
        subprocess.run(cmd, shell=True, check=True)

    # create index of vcf for gatk recalibration
    if not os.path.isfile(vcf + '.idx'):
        cmd = "gatk IndexFeatureFile"
        cmd += " -I " + vcf
        if not quiet:
            print("running cmd: " + cmd)
        subprocess.run(cmd, shell=True, check=True)

    # recalibrate base calls
    prefix = os.path.basename(bam).split('.')[0]
    tbl = outputdir + prefix + '.table'
    bam_out = outputdir + prefix + '_recall.bam'

    cmd = 'gatk BaseRecalibrator'
    cmd += ' -R ' + ref
    cmd += ' -I ' + bam
    cmd += ' --known-sites ' + vcf
    cmd += ' --output ' + tbl

    if not quiet:
        print("running cmd: " + cmd)
    subprocess.run(cmd, shell=True, check=True)
    
    cmd = 'gatk ApplyBQSR'
    cmd += ' -R ' + ref
    cmd += ' -I ' + bam
    cmd += ' -bqsr ' + tbl
    cmd += ' --output ' + bam_out

    if not quiet:
        print("running cmd: " + cmd)
    subprocess.run(cmd, shell=True, check=True)

    return(bam_out)

def main():
    args = get_options()

    # set output and skip if already run
    prefix = os.path.basename(args.read1).replace('_1.fastq.gz', '')
    refix = os.path.splitext(os.path.basename(args.ref))[0]
    if os.path.isfile(args.output_dir + refix + '-' + prefix + '.cram'):
        print('skipping as final output file already exists')
        return

    # make sure trailing forward slash is present
    args.output_dir = os.path.abspath(args.output_dir)
    args.output_dir = os.path.join(args.output_dir, "")
    # make sure paths are full
    args.ref = os.path.abspath(args.ref)
    args.read1 = os.path.abspath(args.read1)
    args.read2 = os.path.abspath(args.read2)

    # generate temporary directory
    temp_dir = os.path.join(os.path.abspath(tempfile.mkdtemp(dir=args.output_dir)), "")

    # change working directory
    prev_dir = os.getcwd()
    os.chdir(temp_dir)

    bam_out, cleaned_cram, out_pileup = run_bwa(args.read1, args.read2, args.ref, temp_dir, args.n_cpu)

    vcf_out = run_lofreq(bam_out, args.ref, temp_dir, args.n_cpu, minCov=10)

    bam_out = run_base_recab(vcf_out, bam_out, args.ref, temp_dir, args.n_cpu)

    vcf_out = run_lofreq(bam_out, args.ref, temp_dir, args.n_cpu, minCov=3, indels=True)

    # os.rename(vcf_out, args.output_dir + refix + '-' + prefix + '.vcf')
    os.rename(cleaned_cram, args.output_dir + refix + '-' + prefix + '.cram')
    os.rename(out_pileup, args.output_dir + refix + '-' + prefix + '.pileup')

    # remove temporary directory
    os.chdir(prev_dir)
    shutil.rmtree(temp_dir)

    return



if __name__ == '__main__':
    main()


