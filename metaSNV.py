#!/usr/bin/env python
# This code is part of the metagenomic SNV calling pipeline (metaSNV)
# Helper script to initiate a new project file structure.
import argparse
from sys import stderr, exit
from os import path
from glob import glob
import os
import shutil
import subprocess
import multiprocessing

basedir = os.path.dirname(os.path.abspath(__file__))


def mkdir_p(dirname):
    '''Equivalent to 'mkdir -p' on the command line'''
    from os import makedirs
    try:
        makedirs(dirname)
    except OSError:
        pass


def create_directories(basedir):
    mkdir_p(basedir)
    for sub in ['cov',
                'bestsplits',
                'snpCaller',
                'filtered',
                'filtered/pop',
                'filtered/ind',
                'distances']:
        mkdir_p(path.join(basedir, sub))



def exit_worker(signum, frame):
    raise RuntimeError("Keyboard Interrupt")
def init_worker():
    import signal
    signal.signal(signal.SIGINT, exit_worker)

def run_sample(sample, command):
    '''Simply wraps subprocess.call and returns contextual information in order
    to provide a good error message (if necessary)'''

    ret = subprocess.call(command)
    return sample, command, ret


def compute_opt(args):
    out_dir = path.join(args.project_dir, 'cov')
    mkdir_p(out_dir)
    p = multiprocessing.Pool(args.threads, init_worker)
    results = []
    for line in open(args.all_samples):
        line = line.rstrip()
        name = path.basename(line)
        cmd = ['{}/src/qaTools/qaCompute'.format(basedir)
                ,'-c' ,'10', '-d',
                '-i', line, '{}/{}.cov'.format(out_dir, name)]
        if args.print_commands:
            print(" ".join(cmd))
        else:
            results.append(p.apply_async(run_sample, (name, cmd)))

    p.close()
    p.join()
    for r in results:
        sample, command, ret = r.get()
        if ret:
            stderr.write("Failure in sample {}".format(sample))
            stderr.write("Call to {} failed.".format(' '.join(cmd)))
            exit(1)

def get_header(args):
    use = open(args.all_samples).readline().rstrip()
    o = subprocess.check_output(["samtools","view","-H",use]).decode("utf-8")
    f = open(args.project_dir + '/bed_header','w')
    for line in o.split('\n')[1:]:
        line = line.rstrip().split('\t')
        if len(line) != 3:
            continue
        if line[0] == "@SQ":
            line[1] = line[1].replace('SN:','')
            line[2] = line[2].replace('LN:','')
            f.write(line[1]+'\t1\t'+line[2]+'\n')
    f.close()
    args.ctg_len = args.project_dir + '/bed_header'

def compute_summary(args):
    '''This information is required by metaSNV_post.py'''

    project_name = path.basename(args.project_dir)
    cov_dir = path.join(args.project_dir, 'cov')
    cov_files = glob(cov_dir + '/*.cov')

    if not cov_files:
        if not args.print_commands:
            stderr.write("Coverage files not found.\n")
        else:
            stderr.write("Coverage files not found.\nFinish running the commands printed above and then run this command again.\n")
        exit(1)
    for f in cov_files:
        cmd = ['python',
                path.join(basedir, 'src/computeGenomeCoverage.py'),
                f,
                f + '.detail',
                f + '.summary']
        subprocess.call(cmd)
    print("\nCoverage summary here: {}".format(args.project_dir))
    print("	Average vertical genome coverage: '{}/{}.all_cov.tab'".format(args.project_dir, project_name))
    print("	Horizontal genome coverage (1X): '{}/{}.all_perc.tab'".format(args.project_dir, project_name))
    print("")
    cmd = ['python',
            '{}/src/collapse_coverages.py'.format(basedir),
            args.project_dir]
    subprocess.call(cmd)

def split_opt(args):

    if args.n_splits > 100:
        stderr.write("Maximum number of splits is 100.\n")
        args.n_splits = 100

    older_files = glob(args.project_dir + '/bestsplits/*')
    if older_files:
        stderr.write("\nremoving old splits.\n")
        for f in older_files:
            os.unlink(f)

    project_name = path.basename(args.project_dir)

    print("\nCalculating best database split:")
    # usage createOptimumSplit.sh <all_cov.tab> <all_perc.tab> <geneDefinitions> <INT_NrSplits> <.outfile>
    cmd = ['python',
            '{}/src/createOptimumSplit.py'.format(basedir),
            "{}/{}.all_cov.tab".format(args.project_dir, project_name),
            "{}/{}.all_perc.tab".format(args.project_dir, project_name),
            args.ctg_len,
            str(args.n_splits),
            path.join(args.project_dir, "bestsplits", "best_split")]
    subprocess.call(cmd)


def execute_snp_call(args, snpCaller, ifile, ofile, split):
    db_ann_args = []
    if args.db_ann != '':
        db_ann_args = ['-g', args.db_ann]
    split_args = []
    if split:
        split_args = ['-l', str(split)]
    samtools_cmd = ['samtools',
                'mpileup',
                '-f', args.ref_db
                ] + split_args + [
                '-B',
                '-b', args.all_samples]
    snpcaller_cmd = [
                snpCaller, '-f', args.ref_db] + db_ann_args + [
                '-i', ifile]
    if args.print_commands:
        print(" ".join(samtools_cmd + ['|'] + snpcaller_cmd + ['>', ofile]))
    else:
        with open(ofile, 'wt') as ofile:
            samtools_call = subprocess.Popen(samtools_cmd, stdout=subprocess.PIPE)
            snpcaller_call = subprocess.Popen(snpcaller_cmd, stdin=samtools_call.stdout, stdout=ofile)
            samtools_call.stdout.close()
            return snpcaller_call.wait()


def snp_call(args):
    out_dir = path.join(args.project_dir, 'snpCaller')
    try:
        os.makedirs(out_dir)
    except:
        pass

    shutil.copy(args.all_samples,args.project_dir+'/all_samples')

    snpCaller = basedir + "/src/snpCaller/snpCall"

    indiv_out = path.join(out_dir, "indiv_called")
    called_SNP = path.join(out_dir, "called_SNPs")


## ACTUAL COMMANDLINE
# 	Note: Due to a bug in samtools v0.1.18, -Q 20 might be erroneous to use.
#	Note: Different phred score scales might be disregarded.
#	Note: If samtools > v0.1.18 is used -Q 20 filtering is highly recommended.

    threads = (args.threads if not args.print_commands else 1)
    p = multiprocessing.Pool(threads, init_worker)
    results = []
    if args.n_splits > 1:
        splits = glob('{}/bestsplits/best_split_*'.format(args.project_dir))
        for split in splits:
            results.append(p.apply_async(execute_snp_call,
                                    (args,
                                    snpCaller,
                                    '{}.{}'.format(indiv_out, path.basename(split)),
                                    '{}.{}'.format(called_SNP, path.basename(split)),
                                    split)))
        p.close()
        p.join()
        for r in results:
            v = r.wait()
            if v > 0:
                stderr.write("SNV calling failed")
                exit(1)
    else:
        v = execute_snp_call(args, snpCaller, indiv_out, called_SNP, None)
        if v > 0:
            stderr.write("SNV calling failed")
            exit(1)



def main():
    parser = argparse.ArgumentParser(description='Compute SNV profiles')
    parser.add_argument('project_dir', metavar='DIR',
                        help='A metaSNP initialized project directory')
    parser.add_argument('all_samples', metavar='FILE',
                        help='File with an input list of bam files, one file per line')
    parser.add_argument("ref_db", metavar='REF_DB_FILE',
                        help='reference multi-sequence FASTA file used for the alignments.')
    parser.add_argument('--db_ann', metavar='DB_ANN_FILE',default='',
                        help='Database gene annotation.')
    parser.add_argument('--print-commands', default=False, action='store_true',
                        help='Instead of executing the commands, simply print them out')
    parser.add_argument('--threads', metavar='INT', default=1,type=int,
                        help='Number of jobs to run simmultaneously. Will create same number of splits, unless n_splits set differently.')
    parser.add_argument('--n_splits', metavar='INT', default=1,type=int,
                        help='Number of bins to split ref into')

    args = parser.parse_args()
    args.project_dir = args.project_dir.rstrip('/')
    if not path.isfile(args.ref_db):
        stderr.write('''
ERROR:	No reference database or annotation file found!"
ERROR:	'{}' is not a file."

SOLUTION: run getRefDB.sh or set up a custom database before running metaSNP caller
        '''.format(args.ref_db))
        parser.print_help()
        exit(1)

    if not path.isfile(basedir+"/src/qaTools/qaCompute") or not path.isfile(basedir+"/src/snpCaller/snpCall"):
        stderr.write('''
ERROR:  No binaries found

SOLUTION: make\n\n'''.format(basedir))
        exit(1)

    if args.threads > 1 and args.n_splits==1:
        args.n_splits=args.threads

    if path.exists(args.project_dir) and not args.print_commands:
        stderr.write("Project directory '{}' already exists\n\n\n".format(args.project_dir))
        exit(1)

    create_directories(args.project_dir)
    compute_opt(args)
    compute_summary(args)
    get_header(args)

    if args.n_splits > 1:
        split_opt(args)
    snp_call(args)

if __name__ == '__main__':
    main()
