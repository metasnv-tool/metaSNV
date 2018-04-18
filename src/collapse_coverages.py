#!/usr/bin/env python
import sys
from glob import glob
from collections import defaultdict
from os import path
avg_cov = defaultdict(lambda : {})
per_cov = defaultdict(lambda : {})

project_dir = sys.argv[1]
project_name = path.basename(project_dir)

bamfiles = []
for f in sorted(glob(project_dir + '/cov/*.summary')):
    bamfile = path.basename(f)[:-len('.cov.summary')]
    for i,line in enumerate(open(f)):
        if i == 0:
            continue
        tokens = line.rstrip().split()
        avg_cov[tokens[0]][bamfile] = tokens[1]
        per_cov[tokens[0]][bamfile] = tokens[2]
    bamfiles.append(bamfile)

def write_matrix(cov, header, ofile):
    out = open(ofile, 'wt')
    out.write('\t')
    out.write('\t'.join(bamfiles))
    out.write('\n')
    out.write('TaxId\t')
    out.write('\t'.join([header for _ in bamfiles]))
    out.write('\n')
    for taxid in sorted(avg_cov.keys()):
        c = cov[taxid]
        out.write('{}\t'.format(taxid))
        out.write('\t'.join([c[bf] for bf in bamfiles]))
        out.write('\n')
    out.close()

write_matrix(avg_cov, 'Average_cov', path.join(project_dir, '{}.all_cov.tab'.format(project_name)))
write_matrix(per_cov, 'Percentage_1x', path.join(project_dir, '{}.all_perc.tab'.format(project_name)))
