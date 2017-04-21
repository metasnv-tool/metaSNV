import sys
import operator
from collections import defaultdict

def argmin(xs):
    return xs.index(min(xs))

cov = open(sys.argv[1],'r')
perc = open(sys.argv[2],'r')
genomes = open(sys.argv[3],'r')
nrToSplit = int(sys.argv[4])
outf = sys.argv[5]

genomeLen = defaultdict(int)
genomeContigs = defaultdict(list)

for line in genomes:
    genome = line.split('\t')[0].split('.')[0]
    leng = int(line.rstrip().split('\t')[2])
    genomeLen[genome] += leng
    genomeContigs[genome].append(line)

print('Found {0} genomes'.format(len(genomeLen)))

#drop header
cov.readline()
cov.readline()

coverage = dict()

for line in cov:
    #Get the sum coverage
    s = 0.0
    l = line.rstrip().split('\t')
    for i in range(1,len(l)):
        s += float(l[i])
    coverage[l[0]] = s

perc.readline()
perc.readline()


table = []
# Get an approximation of how many reads hit each genome. This is as close as
# you will get to figuring out how long running it is going to take:
for k in genomeLen.keys():
    read = genomeLen[k]*coverage[k]
    table.append((read, k))


# We greedily assign each genome in descending order to the less heavily used
# bin:
print(nrToSplit)
outputs = [open('{}_{}'.format(outf, i), 'w') for i in range(nrToSplit)]
weight = [0 for _ in range(nrToSplit)]
for (w, g) in sorted(table, reverse=True):
    pos = argmin(weight)
    weight[pos] += w
    for c in genomeContigs[g]:
        outputs[pos].write(c)
for o in outputs:
    o.close()
