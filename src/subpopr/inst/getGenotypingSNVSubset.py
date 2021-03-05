import sys
import glob
import os.path

hapDir = sys.argv[1] # '../costea2017_data/extra/'  *hap_positions.tab
metaSNVdir = sys.argv[2] # '../../SNP_calling/SNPs_best_split_?'

print("Getting subspecies genotyping info from: "+hapDir+'/*hap_positions.tab')
print("Getting SNV for all data from raw SNV calls: "+metaSNVdir+'/snpCaller/called_SNPs*')

fileDictionary = dict()
positionDictionary = dict()

if(len(glob.glob(hapDir+'/*hap_positions.tab')) < 1):
  sys.exit("Error: no *hap_positions.tab files")
if(len(glob.glob(metaSNVdir+'/snpCaller/called_SNPs*')) < 1):
  sys.exit("Error: no /snpCaller/called_SNPs* files in metaSNV output directory")

for f in glob.glob(hapDir+'/*hap_positions.tab'):
    spec = os.path.basename(f).replace('_hap_positions.tab','')
    #spec = '_'.join(os.path.basename(f).split('_')[0:2]) # fails if species name has '_' in it
    #print(spec)
    if spec+'.pos' not in fileDictionary:
        fileDictionary[spec+'.pos'] = open(hapDir+"/"+spec+'.pos','w')
    inf = open(f)
    inf.readline()
    for line in inf:
        l = line.rstrip().split('\t')
        c = l[1].split(':')
        #code = c[0]+':'+c[1]+':'+c[3]
        code = c[0]+':'+c[2] # ref seq ID : position
        if code not in positionDictionary:
            positionDictionary[code] = []
        if fileDictionary[spec+'.pos'] not in positionDictionary[code]:
            positionDictionary[code].append(fileDictionary[spec+'.pos'])

if len(positionDictionary) == 0:
  sys.exit("Error: no parse-able data in "+hapDir+"/*hap_positions.tab files")

#Now go through all of the SNPs and get these lines out
for a in glob.glob(metaSNVdir+'/snpCaller/called_SNPs*'):
#    print(a)
    for line in open(a):
        l = line.split('\t')
        code = l[0]+':'+l[2] # reference seq ID : postion 
        if code in positionDictionary:
            for ff in positionDictionary[code]:
                ff.write(line)






