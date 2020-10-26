import sys

inf = open(sys.argv[1],'r') # .pos file
outFile = open(sys.argv[1]+".freq",'w')
minDepth = int(sys.argv[2]) # was hard coded as 5

for line in inf:
    c = line.rstrip().split('\t')
    Id = c[0] + ':' +c[1]+':'+ c[2] # ref seq ID : - : position
    cov = c[4].split('|') # this is the first set of "|"-delimited numbers
    # which are the overall depth of coverage per position
    for snp in c[5].split(','): # for each variant allele 
        s = snp.split('|')
        base = s[1] # variant allele
        line = Id + ':' + base
        for i in range(3,len(s)): # samples' data start at position 3
            freq = -5
            if int(cov[i-3]) < minDepth: # if vertical coverage at this position is less than x in this sample, set the SNV freq to -1
                freq = -1
            else:
                freq = float(s[i])/int(cov[i-3]) * 100
            line += '\t' + str(freq)
        outFile.write(line+"\n")
#        print(line, file=outFile)
#        print(line)
inf.close()
outFile.close()
