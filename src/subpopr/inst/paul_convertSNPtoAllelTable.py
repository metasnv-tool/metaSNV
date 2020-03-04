import sys

inf = open(sys.argv[1],'r')
outFile = open(sys.argv[1]+".freq",'w')
minDepth = int(sys.argv[2]) # was hard coded as 5

for line in inf:
    c = line.rstrip().split('\t')
    Id = c[0] + ':' +c[1]+':'+ c[2]
    cov = c[4].split('|') # this is the first set of "|"-delimited numbers
    for snp in c[5].split(','):
        s = snp.split('|')
        base = s[1]
        line = Id + ':' + base
        for i in range(3,len(s)):
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
