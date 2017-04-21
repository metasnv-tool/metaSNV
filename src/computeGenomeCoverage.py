import sys

# We assume that contings are sorted within the genome, so that a contig won't
# show up randomly in the file, but together with its friends.

cov = open(sys.argv[1],'r')
xcov = open(sys.argv[2],'r')

genomeMap = dict()
#Drop cov header
cov.readline()

while True:
    covL = cov.readline()
    xcovL = xcov.readline()
    if not xcovL:
        break
    covL = covL.split('\t')
    xcovL = xcovL.split('\t')
    name = covL[0]
    namex = xcovL[0]
    if name != namex:
        print("Mismatch in names {} != {}".format(name, namex))

    taxId = name.split('.')[0]
    if taxId not in genomeMap:
        genomeMap[taxId] = [0.0, 0.0, 0.0, 0.0]

    #Add the length
    genomeMap[taxId][0] += int(covL[1])

    #Add the average coverage weighted with the length
    genomeMap[taxId][1] += float(covL[2]) * int(covL[1])

    #Add the number of bases covered at 1x
    genomeMap[taxId][2] += int(xcovL[2])

    #Add the number of bases covered at 2x
    genomeMap[taxId][3] += int(xcovL[3])
cov.close()
xcov.close()


with open(sys.argv[3],'w') as outf:
    #Write header
    outf.write('TaxId\tAverage_cov\tPercentage_1x\tPercentage_2x\n')
    #Print the average coverage over the taxId's
    for k in genomeMap:
        outf.write('%s\t%f\t%f\t%f\n'%(k,
                                        genomeMap[k][1]/genomeMap[k][0],
                                        genomeMap[k][2]/genomeMap[k][0]*100,
                                        genomeMap[k][3]/genomeMap[k][0]*100))

