# This script transforms an annotation file of the .gff format into 
# an annotation file that fulfils the requirements for metaSNV.py.

import pandas as pd
import sys

# these paths need to to be adjusted
input_gff            = r"path\annotations.gff"
output_folder        = r"path\folder\\"

# optional, set rename_contigs = True if renaming contigs is necessary
rename_contigs       = False
contig_keys_path     = r"path\ids_key.txt"


# read in the .gff annotation file and create the intermediate file
check = 0
with open(input_gff,"r") as gff:
    with open(output_folder + "intermediate.txt","w+") as out:
        for line in gff:
            if check == 0:
                if line.startswith("##") == False:
                    out.write(line)
                    check = 1
                else: pass
            elif check == 1:
                if line.startswith("##") == False:
                    out.write(line)
                else: break
                    
# read out the relevant annotation data from the intermediate file
gff_data = pd.read_csv(output_folder + r"\intermediate.txt", header = None, sep="\t")

# use the contig keys file to rename the contigs
if rename_contigs:
    contig_keys = pd.read_csv(contig_keys_path, header = None, sep="\t")
    contig_keys.rename(columns ={0:"contig_name",1:"prokka_id"}, inplace=True)
    merged = pd.merge(gff_data, contig_keys, left_on=0, right_on="prokka_id", how="left")
    if len(merged) != len(gff_data): sys.exit("Unique merge is necessary.")

# create a dataframe for the metaSNV annotation format
metaSNV = pd.DataFrame(columns = ["gene_id","external_id","sequence_id",
                                  "type","gene_info","length","start","end",
                                 "strand","start_codon","stop_codon","gc"])
# transform the .gff data to fit the metaSNV annotation format
count_seqid = {}
for i in range(0,len(gff_data)):
    metaSNV.loc[i,"sequence_id"]  = gff_data.loc[i,0]
    if rename_contigs:
        metaSNV.loc[i,"sequence_id"]  = merged.loc[i,"contig_name"]
    metaSNV.loc[i,"type"]         = gff_data.loc[i,2]
    metaSNV.loc[i,"gene_info"]    = "<annotation {}>".format(gff_data.loc[i,8])
    metaSNV.loc[i,"start"]        = gff_data.loc[i,3]
    metaSNV.loc[i,"end"]          = gff_data.loc[i,4]
    metaSNV.loc[i,"length"]       = (metaSNV.loc[i,"end"] - metaSNV.loc[i,"start"]) + 1
    metaSNV.loc[i,"strand"]       = gff_data.loc[i,6]
    metaSNV.loc[i,"start_codon"]  = ""
    metaSNV.loc[i,"stop_codon"]   = ""
    metaSNV.loc[i,"gc"]           = ""
    
metaSNV = metaSNV[metaSNV.type == "CDS"] #remove non-CDS entries
metaSNV.reset_index(drop=True, inplace= True)
for i in range(0,len(metaSNV)):
    if metaSNV.loc[i,"sequence_id"] in count_seqid:
        count_seqid[metaSNV.loc[i,"sequence_id"]] += 1
    else:
        count_seqid[metaSNV.loc[i,"sequence_id"]] = 1
    metaSNV.loc[i,"external_id"]  = "{}.{}".format(metaSNV.loc[i,"sequence_id"], count_seqid[metaSNV.loc[i,"sequence_id"]])
    metaSNV.loc[i,"gene_id"]      = i+1

#output the transformed .gff data into a metaSNV annotation file 
metaSNV.to_csv(output_folder + r"\metaSNV_anntotations.txt", mode="w+", sep="\t", index = False)