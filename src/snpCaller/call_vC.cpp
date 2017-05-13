#include <algorithm>
#include <iostream>
#include <sstream>
#include <unistd.h>
//getopt
//#include <ctype.h>
//#include <unistd.h>
//Other
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <sstream>
#include <vector>
#include <list>
#include <map>
#include "gene.h"

//Get the boost interval support
#include <boost/icl/discrete_interval.hpp>
#include <boost/icl/closed_interval.hpp>
#include <boost/icl/split_interval_map.hpp>
#include <boost/lexical_cast.hpp>
using namespace boost::icl;

struct SNPCallOptions {
    SNPCallOptions()
        :min_coverage(4)
        ,calling_threshold(4)
        ,calling_min_fraction(0.01)
    { }

    int min_coverage;
    int calling_threshold;
    double calling_min_fraction;
};

#define READ_BUFFER_SIZE 10000000


//#define DEBUG 1

//TODO: make nice DEBUG switch (-d flag)
void use_debug() {
//	#define DEBUG 1
}

struct filePosition {
  unsigned long start;
  unsigned long length;
  unsigned long lineCount;
};


/*std::string usage() {
  return "Run as:\nsnpCall <reference_file> <geneAnnotation> individual_SNP_output";
}
*/
static int print_usage()
{
    fprintf(stderr, "\n");
    fprintf(stderr, "metaSNV --- metagenomic SNV caller\n\n");
    fprintf(stderr, "Version: 1.0\n");
    fprintf(stderr, "Contact: Paul Costea <costea@embl.de>,\n\t Robin Muench <rmuench@embl.de>\n\n");
    fprintf(stderr, "Usage:   snpCall [options] <stdin.mpileup> \n");
    fprintf(stderr, "Options: \n");
    fprintf(stderr, "     -f,     faidx indexed reference metagenome \n ");
    fprintf(stderr, "    -g,     gene annotation file [NULL].\n");
    fprintf(stderr, "     -i,     individual SNPs output file [NULL].\n\n");
    fprintf(stderr, "SNP definition: \n");
    fprintf(stderr, "     -c,     minimum coverage (mapped reads) per position [4]\n ");
    fprintf(stderr, "    -p,     minimum non-reference nucleotide allele frequency per position [0.01].\n");
    fprintf(stderr, "     -t,     minimum number of non-reference nucleotides per position [4].\n\n");
    //  fprintf(stderr, "	    -d,-a,	for debugging only \n\n");
    fprintf(stderr, "Note: Expecting samtools mpileup as standard input\n\n");
    return 1;
}


std::map<std::string, Genome*> mapGenomes;
std::map<std::string, filePosition> mapGenes;
//Here's where we save the genes
split_interval_map<long,GeneDef> geneIntervals;

/**
 * @brief Fast thread-safe string tokenizer.
 *
 * It ignores any spaces at the beginning of the string.
 *
 * Returns position after the token separator
 */
const char *toksplit(const char *src, /* Source of tokens */
        char tokchar, /* token delimiting char */
        char *token, /* receiver of parsed token */
        size_t lgh) /* length token can receive */
    /* not including final '\0' */
{
    if (src) {
        while (' ' == *src) src++;
        while (*src && (tokchar != *src)) {
            if (lgh) {
                *token++ = *src;
                --lgh;
            }
            src++;
        }
        if (*src && (tokchar == *src)) src++;
    }
    *token = '\0';
    return src;
}

/**
 * @brief Get a decent index of the genes for fast access and future annotation of SNPs.
 */
bool indexGenomeAndGenes(FILE* refGenome, FILE* refGenes) {
    char line[10000];
    char *tok = new char[10000];
    /*===================================
      We have the reference genomes indexed. Let's get the same for the gene definitions.
      Would be nice if they also has a sort of fai indexing, but that's not in the plan right now.
      ===================================*/

    //Get the header and remember that length
    unsigned long filePosStart = 0;
    unsigned long fileConsumed = 0;
    unsigned long lineCount = 0;
    std::string name = "";
    fgets(line,10000,refGenes);
    filePosStart = strlen(line);//This includes the \n
    filePosition p1;
    while (fgets(line,10000,refGenes)) {
        int pos = 0;
        const char* rest = toksplit(line,'\t',tok,10000);
        while (*rest) {
            if (pos == 2) {//This is the name!
                if (name.compare("") == 0) {//This is the first one
                    name = tok;
                } else if (name.compare(tok) != 0) {//Boom, this is a new one
                    p1.lineCount = lineCount;
                    lineCount = 0;
                    p1.start = filePosStart;
                    filePosStart += fileConsumed;
                    fileConsumed = 0;
                    mapGenes[name] = p1;
                    name = tok;
                    //fprintf(stderr,"%s ",name.c_str());
                }
                break;
            }
            ++pos;
            rest = toksplit(rest,'\t',tok,10000);
        }
	fileConsumed += strlen(line);
        lineCount += 1;
    }
    //Add the last one!
    p1.start = filePosStart;
    p1.lineCount = lineCount;
    mapGenes[name] = p1;

    //==================================================================
    //==================================================================

    name = "";
    std::string genome = "";
    bool skip = false;
    //Make the map
    while (fgets(line,10000,refGenome)) {//Get name, start, length
        line[strlen(line)-1] = '\0';
        if (line[0] == '>') {//New genome
            //Do we already have one?
            if (genome.length() > 0 && !skip) {//Save it!
	        mapGenomes[name] = new Genome(genome);
                //fprintf(stderr,"Loaded: %s|\n",name.c_str());
                genome = "";
            }
            name = line;
            name = name.erase(0,1);//Remove the >
            if (mapGenes.find(name) == mapGenes.end()){//This genome has no gene definitions, don't load it!
                skip = true;
            } else {
                skip = false;
            }
        } else {
            if (skip)
                continue;
            //Add to genome
            genome += line;
        }
    }
    //Add last genome
    mapGenomes[name] = new Genome(genome);

    std::cerr << "Genomes loaded!\n";

    delete[] tok;
    return true;
}

/**
 *  @brief Load and encode the genomes to save some space. Probably overkill.
 */

bool loadGenome(std::string gName, FILE* refGenes, bool* hasGenes) {
    //Drop old one
    geneIntervals.clear();

    if (mapGenes.find(gName) == mapGenes.end()) {//We don't have genes in this genome?
        *hasGenes = false;
        //And, just return
        return true;
    }
    *hasGenes = true;

    char line[10000];
    char *tok = new char[10000];

    if (mapGenomes.find(gName) == mapGenomes.end()) {
        fprintf(stderr,"Weird...%s\n",gName.c_str());
    }

    //Now, most of the read genome is useless
    //Read the genes and discard "junk" dna
    long start=0;
    long end=0;
    char strand='x';

    filePosition p = mapGenes[gName];
    //Scroll file to "start" of gene definitions
    if (fseek(refGenes,p.start,SEEK_SET) != 0) {
        fprintf(stderr,"File seek failed %ld \n",p.start);
        return false;
    }
    unsigned long lc = 0;
    std::string geneName = "";
    while (lc < p.lineCount) {
        if (!fgets(line,10000,refGenes)) {
            fprintf(stderr,"Read failed. Seeking file to %ld\n",p.start);
        }

        ++lc;
        const char* rest = toksplit(line,'\t',tok,10000);
        int pos = 0;
        while (tok) {
            if (pos == 1) {
                geneName = tok;
            }
            if (pos == 2) {//This is the name!
                if (gName.compare(tok) != 0) {
                    fprintf(stderr,"Reading wrong gene defintion for %s\n. Scafold supposed to be %s, but is %s\n.",geneName.c_str(),tok,gName.c_str());
                    break;
                }
            } else if (pos == 6) {//Start
                start = atol(tok)-1;//These are 1 based!!!!

            } else if (pos == 7) {//End
                end = atol(tok)-1;//These are 1 based!!!!

            } else if (pos == 8) {//Strand
                strand = tok[0];
                break;
            }
            ++pos;
            rest = toksplit(rest,'\t',tok,10000);
        }
        if (pos < 8) {
            // Error?
        }

        Gene g(start, end, geneName, strand);
        //Copy the seq into gene and add to genome
        if (start > end) {//goes around!
            fprintf(stderr,"This gene goes around :(.\nPretending we didn't see it.\n");
        } else {
            discrete_interval<long> gene_interval = construct<discrete_interval<long> >(start,end,interval_bounds::closed());
            GeneDef gD(g);
            geneIntervals += make_pair(gene_interval,gD);
        }
    }

    delete[] tok;
    return true;
}

typedef std::map<char, std::vector<long>> base_count_map_t;
int getSum(const base_count_map_t& bpCounts, std::string c, int sample=0) {
    int sum = 0;
    for (int i=0; i<c.length(); ++i) {
        sum += bpCounts.find(c[i])->second[sample];
    }
    return sum;
}

/**
 * @brief Return a reverse complemented copy of the input
 */

std::string revComplement(std::string codon) {
    std::string res;
    res.reserve(codon.length());
    for (int i = codon.length()-1; i>=0; --i) {
        if (codon[i] == 'A') {
            res += 'T';
        } else if (codon[i] == 'T') {
            res += 'A';
        } else if (codon[i] == 'C') {
            res += 'G';
        } else if (codon[i] == 'G') {
            res += 'C';
        }
    }
    return res;
}

std::string getCoverageString(const base_count_map_t& bpCounts, int nrSamples, std::string c) {
    std::ostringstream out;
    //Get sum for each sample
    for (int i=1; i<=nrSamples; ++i) {
        const int cov = getSum(bpCounts, c, i);
        out << cov << '|';
    }
    std::string res = out.str();
    return res.substr(0, res.size() - 1); // remove final '|'
}

/**
 * @brief Main function
 */
int main(int argc, char** argv) {

    FILE* genomes = NULL;
    FILE* genes = NULL;
    FILE* individualFile = NULL;

    int index;
    int c;

    opterr = 0;

    char* line = new char[READ_BUFFER_SIZE];

    base_count_map_t bpCounts;
    SNPCallOptions options;

    while ((c = getopt (argc, argv, "hdab:f:g:i:c:p:t:")) != -1)
        switch (c)
        {
            case 'h':		// help message
                print_usage();
                return -1;
            case 'a':		// unused flag
                break;
            case 'd':         // Debug switch not working!
                if (optopt == 'd')
                    use_debug();
                break;
            case 'b':		// TODO: list of bam files for header generation [optional]
                break;
            case 'f':		//reference fasta file [required]
                genomes = fopen(optarg, "r");
                if (genomes == NULL) {
                    fprintf(stderr,"Cannot open %s\n",optarg);
                    return -1;
                }
                break;
            case 'g':		//gene annotation file [optional] (TODO format?)
                genes = fopen(optarg, "r");
                if (genes == NULL) {
                    fprintf(stderr,"Cannot open %s\n",optarg);
                    return -1;
                }
                break;
            case 'i':		// output filename [required]
                individualFile = fopen(optarg, "w");
                if (individualFile == NULL) {
                    fprintf(stderr,"Cannot open %s\n",optarg);
                    return -1;
                }
                //        printf("Writing individual SNPs to %s \n", optarg);
                break;
            case 'c':
                options.min_coverage = atol(optarg);
                break;
            case 'p':
                options.calling_min_fraction = atof(optarg);
                break;
            case 't':
                options.calling_threshold = atol(optarg);
                break;
            case '?':
                if (( optopt == 'f') || ( optopt == 'g') || ( optopt == 'i') ){
                    if ( optopt == 'f'){
                        fprintf (stderr, "Option -%c requires a reference file.\n", optopt);
                    }
                    if ( optopt == 'g'){
                        fprintf (stderr, "Option -%c requires an annotation file.\n", optopt);
                    }
                    if ( optopt == 'i'){
                        fprintf (stderr, "Option -%c requires an output filename.\n", optopt);
                    }
                } else if (isprint (optopt)){
                    fprintf (stderr, "Unknown option `-%c'.\n", optopt);
                } else {
                    fprintf (stderr, "Unknown option character `\\x%x'.\n", optopt); //e.g. umlauts
                    return 1;
                }
            default:
                abort ();
        }

    // Fish out non option arguments
    for (index = optind; index < argc; index++){
        printf ("Non-option argument %s\n", argv[index]);
        return 0;
    }



    //TODO: do the counting and generate header from bam file
    //fgets(line,READ_BUFFER_SIZE,stdin) //gets the very first line of the pileup (drops it!)

    fgets(line,READ_BUFFER_SIZE,stdin);
    unsigned int number_of_tabs = 0;
    for (int i = 0; i < strlen(line); i++){
        if ('\t' == line[i]){
            ++number_of_tabs;
        }
    }
    //+1 (last is newline) -3 (Omit first three: geneID \t SEQ \t REF \t START)
    int nrSamples = int(number_of_tabs+1-3)/3; // replacing atoi(argv[3]);


    fprintf(stderr,"Identified %d samples\n",nrSamples);
    bpCounts['.'] = std::vector<long>(nrSamples + 1);
    bpCounts[','] = std::vector<long>(nrSamples + 1);
    bpCounts['a'] = std::vector<long>(nrSamples + 1);
    bpCounts['c'] = std::vector<long>(nrSamples + 1);
    bpCounts['t'] = std::vector<long>(nrSamples + 1);
    bpCounts['g'] = std::vector<long>(nrSamples + 1);
    bpCounts['A'] = std::vector<long>(nrSamples + 1);
    bpCounts['C'] = std::vector<long>(nrSamples + 1);
    bpCounts['T'] = std::vector<long>(nrSamples + 1);
    bpCounts['G'] = std::vector<long>(nrSamples + 1);


    //Now read and index genomes and genes if reference genome and annotationfile are given
    if ( (genomes != NULL) && (genes != NULL) ) {
        fprintf(stderr,"Found reference genomes and annotation file.\nLoading Genomes...\n");
        indexGenomeAndGenes(genomes,genes);
        fclose(genomes);
    }


    std::map<std::string, char> CodonMap(codons, codons + sizeof codons / sizeof codons[0]);


    //Reset file
    //fseek(in,0,SEEK_SET);
    std::string name = "";
    int cov = 0;
    std::string num;
    int skip = 0;
    int pos = 0;
    bool genomeLoaded = false;
    while (fgets(line,READ_BUFFER_SIZE,stdin)) {
#ifdef DEBUG
        fprintf(stderr,"\nLINE READ:\n%s\n",line);
#endif
        int lLen = strlen(line);
        if (lLen == READ_BUFFER_SIZE-1) {//This line is longer than buffer!
            fprintf(stderr,"You need a bigger buffer!!!\n%d Is too little\n",lLen);
        }

        line[--lLen]='\0';
        //Reset map
        for (auto it = bpCounts.begin(); it != bpCounts.end(); ++it) {
            std::fill(it->second.begin(), it->second.end(), 0);
        }
        pos = 0;

        char *tok = new char[10000];
        const char *rest = toksplit(line,'\t',tok,10000);
#ifdef DEBUG
        fprintf(stderr,"Splitting string\n");
        fprintf(stderr,"Split1: %s\n",tok);
#endif
        int lP = 0;
        char base;
        while (*rest) {
            if (pos == 0) {//genome name
                if (name.compare(tok)!=0) {//New genome. Load it!
                    name = tok;
                    genomeLoaded = false;
                }
            } else if (pos==1) {
                //==================================================================
                //====================== IMPORTANT!!! ==============================
                lP = atol(tok)-1;//Because the positions in the pileup are 1 based!!!
                //==================================================================
            } else if (pos==2) {
                base = tok[0];
            } else if ((pos > 3) && (pos % 3 == 1)){//These are the base calls
#ifdef DEBUG
                fprintf(stderr,"\n====\n%s\n====\n",tok);
#endif
                int i = 0;
                int len = strlen(tok);
                while (i < len) {
                    switch (tok[i]) {
                        case '^':
                            //Skip next
                            ++i;
                            break;
                        case '+':
                        case '-':
                            //Skip x number of bases, as given by number after this character
                            num = "";
                            while (isdigit(tok[++i])) num += tok[i];
                            skip = atoi(num.c_str());
                            i += skip-1;
                            break;
                        case '*':
                        case '$':
                        case 'N':
                        case 'n':
                            break;
                        default:
                            ++bpCounts[tok[i]][0];
                            ++bpCounts[tok[i]][pos/3];
                            break;
                    }
                    ++i;
                }
            }
            ++pos;
#ifdef DEBUG
            fprintf(stderr,"Pos now: %d\n",pos);
#endif
            rest = toksplit(rest,'\t',tok,10000);
        }

        delete[] tok;

        cov = getSum(bpCounts, "actgACTG,.");

        if (cov < options.min_coverage) {//We don't have enough coverage
            continue;
        }
        if (getSum(bpCounts, "actgACTG") < options.calling_threshold) {//We don't have enough SNP calls
            continue;
        }

        bool hasGenes;
        //Is genome loaded?
        if (!genomeLoaded) {
            loadGenome(name,genes,&hasGenes);
            genomeLoaded = true;
        }

        std::string snps = "actg";
        long snpCount = 0;
        std::string s = "";
        std::string indiv = "";
        //Is it on a gene?
        std::string geneName = "-";
        GeneDef def = geneIntervals(lP);
        Gene g(0,0,"",'*');
        bool isInGene = false;
        if (def.hasGene()) {
            g = def.getGene();
            geneName = g.name;
            isInGene = true;
        }
        std::string oldCodon = "";
        bool write = false;
        for (int i=0; i<snps.length(); ++i) {

            //Skip same base
            if (snps[i] == base) continue;
            std::ostringstream internal;
            char check[3];
            sprintf(check,"%c%c",snps[i],toupper(snps[i]));
            snpCount = getSum(bpCounts, check);//How many times do we see this position as a variant in the entire population?
            bool writeThis = false;
            std::string* sToWrite = NULL;

            if ((snpCount >= options.calling_threshold) && (snpCount >= cov*options.calling_min_fraction)) { // This is a "common" variant
                write = true;
                writeThis = true;
                sToWrite = &s;
            } else { //May be a individual variant
                for (int i=1; i<=nrSamples; ++i) {
                    int s = getSum(bpCounts, check, i);
                    if (s >= options.calling_threshold) {//We observe the variant "options.calling_threshold" times in one sample.
                        writeThis = true;
                        sToWrite = &indiv;
                        break;
                    }
                }
            }

            if (writeThis) {
                if ((hasGenes) && (isInGene)) {
                    //TODO: consider multiple genes!

                    //Get codon spanning this base
                    long codonStart;
                    int codonPosition;//i.e., the position of the base within the codon
                    if (g.start < g.end) {//Normal gene
                        codonPosition = (lP - g.start)%3;
                        codonStart = lP - codonPosition;
                        oldCodon = mapGenomes.find(name)->second->getSequence(codonStart,codonStart+2);
                    } else {//This gene is "circular"
                        fprintf(stderr,"Will not handle circular genes\n");
                        continue;
                    }
                    std::string newCodon = oldCodon;
                    newCodon[codonPosition] = toupper(snps[i]);
                    //What strand are we on?
                    if (g.strand == '-') {//Reverse complement
                        oldCodon = revComplement(oldCodon);
                        newCodon = revComplement(newCodon);
                    }
                    internal << snpCount << '|' << check[1] << '|';
                    //Check if syonymous
                    if (CodonMap[newCodon] == CodonMap[oldCodon]) {
                        internal << "S";
                    } else {
                        internal << "N";
                    }
                    internal << "[" << oldCodon << "-" << newCodon + "]|" << getCoverageString(bpCounts, nrSamples, check);
                    *sToWrite += ',' + internal.str();
                } else {//Just get the small string
                    internal << snpCount << '|' << check[1] << "|.|" << getCoverageString(bpCounts, nrSamples, check);
                    *sToWrite += ',' + internal.str();
                }
            }
        }
        //Writing OutputFile
        if (write) {
            s.erase(0,1);//There's an extra comma here

            //We now have information about all of the snp's. If multiple, sperate them by ,
            fprintf(stdout,"%s\t%s\t%ld\t%c\t%s\t%s\n",
                    name.c_str(),
                    geneName.c_str(),
                    long(lP+1),//again, this one gets converted to 1 based representation
                    base,
                    getCoverageString(bpCounts, nrSamples, "actgACTG,.").c_str(),
                    s.c_str());
        }
        if (indiv.length() != 0) {
            indiv.erase(0,1);
            if (individualFile == NULL) {
                fprintf(stderr, "Individual SNPs detected, but no individual output file specified (-i option).\n");
                individualFile = fopen("/dev/null", "w"); // This way, in the next iteration, warning above will not be printed
            } else {
                fprintf(individualFile,"%s\t%s\t%ld\t%c\t%s\t%s\n",
                        name.c_str(),
                        geneName.c_str(),
                        long(lP+1),//again, this one gets converted to 1 based representation
                        base,
                        getCoverageString(bpCounts, nrSamples, "actgACTG,.").c_str(),
                        indiv.c_str());
            }
        }
    }


    if (genes != NULL){
        fclose(genes);
    }
    if (individualFile != NULL){
        fclose(individualFile);
    }

    return 0;
}
