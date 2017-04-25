/*
    qaTools - Just more qa tools.
    Copyright (C) 2011  P. Costea(paul.igor.costea@scilifelab.se)

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <time.h>
#include <string>
#include <map>
#include <list>
#include <getopt.h>

#include <htslib/sam.h>
#include <htslib/khash.h>

#include "radix.h"

KHASH_SET_INIT_STR(rg)
typedef khash_t(rg) *rghash_t;

typedef struct
{
  int doMedian,maxCoverage,minQual,maxInsert,windowSize;
  bool spanCov,silent;
  uint32_t subsam_seed;
  double subsample;
  FILE* detailed,*profile,*specific,*regionDef;
}Options;

typedef struct
{
  int start,end;
  std::string alias;
}Interval;


#define MIN(x,y) \
  ((x) < (y)) ? (x) : (y)

#define EXIT_IF_NULL(P) \
  if (P == NULL) {\
  fprintf(stderr,"NULL pointer error in line %d, file(%s)\n",__LINE__,__FILE__);\
    return 1;\
  }

std::map<std::string,std::list<Interval> > iMap;

/**
 * Check if read is properly mapped
 * @return true if read mapped, false otherwise
 */
static bool is_mapped(const bam1_core_t *core)
{

  if (core->flag&BAM_FUNMAP) {
    return false;
  }

  return true;
}

/**
 * Print usage instructions
 */
static int print_usage()
{
  fprintf(stderr, "\n");
  fprintf(stderr, "Version: 1.5\n");
  fprintf(stderr, "Contact: Paul Costea <paul.igor.costea@scilifelab.se>\n\n");
  fprintf(stderr, "Usage:   qaCompute [options] <in.bam/sam> <output.out>\n");
  fprintf(stderr, "Options: \n");
  fprintf(stderr, "         -m            Also compute median coverage\n");
  fprintf(stderr, "         -q            Quality threshold. (min quality to consider) [1].\n");
  fprintf(stderr, "         -d            Print per-chromosome histogram [<output.out>.detail]\n");            
  fprintf(stderr, "         -p [INT]      Print coverage profile at INT window size to bed file [<output.out>.profile] [50000]\n");
  fprintf(stderr, "         -x [STR]      Print coverage over the features in STR. Should be of format: Chr<TAB>Start<TAB>End<TAB>Alias\n");
  fprintf(stderr, "         -a [FLOAT]    Subsample reads with probability [FLOAT]. [1.0]\n");
  fprintf(stderr, "         -i            Silent.Don't print too much stuff!\n");
  fprintf(stderr, "         -s [INT]      Compute 'span coverage' rather than base coverage, limiting insert size to INT. -1 -> consider all!\n");
  fprintf(stderr, "         -c [INT]      Maximum coverage to consider in histogram [30]\n");
  fprintf(stderr, "         -h [STR]      Use header from specified file. .sam OR .bam (Header must match the info in your input file. Otherwise, output is meaningless!)\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Note: Input file should be sorted\n\n");
  return 1;
}

static void specific_print_cov(FILE* outputFile, int* data, char* name, const uint32_t chrSize)
{
  //Are we at all interested in this contig?
  std::map<std::string,std::list<Interval> >::iterator it = iMap.find(std::string(name));
  if (it != iMap.end()) {//We apparently care
	//For each thing in this list, compute the average
	std::list<Interval>::iterator i = it->second.begin();
	uint64_t covSum = 0;
	uint32_t c;
	while (i != it->second.end()) {
	  covSum = 0;
	  if (i->end > chrSize) {
	    fprintf(stderr,"Your intervals seem to hang over the edge of the contigs. will be segfaulting now\n");
	  }
	  for (c=i->start; c<=i->end; ++c) {
	    covSum += data[c];
	  }
	  fprintf(outputFile,"%s\t%4.5f\n",i->alias.c_str(),(double)covSum/(i->end-i->start+1));
	  ++i;
	}
	//Now that we're done with this, we should really drop it
	iMap.erase(it);
  }
}

static void compute_print_cov(FILE* outputFile, Options userOpt, int* data, char* name,const uint32_t chrSize, uint64_t* coverageHist,const int currentTid)
{
  //clock_t start = clock();
  int32_t covVal = 0;
  uint64_t covSum = 0;
  uint32_t i;
  uint64_t wSum = 0;
  //Histogram vector.
  uint64_t* localCoverageHist = NULL;
  if (userOpt.detailed) {
    //Allocate
    localCoverageHist = new uint64_t[userOpt.maxCoverage+1];
    //Clear
    memset( localCoverageHist, 0, (userOpt.maxCoverage+1)*sizeof(uint64_t));
  }

  //Go through chromosome and count avarage covarage.
  for (i=0; i<chrSize; ++i){
    covVal += data[i];
    if (covVal < 0) {//int overrun!? Silly check!
      fprintf(stderr,"Probably really good coverage, since variables overrun!\n");
    }
    //This will be sorted later.
    data[i] = covVal;
    uint64_t prev = covSum;
    covSum += covVal;
    if (prev > covSum) {
      fprintf(stderr,"That's a big sum: %lu...\n",prev);
    }
    //Add value to histogram
    if (covVal > userOpt.maxCoverage) {
      ++coverageHist[userOpt.maxCoverage];
      if (localCoverageHist)
	++localCoverageHist[userOpt.maxCoverage];
    } else {
      ++coverageHist[covVal];
      if (localCoverageHist)
	++localCoverageHist[covVal];
    }

  }

  //Do we want specific intervals?
  if (userOpt.specific != NULL) {
    specific_print_cov(userOpt.specific,data,name,chrSize);
  }

  //Print coverage profile?
  if (userOpt.profile != NULL) {
    wSum = data[0];
    for (i=1; i<chrSize; ++i) {
      wSum += data[i];
      if (i % userOpt.windowSize == 0) {
	//Print to detailed file
	fprintf(userOpt.profile,"%s\t%d\t%d\t%4.5f\n",name,i-userOpt.windowSize+1,i,(double)wSum/userOpt.windowSize);
	wSum = 0;
      }
    }
    if ((i-1) % userOpt.windowSize != 0) {//Write last interval!
      fprintf(userOpt.profile,"%s\t%d\t%d\t%4.5f\n",name,i-(i%userOpt.windowSize)+1,i,(double)wSum/(i%userOpt.windowSize));
    }
  }

  if (userOpt.doMedian)
    //Sort entireChr
    radix_sort(data, chrSize);

  if (localCoverageHist) {//Print details!
    fprintf(userOpt.detailed,"%s\t%d\t",name,chrSize);
    for (int i=1; i<=userOpt.maxCoverage; ++i) {
        uint64_t coverage = 0;
        //All that has been covered i, had been covered i+1, i+2 and so on times. Thus, do this addition                                              
	for (int x = i; x<=userOpt.maxCoverage; ++x) coverage += localCoverageHist[x];
        fprintf(userOpt.detailed, "%d\t", int(coverage));
    }
    fprintf(userOpt.detailed,"\n");
    fflush(userOpt.detailed);
    //Clean histogram!
    delete[] localCoverageHist;
    localCoverageHist = NULL;
  }

  //Printout avarage coverage over this chrom
  if (!userOpt.silent) {
    fprintf(stdout,"Coverage sum %lu ! \n", covSum);
    fprintf(stdout,"Average coverage over %s : %3.2f\n", name, (double)covSum / chrSize);
    if (userOpt.doMedian)
      fprintf(stdout,"Median coverage over %s : %d\n", name, data[chrSize/2]);
  }
  if (userOpt.doMedian)
    fprintf(outputFile, "%s\t%d\t%3.5f\t%d\n", name, chrSize, (double)covSum / chrSize, data[chrSize/2]);
  else
    fprintf(outputFile, "%s\t%d\t%3.5f\n", name, chrSize, (double)covSum / chrSize);

  //clock_t end = clock();
  //printf("time to compute this: %3.2f\n", (end-start)/CLOCKS_PER_SEC);
}

/**
 * Print 0 coverage over some id's
 */
void printSkipped(FILE* outputFile, Options userOpt, bam_hdr_t* head, int start, int end)
{
  for (int i = start; i < end; ++i) {
	  if (!userOpt.silent) {
		  printf("Computing %s of size %u... \n",head->target_name[i],head->target_len[i]);
		  printf("Coverage sum %d ! \n", 0);
		  printf("Average coverage over %s : %3.5f\n", head->target_name[i], 0.0);
		  if (userOpt.doMedian)
			  printf("Median coverage over %s : %d\n", head->target_name[i], 0);
	  }
	  if (userOpt.doMedian)
		  fprintf(outputFile, "%s\t%d\t%3.5f\t%d\n", head->target_name[i],head->target_len[i], 0.0, 0);
	  else
		  fprintf(outputFile, "%s\t%d\t%3.5f\n", head->target_name[i],head->target_len[i], 0.0);

	  if (userOpt.detailed) {//Print details!
	    fprintf(userOpt.detailed,"%s\t%d\t",head->target_name[i],head->target_len[i]);
	    for (int k=1; k<=userOpt.maxCoverage; ++k) {
	      fprintf(userOpt.detailed,"%d\t",0);
	    }
	    fprintf(userOpt.detailed,"\n");
	  }

	  //Print coverage profile?
	  if (userOpt.profile != NULL) {
	    int k;
	    for (k = 1; k < head->target_len[i]; ++k) {
	      if (k % userOpt.windowSize == 0) {
		//Print to detailed file
		fprintf(userOpt.profile,"%s\t%d\t%d\t%4.5f\n",head->target_name[i],k-userOpt.windowSize+1,k,0.0);
	      }
	    }
	    if ((k-1) % userOpt.windowSize != 0) {//Write last interval!
	      fprintf(userOpt.profile,"%s\t%d\t%d\t%4.5f\n",head->target_name[i],k-(k % userOpt.windowSize)+1,k,0.0);
	    }
	  }
  }
}

/**
 * Open a .sam/.bam file. 
 * @returns NULL is open failed.
 */
htsFile* open_alignment_file(std::string path)
{
  std::string flag = "r";
  //  if (path.substr(path.size()-3).compare("bam") == 0) {                                                                                                                                               
    //BAM file!                                                                                                                                    
    flag += "b";                                                                                                                                                                                                             
    //}
  htsFile* fp = sam_open(path.c_str(), flag.c_str());
  if (!fp) {
    fprintf(stderr, "qaCompute: Failed to open file %s\n", path.c_str());
  }
  return fp;
}

/**
 * Main of app
 */
int main(int argc, char *argv[])
{
  FILE *outputFile;
  Options userOpt;
  std::string headerFile = "";
  userOpt.doMedian = 0;
  userOpt.maxCoverage = 30;
  userOpt.windowSize = 50000;
  userOpt.profile = NULL;
  //Pointer to file for storing specific interval coverages
  userOpt.specific = NULL;
  //File that defines the regions to be computed
  userOpt.regionDef = NULL;
  userOpt.spanCov = false;
  userOpt.silent = false;
  userOpt.detailed = NULL;
  userOpt.minQual = 1;
  userOpt.maxInsert = -1;
  userOpt.subsample = -1.0;
  userOpt.subsam_seed = 0.0;
  bool doDetail = false;
  bool doProfile = false;
  std::map<std::string,std::list<Interval> > intervalMap;
  int arg;
  char *q;
  //Get args                                                                                                                                               
  while ((arg = getopt(argc, argv, "mdip:s:q:c:h:x:a:")) >= 0) {
    switch (arg) {
    case 'm': userOpt.doMedian = 1; break;
    case 'd': doDetail = true; break;
    case 'i': userOpt.silent = true; break;
    case 'q': userOpt.minQual = atoi(optarg); break;
    case 'c': userOpt.maxCoverage = atoi(optarg); break;
    case 'a': 
      if ((userOpt.subsam_seed = strtol(optarg, &q, 10)) != 0) {
	srand(userOpt.subsam_seed);
	userOpt.subsam_seed = rand();
      }
      userOpt.subsample = strtod(q, &q);
      break;
    case 'h': headerFile = optarg;
      fprintf(stdout,"Using header from %s\n",optarg);
      break;
    case 'p': userOpt.windowSize = atoi(optarg);
      doProfile = true;
      break;
    case 's': userOpt.spanCov = true;
      userOpt.maxInsert = atoi(optarg);
      fprintf(stdout,"Max insert size %d\n",userOpt.maxInsert);
      break;
    case 'x': userOpt.regionDef = fopen(optarg,"r");
	if (userOpt.regionDef == 0) {
	   fprintf(stderr,"Unable to open region definition file %s\n",optarg);
	   return -1;
	}
	//Read all of this in
	int s,e;
	char n[10000],a[10000];
	while (fscanf(userOpt.regionDef,"%s\t%d\t%d\t%s",&n,&s,&e,&a) != EOF) {
	    Interval i;
	    i.start = s; i.end=e; i.alias = a;
	    iMap[std::string(n)].push_back(i);
	}
	break;
    default:
      fprintf(stderr,"Read wrong argument %d with value %s\n",arg,optarg);
      return -1;
    }
  }

  if (argc-optind != 2) {
    print_usage();
    return 1;
  }

  bool outsideHeader = !headerFile.empty();
  htsFile * headerF = NULL;

  std::string alignFile(argv[optind]);
  htsFile *fp = open_alignment_file(alignFile);
  EXIT_IF_NULL(fp);
  bam_hdr_t* head = sam_hdr_read(fp);
  EXIT_IF_NULL(head);

  if (outsideHeader) {
    headerF = open_alignment_file(headerFile);
    ::free(head);
    head = sam_hdr_read(headerF);
  }

  if ((outputFile = fopen(argv[optind+1], "wt")) == 0) {
    fprintf(stderr, "qaCompute: Filed to create output file %s\n", argv[optind+1]);
    return 1;
  }
  if (doDetail) {//Create detailed output file.
    std::string fName = argv[optind+1];
    fName += ".detail";
    userOpt.detailed = fopen(fName.c_str(),"wt");
    if (userOpt.detailed == NULL) {
      fprintf(stderr,"qaCompute: Unable to create detailed output file %s. No details will be printed!\n",fName.c_str());
    }
    fprintf(stdout,"Printing details in %s!\n",fName.c_str());
  }
  if (doProfile) {//Create profile output file.
    std::string fName = argv[optind+1];
    fName += ".profile";
    userOpt.profile = fopen(fName.c_str(),"wt");
    if (userOpt.profile == NULL) {
      fprintf(stderr,"qaCompute: Unable to create profile output file %s. Profile will not be printed!\n",fName.c_str());
    }
  }
  if (userOpt.regionDef != NULL) {//Create specific output file
    std::string fName = argv[optind+1];
    fName += ".specific";
    userOpt.specific = fopen(fName.c_str(),"wt");
    if (userOpt.specific == NULL) {
      fprintf(stderr,"qaCompute: Unable to create specific output file %s. Specific coverage will not be printed!\n",fName.c_str());
    }
    fclose(userOpt.regionDef);
  }

    //Initialize bam entity
    bam1_t *b = bam_init1();

    //All var declarations
    uint64_t totalGenomeLength = 0;
    uint32_t unmappedReads = 0;
    uint32_t zeroQualityReads = 0;
    uint32_t totalNumberOfReads = 0;
    uint32_t totalProperPaires = 0;
    uint32_t chrSize = 0;
    uint32_t interChr = 0;
    uint32_t duplicates = 0;
    uint32_t usedReads = 0;
 
    int *entireChr = NULL;
    //Keep header for further reference
    
    //Compute genome length
    for (int i=0; i<head->n_targets; ++i) {
    	totalGenomeLength += head->target_len[i];
    }

    int32_t currentTid = -1;

    //Create "map" vector for histogram
    uint64_t* coverageHist= (uint64_t*)malloc((userOpt.maxCoverage+1)*sizeof(uint64_t)); 
    memset( coverageHist, 0, (userOpt.maxCoverage+1)*sizeof(uint64_t));

    //Write file table header
    if (userOpt.doMedian == 1)
      fprintf(outputFile, "Chromosome\tSeq_len\tAvg_Cov\tMedian_Cov\n");
    else
      fprintf(outputFile, "Chromosome\tSeq_lem\tAvg_Cov\n");

    while (sam_read1(fp, head, b) >= 0) {
      
      //uint32_t* cigar = bam1_cigar(b);
      //Get bam core.
      const bam1_core_t *core = &b->core;

      if (core == NULL) {
	//There is something wrong with the read/file
	printf("Input file is corrupt!");
	//Leak everything and exit!
	return -1;
      }

      //Are we subsampling?
      if (userOpt.subsample > 0.) {
        uint32_t k = __ac_Wang_hash(__ac_X31_hash_string(bam_get_qname(b)) ^ userOpt.subsam_seed);
        if ((double)(k&0xffffff) / 0x1000000 >= userOpt.subsample) continue;
      }

      //BAM block has been read
      if (!is_mapped(core))
	++unmappedReads;
      else {

	if (core->tid != currentTid) {
	  
	  if (core->tid == -1) {//This read is not actually mapped..ufff
            fprintf(stderr, "Read a read that has mapped flags, but isn't actually mapped: %s\nTrying to recover\n", bam_get_qname(b));
	    ++unmappedReads;
	    ++totalNumberOfReads;
            //Skip 
            continue;
          }

	  //Count coverage!
	  if (currentTid != -1) {
	    if (!userOpt.silent)
	      fprintf(stdout,"Basing coverage on %u reads\n",usedReads);
	    usedReads = 0;
	    compute_print_cov(outputFile, userOpt, entireChr, head->target_name[currentTid], chrSize, coverageHist, currentTid);
	  }

	  //Get length of next section                                                                                       
          chrSize = head->target_len[core->tid];
	  if (chrSize < 1) {//We can't have such sizes! this can't be right
	    fprintf(stderr,"%s has size %d, which can't be right!\nCheck bam header!",head->target_name[core->tid],chrSize);
	  }

	  //Done with current section.
	  //Allocate memory
	  entireChr = (int*)realloc(entireChr, (chrSize+1)*sizeof(int));
	  
	  if (entireChr == NULL) {
	    fprintf(stderr,"Allocation failed! \n");
	    return -1;
	  }
	  memset(entireChr, 0, (chrSize+1)*sizeof(int));

	  //Have we skipped some contigs/chr from header? Should print 0 coverage on them
	  if ((currentTid + 1 != core->tid) && (currentTid != -1)) {//Since this is a sorted file!
		  printSkipped(outputFile,userOpt,head,currentTid+1,core->tid);
		  //Add entire length to 0-coverage bases!
		  coverageHist[0] += head->target_len[core->tid];
	  }

	  if (currentTid == -1) {
		  currentTid = core->tid;
		  printSkipped(outputFile,userOpt,head,0,currentTid);
	  } else {
		  currentTid = core->tid;
	  }
	  if (!userOpt.silent)
		  printf("Computing %s of size %u... \n",head->target_name[core->tid],chrSize);
	
	}
	
	//If read has quality == 0, we won't count it as mapped
	if (core->qual >= userOpt.minQual) {
	 if (core->flag&BAM_FPROPER_PAIR) {
	    //Is part of a proper pair
	    ++totalProperPaires;
	  }

	 if (core->flag&BAM_FDUP) {
	   //This is a duplicate. Don't count it!.
	   ++duplicates;
	 } else {
	   if (!userOpt.spanCov) {
	     //All entries in SAM file are represented on the forward strand! (See specs of SAM format for details)    
	     ++entireChr[core->pos];
	     ++usedReads;
	     if ((uint32_t)(core->pos+core->l_qseq) >= chrSize)
	       --entireChr[chrSize-1];
	     else
	       --entireChr[core->pos+core->l_qseq];
	   } else {
	     //Computing span coverage. 
	     //Only consider first read in pair! and extend a bit to the end of the insert
	     if ((core->flag&BAM_FREAD1) //First in pair
		 && !(core->flag&BAM_FMUNMAP) /*Mate is also mapped!*/
		 && (core->tid == core->mtid) /*Mate on the same chromosome*/
		 ) {
	       int32_t start = MIN(core->pos,core->mpos);
	       int32_t end = start+abs(core->isize);
	       int32_t iSize = end-start;
	       if ((userOpt.maxInsert == -1) || (iSize <= userOpt.maxInsert)) {
		 ++entireChr[start];
		 if ((uint32_t)end >= chrSize)
		   --entireChr[chrSize-1];
		 else
		   --entireChr[end];
		 ++usedReads;
	       }
	     } else if (core->tid != core->mtid) {
	       //Count inter-chrom mates
	       ++interChr;
	     }
	   }
	 }

	} else {
	  //Count is as unmapped?
	  ++zeroQualityReads;
	}
      }

      ++totalNumberOfReads;
      
    }

    //Compute coverage for the last "chromosome"
    compute_print_cov(outputFile, userOpt, entireChr, head->target_name[currentTid], chrSize, coverageHist, currentTid);

    //Is this the last!???
    //Print all other contings with coverage 0!
    if (currentTid != head->n_targets) {//Print to that.
    	printSkipped(outputFile,userOpt,head,currentTid+1,head->n_targets);
	}

    //Print all the stuff left in the iMap if there were any specifics to be printed
    if (userOpt.specific != NULL) {
	std::map<std::string,std::list<Interval> >::iterator it = iMap.begin();
	while (it != iMap.end()) {
	    std::list<Interval>::iterator i = it->second.begin();
	    while (i != it->second.end()) {
		fprintf(userOpt.specific,"%s\t%4.5f\n",i->alias.c_str(),0.0);
		++i;
            }
	    ++it;
	}
    }

    bam_destroy1(b);
    free(entireChr);

    //fprintf(stdout,"\nDuplicates:%u \n", duplicates);

    //Print header for next table in output file
    fprintf(outputFile,"\nCov*X\tPercentage\tNr. of bases\n");

    //fprintf(stdout,"Total genome lenght %lu \n", totalGenomeLength);
    //Compute procentages of genome cover!.
    int i;
    for (i=0; i<=userOpt.maxCoverage; ++i) {
      if (i == 0) {
	//Non-covered!
	//fprintf(stdout,"%3.2f of genome has not been covered\n", (double)(coverageHist[i])/totalGenomeLength*100);
      } else {
	uint64_t coverage = 0;
	//All that has been covered i, had been covered i+1, i+2 and so on times. Thus, do this addition
	for (int x = i; x<=userOpt.maxCoverage; ++x) coverage += coverageHist[x];
	//fprintf(stdout,"%3.2f of genome has been covered at least %dX \n", (double)(coverage)/totalGenomeLength*100, i);
	fprintf(outputFile,"%d\t%3.5f\t%lu\n",i, (double)(coverage)/totalGenomeLength*100, coverage);
      }

    }

    fprintf(outputFile,"\nOther\n");

    //Printout procentage of mapped/unmapped reads                                                                                                     
    double procentageOfUnmapped = 100*((double)unmappedReads/totalNumberOfReads);
    double procentageOfZeroQuality = 100*((double)zeroQualityReads/totalNumberOfReads);
    fprintf(outputFile,"Total number of reads: %u\n", totalNumberOfReads);
    fprintf(outputFile,"Total number of duplicates found and ignored: %u\n", duplicates);
    fprintf(outputFile,"Percentage of unmapped reads: %3.5f\n", procentageOfUnmapped);
    fprintf(outputFile,"Percentage of sub-par quality mappings: %3.5f\n", procentageOfZeroQuality);
    int32_t nrOfPaires = totalNumberOfReads/2;
    double procOfProperPaires = (double)(100*(double)totalProperPaires/2)/nrOfPaires;
    fprintf(outputFile,"Number of proper paired reads: %u\n", totalProperPaires);
    fprintf(outputFile,"Percentage of proper pairs: %3.5f\n", procOfProperPaires);
    if (userOpt.spanCov) {
      fprintf(outputFile, "Number of interchromosomal pairs: %u\n",interChr);
    }

    //printf("Out of %u reads, you have %3.5f unmapped reads\n and %3.5f sub-par quality mappings\n", totalNumberOfReads ,procentageOfUnmapped, procentageOfZeroQuality);
    

    free(coverageHist);

  
    fclose(outputFile);
    if (outsideHeader) {
      hts_close(headerF);
    }
    hts_close(fp);
    
    if (userOpt.detailed)
      fclose(userOpt.detailed);

    if (userOpt.profile)
      fclose(userOpt.profile);
	
    if (userOpt.specific)
      fclose(userOpt.specific);
  
  return 0;
}
