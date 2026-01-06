# Author : Carolane Charest

# Script to extract repeat quantifications from sample produced by TRGT

import os
import sys
import gzip
import re

def extractRepeatCount(inputDir, TRID, TRMOTIF, outputFile) :

    files = os.listdir(inputDir)
    files.sort()
    repeatCountContent = []
    tmpContent = []
    outfile = open(outputFile, "wt")


    for fileCount, file in enumerate(files) :

        with gzip.open(inputDir+file, "rt", encoding="latin-1") as vcf :
            for line in vcf :

                lineContent = line.split()

                if line.startswith("chr") : 

                    if lineContent[7].startswith('TRID=' + TRID) : 

                        motifs = []
                        motifsCount = []

                        motifs = re.search("MOTIFS=(.*);", lineContent[7]).group(1).split(",")
                        motifsCount = re.search(".*:.*:.*:(.*):.*:.*:.*", lineContent[9]).group(1).split(",")

                        for i, count in enumerate(motifsCount) : motifsCount[i] = count.split("_")
					
                        for i, motif in enumerate(motifs) :

                            if TRMOTIF == motif :
                                repeatCountContent.append([file.split("_")[0], motifsCount[0][i], motifsCount[1][i]])
        
    for item in repeatCountContent:
        outfile.write("\t".join(item) + "\n")
    outfile.close()


extractRepeatCount(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])




