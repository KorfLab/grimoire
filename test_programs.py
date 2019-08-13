#!/usr/bin/env python3

import sys
import os

# kalki - convert gff3, gtf, bed12
#assert(os.system('gunzip -c data/C.elegans.1percent.gff3.gz | python3 bin/kalki --source wb > tutorial/ce1.gff') == 0)
#assert(os.system('gunzip -c data/C.elegans.1percent.gtf.gz | python3 bin/kalki --source gtf > tutorial/ce1.gtf') == 0)

# haman - split data into training and testing sets

#assert(os.system('python3 bin/haman --fasta data/A.thaliana.1percent.fasta.gz --gff data/A.thaliana.1percent.gff3.gz --out data/set --segment gene --split 2 --padding 100') == 0)


# make sure the annotation from converting gff and gtf is the same
# make sure that araport == tair10 too
