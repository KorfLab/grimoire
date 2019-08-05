#!/usr/bin/env python3

import sys
import os

# kalki - convert gff3, gtf, bed12


# haman - split data into training and testing sets

assert(os.system('python3 bin/haman --fasta data/A.thaliana.1percent.fasta.gz --gff data/A.thaliana.1percent.gff3.gz --out data/set --segment gene --split 2 --padding 100') == 0)

