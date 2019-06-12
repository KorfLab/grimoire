#!/usr/bin/python3

"""test harness for grimoire"""

import sys
import json

### FASTA ###

from toolbox.fasta import Fasta
ff = Fasta('data/TAIR10_1.fasta')
for id in ff.ids:
	entry = ff.get(id)
	print(entry.id, entry.desc, entry.seq[0:50])

### GFF ###

from toolbox.gff import Gff
gf = Gff('data/TAIR10_1.gff3')
print(gf.types)
print(gf.chroms)
for f in gf.get(chrom='Chr1', beg=7500, end=8500):
	print(f.beg, f.end, f.type, f.strand)

