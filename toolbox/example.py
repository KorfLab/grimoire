#!/usr/bin/python3

"""test harness for stuff in the toolbox"""

import sys
import json


from fasta import Fasta

ff = Fasta(sys.argv[1])
for id in ff.ids:
	entry = ff.get(id)
	print(entry.seq)

from gff import Gff

gf = Gff(sys.argv[2])

print(json.dumps(gf.get(range='Chr1:7500-8500', type='CDS'), indent = 4))

