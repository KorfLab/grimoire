#!/usr/bin/python3

"""test harness for stuff in the toolbox"""

import sys


from fasta import Fasta

ff = Fasta(sys.argv[1])
for id in ff.ids:
	entry = ff.get(id)
	print(entry.seq)