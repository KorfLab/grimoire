#!/usr/bin/env python3

import copy
from subprocess import run
import re

CMD = [
	'bin/forge',
	'--fasta',  'data/C.elegans.1percent.fasta.gz',
	'--gff3',  'data/C.elegans.1percent.gff3.gz',
	'--model', 'don',
	'--replicant']

print('\t'.join(['dlen', 'dctx', 'ictx', 'ectx', 'perf']))
for dl in range(8, 9):
	for dc in range(0, 1):
		for ic in range(0, 2):
			for ec in range(0, 2):
				cmd = copy.copy(CMD)						
				cmd.append('--don_len')
				cmd.append(str(dl))
				cmd.append('--don_ctx')
				cmd.append(str(dc))
				cmd.append('--int_ctx')
				cmd.append(str(ic))
				cmd.append('--exon_ctx')
				cmd.append(str(ec))
				p = run(cmd, capture_output=True)
				out = str(p.stdout, 'utf-8')
				m = re.search('Accuracy: (\S+)', out)
				print('\t'.join([str(dl), str(dc), str(ic), str(ec), m[1]]))

