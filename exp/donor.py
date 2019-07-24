#!/usr/bin/env python3

import copy
import subprocess
import re

CMD = [
	'bin/forge',
	'--fasta',  'data/C.elegans.1percent.fasta.gz',
	'--gff3',  'data/C.elegans.1percent.gff3.gz',
	'--model', 'don',
	'--replicant']

print('\t'.join(['dlen', 'dctx', 'ictx', 'ectx', 'perf']))
for dl in range(2, 9):
	for dc in range(0, 2):
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
				p = subprocess.run(cmd, stdout=subprocess.PIPE)
				out = str(p.stdout, 'utf-8')
				m = re.search('Accuracy: (\S+)', out)
				print('\t'.join([str(dl), str(dc), str(ic), str(ec), m[1]]))

