#!/usr/bin/env python3

import copy
from subprocess import run
import re

CMD = [
	'bin/forge',
	'--fasta',  'data/C.elegans.1percent.fasta.gz',
	'--gff3',  'data/C.elegans.1percent.gff3.gz',
	'--model', 'acc',
	'--replicant']

print('\t'.join(['alen', 'actx', 'ictx', 'ectx', 'perf']))
for al in range(2, 10):
	for ac in range(0, 2):
		for ic in range(0, 2):
			for ec in range(0, 2):
				cmd = copy.copy(CMD)						
				cmd.append('--acc_len')
				cmd.append(str(al))
				cmd.append('--acc_ctx')
				cmd.append(str(ac))
				cmd.append('--int_ctx')
				cmd.append(str(ic))
				cmd.append('--exon_ctx')
				cmd.append(str(ec))
				p = run(cmd, capture_output=True)
				out = str(p.stdout, 'utf-8')
				m = re.search('Accuracy: (\S+)', out)
				print('\t'.join([str(al), str(ac), str(ic), str(ec), m[1]]))

