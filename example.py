#!/usr/bin/python3

"""test harness for grimoire"""


import toolbox
import hmm
import genome

### FASTA ###

ff = toolbox.FASTA_file('data/TAIR10_1.fasta')
for id in ff.ids:
 	entry = ff.get(id)

### GFF ###

gf = toolbox.GFF_file('data/TAIR10_1.gff3')
for f in gf.get(chrom='Chr1', beg=7500, end=8500):
	type = f.type

## DNA ##

seq = 'ACCCCGAGGAGAGGACCCATAGGC'
rev = toolbox.revcomp(seq)
pro = toolbox.translate(seq)

## HMM ##
model = hmm.HMM.read('toy.hmm')
ff = toolbox.FASTA_file('toy.fasta')
for id in ff.ids:
	seq = ff.get(id).seq
	print(seq)
# 	path, score = viterbi.decode(model=hmm, seq=seq)
# 	print(path)
# 	print(score)
	paths, scores = hmm.stochastic(model=model, seq=seq, n=10)
	for i in range(len(paths)):
		print(paths[i])
		print(scores[i])
	break