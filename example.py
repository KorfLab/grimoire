#!/usr/bin/python3

"""test harness for grimoire"""

import sys
import math

import sequence
import toolbox
import hmm
import genome
import decode

model = hmm.HMM.read('toy.hmm')
fasta = toolbox.FASTA_stream('toy.fasta')
for entry in fasta:
	dna = sequence.DNA(name=entry.id, seq=entry.seq)
	dna.check_alphabet()
	
	v1 = decode.Viterbi(model=model, dna=dna, log=False)
	print(v1.score, v1.path)
	decode.inspect_matrix(v1, 0, 9, 'score')

	v2 = decode.StochasticViterbi(model=model, dna=dna, log=False)
	decode.inspect_matrix(v2, 0, 9, 'score')
	
	v3 = decode.Posterior(model=model, dna=dna, log=True)
	decode.inspect_matrix(v3, 0, 9, 'fwd')
	
	sys.exit(1)


sys.exit(1)

## Sequence ##
s1 = sequence.DNA(name='foo', seq='ACGTAAACCCGGGTTT')
s2 = s1.revcomp()
print(s1.fasta(wrap=10))
print(s2.fasta(wrap=10))
p1 = s1.translate()
print(p1.seq)
k2 = sequence.generate_kmers(k=2)
print(k2)

## Toolbox ##
ff = toolbox.FASTA_file('toy.fasta')
for id in ff.ids:
	entry = ff.get(id)
	print(entry.id, entry.seq)

gf = toolbox.GFF_file('data/TAIR10_1.gff3')
for f in gf.get(chrom='Chr1', beg=7500, end=8500):
	print(f.type)

## HMM ##
model = hmm.HMM.read('toy.hmm')
ff = toolbox.FASTA_file('toy.fasta')

## Genome ##
gen = genome.Genome(fasta='data/TAIR10_1.fasta', gff3='data/TAIR10_1.gff3', check_alphabet=False)
