#!/usr/bin/python3

"""test harness for grimoire"""

import sys
import math
import copy

import sequence
import toolbox
import hmm
import genome
import decode

modelP = hmm.HMM.read('toy.hmm')
modelL = copy.deepcopy(modelP)
modelL.convert2log()
fasta = toolbox.FASTA_stream('toy.fasta')


for entry in fasta:
	dna = sequence.DNA(name=entry.id, seq=entry.seq)
	dna.check_alphabet()
	
	v1P = decode.Viterbi(model=modelP, dna=dna)
	v1P.inspect(0, len(dna.seq), 'score')
	parse = v1P.generate_path()
	print(parse.score, parse.path)
	
	v1L = decode.Viterbi(model=modelL, dna=dna)
	v1L.inspect(0, len(dna.seq), 'score')
	parse = v1L.generate_path()
	print(parse.score, parse.path)
	
	v2P = decode.StochasticViterbi(model=modelP, dna=dna, seed=1)
	parse_list = v2P.generate_paths(10)
	v2P.inspect(0, len(dna.seq), 'score')
	for parse in parse_list:
		print(parse.score, parse.path, parse.freq)
	
	v2L = decode.StochasticViterbi(model=modelL, dna=dna, seed=1)
	parse_list = v2L.generate_paths(10)
	v2L.inspect(0, len(dna.seq), 'score')
	for parse in parse_list:
		print(parse.score, parse.path, parse.freq)
		
		#for f in parse.features(dna=dna, labels=['S1', 'S2']):
		#	print(f.gff())
	
	
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
