import argparse
import sys

import grimoire.toolbox as toolbox
import grimoire.sequence as sequence
import grimoire.decode as decode
import grimoire.hmm as hmm

parser = argparse.ArgumentParser(description='2 phase gene finder')
parser.add_argument('--fasta', required=True, type=str,
	metavar='<file>', help='fasta file to decode')
parser.add_argument('--gene', required=True, type=str,
	metavar='<file>', help='gene hmm')
parser.add_argument('--mRNA', required=True, type=str,
	metavar='<file>', help='mRNA hmm')
parser.add_argument('--s1', required=True, type=int,
	metavar='<file>', help='samples for first pass')
parser.add_argument('--s2', required=True, type=int,
	metavar='<file>', help='samples for second pass')
arg = parser.parse_args()

m1 = hmm.HMM.read(arg.gene)
m1.convert2log()
m2 = hmm.HMM.read(arg.mRNA)
m2.convert2log()
fasta = toolbox.FASTA_stream(arg.fasta)
for entry in fasta:
	dna = sequence.DNA(seq=entry.seq, name=entry.id)
	v = decode.Viterbi(model=m1, dna=dna)
	p = v.generate_path()
	seq = []
	for f in p.features():
		if f.type == 'EXON':
			seq.append(f.seq_str())
	if len(seq) > 0:
		tx = sequence.DNA(seq=''.join(seq), name=dna.name)
		v = decode.Viterbi(model=m2, dna=tx)
		p = v.generate_path()
		for f in p.features():
			print(f)
	
