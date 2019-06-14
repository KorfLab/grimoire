#!/usr/bin/python3

import argparse
import sys
import json

import toolbox
import hmm
from hmm import HMM

## Command line stuff ##

parser = argparse.ArgumentParser(description='HMM trainer for genes')
parser.add_argument('--fasta', required=True, type=str,
	metavar='<path>', help='path to fasta file (%(type)s)')
parser.add_argument('--gff', required=True, type=str,
	metavar='<path>', help='path to GFF file (%(type)s)')
parser.add_argument('--out', required=True, type=str,
	metavar='<path>', help='path to output file (%(type)s)')
parser.add_argument('--model', required=True, type=str,
	metavar='<model>', help='exon|cds')
parser.add_argument('--acc', required=False, type=int, default=0,
	metavar='<int>', help='acceptor context [%(default)d]')
parser.add_argument('--don', required=False, type=int, default=0,
	metavar='<int>', help='donor context [%(default)d]')
parser.add_argument('--exon', required=False, type=int, default=0,
	metavar='<int>', help='exon context [%(default)d]')
parser.add_argument('--gen', required=False, type=int, default=0,
	metavar='<int>', help='genomic context [%(default)d]')
parser.add_argument('--int', required=False, type=int, default=0,
	metavar='<int>', help='intron context [%(default)d]')
parser.add_argument('--test', action='store_true')
arg = parser.parse_args()

ff = toolbox.fasta.Fasta('data/TAIR10_1.fasta')
gf = toolbox.gff.Gff('data/TAIR10_1.gff3')

if arg.model == 'exon':
	acc_seqs = []
	don_seqs = []
	exon_seqs = []
	int_seqs = []
	chroms = gf.chroms
	for c in chroms:
		seq = ff.get(c).seq
		# exons, acceptors, donors
		for f in gf.get(type='exon', chrom=c):
			if (f.strand == '+'): # for now, just doing + strand
				acc = seq[f.beg-3:f.beg-1]
				don = seq[f.end:f.end+2]
				if acc == 'AG' and don == 'GT':
					acc_seqs.append(seq[f.beg-11:f.beg-1])
					don_seqs.append(seq[f.end:f.end+9])
					exon_seqs.append(seq[f.beg:f.end-1])
		# intergenic
		genes = []
		for f in gf.get(type='gene', chrom=c):
			genes.append({'beg':f.beg, 'end':f.end})
		genes = sorted(genes, key=lambda i: i['beg'])
		for i in range(len(genes) -1):
			if arg.test and i > 100: break # too much intergenic
			beg = genes[i]['end'] + 1
			end = genes[i+1]['beg'] -1
			int_seqs.append(seq[beg:end])
			int_seqs.append(toolbox.dna.revcomp(seq[beg:end]))
	
	acc_emits = hmm.state.train_emissions(acc_seqs, context=arg.acc)
	don_emits = hmm.state.train_emissions(don_seqs, context=arg.don)
	exon_emits = hmm.state.train_emission(exon_seqs, context=arg.exon)
	int_emits = hmm.state.train_emission(int_seqs, context=arg.gen)
	
	acc_states = hmm.state.state_factory('ACC', acc_emits)
	don_states = hmm.state.state_factory('DON', don_emits)
	exon_state = hmm.state.State(name='EXON', context=arg.exon, emits=exon_emits)
	int_state = hmm.state.State(name='GEN', context=arg.gen, emits=int_emits)

	hmm.connect_all(acc_states)
	hmm.connect2(acc_states[-1], exon_state, 1)
	hmm.connect2(exon_state, exon_state, 0.99)
	hmm.connect2(exon_state, don_states[0], 0.01)
	hmm.connect_all(don_states)
	hmm.connect2(don_states[-1], int_state, 1)
	hmm.connect2(int_state, int_state, 0.99)
	hmm.connect2(int_state, acc_states[0], 0.01)
	
	model = HMM(name='simple_exon',
		states=[int_state] + acc_states + [exon_state] + don_states)
	model.write(arg.out)
	
elif arg.model == 'cds':
	pass

elif arg.model == 'prok':
	pass

elif arg.model == 'euk':
	pass

