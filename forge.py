#!/usr/bin/python3

import argparse
import sys
import json

import toolbox
import hmm

## Command line stuff ##

parser = argparse.ArgumentParser(description='Simple exon trainer (for now)')
parser.add_argument('--fasta', required=True, type=str,
	metavar='<path>', help='path to fasta file (%(type)s)')
parser.add_argument('--gff', required=True, type=str,
	metavar='<path>', help='path to GFF file (%(type)s)')
parser.add_argument('--model', required=True, type=str,
	metavar='<model>', help='exon|cds')
args = parser.parse_args()

ff = toolbox.fasta.Fasta('data/TAIR10_1.fasta')
gf = toolbox.gff.Gff('data/TAIR10_1.gff3')



if args.model == 'exon':
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
			beg = genes[i]['end'] + 1
			end = genes[i+1]['beg'] -1
			int_seqs.append(seq[beg:end])
			int_seqs.append(toolbox.dna.revcomp(seq[beg:end]))
	
	model = hmm.create_exon_hmm(
		exon_seqs = exon_seqs, exon_context=3,
		acc_seqs = acc_seqs, acc_context=1,
		don_seqs = don_seqs, don_context=1,
		int_seqs = int_seqs, int_context=3)
	
	model.dump()
	
