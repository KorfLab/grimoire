#!/usr/bin/python3

import argparse
import sys
import re
import json
import numpy
import matplotlib.pyplot as plt

import toolbox
import genome

## Command line stuff ##

parser = argparse.ArgumentParser(description='Genome annotation tool')
parser.add_argument('--fasta', required=True, type=str,
	metavar='<path>', help='path to fasta file (%(type)s)')
parser.add_argument('--gff3', required=True, type=str,
	metavar='<path>', help='path to GFF file (%(type)s)')
parser.add_argument('--stats', action='store_true')
parser.add_argument('--split', type=int, metavar='<int>',
	help='number of discrete parts to split into')
arg = parser.parse_args()

if arg.stats:
	tpg = [] # transcript per gene
	ept = [] # exons per transcript
	ipt = [] # introns per transcript
	u5pt = [] # 5'utrs per transcript
	u3pt = [] # 3'utrs per transcript
	ilen = [] # intron lengths
	elen = [] # exon lengths
	u5len = [] # 5'utr lengths
	u3len = [] # 3'utr lengths
	clen = [] # cds lengths (complete)

	gen = genome.Genome(fasta=arg.fasta, gff3=arg.gff3)
	for chr in gen.chromosomes:
		for gene in chr.features:
			#tpg.append(len(gene.children))
			for tx in gene.children:
				print(gene.id, tx.id)
				
				"""
				## count distributions
				ept.append(len(tx.exons))
				ipt.append(len(tx.introns))
				u5pt.append(len(tx.utr5s))
				u3pt.append(len(tx.utr3s))
								
				#  lengths
				for e in tx.exons: elen.append(e.end - e.beg + 1)
				for u in tx.utr5s: u5len.append(u.end - u.beg + 1)
				for u in tx.utr3s: u3len.append(u.end - u.beg + 1)
				for i in tx.introns: ilen.append(i.end - i.beg + 1)
				clen.append(len(tx.cds_str()))
				"""
	
	
#	plt.hist(clen, bins='auto')
#	plt.title('length')
#	plt.savefig('plot.png')
	

