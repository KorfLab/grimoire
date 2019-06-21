#!/usr/bin/python3

import argparse
import sys
import re
import json

import toolbox
import genome

## Command line stuff ##

parser = argparse.ArgumentParser(description='Genome annotation tool')
parser.add_argument('--fasta', required=True, type=str,
	metavar='<path>', help='path to fasta file (%(type)s)')
parser.add_argument('--gff', required=True, type=str,
	metavar='<path>', help='path to GFF file (%(type)s)')
parser.add_argument('--validate', action='store_true')
parser.add_argument('--stats', action='store_true')
arg = parser.parse_args()

if arg.stats:
	tpg = [] # transcript per gene
	ept = [] # exons per transcript
	ipt = [] # introns per transcript
	u5pt = [] # 5'utrs per transcript
	u3pt = [] # 3'utrs per transcript
	start = {} # start sequences
	stop = {} # stop sequences
	ifs = [] # names of genes with in frame stops
	nsc = [] # names of genes with no stop codon
	splice = {} # splice sequences
	ilen = [] # intron lengths
	elen = [] # exon lengths
	u5len = [] # 5'utr lengths
	u3len = [] # 3'utr lengths

	gen = genome.Genome(fasta=arg.fasta, gff=arg.gff)
	for chr in gen.chromosomes:
		for gene in chr.genes:
			tpg.append(len(gene.transcripts))
			for tx in gene.transcripts:
				## count distributions
				ept.append(len(tx.exons))
				ipt.append(len(tx.introns))
				u5pt.append(len(tx.utr5s))
				u3pt.append(len(tx.utr3s))
				
				## start and stop
				cds = tx.cds()
				icodon = cds[0:3]
				if icodon not in start: start[icodon] = 0
				start[icodon] += 1
				tcodon = cds[len(cds)-3:len(cds)]
				if tcodon not in stop: stop[tcodon] = 0
				stop[tcodon] += 1
				
				## coding sequence
				pro = tx.protein()
				if '*' in pro[0: -1]: ifs.append(tx.id)
				if pro[-1:len(pro)] != '*': nsc.append(tx.id)
				
				## introns
				for i in tx.introns:
					iseq = i.seq()
					sd = iseq[0:2]
					sa = iseq[-2:len(iseq)]
					ss = sd + '-' + sa
					if ss not in splice: splice[ss] = 0
					splice[ss] += 1
				
				#  lengths
				for e in tx.exons: elen.append(e.end - e.beg + 1)
				for u in tx.utr5s: u5len.append(u.end - u.beg + 1)
				for u in tx.utr3s: u3len.append(u.end - u.beg + 1)
				for i in tx.introns: ilen.append(i.end - i.beg + 1)
	
	
	# do stuff...
	print(splice)