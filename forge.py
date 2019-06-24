#!/usr/bin/python3

import argparse
import sys
import json

import toolbox
import hmm
import genome
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
parser.add_argument('--acc_len', required=False, type=int, default=5,
	metavar='<int>', help='acceptor length [%(default)d]')
parser.add_argument('--acc_ctx', required=False, type=int, default=0,
	metavar='<int>', help='acceptor context [%(default)d]')
parser.add_argument('--don_len', required=False, type=int, default=5,
	metavar='<int>', help='donor length [%(default)d]')
parser.add_argument('--don_ctx', required=False, type=int, default=0,
	metavar='<int>', help='donor context [%(default)d]')
parser.add_argument('--exon_ctx', required=False, type=int, default=0,
	metavar='<int>', help='exon context [%(default)d]')
parser.add_argument('--gen_ctx', required=False, type=int, default=0,
	metavar='<int>', help='genomic context [%(default)d]')
parser.add_argument('--int_ctx', required=False, type=int, default=0,
	metavar='<int>', help='intron context [%(default)d]')
parser.add_argument('--test', action='store_true')
arg = parser.parse_args()

ff = toolbox.fasta.FastaFile('data/TAIR10_1.fasta')
gf = toolbox.gff.Gff('data/TAIR10_1.gff3')


	
if arg.model == 'internal_exon':
	# skippng
	#	genes with issues
	#	genes with multiple transcripts
	#	genes with fewer than 3 exons
	gen = genome.Genome(gff=arg.gff, fasta=arg.fasta)
	acc_seqs = []
	don_seqs = []
	exon_seqs = []
	for chr in gen.chromosomes:
		for gene in chr.genes:
			tx = gene.transcripts[0]
			if gene.issues: continue
			if len(gene.transcripts) > 1: continue
			if len(tx.exons) < 3: continue
			for i in range(1, len(tx.exons) -1):
				exon_seqs.append(tx.exons[i].seq())
			for intron in tx.introns:
				iseq = intron.seq()
				acc_seqs.append(iseq[-arg.acc_len:len(iseq)])
				don_seqs.append(iseq[0:arg.don_len])
	
	acc_emits = hmm.state.train_emissions(acc_seqs, context=arg.acc_ctx)
	don_emits = hmm.state.train_emissions(don_seqs, context=arg.don_ctx)
	exon_emits = hmm.state.train_emission(exon_seqs, context=arg.exon_ctx)
	
	acc_states = hmm.state.state_factory('ACC', acc_emits)
	don_states = hmm.state.state_factory('DON', don_emits)
	exon_state = hmm.state.State(name='EXON', context=arg.exon_ctx, emits=exon_emits)
	acc_states[0].init = 1
	don_states[arg.don_len-1].term = 1

	hmm.connect_all(acc_states)
	hmm.connect2(acc_states[-1], exon_state, 1)
	hmm.connect2(exon_state, exon_state, 0.99)
	hmm.connect2(exon_state, don_states[0], 0.01)
	hmm.connect_all(don_states)
	
	model = HMM(name=arg.out, states=acc_states + [exon_state] + don_states)
	model.write(arg.out)


elif arg.model == 'splicing':
	gen = genome.Genome(gff=arg.gff, fasta=arg.fasta)
	ep_seqs = []
	en_seqs = []
	don_seqs = []
	acc_seqs = []
	intron_seqs = []
	for chr in gen.chromosomes:
		for gene in chr.genes:
			tx = gene.transcripts[0]
			if gene.issues: continue
			if len(gene.transcripts) > 1: continue
			if len(tx.exons) < 2: continue
			for i in range(len(tx.exons) - 1):
				ep_seqs.append(tx.exons[i].seq())
				en_seqs.append(tx.exons[i + 1].seq())
			for intron in tx.introns:
				iseq = intron.seq()
				don_seqs.append(iseq[0:arg.don_len])
				acc_seqs.append(iseq[-arg.acc_len:len(iseq)])
				intron_seqs.append(iseq[arg.don_len:-arg.acc_len])
	
	ep_emits = hmm.state.train_emission(ep_seqs, context=arg.exon_ctx)
	don_emits = hmm.state.train_emissions(don_seqs, context=arg.don_ctx)
	intron_emits = hmm.state.train_emission(intron_seqs, context=arg.int_ctx)
	acc_emits = hmm.state.train_emissions(acc_seqs, context=arg.acc_ctx)
	en_emits = hmm.state.train_emission(en_seqs, context=arg.exon_ctx)
	
	ep_state = hmm.state.State(name='EXP', context=arg.exon_ctx, emits=ep_emits)
	ep_state.init = 1
	don_states = hmm.state.state_factory('DON', acc_emits)
	intron_state = hmm.state.State(name='INT', context=arg.int_ctx, emits=intron_emits)
	acc_states = hmm.state.state_factory('ACC', don_emits)
	en_state = hmm.state.State(name='EXN', context=arg.exon_ctx, emits=en_emits)
	en_state.term = 1
	
	hmm.connect2(ep_state, ep_state, 0.99)
	hmm.connect2(ep_state, don_states[0], 0.01)
	hmm.connect_all(don_states)
	hmm.connect2(don_states[arg.don_len-1], intron_state, 1)
	hmm.connect2(intron_state, intron_state, 0.98)
	hmm.connect2(intron_state, acc_states[0], 0.02)
	hmm.connect_all(acc_states)
	hmm.connect2(acc_states[arg.acc_len-1], en_state, 1)
	hmm.connect2(en_state, en_state, 1)
	
	model = HMM(name=arg.out, states=[ep_state] + don_states + [intron_state] + acc_states + [en_state])
	model.write(arg.out)

elif arg.model == 'mRNA':
	gen = genome.Genome(gff=arg.gff, fasta=arg.fasta)
	u5_seqs = []
	koz_seqs = []
	atg_seqs = []
	cds_seqs = []
	u3_seqs = []
	for chr in gen.chromosomes:
		for gene in chr.genes:
			tx = gene.transcripts[0]
			if gene.issues: continue
			if len(gene.transcripts) > 1: continue
			cds = tx.cds()
			ptx = tx.primary_tx()
			start = ptx.find(cds)
			end = start + len(cds)
			if start > arg.u5_ctx + arg.koz_len:
				koz_seqs.append(ptx[start-arg.koz_len:start])
				u5_seqs.append(ptx[0 : start - arg.koz_len])
			if len(ptx) - end > arg.u3_ctx:
				u3_seqs.append(ptx[end:len(ptx)])
			atg_seqs.append(cds[0:3])
			cds_seqs.append(cds[3:len(cds)])
			
	u5_emits = hmm.state.train_emission(u5_seqs, context=arg.u5_ctx)
	koz_emits = hmm.state.train_emissions(koz_seqs, context=arg.koz_ctx)
	atg_emits = hmm.state.train_emissions(atg_seqs, context=arg.atg_ctx)
	#cds_emits = hmm.state.train_cds(cds_seqs, context=arg.cds_ctx)
	u3_emits = hmm.state.train_emission(u3_seqs, context=arg.u3_ctx)
	
	"""
	acc_states = hmm.state.state_factory('ACC', acc_emits)
	don_states = hmm.state.state_factory('DON', don_emits)
	exon_state = hmm.state.State(name='EXON', context=arg.exon_ctx, emits=exon_emits)
	acc_states[0].init = 1
	don_states[arg.don_len-1].term = 1

	hmm.connect_all(acc_states)
	hmm.connect2(acc_states[-1], exon_state, 1)
	hmm.connect2(exon_state, exon_state, 0.99)
	hmm.connect2(exon_state, don_states[0], 0.01)
	hmm.connect_all(don_states)
	
	model = HMM(name=arg.out, states=acc_states + [exon_state] + don_states)
	model.write(arg.out)
	"""

