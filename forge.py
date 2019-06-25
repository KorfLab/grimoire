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
	metavar='<path>', help='path to GFF3 file (%(type)s)')
parser.add_argument('--hmm', required=True, type=str,
	metavar='<path>', help='path to output HMM file (%(type)s)')
parser.add_argument('--sources', required=False, type=str,
	metavar='<path>', help='path to sources GFF (%(type)s)')
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
parser.add_argument('--u5_len', required=False, type=int, default=5,
	metavar='<int>', help='UTR5 length [%(default)d]')
parser.add_argument('--u5_ctx', required=False, type=int, default=0,
	metavar='<int>', help='UTR5 context [%(default)d]')
parser.add_argument('--u3_len', required=False, type=int, default=5,
	metavar='<int>', help='UTR3 length [%(default)d]')
parser.add_argument('--u3_ctx', required=False, type=int, default=0,
	metavar='<int>', help='UTR3 context [%(default)d]')
parser.add_argument('--koz_len', required=False, type=int, default=5,
	metavar='<int>', help='Kozak length [%(default)d]')
parser.add_argument('--koz_ctx', required=False, type=int, default=0,
	metavar='<int>', help='Kozak context [%(default)d]')
parser.add_argument('--atg_ctx', required=False, type=int, default=0,
	metavar='<int>', help='ATG context [%(default)d]')
parser.add_argument('--cds_ctx', required=False, type=int, default=0,
	metavar='<int>', help='CDS context [%(default)d]')
parser.add_argument('--test', action='store_true')
arg = parser.parse_args()

ff = toolbox.fasta.FastaFile('data/TAIR10_1.fasta')
gf = toolbox.gff.Gff('data/TAIR10_1.gff3')


	
if arg.model == 'internal_exon':
	gen = genome.Genome(gff=arg.gff, fasta=arg.fasta)
	acc_seqs = []
	don_seqs = []
	exon_seqs = []
	exon_len = 0
	splices = 0
	for chr in gen.chromosomes:
		for gene in chr.genes:
			tx = gene.transcripts[0]
			if gene.issues: continue
			if len(gene.transcripts) > 1: continue
			if len(tx.exons) < 3: continue
			for i in range(1, len(tx.exons) -1):
				exon_seqs.append(tx.exons[i].seq())
				exon_len += tx.exons[i].end - tx.exons[i].beg + 1
				splices += 1
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
	hmm.connect2(exon_state, exon_state, 1 - splices/exon_len)
	hmm.connect2(exon_state, don_states[0], splices/exon_len)
	hmm.connect_all(don_states)
	
	model = HMM(name=arg.hmm, states=acc_states + [exon_state] + don_states)
	model.write(arg.hmm)


elif arg.model == 'splicing':
	gen = genome.Genome(gff=arg.gff, fasta=arg.fasta)
	ep_seqs = []
	en_seqs = []
	don_seqs = []
	acc_seqs = []
	intron_seqs = []
	exon_len = 0
	splices = 0
	intron_len = 0
	for chr in gen.chromosomes:
		for gene in chr.genes:
			tx = gene.transcripts[0]
			if gene.issues: continue
			if len(gene.transcripts) > 1: continue
			if len(tx.exons) < 2: continue
			for i in range(len(tx.exons) - 1):
				ep_seqs.append(tx.exons[i].seq())
				en_seqs.append(tx.exons[i + 1].seq())
				exon_len += tx.exons[i].end - tx.exons[i].beg + 1
				splices += 1
			for intron in tx.introns:
				iseq = intron.seq()
				intron_len += len(iseq) - arg.don_len - arg.acc_len
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
	
	hmm.connect2(ep_state, ep_state, 1 - splices/exon_len)
	hmm.connect2(ep_state, don_states[0], splices/exon_len)
	hmm.connect_all(don_states)
	hmm.connect2(don_states[arg.don_len-1], intron_state, 1)
	hmm.connect2(intron_state, intron_state, 1 - splices/intron_len)
	hmm.connect2(intron_state, acc_states[0], splices/intron_len)
	hmm.connect_all(acc_states)
	hmm.connect2(acc_states[arg.acc_len-1], en_state, 1)
	hmm.connect2(en_state, en_state, 1)
	
	model = HMM(name=arg.hmm, states=[ep_state] + don_states + [intron_state] + acc_states + [en_state])
	model.write(arg.hmm)

elif arg.model == 'mRNA':
	gen = genome.Genome(gff=arg.gff, fasta=arg.fasta)
	u5_seqs = []
	koz_seqs = []
	atg_seqs = []
	cds_seqs = []
	u3_seqs = []
	u5_len = 0
	mRNAs = 0
	cds_len = 0
	for chr in gen.chromosomes:
		for gene in chr.genes:
			tx = gene.transcripts[0]
			if gene.issues: continue
			if len(gene.transcripts) > 1: continue
			cds = tx.cds()
			ptx = tx.primary_tx()
			beg = ptx.find(cds)
			end = beg + len(cds)
			if beg < arg.u5_ctx + arg.koz_len: continue
			if len(ptx) - end < arg.u3_ctx: continue
			
			u5_seqs.append(ptx[0 : beg - arg.koz_len])
			koz_seqs.append(ptx[beg-arg.koz_len:beg])
			u3_seqs.append(ptx[end:len(ptx)])
			atg_seqs.append(cds[0:3])
			cds_seqs.append(cds[3:len(cds)])
			
			u5_len += beg - arg.koz_len
			cds_len += len(cds) - 3
			mRNAs += 1
						
	u5_emits = hmm.state.train_emission(u5_seqs, context=arg.u5_ctx)
	koz_emits = hmm.state.train_emissions(koz_seqs, context=arg.koz_ctx)
	atg_emits = hmm.state.train_emissions(atg_seqs, context=arg.atg_ctx)
	cds_emits = hmm.state.train_cds(cds_seqs, context=arg.cds_ctx)
	u3_emits = hmm.state.train_emission(u3_seqs, context=arg.u3_ctx)
	
	
	u5_state = hmm.state.State(name='UTR5', context=arg.u5_ctx, emits=u5_emits)
	koz_states = hmm.state.state_factory('KOZ', koz_emits)
	atg_states = hmm.state.state_factory('ATG', atg_emits)
	cds_states = hmm.state.state_factory('CDS', cds_emits)
	u3_state = hmm.state.State(name='UTR3', context=arg.u3_ctx, emits=u3_emits)
	u5_state.init = 1
	u3_state.term = 1

	hmm.connect2(u5_state, u5_state, 1 - mRNAs/u5_len)
	hmm.connect2(u5_state, koz_states[0], mRNAs/u5_len)
	hmm.connect_all(koz_states)
	hmm.connect2(koz_states[-1], atg_states[0], 1)
	hmm.connect_all(atg_states)
	hmm.connect2(atg_states[-1], cds_states[0], 1)
	hmm.connect_all(cds_states)
	hmm.connect2(cds_states[2], cds_states[0], 1 - mRNAs / (cds_len/3))
	hmm.connect2(cds_states[2], u3_state, mRNAs / (cds_len/3))
	
	model = HMM(name=arg.hmm, states=[u5_state] + koz_states + atg_states
		+ cds_states + [u3_state])
	model.write(arg.hmm)

