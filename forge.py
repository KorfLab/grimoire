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
parser.add_argument('--gff3', required=True, type=str,
	metavar='<path>', help='path to GFF3 file (%(type)s)')
parser.add_argument('--hmm', required=True, type=str,
	metavar='<path>', help='path to output HMM file (%(type)s)')
parser.add_argument('--sources', required=False, type=str,
	metavar='<path>', help='path to sources GFF (%(type)s)')
parser.add_argument('--replicant', required=False, type=str,
	metavar='<path>', help='path to replicant report (%(type)s), req --sources')
parser.add_argument('--model', required=True, type=str,
	metavar='<model>', help='exon|cds')
parser.add_argument('--null_ctx', required=False, type=int, default=0,
	metavar='<int>', help='null model context [%(default)d]')
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

ff = toolbox.FASTA_file('data/TAIR10_1.fasta')
gf = toolbox.GFF_file('data/TAIR10_1.gff3')

def output_sources(features, filename):
	fp = open(filename, 'w+')
	for feature in features:
		fp.write(feature.gff())

def output_replicant(model, features, filename):
#	perf = hmm.Performance()
	for feature in features:
		pass
#		decode = model.decode(algorithm='viterbi', transform='log',
#			sequence=feature.seq())
#		perf.add_run(decode)


if arg.model == 'internal_exon':
	gen = genome.Genome(gff3=arg.gff3, fasta=arg.fasta)
	acc_seqs = []
	don_seqs = []
	exon_seqs = []
	txa = []
	exon_len = 0
	splices = 0
	for chr in gen.chromosomes:
		for gene in chr.genes:
			tx = gene.transcripts[0]
			if gene.issues: continue
			if len(gene.transcripts) > 1: continue
			if len(tx.exons) < 3: continue
			for i in range(1, len(tx.exons) -1):
				iprev = tx.introns[i-1].seq()
				acc_seqs.append(iprev[-arg.acc_len:len(iprev)])
				exon_seqs.append(tx.exons[i].seq())
				exon_len += tx.exons[i].end - tx.exons[i].beg + 1
				inext = tx.introns[i].seq()
				don_seqs.append(inext[0:arg.don_len])
				splices += 1
				if arg.sources or arg.replicant:
					txid = 'model-' + str(splices)
					parent = genome.Feature(chr, tx.introns[i-1].end - arg.acc_len + 1,
						tx.introns[i].beg + arg.don_len - 1, tx.strand, 'model', id=txid)
					parent.add_child(genome.Feature(chr, tx.introns[i-1].end - arg.acc_len + 1,
						tx.introns[i-1].end, tx.introns[i-1].strand, 'ACC', parent=txid))
					parent.add_child(genome.Feature(chr, tx.exons[i].beg, tx.exons[i].end,
						tx.exons[i].strand, 'E', parent=txid))
					parent.add_child(genome.Feature(chr, tx.introns[i].beg,
						tx.introns[i].beg + arg.don_len - 1, tx.introns[i].strand, 'DON',
						parent=txid))
					txa.append(parent)
	
	acc_emits = hmm.train_emissions(acc_seqs, context=arg.acc_ctx)
	don_emits = hmm.train_emissions(don_seqs, context=arg.don_ctx)
	exon_emits = hmm.train_emission(exon_seqs, context=arg.exon_ctx)
	
	acc_states = hmm.state_factory('ACC', acc_emits)
	don_states = hmm.state_factory('DON', don_emits)
	exon_state = hmm.State(name='EXON', context=arg.exon_ctx, emits=exon_emits)
	acc_states[0].init = 1
	don_states[arg.don_len-1].term = 1

	hmm.connect_all(acc_states)
	hmm.connect2(acc_states[-1], exon_state, 1)
	hmm.connect2(exon_state, exon_state, 1 - splices/exon_len)
	hmm.connect2(exon_state, don_states[0], splices/exon_len)
	hmm.connect_all(don_states)
	
	null_emits = hmm.train_emission(chr.seq, context=arg.null_ctx)	
	null_state = hmm.State(name='NULL', context=arg.null_ctx, emits=null_emits)
	hmm.connect2(null_state, null_state, 1)
	
	model = HMM(name=arg.hmm, states=acc_states + [exon_state] + don_states,
		null = null_state)
	model.write(arg.hmm)
	
	if arg.sources: output_sources(txa, arg.sources)
	if arg.replicant: output_replicant(model, txa, arg.replicant)
		

elif arg.model == 'splicing':
	gen = genome.Genome(gff3=arg.gff3, fasta=arg.fasta)
	ep_seqs = []
	en_seqs = []
	don_seqs = []
	acc_seqs = []
	intron_seqs = []
	txa = []
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
				iseq = tx.introns[i].seq()
				don_seqs.append(iseq[0:arg.don_len])
				intron_seqs.append(iseq[arg.don_len:-arg.acc_len])
				intron_len += len(iseq) - arg.don_len - arg.acc_len
				en_seqs.append(tx.exons[i + 1].seq())
				acc_seqs.append(iseq[-arg.acc_len:len(iseq)])
				exon_len += tx.exons[i].end - tx.exons[i].beg + 1
				splices += 1
				txid = 'model-' + str(splices)
				
				if arg.sources or arg.replicant:
					# create feature object for source sequence
					parent = genome.Feature(chr, tx.exons[i].beg, tx.exons[i + 1].end, tx.strand,
							'model', id=txid)
					exp = genome.Feature(chr, tx.exons[i].beg, tx.exons[i].end,
							tx.exons[i].strand, 'EXP', parent=txid)
					parent.add_child(exp)
					don = genome.Feature(chr, tx.introns[i].beg, tx.introns[i].beg + arg.don_len - 1,
							tx.introns[i].strand, 'DON', parent=txid)
					parent.add_child(don)
					int = genome.Feature(chr, tx.introns[i].beg + arg.don_len,
							tx.introns[i].end - arg.acc_len, tx.introns[i].strand,
							'INT', parent=txid)
					parent.add_child(int)
					acc = genome.Feature(chr, tx.introns[i].end - arg.acc_len + 1,
							tx.introns[i].end, tx.introns[i].strand, 'ACC', parent=txid)
					parent.add_child(acc)
					exn = genome.Feature(chr, tx.exons[i + 1].beg, tx.exons[i + 1].end,
							tx.exons[i + 1].strand, 'EXN', parent=txid)
					parent.add_child(exn)
				
					txa.append(parent)				
	
	ep_emits = hmm.train_emission(ep_seqs, context=arg.exon_ctx)
	don_emits = hmm.train_emissions(don_seqs, context=arg.don_ctx)
	intron_emits = hmm.train_emission(intron_seqs, context=arg.int_ctx)
	acc_emits = hmm.train_emissions(acc_seqs, context=arg.acc_ctx)
	en_emits = hmm.train_emission(en_seqs, context=arg.exon_ctx)
	
	ep_state = hmm.State(name='EXP', context=arg.exon_ctx, emits=ep_emits)
	ep_state.init = 1
	don_states = hmm.state_factory('DON', acc_emits)
	intron_state = hmm.State(name='INT', context=arg.int_ctx, emits=intron_emits)
	acc_states = hmm.state_factory('ACC', don_emits)
	en_state = hmm.State(name='EXN', context=arg.exon_ctx, emits=en_emits)
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
	
	null_emits = hmm.train_emission(chr.seq, context=arg.null_ctx)	
	null_state = hmm.State(name='NULL', context=arg.null_ctx, emits=null_emits)
	hmm.connect2(null_state, null_state, 1)
	
	model = HMM(name=arg.hmm, null=null_state,
		states=[ep_state] + don_states + [intron_state] + acc_states + [en_state])
	model.write(arg.hmm)
	if arg.sources: output_sources(txa, arg.sources)
	if arg.replicant: output_replicant(model, txa, arg.replicant)
	

elif arg.model == 'mRNA':
	gen = genome.Genome(gff3=arg.gff3, fasta=arg.fasta)
	u5_seqs = []
	koz_seqs = []
	atg_seqs = []
	cds_seqs = []
	u3_seqs = []
	txa = []
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
			
			if arg.sources or arg.replicant:
				txid = 'model-' + str(mRNAs)
				parent = genome.Feature(chr, tx.beg, tx.end, tx.strand, 'model', id=txid)
				parent.add_child(genome.Feature(chr, tx.beg,
					tx.beg + (beg - arg.koz_len) - 1, tx.strand, 'UTR5', parent=txid))
				parent.add_child(genome.Feature(chr, tx.beg + (beg - arg.koz_len),
					tx.beg + beg - 1, tx.strand, 'KOZ', parent=txid))
				parent.add_child(genome.Feature(chr, tx.beg + beg, tx.beg + beg + 2,
					tx.strand, 'ATG', parent=txid))
				parent.add_child(genome.Feature(chr, tx.beg + beg + 3, tx.beg + end, tx.strand, 'CDS',
					parent=txid))
				parent.add_child(genome.Feature(chr, tx.beg + end + 1, tx.end, tx.strand, 'UTR3',
					parent=txid))
				txa.append(parent)
					
	u5_emits = hmm.train_emission(u5_seqs, context=arg.u5_ctx)
	koz_emits = hmm.train_emissions(koz_seqs, context=arg.koz_ctx)
	atg_emits = hmm.train_emissions(atg_seqs, context=arg.atg_ctx)
	cds_emits = hmm.train_cds(cds_seqs, context=arg.cds_ctx)
	u3_emits = hmm.train_emission(u3_seqs, context=arg.u3_ctx)
	
	u5_state = hmm.State(name='UTR5', context=arg.u5_ctx, emits=u5_emits)
	koz_states = hmm.state_factory('KOZ', koz_emits)
	atg_states = hmm.state_factory('ATG', atg_emits)
	cds_states = hmm.state_factory('CDS', cds_emits)
	u3_state = hmm.State(name='UTR3', context=arg.u3_ctx, emits=u3_emits)
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
	
	null_emits = hmm.train_emission(chr.seq, context=arg.null_ctx)	
	null_state = hmm.State(name='NULL', context=arg.null_ctx, emits=null_emits)
	hmm.connect2(null_state, null_state, 1)
	
	model = HMM(name=arg.hmm, states=[u5_state] + koz_states + atg_states
		+ cds_states + [u3_state], null=null_state)
	model.write(arg.hmm)
	
	if arg.sources: output_sources(txa, arg.sources)
	if arg.replicant: output_replicant(model, txa, arg.replicant)

else:
	print('unrecognized model name')
	sys.exit(1)