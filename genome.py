
import re
import json
import sys
import operator

import toolbox

class AnnotationError(Exception):
	pass

class Feature:

	def __init__(self, chrom, beg, end, strand, type, id=None):
		self.chrom = chrom
		self.beg = beg
		self.end = end
		self.strand = strand
		self.type = type
		self.id = id

	def seq(self):
		seq = self.chrom.seq[self.beg-1:self.end]
		if self.strand == '-': seq = toolbox.dna.revcomp(seq)
		return seq

class mRNA:

	def __init__(self, id=None, features=None):
			
		## sanity checks for chrom, beg, end, strand
		chrom = features[0].chrom
		strand = features[0].strand
		for f in features:
			if f.chrom != chrom: raise AnnotationError('mixed chroms')
			if f.strand != strand: raise AnnotationError('mixed strands')
			if f.beg > f.end: raise AnnotationError('beg > end')
		self.chrom = chrom
		self.strand = strand
		
		
		## build mRNA from features
		self.exons = []
		self.cdss = []
		self.utr5s = []
		self.utr3s = []
		self.id = id
		for f in features:
			if f.type == 'exon': self.exons.append(
				Feature(f.chrom, f.beg, f.end, f.strand, f.type))
			elif f.type == 'CDS': self.cdss.append(
				Feature(f.chrom, f.beg, f.end, f.strand, f.type))
			elif f.type == 'five_prime_UTR': self.utr5s.append(
				Feature(f.chrom, f.beg, f.end, f.strand, f.type))
			elif f.type == 'three_prime_UTR': self.utr3s.append(
				Feature(f.chrom, f.beg, f.end, f.strand, f.type))
			else:
				raise AnnotationError('unknown type in feature table')
				
		# sort lists
		self.exons.sort(key = operator.attrgetter('beg'))
		self.cdss.sort(key = operator.attrgetter('beg'))
		self.utr5s.sort(key = operator.attrgetter('beg'))
		self.utr3s.sort(key = operator.attrgetter('beg'))

		self.beg = self.exons[0].beg
		self.end = self.exons[-1].end
		
		# create introns
		self.introns = []
		for i in range(len(self.exons)-1):
			beg = self.exons[i].end +1
			end = self.exons[i+1].beg -1
			self.introns.append(
				Feature(self.chrom, beg, end, self.strand, 'intron'))

	def cds(self):
		seq = []
		for exon in self.cdss: seq.append(exon.seq())
		if self.strand == '-': seq.reverse()
		return ''.join(seq)

	def protein(self):
		return toolbox.dna.translate(self.cds())


class Gene:

	def __init__(self, id=None, transcripts=None):
		self.id = id
		self.beg = None
		self.end = None
		self.strand = None
		self.transcripts = transcripts

class Chromosome:
	
	def __init__(self, fasta=None, gff=None):
		self.id = fasta.id
		self.seq = fasta.seq
		self.genes = []
		ft = {}
		parts = ['CDS', 'exon', 'five_prime_UTR', 'three_prime_UTR']
		for part in parts:
			for entry in gff.get(type=part, chrom=self.id):
				pid = re.search('Parent=([\w\.]+)', entry.attr)[1]
				if pid not in ft: ft[pid] = []
				ft[pid].append(Feature(self, entry.beg, entry.end,
					entry.strand, entry.type))
		gt = {}
		for entry in gff.get(type='mRNA', chrom=self.id):
			tid = re.search('ID=([\w\.]+)', entry.attr)[1]
			gid = re.search('Parent=([\w\.]+)', entry.attr)[1]
			if gid not in gt: gt[gid] = {}
			gt[gid][tid] = entry
				
		for gid in gt:
			txs = []
			for tid in gt[gid]:
				tx = mRNA(id=tid, features=ft[tid])
				txs.append(tx)
			gene = Gene(id=gid, transcripts=txs)
			self.genes.append(gene)


class Genome:

	def __init__(self, species=None, fasta=None, gff=None):
		self.species = species
		self.chromosomes = []
		gf = toolbox.gff.Gff(gff)
		ff = toolbox.fasta.FastaFile(fasta)
		for id in ff.ids:
			entry = ff.get(id)
			self.chromosomes.append(Chromosome(fasta=entry, gff=gf))
		
