
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
		self.issues = []
		if beg < 0: self.issues.append('beg < 0')
		if beg > end: self.issues.append('beg > end')
		if end > chrom.length: self.issues.append('end out of range')

	def seq(self):
		seq = self.chrom.seq[self.beg-1:self.end]
		if self.strand == '-': seq = toolbox.dna.revcomp(seq)
		return seq

class mRNA:

	def __init__(self, id=None, features=None, rules='std'):
			
		## sanity checks for chrom, beg, end, strand
		self.issues = []
		chrom = features[0].chrom
		strand = features[0].strand
		for f in features:
			if f.issues: self.issues.append('feature issue')
			if f.chrom != chrom: self.issues.append('mixed chroms')
			if f.strand != strand: self.issues.append('mixed strands')
			if f.beg > f.end: self.issues.append('beg > end')
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
		
		# check for overlapping features
		if self.overlaps(self.exons): self.issues.append('exon overlaps')
		if self.overlaps(self.cdss): self.issues.append('cds overlaps')
		if self.overlaps(self.utr5s): self.issues.append('utr5 overlaps')
		if self.overlaps(self.utr3s): self.issues.append('utr3 overlaps')

		# assignments
		self.beg = self.exons[0].beg
		self.end = self.exons[-1].end
		
		# create introns
		self.introns = []
		for i in range(len(self.exons)-1):
			beg = self.exons[i].end +1
			end = self.exons[i+1].beg -1
			self.introns.append(
				Feature(self.chrom, beg, end, self.strand, 'intron'))

		# gene structure warnings
		min_intron = 30
		max_intron = 10000
		dons = {'GT':True}
		accs = {'AG':True}
		starts = {'ATG':True}
		stops = {'TAA':True, 'TGA':True, 'TAG':True}
		if rules == 'std':
			pass # don't override parameters above
		elif rules == 'mammal':
			min_intron = 50
			max_intron = 50000
		else:
			pass # there ought to be other rule sets...
			
	
		# canonical splicing and intron length
		for intron in self.introns:
			s = intron.seq()
			don = s[0:2]
			acc = s[-2:len(s)]
			if don not in dons: self.issues.append('don:' + don)
			if acc not in accs: self.issues.append('acc:' + acc)
			if len(s) < min_intron or len(s) > max_intron:
				self.issues.append('int_len:' + str(len(s)))
			
		# translation checks
		cds = self.cds()
		pro = toolbox.dna.translate(cds)
		start = cds[0:3]
		stop = cds[-3:len(cds)]
		if start not in starts: self.issues.append('start:' + start)
		if stop not in stops: self.issues.append('stop:' + stop)
		for i in range(len(pro) - 1):
			if pro[i:i+1] == '*': self.issues.append('ptc:' + str(i))
	
	def cds(self):
		seq = []
		for exon in self.cdss: seq.append(exon.seq())
		if self.strand == '-': seq.reverse()
		return ''.join(seq)

	def protein(self):
		return toolbox.dna.translate(self.cds())

	def overlaps(self, list):
		if len(list) > 1:
			for i in range(len(list)- 1):
				if list[i].end >= list[i+1].beg: return True
		return False
			

class Gene:

	def __init__(self, id=None, transcripts=None):
		self.id = id
		self.beg = None
		self.end = None
		self.strand = None
		self.transcripts = transcripts
		self.issues = []
		
		self.beg = transcripts[0].beg
		self.end = transcripts[0].end
		self.strand = transcripts[0].strand
		for tx in transcripts:
			if tx.beg < self.beg: self.beg = tx.beg
			if tx.end > self.end: self.end = tx.end
			if tx.strand != self.strand: self.issues.append('mixed strands')
			if tx.issues: self.issues.append('tx issue')

class Chromosome:
	
	def __init__(self, fasta=None, gff=None):
		self.id = fasta.id
		self.seq = fasta.seq
		self.genes = []
		self.length = len(fasta.seq)
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
	
	