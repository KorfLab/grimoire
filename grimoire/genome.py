"""
Genome

This modeule contains classes and functions used for the representation and
validation of genomes.

"""

import re
import json
import sys
import operator
import gzip

import grimoire.io as io
import grimoire.sequence as sequence

class GenomeError(Exception):
	pass


class Genome:
	"""Class representing a genome, which has chromosomes (DNA objects)"""

	def __init__(self, species=None, fasta=None, gff3=None,
			check_alphabet=False):
		"""
		Parameters
		----------
		species: str
		   Species of genome, not really used yet
		fasta: file
			Path to fasta file, may be compressed
		gff3: file
			Path to gff3 file, may be compressed
		check_alphabet: bool
			Check whether it is the correct alphabet (good idea)
		"""

		self.species = species
		self.chromosomes = []
		ff = io.FASTA_stream(fasta)
		gf = io.GFF_file(gff3)

		mRNA_parts = ['CDS', 'exon']
		for entry in ff:
			chrom = sequence.DNA(name=entry.id, seq=entry.seq)
			if check_alphabet: chrom.check_alphabet()

			# convert protein-coding gene-based GFF to Features
			genes = {}
			mRNAs = {}
			parts = []
			for g in gf.get(chrom=chrom.name):
				id, pid = None, None
				im = re.search('ID=([\w\.\:]+)', g.attr)
				pm = re.search('Parent=([\w\.:]+)', g.attr)
				if im: id = im[1]
				if pm: pid = pm[1]
				if g.type == 'gene':
					genes[id] = ProteinCodingGene(chrom, g.beg, g.end, g.strand, g.type, id=id, parent_id=pid)
				elif g.type == 'mRNA':
					mRNAs[id] = mRNA(chrom, g.beg, g.end, g.strand, g.type, id=id, parent_id=pid)
				elif g.type in mRNA_parts:
					parts.append(Feature(chrom, g.beg, g.end, g.strand, g.type, id=id, parent_id=pid))

			# add parts to mRNAs
			for f in parts:
				f.validate()
				if f.parent_id in mRNAs:
					mRNAs[f.parent_id].add_child(f)

			# add mRNAs to genes
			for txid in mRNAs:
				f = mRNAs[txid]
				if f.parent_id in genes:
					genes[f.parent_id].add_child(f)

			# add genes to chromosome
			for gid in genes:
				chrom.features.append(genes[gid])

			# add chromosome to genome
			self.chromosomes.append(chrom)

class Genomic:
	"""Class for iterating through DNA objects with attached feature tables"""

	def __init__(self, fasta=None, gff=None, check_alphabet=None):
		"""
		Use path to file or file pointer.

		Parameters
		----------
		fasta: str
			Path to fasta file (may be compressed)
		gff: str
			Path to gff file (may be compressed)
		"""

		self.fp = None
		self.gz = False
		self.gff = io.GFF_file(gff)
		self.check = check_alphabet

		if re.search('\.gz$', fasta):
			self.fp = gzip.open(fasta)
			self.gz = True
		else:
			self.fp = open(fasta, 'r')
		self.lastline = ''
		self.done = False

	def __iter__(self):
		return self

	def __next__(self):
		return self.next()

	def next(self):
		"""
		Retrieves the next entry of the FASTA file as DNA with features
		"""

		if self.done: raise StopIteration()
		header = None
		if self.lastline[0:1] == '>':
			header = self.lastline
		else:
			header = self.fp.readline()
			if self.gz: header = str(header, 'utf-8')

		m = re.search('>\s*(\S+)\s*(.*)', header)
		id = m[1]
		desc = m[2]
		seq = []

		while (True):
			line = self.fp.readline()
			if self.gz: line = str(line, 'utf-8')
			if line[0:1] == '>':
				self.lastline = line
				break
			if line == '':
				self.done = True
				self.fp.close()
				break

			line = line.replace(' ', '')
			seq.append(line.strip())
		
		dna = sequence.DNA(name = id, seq=''.join(seq))
		if self.check: dna.check_alphabet()

		# add features
		for g in self.gff.get(chrom=dna.name):
			id, pid = None, None
			im = re.search('ID=([\w\.\:]+)', g.attr)
			pm = re.search('Parent=([\w\.:]+)', g.attr)
			if im: id = im[1]
			if pm: pid = pm[1]
			dna.features.append(Feature(dna, g.beg, g.end, g.strand,
				g.type, source=g.source, score=g.score,
				id=id, parent_id=pid))
				
		# add genes (intentionally redudant with above: separate addresses)
		mRNA_parts = ['CDS', 'exon']
		genes = {}
		mRNAs = {}
		parts = []
		for g in self.gff.get(chrom=dna.name):
			id, pid = None, None
			im = re.search('ID=([\w\.\:]+)', g.attr)
			pm = re.search('Parent=([\w\.:]+)', g.attr)
			if im: id = im[1]
			if pm: pid = pm[1]

			if g.type == 'gene':
				genes[id] = ProteinCodingGene(dna, g.beg, g.end, 
					g.strand, g.type, id=id, parent_id=pid)
			elif g.type == 'mRNA':
				mRNAs[id] = mRNA(dna, g.beg, g.end, g.strand, g.type,
					id=id, parent_id=pid)
			elif g.type in mRNA_parts:
				parts.append(Feature(dna, g.beg, g.end, g.strand, g.type,
					id=id, parent_id=pid))
		
			# add parts to mRNAs
			for f in parts:
				f.validate()
				if f.parent_id in mRNAs:
					mRNAs[f.parent_id].add_child(f)

			# add mRNAs to genes
			for txid in mRNAs:
				f = mRNAs[txid]
				if f.parent_id in genes:
					genes[f.parent_id].add_child(f)

			# add genes to dna
			for gid in genes:
				dna.genes.append(genes[gid])
			
		return dna