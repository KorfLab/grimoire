import unittest
import os

import grimoire.io as io

class TestIO(unittest.TestCase):

	def test_FASTA_entry(self):
		entry = io.FASTA_entry('foo', 'bar', 'ACGT')
		self.assertIsInstance(entry, io.FASTA_entry)
		s = entry.string(wrap=2)
		self.assertEqual(s, '>foo bar\nAC\nGT\n')
	
	def test_FASTA_file(self):
		if os.system('gunzip -c data/C.elegans.1percent.fasta > data/test.fasta'):
			raise 'unable to create test file'
		ff = io.FASTA_file('data/test.fasta')
		self.assertIsInstance(ff, io.FASTA_file)
		c4 = ff.get('IV')
		self.assertEqual(len(c4.seq), 174938)
		os.system('rm data/test.fasta')
		with self.assertRaises(io.FASTA_error):
			ff = io.FASTA_file('data/README.md')

	def test_FASTA_stream(self):
		ff = io.FASTA_stream('data/A.thaliana.1percent.fasta.gz')
		text = ''
		for e in ff: text += e.id
		self.assertEqual(text, 'Chr1Chr2Chr3Chr4Chr5ChrMChrC')
		ff = io.FASTA_stream('data/A.thaliana.1percent.gff3.gz')
		with self.assertRaises(io.FASTA_error):
			for e in ff: print(e.id)
		
	def test_GFF_entry(self):
		e = io.GFF_entry('I	WormBase	gene	3747	3909	.	-	.	ID=Gene:WBGene00023193;Name=WBGene00023193;interpolated_map_position=-21.9001;sequence_name=Y74C9A.6;biotype=snoRNA;so_term_name=snoRNA_gene;curie=WB:WBGene00023193;Alias=Y74C9A.6')
		self.assertIsInstance(e, io.GFF_entry)
		self.assertEqual(e.source, 'WormBase')
		with self.assertRaises(io.GFF_skip):
			c = io.GFF_entry('#comment line')
		with self.assertRaises(io.GFF_skip):
			b = io.GFF_entry('\n')
		with self.assertRaises(io.GFF_error):
			f = io.GFF_entry('I	WormBase	gene	3747	3909')
		
	def test_GFF_file(self):
		gff = io.GFF_file('data/A.thaliana.1percent.gff3.gz')
		self.assertIsInstance(gff, io.GFF_file)
		stuff = gff.get(chrom='Chr1', type='exon', beg=5000, end=9999)
		self.assertEqual(len(stuff), 20)

	def test_GFF_stream(self):
		gff = io.GFF_stream('data/C.elegans.1percent.gff3.gz')
		self.assertIsInstance(gff, io.GFF_stream)
		count = 0
		for e in gff:
			if e.type == 'exon': count += 1
		self.assertEqual(count, 2563)

if __name__ == '__main__':
	unittest.main(verbosity=2)
