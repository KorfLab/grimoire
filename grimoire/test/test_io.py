import unittest
import os

import grimoire.io as io

class TestIO(unittest.TestCase):
	
	def test_FASTA_file(self):
		if os.system('gunzip -c data/C.elegans.1percent.fasta > data/test.fasta'):
			raise 'unable to create test file'
		ff = io.FASTA_file('data/test.fasta')
		c4 = ff.get('IV')
		self.assertEqual(len(c4.seq), 174938)
		os.system('rm data/test.fasta')

	def test_FASTA_stream(self):
		ff = io.FASTA_stream('data/A.thaliana.1percent.fasta.gz')
		text = ''
		for e in ff: text += e.id
		self.assertEqual(text, 'Chr1Chr2Chr3Chr4Chr5ChrMChrC')
		
	def test_GFF_file(self):
		gff = io.GFF_file('data/A.thaliana.1percent.gff3.gz')
		stuff = gff.get(chrom='Chr1', type='exon', beg=5000, end=9999)
		self.assertEqual(len(stuff), 20)

	def test_GFF_stream(self):
		gff = io.GFF_stream('data/C.elegans.1percent.gff3.gz')
		count = 0
		for e in gff:
			if e.type == 'exon': count += 1
		self.assertEqual(count, 2563)

if __name__ == '__main__':
	unittest.main(verbosity=2)