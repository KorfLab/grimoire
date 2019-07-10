
import unittest

import grimoire.toolbox as toolbox

class TestToolbox(unittest.TestCase):

	def test_log(self):
		self.assertEqual(toolbox.log(1), 0)
		self.assertEqual(toolbox.log(0), -999)
		
	def test_sumlog(self):
		self.assertEqual(toolbox.sumlog(-1, -1), -0.3068528194400547)

	def test_GFF_file(self):
		gff = toolbox.GFF_file('data/A.thaliana.1percent.gff3')
		stuff = gff.get(chrom='Chr1', type='exon', beg=5000, end=9999)
		self.assertEqual(len(stuff), 20)

	def test_GFF_stream(self):
		gff = toolbox.GFF_stream('data/C.elegans.1percent.gff3')
		count = 0
		for e in gff:
			if e.type == 'exon': count += 1
		self.assertEqual(count, 2563)

	def test_FASTA_file(self):
		ff = toolbox.FASTA_file('data/C.elegans.1percent.fasta')
		c4 = ff.get('IV')
		self.assertEqual(len(c4.seq), 174938)

	def test_FASTA_stream(self):
		ff = toolbox.FASTA_stream('data/A.thaliana.1percent.fasta')
		text = ''
		for e in ff: text += e.id
		self.assertEqual(text, '12345mitochondriachloroplast')


if __name__ == '__main__':
	unittest.main(verbosity=2)


