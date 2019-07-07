
import unittest

import toolbox
import sequence
import genome
import hmm

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

class TestSequence(unittest.TestCase):
	
	def test_generate_kmers(self):
		k = ''.join(sequence.generate_kmers(alphabet='nt', k=2).keys())
		self.assertEqual(k, 'AAACAGATCACCCGCTGAGCGGGTTATCTGTT')
	
	def test_revcomp_str(self):
		s = 'AAAACCCGGT'
		self.assertEqual(sequence.revcomp_str(s), 'ACCGGGTTTT')

	def test_translate_str(self):
		s = 'ATAGCGAAT'
		self.assertEqual(sequence.translate_str(s), 'IAN')

	def test_DNA(self):
		s = sequence.DNA(seq='CATGCCAACAATGCNCAY')
		s.check_alphabet()
		p = s.translate()
		p.check_alphabet()
		self.assertEqual(p.seq, 'HANNAH')

class TestGenome(unittest.TestCase):
	
	def test_Genome(self):
		issues = 0
		gen = genome.Genome(fasta='data/C.elegans.1percent.fasta',
			gff3='data/C.elegans.1percent.gff3')
		for c in gen.chromosomes:
			for g in c.features:
				for m in g.mRNAs():
					if m.issues: issues += 1
		self.assertEqual(issues, 94)

class TestHMM(unittest.TestCase):
	pass


class TestDecode(unittest.TestCase):
	pass
	

	


if __name__ == '__main__':
	unittest.main()