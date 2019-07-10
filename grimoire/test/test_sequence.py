
import unittest

import grimoire.sequence as sequence

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

if __name__ == '__main__':
	unittest.main(verbosity=2)


