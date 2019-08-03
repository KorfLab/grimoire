
import unittest

import grimoire.sequence as sequence

class TestSequence(unittest.TestCase):
	
	def test_DNA(self):
		s = sequence.DNA(seq='CATGCCAACAATGCNCAY')
		s.check_alphabet()
		s.revcomp()
		p = s.translate()
		p.check_alphabet()
		self.assertEqual(p.seq, 'HANNAH')
		
	def test_Protein(self):
		s = sequence.Protein(seq='MDDDIAALVVDNGSGMCKAGFAGDDAP', name='NP_990849.1')
		s.fasta()
		s.check_alphabet()

if __name__ == '__main__':
	unittest.main(verbosity=2)

