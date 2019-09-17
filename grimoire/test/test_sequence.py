
import unittest

import grimoire.sequence as sequence

class TestSequence(unittest.TestCase):
	
	def test_DNA(self):
		s = sequence.DNA(seq='CATGCCAACAATGCNCAY')
		self.assertIsInstance(s, sequence.DNA)
		try:
			s.check_alphabet()
		except sequence.SequenceError:
			self.fail('DNA object failed check_alphabet()')
		p = s.translate()
		self.assertIsInstance(p, sequence.Protein)
		self.assertEqual(p.seq, 'HANNAH')
		s.revcomp()
		self.assertIsInstance(s, sequence.DNA)
		self.assertEqual(s.seq, 'RTGNGCATTGTTGGCATG')
		
		s2 = sequence.DNA(seq='CATGCCAACAATGCNCAQ')
		with self.assertRaises(sequence.SequenceError):
			s2.check_alphabet()
		with self.assertRaises(sequence.SequenceError):
			p2 = s2.translate()
		with self.assertRaises(sequence.SequenceError):
			s2.revcomp()

		
	def test_Protein(self):
		s = sequence.Protein(seq='MDDDIAALVVDNGSGMCKAGFAGDDAP',
			name='NP_990849.1')
		self.assertIsInstance(s, sequence.Protein)
		try:
			s.check_alphabet()
		except sequence.SequenceError:
			self.fail('Protein object failed check_alphabet()')

if __name__ == '__main__':
	unittest.main(verbosity=2)

