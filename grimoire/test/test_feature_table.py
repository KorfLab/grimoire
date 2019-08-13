
import unittest

from grimoire.feature_table import FeatureTable
from grimoire.hmm import HMM
from grimoire.sequence import DNA
from grimoire.feature import Feature

class TestFeatureTable(unittest.TestCase):
	
	def test_build_genes(self):
		pass # this is tested in test_genome
	
	def test_compare(self):
		dna = DNA(seq='AAAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTTTTTTTT')
		ft1 = FeatureTable(dna=dna, features=[
			Feature(dna,  1, 10, '+', 'PolyA'),
			Feature(dna, 11, 20, '+', 'PolyC'),
			Feature(dna, 21, 30, '+', 'PolyG'),
			Feature(dna, 31, 40, '+', 'PolyT')
		])
		ft2 = FeatureTable(dna=dna, features=[
			Feature(dna,  1, 10, '+', 'PolyA'),
			Feature(dna, 11, 20, '-', 'PolyG'),
			Feature(dna, 21, 40, '+', 'PolyK')
		])
		stats = ft1.compare(ft2)
		self.assertEqual(stats['nt_same'], 10)
		self.assertEqual(stats['matrix']['PolyA']['PolyA'], 10)

if __name__ == '__main__':
	unittest.main(verbosity=2)

