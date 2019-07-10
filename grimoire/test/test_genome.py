
import unittest

import grimoire.genome as genome

class TestGenome(unittest.TestCase):
	
	def test_Genome(self):
		issues = 0
		map = {'1':'Chr1', '2':'Chr2', '3':'Chr3', '4':'Chr4', '5':'Chr5',
			'mitochondria':'ChrM', 'chloroplast':'ChrC'}
		gen = genome.Genome(
			fasta='data/A.thaliana.1percent.fasta',
			gff3='data/A.thaliana.1percent.gff3',
			chr_map=map)
		for c in gen.chromosomes:
			for g in c.features:
				for m in g.mRNAs():
					if m.issues: issues += 1
		self.assertEqual(issues, 30)

if __name__ == '__main__':
	unittest.main(verbosity=2)


