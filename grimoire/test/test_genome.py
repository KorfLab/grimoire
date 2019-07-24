
import unittest

import grimoire.genome as genome

class TestGenome(unittest.TestCase):
	
	def test_Genome(self):
		issues = 0
		gen = genome.Genome(
			fasta='data/A.thaliana.1percent.fasta.gz',
			gff3='data/A.thaliana.1percent.gff3.gz')
		for c in gen.chromosomes:
			for g in c.features:
				for m in g.mRNAs():
					if m.issues: issues += 1
		self.assertEqual(issues, 31)

if __name__ == '__main__':
	unittest.main(verbosity=2)


