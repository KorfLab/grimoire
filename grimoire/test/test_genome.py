
import unittest

import grimoire.genome as genome

class TestGenome(unittest.TestCase):
	
	def test_Genome(self):
		gerr = 0
		terr = 0
		ferr = 0
		gen = genome.Reader(
			fasta='data/C.elegans.1percent.fasta.gz',
			gff='data/C.elegans.1percent.gff3.gz')
		for chrom in gen:
			genes = genome.build_genes(chrom)
			for g in genes:
				if g.issues: gerr += 1
				for t in g.transcripts():
					if t.issues: terr += 1
					for f in t.children:
						if f.issues: ferr += 1
		
		self.assertEqual(gerr, 29)
		self.assertEqual(terr, 47)
		self.assertEqual(ferr, 0)

if __name__ == '__main__':
	unittest.main(verbosity=2)

