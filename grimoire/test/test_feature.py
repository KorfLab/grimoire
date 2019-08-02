import unittest

import grimoire.io as io
import grimoire.sequence as sequence
import grimoire.feature as feature

class TestFeature(unittest.TestCase):

	def test_Feature(self):
		ff = io.FASTA_stream('data/A.thaliana.1percent.fasta.gz')
		entry = ff.next()
		dna = sequence.DNA(seq=entry.seq, name=entry.id)
		f = feature.Feature(dna, 50, 500, '+', 'exon', id='foo')
		f.add_child(feature.Feature(dna, 350, 353, '+', 'start'))
		f.validate()
		f.seq_str()
		g = f.gff()
		s = 'Chr1\t.\texon\t50\t500\t.\t+\t.\tID=foo\nChr1\t.\tstart\t350\t353\t.\t+\t.\t\n'
		self.assertEqual(g, s)
		self.assertEqual(f.overlap(feature.Feature(dna, 490, 1000, '-', 'exon')), True)
		self.assertEqual(str(f), g)
	
	def test_mRNA(self):
		pass


if __name__ == '__main__':
	unittest.main(verbosity=2)