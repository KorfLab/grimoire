import unittest

import grimoire.io as io
from grimoire.sequence import DNA
from grimoire.feature import Feature, mRNA, ncRNA, Gene

class TestFeature(unittest.TestCase):

	def setUp(self):
		ff = io.FASTA_stream('data/A.thaliana.1percent.fasta.gz')
		entry = ff.next()
		self.dna = DNA(seq=entry.seq, name=entry.id)

	def test_Feature(self):
		f = Feature(self.dna, 50, 500, '+', 'exon', id='foo')
		f.add_child(Feature(self.dna, 350, 353, '+', 'start'))
		f.validate()
		f.seq_str()
		g = f.gff()
		s = 'Chr1\t.\texon\t50\t500\t.\t+\t.\tID=foo\nChr1\t.\tstart\t350\t353\t.\t+\t.\t\n'
		self.assertEqual(g, s)
		self.assertEqual(f.overlap(Feature(self.dna, 490, 1000, '-', 'exon')), True)
		self.assertEqual(str(f), g)
	
	def test_Gene_and_Transcripts(self):
		
		m = mRNA(self.dna, 11649, 13714, '-', 'mRNA', id='m1')
		m.add_child(Feature(self.dna, 13335, 13714, '-', 'exon'))
		m.add_child(Feature(self.dna, 11864, 12940, '-', 'CDS'))
		m.add_child(Feature(self.dna, 11649, 13173, '-', 'exon'))
		m.validate()
		self.assertEqual(len(m.issues), 0)
		self.assertEqual(len(m.tx_str()), 1905)
		self.assertEqual(len(m.cds_str()), 1077)
		self.assertEqual(len(m.protein_str()), 359)
		
		n = ncRNA(self.dna, 100, 200, '+', 'RNA', id='n1')
		n.add_child(Feature(self.dna, 100, 180, '+', 'exon'))
		n.add_child(Feature(self.dna, 190, 200, '+', 'exon'))
		n.validate()
		self.assertEqual(n.issues['short_exon'], True)
		self.assertEqual(n.issues['short_intron'], True)
		self.assertEqual(n.issues['donor'], True)
		self.assertEqual(n.issues['acceptor'], True)
		
		g = Gene(self.dna, 11649, 13714, '-', 'gene', id='g1')
		g.add_child(m)
		g.add_child(n)
		g.validate()
		self.assertEqual(g.issues['mixed_strands'], True)
		self.assertEqual(g.issues['child.beg<parent.beg'], True)

if __name__ == '__main__':
	unittest.main(verbosity=2)
	