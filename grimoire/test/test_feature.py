import unittest
import json

import grimoire.io as io
from grimoire.sequence import DNA
from grimoire.feature import Feature, mRNA, ncRNA, Gene, FeatureTable, FeatureError

class TestFeature(unittest.TestCase):

	def setUp(self):
		ff = io.FASTA_stream('data/at10.fa.gz')
		for e in ff:
			self.dna = DNA(seq=e.seq, name=e.id)
			break

	def test_Feature(self):
		f = Feature(self.dna, 50, 500, '+', 'exon', id='foo')
		f.validate()
		self.assertTrue(f.validated)
		self.assertEqual(f.issues, {})
		f.add_child(Feature(self.dna, 350, 353, '+', 'start', pid='foo'))
		self.assertFalse(f.validated)
		f.validate()
		self.assertEqual(f.issues, {})
		f.add_parent('bar')
		self.assertEqual(len(f.seq_str()), 451)
		g = f.gff()
		s = 'Chr1\t.\texon\t50\t500\t.\t+\t.\tID=foo;Parent=bar\nChr1\t.\tstart\t350\t353\t.\t+\t.\tParent=foo\n'
		self.assertEqual(g, s)
		self.assertTrue(f.overlap(Feature(self.dna, 490, 1000, '-', 'exon')))
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
		self.assertTrue(n.issues['short_exon'])
		self.assertTrue(n.issues['short_intron'])
		self.assertTrue(n.issues['donor'])
		self.assertTrue(n.issues['acceptor'])
		
		g = Gene(self.dna, 11649, 13714, '-', 'gene', id='g1')
		g.add_child(m)
		g.add_child(n)
		g.validate()
		self.assertTrue(g.issues['mixed_strands'])
		self.assertTrue(g.issues['child.beg<parent.beg'])
		
	def test_FeatureTable(self):
		ft = FeatureTable(dna=self.dna)
		ft.add_feature(Feature(self.dna, 3631, 5899, '+', 'gene', id='AT1G01010'))
		ft.add_feature(Feature(self.dna, 3631, 5899, '+', 'mRNA', id='AT1G01010.1', pid='AT1G01010'))
		ft.add_feature(Feature(self.dna, 3631, 3913, '+', 'exon', pid='AT1G01010.1'))
		ft.add_feature(Feature(self.dna, 3996, 4276, '+', 'exon', pid='AT1G01010.1'))
		ft.add_feature(Feature(self.dna, 4486, 4605, '+', 'exon', pid='AT1G01010.1'))
		ft.add_feature(Feature(self.dna, 4706, 5095, '+', 'exon', pid='AT1G01010.1'))
		ft.add_feature(Feature(self.dna, 5174, 5326, '+', 'exon', pid='AT1G01010.1'))
		ft.add_feature(Feature(self.dna, 5439, 5899, '+', 'exon', pid='AT1G01010.1'))
		g = ft.build_genes()
		self.assertEqual(len(g), 1)
		self.assertIsInstance(g[0], Gene)
		ft2 = FeatureTable(dna=self.dna)
		ft2.add_feature(Feature(self.dna, 3631, 4276, '+', 'gene', id='AT1G01010'))
		ft2.add_feature(Feature(self.dna, 3631, 4276, '+', 'mRNA', id='AT1G01010.1', pid='AT1G01010'))
		ft2.add_feature(Feature(self.dna, 3631, 3913, '+', 'exon', pid='AT1G01010.1'))
		ft2.add_feature(Feature(self.dna, 3996, 4276, '+', 'exon', pid='AT1G01010.1'))
		(same, diff, mat) = ft.nt_compare(ft2)
		self.assertEqual(same, 302653)
		ff = io.FASTA_stream('data/ce270.fa.gz')
		dna = None
		for e in ff:
			dna = DNA(seq=e.seq, name=e.id)
			break
		with self.assertRaises(FeatureError):
			ft.add_feature(Feature(dna, 900, 1200, '+', 'exon'))
		ft3 = FeatureTable(dna=self.dna)
		ft3.add_feature(Feature(self.dna, 50, 500, '+', 'exon'))
		ft3.add_feature(Feature(self.dna, 350, 520, '+', 'exon'))
		ft3.add_feature(Feature(self.dna, 1000, 2000, '+', 'exon'))
		f3 = ft3.fetch(900, 3000)
		self.assertEqual(len(f3), 1)
		

if __name__ == '__main__':
	unittest.main(verbosity=2)	