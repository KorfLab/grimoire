
import unittest
import copy
import os
import json

import grimoire.sequence as sequence
import grimoire.genome as genome
import grimoire.hmm as hmm
from grimoire.feature import Feature

class TestHMM(unittest.TestCase):

	def setUp(self):
		self.hmm = None
		
		donor_seqs = []
		gen_seqs = []
		cds_seqs = []

		gen = genome.Reader(
			gff='data/ce270.gff3.gz',
			fasta='data/ce270.fa.gz')
		for chrom in gen:
			genes = chrom.ftable.build_genes()
			gen_seqs.append(chrom.seq)
			for gene in genes:
				if gene.issues: continue
				for mrna in gene.transcripts():
					if not mrna.is_coding: continue
					cds_seqs.append(mrna.cds_str())
					for intron in mrna.introns:
						donor_seqs.append(intron.seq_str()[0:6])
		
		donor_emits = hmm.train_emissions(donor_seqs, context=0)
		donor_states = hmm.state_factory('DON', donor_emits)
		hmm.connect_all(donor_states)
		
		genome_emits = hmm.train_emission(gen_seqs, context=0)
		genome_state = hmm.State(name='GEN', context=0, emits=genome_emits)
		genome_state.init = 1
		genome_state.term = 1
		hmm.connect2(genome_state, genome_state, 0.99)
		hmm.connect2(genome_state, donor_states[0], 0.01)
		hmm.connect2(donor_states[-1], genome_state, 1)

		null_state = hmm.null_state_factory(file='data/ce270.fa.gz')

		self.hmm = hmm.HMM(name='test', null=null_state,
			states=donor_states + [genome_state])
		self.seqs = donor_seqs
		self.cdss = cds_seqs
		
	def test_emission_model(self):
		em = hmm.emission_model(context=2, alphabet='nt')
		self.assertEqual(len(em), 16)
		self.assertEqual(em['AA']['A'], 0)
		with self.assertRaises(hmm.HMMError):
			em = hmm.emission_model(context=-1, alphabet='nt')
		
	def test_train_emission(self):
		em = hmm.train_emission(self.seqs, context=1)
		self.assertEqual(len(em), 4)
		self.assertEqual(em['A']['A'], 0.2994)
		with self.assertRaises(hmm.HMMError):
			em = hmm.train_emission(self.seqs, context=-1)
		
	def test_train_emissions(self):
		em = hmm.train_emissions(self.seqs, context=1)
		self.assertEqual(len(em), len(self.seqs[0]))
		self.assertEqual(len(em[0]), 4)
		with self.assertRaises(hmm.HMMError):
			em = hmm.train_emission(self.seqs, context=-1)
		
	def test_train_cds(self):
		cds = hmm.train_cds(self.cdss, context=2)
		# there should be very few stop codons in frame
		self.assertTrue(cds[2]['TA']['A'] < cds[2]['TA']['T'])	
		
	def test_state_factory_and_connect_all(self):
		sl = hmm.state_factory('DON', hmm.train_emissions(self.seqs, context=0))
		self.assertEqual(len(sl), 6)
		self.assertEqual(sl[0].name, 'DON-0')
		hmm.connect_all(sl)
		self.assertEqual(sl[0].next['DON-1'], 1)
		
	def test_connect2(self):
		s1 = hmm.State(name='s1')
		s2 = hmm.State(name='s2')
		hmm.connect2(s1, s2, 0.5)
		self.assertEqual(s1.next['s2'], 0.5)
		
	def test_HMM_emit(self):
		self.assertEqual(self.hmm.null.emit['A'], 0.3197)
	
	def test_HMM_io(self):
		p1 = '/tmp/donor.hmm'
		p2 = '/tmp/donor.hmm.gz'
		self.hmm.write(p1)
		self.hmm.write(p2)
		self.assertLess(os.path.getsize(p2), os.path.getsize(p1))
		m1 = hmm.HMM.read(p1)
		self.assertIsInstance(m1, hmm.HMM)
		m2 = hmm.HMM.read(p2)
		self.assertIsInstance(m2, hmm.HMM)
		self.assertEqual(len(m1.states), len(m2.states), 7)
		self.assertFalse(m1.log)
		self.assertEqual(m1.log, m2.log)
		m1.convert2log()
		self.assertTrue(m1.log)
		with self.assertRaises(hmm.HMMError):
			m1.convert2log()

	def test_Viterbi_and_Parse(self):
		dna = sequence.DNA(name='test', seq='AAAAGTAAGTTTTT')
		v = hmm.Viterbi(model=self.hmm, dna=dna)
		p = v.generate_path()
		self.assertIsInstance(p, hmm.Parse)
		self.assertEqual(p.freq, None)
		self.assertEqual(p.path[0], 'GEN')
		self.assertAlmostEqual(p.score, 4.837415437963225)
		ft = p.features()
		self.assertEqual(ft[1].type, 'DON')
		self.assertEqual(ft[1].beg, 5)

	def test_Viterbi_logspace(self):
		dna = sequence.DNA(name='test', seq='AAAAGTAAGTTTTT')
		m = copy.deepcopy(self.hmm)
		m.convert2log()
		v = hmm.Viterbi(model=m, dna=dna)
		p = v.generate_path()
		self.assertAlmostEqual(p.score, 1.5763805776787159)

	def test_StochasticViterbi_and_Parse(self):
		dna = sequence.DNA(name='test', seq='TTTTTGTAAGTAAGTTTTT')
		v = hmm.StochasticViterbi(model=self.hmm, dna=dna, seed=1)
		ps = v.generate_paths(1000)
		self.assertIsInstance(ps[0], hmm.Parse)
		self.assertAlmostEqual(ps[0].score, 4.600333948735676)
		self.assertAlmostEqual(ps[0].freq, 0.459)

	def test_StochasticViterbi_logspace(self):
		dna = sequence.DNA(name='test', seq='TTTTTGTAAGTAAGTTTTT')
		m = copy.deepcopy(self.hmm)
		m.convert2log()
		v = hmm.StochasticViterbi(model=m, dna=dna, seed=1)
		ps = v.generate_paths(1000)
		self.assertAlmostEqual(ps[0].score, 1.5261288984112191)
		self.assertAlmostEqual(ps[0].freq, 0.459)

	def test_ForwardBackward(self):
		dna = sequence.DNA(name='test', seq='AGT')
		p = hmm.ForwardBackward(model=self.hmm, dna=dna)
		self.assertAlmostEqual(p.posterior(state_name='GEN', i=0),
			0.021213450734526202, places=8)
		self.assertAlmostEqual(p.posterior(state_name='DON-1', i=2),
			0.0031969999999999998, places=8)
		self.assertAlmostEqual(p.posterior(state_name='DON-0', i=1),
			0.0031969999999999998, places=8)
		self.assertAlmostEqual(p.posterior(state_name='GEN', i=2),
			0.0180164507345262, places=8)
	
	def test_Transcoder(self):
		dna = sequence.DNA(name='test', seq='AAAGTGAGCCCC')
		dna.ftable.add_feature(Feature(dna, 1, 3, '+', 'GEN'))
		dna.ftable.add_feature(Feature(dna, 4, 9, '+', 'DON'))
		dna.ftable.add_feature(Feature(dna, 10, 12, '+', 'GEN'))
		tc = hmm.Transcoder(model=self.hmm, dna=dna)
		self.assertAlmostEqual(tc.score, 1.250948871669515)
		
	
if __name__ == '__main__':
	unittest.main(verbosity=2)

