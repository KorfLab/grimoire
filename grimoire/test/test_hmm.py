
import unittest
import copy
import os

import grimoire.sequence as sequence
import grimoire.genome as genome
import grimoire.hmm as hmm

class TestHMM(unittest.TestCase):

	def setUp(self):
		self.hmm = None

		gen = genome.Genome(
			gff3='data/C.elegans.1percent.gff3.gz',
			fasta='data/C.elegans.1percent.fasta.gz')
		donor_seqs = []
		null_seqs = []
		for chr in gen.chromosomes:
			null_seqs.append(chr.seq)
			for gene in chr.features:
				if not gene.mRNAs(): continue # some tx are miRNA etc
				if gene.issues: continue
				for mrna in gene.mRNAs():
					for intron in mrna.introns:
						donor_seqs.append(intron.seq_str()[0:6])
		donor_emits = hmm.train_emissions(donor_seqs, context=0)
		donor_states = hmm.state_factory('DON', donor_emits)
		hmm.connect_all(donor_states)

		genome_emits = hmm.train_emission(null_seqs, context=1)
		genome_state = hmm.State(name='GEN', context=1, emits=genome_emits)
		genome_state.init = 1
		genome_state.term = 1
		hmm.connect2(genome_state, genome_state, 0.99)
		hmm.connect2(genome_state, donor_states[0], 0.01)
		hmm.connect2(donor_states[-1], genome_state, 1)

		null_emits = hmm.train_emission(chr.seq, context=0)
		null_state = hmm.State(name='NULL', context=0, emits=null_emits)
		null_state.init = 1
		null_state.term = 1
		hmm.connect2(null_state, null_state, 1)

		self.hmm = hmm.HMM(name='test', null=null_state,
			states=donor_states + [genome_state])

"""
	def test_HMM_emit(self):
		self.assertEqual(self.hmm.null.emit['A'], 0.3504)

	def test_HMM_write(self):
		self.hmm.write('data/donor.hmm')
		self.assertEqual(os.path.getsize('data/donor.hmm'), 3257)
		self.hmm.write('data/donor.hmm')
		self.assertEqual(os.path.getsize('data/donor.hmm.gz'), 489)

	def test_HMM_read(self):
		model = hmm.HMM.read('data/donor.hmm')
		self.assertEqual(len(model.states), 7)
		self.assertEqual(model.name, 'test')
		self.assertEqual(model.log, False)
		model = hmm.HMM.read('data/donor.hmm.gz')
		self.assertEqual(len(model.states), 7)
		self.assertEqual(model.name, 'test')
		self.assertEqual(model.log, False)

	def test_Decode_Viterbi(self):
		dna = sequence.DNA(name='test', seq='AAAAGTAAGTTTTT')
		v = decode.Viterbi(model=self.hmm, dna=dna)
		p = v.generate_path()
		self.assertEqual(p.score, 4.401452575029605)
		ft = p.features()
		self.assertEqual(ft[1].type, 'DON')
		self.assertEqual(ft[1].beg, 5)

	def test_Decode_Viterbi_logspace(self):
		dna = sequence.DNA(name='test', seq='AAAAGTAAGTTTTT')
		m = copy.deepcopy(self.hmm)
		m.convert2log()
		v = decode.Viterbi(model=m, dna=dna)
		p = v.generate_path()
		self.assertEqual(p.score, 1.481934617131616)

	def test_Decode_StochasticViterbi(self):
		dna = sequence.DNA(name='test', seq='TTTTTGTAAGTAAGTTTTT')
		v = decode.StochasticViterbi(model=self.hmm, dna=dna, seed=1)
		ps = v.generate_paths(1000)
		self.assertEqual(ps[0].score, 0.8980462178059174)
		self.assertEqual(ps[0].freq, 0.524)

	def test_Decode_StochasticViterbi_logspace(self):
		dna = sequence.DNA(name='test', seq='TTTTTGTAAGTAAGTTTTT')
		m = copy.deepcopy(self.hmm)
		m.convert2log()
		v = decode.StochasticViterbi(model=m, dna=dna, seed=1)
		ps = v.generate_paths(1000)
		self.assertEqual(ps[0].score, -0.10753374451446263)
		self.assertEqual(ps[0].freq, 0.524)

	def test_Decode_ForwardBackward(self):
		dna = sequence.DNA(name='test', seq='AGT')
		p = decode.ForwardBackward(model=self.hmm, dna=dna)
		self.assertAlmostEqual(p.posterior(state_name='GEN', i=0),
			0.01224654564, places=8)
		self.assertAlmostEqual(p.posterior(state_name='DON-1', i=2),
			0.0025, places=8)
		self.assertAlmostEqual(p.posterior(state_name='DON-0', i=1),
			0.0025, places=8)
		self.assertAlmostEqual(p.posterior(state_name='GEN', i=2),
			0.009746545644, places=8)
"""

if __name__ == '__main__':
	unittest.main(verbosity=2)
