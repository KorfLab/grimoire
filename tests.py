
import unittest

import toolbox
import sequence
import genome
import hmm
import decode

class TestToolbox(unittest.TestCase):

	def test_log(self):
		self.assertEqual(toolbox.log(1), 0)
		self.assertEqual(toolbox.log(0), -999)
		
	def test_sumlog(self):
		self.assertEqual(toolbox.sumlog(-1, -1), -0.3068528194400547)

	def test_GFF_file(self):
		gff = toolbox.GFF_file('data/A.thaliana.1percent.gff3')
		stuff = gff.get(chrom='Chr1', type='exon', beg=5000, end=9999)
		self.assertEqual(len(stuff), 20)

	def test_GFF_stream(self):
		gff = toolbox.GFF_stream('data/C.elegans.1percent.gff3')
		count = 0
		for e in gff:
			if e.type == 'exon': count += 1
		self.assertEqual(count, 2563)

	def test_FASTA_file(self):
		ff = toolbox.FASTA_file('data/C.elegans.1percent.fasta')
		c4 = ff.get('IV')
		self.assertEqual(len(c4.seq), 174938)

	def test_FASTA_stream(self):
		ff = toolbox.FASTA_stream('data/A.thaliana.1percent.fasta')
		text = ''
		for e in ff: text += e.id
		self.assertEqual(text, '12345mitochondriachloroplast')

class TestSequence(unittest.TestCase):
	
	def test_generate_kmers(self):
		k = ''.join(sequence.generate_kmers(alphabet='nt', k=2).keys())
		self.assertEqual(k, 'AAACAGATCACCCGCTGAGCGGGTTATCTGTT')
	
	def test_revcomp_str(self):
		s = 'AAAACCCGGT'
		self.assertEqual(sequence.revcomp_str(s), 'ACCGGGTTTT')

	def test_translate_str(self):
		s = 'ATAGCGAAT'
		self.assertEqual(sequence.translate_str(s), 'IAN')

	def test_DNA(self):
		s = sequence.DNA(seq='CATGCCAACAATGCNCAY')
		s.check_alphabet()
		p = s.translate()
		p.check_alphabet()
		self.assertEqual(p.seq, 'HANNAH')

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

class TestHMM(unittest.TestCase):
	
	def setUp(self):
		self.hmm = None
	
		gen = genome.Genome(
			gff3='data/C.elegans.1percent.gff3',
			fasta='data/C.elegans.1percent.fasta')
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

	def test_HMM_emit(self):
		self.assertEqual(self.hmm.null.emit['A'], 0.3504)

	def test_Decode_Viterbi(self):
		dna = sequence.DNA(name='test', seq='AAAAGTAAGTTTTT')
		v = decode.Viterbi(model=self.hmm, dna=dna)
		p = v.generate_path()
		ft = p.features(labels=['DON', 'GEN'], dna=dna)
		self.assertEqual(ft[1].type, 'DON')
		self.assertEqual(ft[1].beg, 5)



if __name__ == '__main__':
	unittest.main()


