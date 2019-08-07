
import unittest

import grimoire.toolbox as toolbox

class TestToolbox(unittest.TestCase):

	def test_prod(self):
		l = [2, 3.5]
		self.assertAlmostEqual(toolbox.prod(l), 7)
		
	def test_log(self):
		self.assertAlmostEqual(toolbox.log(1), 0)
		self.assertAlmostEqual(toolbox.log(0), -999)
		
	def test_sumlog(self):
		self.assertAlmostEqual(toolbox.sumlog(-1, -1), -0.3068528194400547)

	def test_generate_kmers(self):
		k = ''.join(toolbox.generate_kmers(alphabet='nt', k=2).keys())
		self.assertEqual(k, 'AAACAGATCACCCGCTGAGCGGGTTATCTGTT')
	
	def test_revcomp_str(self):
		s = 'AAAACCCGGT'
		self.assertEqual(toolbox.revcomp_str(s), 'ACCGGGTTTT')

	def test_translate_str(self):
		s = 'ATAGCGAAT'
		self.assertEqual(toolbox.translate_str(s), 'IAN')
				
if __name__ == '__main__':
	unittest.main(verbosity=2)
