
import unittest
import os
import subprocess
import re

class Test_Programs(unittest.TestCase):

	def test_tutorial(self):
		
		assert(os.system('calfo --fasta data/ce270.fa.gz --gff data/ce270.gff3.gz --title ce270gff --html /tmp/ce270gff3.html') == 0)
		assert(os.path.getsize('/tmp/ce270gff3.html') == 260615)

		assert(os.system('haman --fasta data/ce270.fa.gz --gff data/ce270.gff3.gz --segment gene --split 2 --out /tmp/set') == 0)
		assert(os.path.getsize('/tmp/set-0.gff3') == 98119)
		
		assert(os.system('milwa --fasta /tmp/set-0.fa --gff /tmp/set-0.gff3 --model don --canonical --first --hmm /tmp/donor.hmm') == 0)
		assert(os.path.getsize('/tmp/donor.hmm') == 2683)
		
		assert(os.system('dumapic --hmm /tmp/donor.hmm --svg /tmp/donor.svg') == 0)
		assert(os.path.getsize('/tmp/donor.svg') == 5623)
		
		assert(os.system('mogref --hmm /tmp/donor.hmm --fasta /tmp/fake.fa --gff /tmp/fake.gff --count 10 --length 200 --seed 1') == 0)
		assert(os.path.getsize('/tmp/fake.gff') == 770)
		
		assert(os.system('milwa --fasta /tmp/set-1.fa --gff /tmp/set-1.gff3 --model don --canonical --first --source /tmp/donors') == 0)
		assert(os.path.getsize('/tmp/donors.fa') == 169307)

		assert(os.system('halito --fasta /tmp/donors.fa --hmm /tmp/donor.hmm > /tmp/out.gff') == 0)
		assert(os.path.getsize('/tmp/out.gff') == 42454)

