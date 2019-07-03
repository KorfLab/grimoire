
import sys
import math
import random
import re

import decode
import genome
import toolbox
import sequence



## build up genes with new methods



#sys.exit(1)



## convert bed12 to GFF3

fp = open('araport.bed', 'r')
source = 'ARAPORT11'
# genes can't be created until all transcripts are created
genes = {}
while (1):
	line = fp.readline()
	if line == '': break
	col = line.split('\t')
	chr_id = col[0]
	chr_beg = int(col[1]) + 1
	chr_end = int(col[2])
	txid = col[3]
	score = col[4]
	strand = col[5]
	cds_beg = int(col[6]) + 1
	cds_end = int(col[7])
	rgb = col[8]
	n = int(col[9])
	sizes = col[10].split(',')
	starts = col[11].split(',')
	gid = re.search('(\w+)\.\d+', txid)[1]
	attr = 'ID=' + txid + ';Parent=' + gid
	
	if gid not in genes:
		print('\t'.join([chr_id, source, 'gene', str(chr_beg), str(chr_end),
			score, strand, 'ID='+gid]))
		genes[gid] = True
	
	print('\t'.join([chr_id, source, 'mRNA', str(chr_beg), str(chr_end),
		score, strand, attr]))

	for i in range(n):
		beg = chr_beg + int(starts[i])
		end = beg + int(sizes[i]) -1
		attr = 'Parent=' + txid
		print('\t'.join([chr_id, source, 'exon', str(beg), str(end),
			score, strand, attr]))

sys.exit(1)


## build from GFF3
gen = genome.Genome(fasta='data/TAIR10_1.fasta', gff3='data/TAIR10_1.gff3')
for chr in gen.chromosomes:
	for gene in chr.genes:
		if gene.issues:
#			print(gene.id, 'has issues')
			for tx in gene.transcripts:
				if tx.issues:
#					print('\t', tx.id, 'has issues')
					for issue in tx.issues:
#						print('\t\t', issue)
						pass


