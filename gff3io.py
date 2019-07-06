
import argparse
import re

parser = argparse.ArgumentParser(description='GFF3 gene feature selector')
parser.add_argument('--file', required=True, type=str,
	metavar='<path>', help='path to annotation file (%(type)s)')
parser.add_argument('--source', required=True, type=str,
	metavar='<string>', help='wb|fb|tair|ap (%(type)s)')
arg = parser.parse_args()

types = {
	'gene' : True,
	'mRNA' : True,
	'exon' : True,
	'CDS' : True,
	'five_prime_utr' : True,
	'three_prime_utr' : True,
}

if arg.source == 'wb':
	with open(arg.file) as file:
		for line in file:
			if line[0:1] == '#': continue
			col = line.rstrip().split('\t')
			if col[1] == 'WormBase' and col[2] in types:
				m = re.search('pseudogene', col[8], re.IGNORECASE)
				if m: continue
				print(line, end='')

if arg.source == 'fb':
	with open(arg.file) as file:
		for line in file:
			if line[0:1] == '#':
				if line[0:7] == '##FASTA':break
				else: continue
			col = line.rstrip().split('\t')
			if col[1] == 'FlyBase' and col[2] in types:
				# SO:0000043 is the SO tag for 'processed_psedudogene'
				m = re.search('SO:0000043', col[8])
				if m: continue
				print(line, end='')

if arg.source == 'tair':
	with open(arg.file) as file:
		for line in file:
			if line[0:1] == '#': continue
			col = line.rstrip().split('\t')
			if col[2] in types:
				m = re.search('pseud', col[2])
				if m: continue
				print(line, end='')


if arg.source == 'ap':
	with open(arg.file) as file:
		pass

"""
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
"""