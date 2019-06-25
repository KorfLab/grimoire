
import genome

gen = genome.Genome(gff='data/TAIR10_1.gff3', fasta='data/TAIR10_1.fasta')

txn = 0
txa = []
for chr in gen.chromosomes:
	for gene in chr.genes:
		if gene.issues:
			print(gene.id, 'has issues')
			for tx in gene.transcripts:
				if tx.issues:
					print('\t', tx.id, 'has issues')
					for issue in tx.issues:
						print('\t\t', issue)
				else:
					txn += 1
					txid = 'tx-' + str(txn)
					parent = genome.Feature(chr, tx.beg, tx.end, tx.strand,
						'FOO', id=txid)
					for exon in tx.exons:
						child = genome.Feature(chr, exon.beg, exon.end,
							exon.strand, 'BAR', parent=txid)
						parent.add_child(child)
					txa.append(parent)

for tx in txa:
	print(tx.gff())