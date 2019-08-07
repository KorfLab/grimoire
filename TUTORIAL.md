Tutorial: grimoire
==================

In this tutorial, you will learn how to train HMMs, decode sequences,
and compare annotations. All the programs you will need are in
grimoire/bin and all the data files are in grimoire/data. Here's the
road map.

1. Sanitizing annotation files with `kalki`
2. Examining genome annotation with `calfo`
3. Creating training and testing sets with `haman`
4. Building an HMM with `milwa`
5. Decoding sequences with `halito`
6. Comparing annotations with `latumapic`
7. Tuning models to improve accuracy

Before we get to the walk-through, lets take a brief moment to review
who should and should not use this tutorial. This is a command-line
tutorial. We make the following assumptions:

+ You have enough Unix knowledge to modify your PATH and permissions
+ You know how to set enviornment variables
+ You have enough Python knowledge to install grimoire
+ You will not email the authors asking for Unix/Python help

It is further assumed that you have an interest in _genomics_. If you
don't know what that means, you probably downloaded the wrong software.
Here are some expected users of grimoire.

+ You want to build a gene prediction program 
+ You want to identify potential errors in genome annotation
+ You want to model some kind of sequence feature
+ You want to develop something new using grimoire

It is further assumed that you have enough bioinformatics knowledge to
know that there are a variety of standard formats and that standard
formats are sometimes interpreted loosely. You're also okay working
around these problems.

+ You have worked with FASTA files before
+ You have used BLAST or similar programs on the command line
+ You can write scripts in Python/Perl to munge files

OK, if none of this scared you off, it's time to start the tutorial.

---

1. Sanitizing annotation files with `kalki`
-------------------------------------------

While we may wish all genome annotation was provided in one standard
format, this isn't the case. Even where such standards exist, such as
GFF3, it's unwise to assume that all GFF3 is formatted properly. Before
we can do anything useful with the grimoire tools, we have to put the
annotation file(s) through a sanitizer that will format the downloaded
file into the preferred format, which for grimoire happens to be GFF3.

In the `grimoire/data` directory you will find the 6 files we will be
using today. These represent 1% of the genome from Arabidopsis thaliana
and Caenorhabditis elegans. The fasta files contain the chromosomes
while the gff, bed12, and gtf files contain the annotations.

	A.thaliana.1percent.fasta.gz
	A.thaliana.1percent.gff3.gz
	araport11.1percent.bed12
	C.elegans.1percent.fasta.gz
	C.elegans.1percent.gff3.gz
	C.elegans.1percent.gtf.gz

`A.thaliana.1percent.gff3.gz` comes from TAIR10 while
`araport11.1percent.bed12` comes from Araport, which is effectively
TAIR11. The file format is completely different and the annotations may
also be different (see below). Importantly, the sequence did not change
from TAIR10 to TAIR11 so we can use the same sequence with both
annotations. A similar situation exists for the C.elegans files, where
the sequence is the same, but the annotation has different versions and
formats (gff3 is WS270 from WormBase, gtf is WS271 from Ensembl).

Before we start running programs, let's create a directory where we can
do all the work. To make life easy, first create a GRIMOIRE environment
variable and set it to the location of the grimoire directory (this
directory contains `bin` and `grimoire` among other things. Depending on
how you installed grimoire, you may want to append $GRIMOIRE to your
$PYTHONPATH and $GRIMOIRE/bin to your $PATH. Again, if any of this is
even remotely abstruse to you, you may want to improve your Unix skills
before proceeding.

Let's create a few symlinks so that the command lines in this tutorial
aren't so long.

	ln -s $GRIMOIRE/data/A.thaliana.1percent.fasta.gz at1.fa.gz
	ln -s $GRIMOIRE/data/A.thaliana.1percent.gff3.gz at1.gf.gz
	ln -s $GRIMOIRE/data/araport11.1percent.bed12 at1.bed

	ln -s $GRIMOIRE/data/C.elegans.1percent.fasta.gz ce1.fa.gz
	ln -s $GRIMOIRE/data/C.elegans.1percent.gff3.gz ce1.gf.gz
	ln -s $GRIMOIRE/data/C.elegans.1percent.gtf.gz ce1.gt.gz
	
The first task is to sanitize the annotation files with `kalki` to
ensure that they will operate with the downstream tools. In some cases,
the output may be identical to the input.

	gunzip -c at1.gf.gz | kalki --source tair > at1.gff3
	cat at1.bed | kalki --source ap > at2.gff3
	
	gunzip -c ce1.gf.gz | kalki --source wb > ce1.gff3
	gunzip -c ce1.gt.gz | kalki --source gtf > ce2.gff3

Since we've created 2 different annotation files for each genome, we
might wonder if they are in fact different? Did any genes change from
TAIR10 to TAIR11 or from WS270 to WS271. Well, we only have 1% of each
genome, so we can't answer that fully, but we can check that 1% with
`latumapic`.






	latumapic --fasta ce1.fa.gz --file1 ce1.gf.gz --file2 ce1.gff3



2. Examining genome annotation with `calfo`
-------------------------------------------



3. Creating training and testing sets with `haman`
--------------------------------------------------

4. Building an HMM with `milwa`
-------------------------------

For debugging purposes, you should run `dumapic`.
For fun, you can generate sequences with `morlis`.

5. Decoding sequences with `halito`
-----------------------------------

6. Comparing annotations with `latumapic`
-----------------------------------------

7. Tuning models to improve accuracy
--------------------------------

+ Lexicalized emissions
+ Optimal lengths of features


