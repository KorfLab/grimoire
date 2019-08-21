Tutorial: grimoire
==================

In this tutorial, you will learn how to train HMMs, decode sequences,
and compare annotations. All the programs you will need are in
`grimoire/bin` and all the data files are in `grimoire/data`. Here's the
road map for this tutorial.

1. Setting up
2. Comparing genome annotations with `latumapic`
3. Examining genome annotation with `calfo`
4. Visualizing annotation files with `kandi`
5. Creating training and testing sets with `haman`
6. Building an HMM with `milwa`
7. Decoding sequences with `halito`
8. Comparing predictions with `latumapic`
9. Tuning models to improve accuracy

## Disclaimers ##

Before we get to the walk-through, lets take a brief moment to review
who should and should not use this tutorial. This is a Linux
_command-line_ tutorial. We make the following assumptions about your
general Linux skills:

+ You can modify file permissions and set environment variables
+ You can install python modules (probably in virtual environments)
+ You can clone Git repositories

It is further assumed that you have an interest in _genomics_. If you
don't know what that means, you probably downloaded the wrong software.
Here are some expected users of grimoire.

+ You want to run a gene prediction program 
+ You want to identify potential errors in genome annotation
+ You want to model some kind of sequence feature
+ You want to develop a new algorithm using grimoire

It is further assumed that you have enough bioinformatics knowledge to
know that there are a variety of standard formats and that standard
formats are sometimes interpreted loosely. You're also okay working
around these problems. Here are a few more assumptions:

+ You have worked with FASTA files before
+ You have worked with GFF files before
+ You have used BLAST or similar programs on the command line
+ You can write scripts in Python/Perl to munge files

Finally, while we are excited when other people use our software, we are
an academic lab with limited time and funding. We cannot offer much
support. We appreciate bug reports and feature requests, but we don't
intend to offer general Linux, Python, or Git help.

OK, if you're still around, let's start the tutorial.

## 1. Setup ##

pip3 install matplotlib, graphviz
Mac: brew install graphviz

It is assumed that you either `pip3` installed the grimoire package
(probably in a virtual environment) or you are working with a cloned git
repo and have set your `PYTHONPATH` to include `grimoire` and your
`PATH` to include `grimoire/bin`.

Grimoire depends on the matplotlib and graphviz packages. Make sure
those are installed (possibly in a virtual environment)

	pip3 install matplotlib
	pip3 install graphviz

You also need to install the graphviz executables and have them
available in your `PATH`.

	sudo apt install graphviz # Ubuntu, Debian
	brew install graphviz     # Mac via homebrew

To make sure the `grimoire` pacakage and its dependencies are installed
correctly, run the unit tests.

	python3 setup.py test

If any of the tests fail, it may be because one of the python library is
missing or the graphviz executables are not in your `PATH`. Fix those
before proceeding. Don't continue the tutorial unless all the unit tests
pass cleanly.

In the `grimoire/data` directory you will find the files we will be
using in this tutorial. These sequence and annotation files represent 1%
of the C. elegans genome. The GFF3 and GTF files contain mostly the same
information, but not exactly. The GFF3 has some RNA genes while the GTF
has only coding genes.

	+ `ce270.fa.gz` from WS270
	+ `ce270.gff3.gz` from WS270
	+ `ce270.gtf.gz` from WS270

To keep things tidy, create a working diretory somewhere and make
symlinks from the files above to your working directory (you don't need
to have set `GRIMOIRE` as shown below, so if you didn't, use the path to
wherever the `grimoire` root directory is).

	ln -s $GRIMOIRE/ce*.gz .

## 2. Comparing genome annotations with `latumapic` ##

One of the first problems any bioinformatician faces is how to deal with
"other peoples' data". In this case, we're talking about genome
annotation in some standard-ish format. WormBase is the primary source
for C. elegans files. But you can also pick them up from Ensemble, for
example. Even at WormBase you can find more than one format for the
annotation. The GFF3 file contains _everything_, which may be more than
you want while the GTF file contains just the genes. You might wonder
how different they are. For such a task, we can use `latumapic`.

	latumapic --fasta ce270.fa.gz --file1 ce270.gff3.gz --file2 ce270.gtf.gz --feature gene

The output of this command is the following:

	112942 64247 177189

This indicates that among the 177189 nucleotides in the genome, 112942
are labeled the same, while 64247 are labeled differently. Why would
there be differences in the genes when both versions of the genome are
WS270? Let's take a closer look at the components of genes.

	latumapic --fasta ce270.fa.gz --file1 ce270.gff3.gz --file2 ce270.gtf.gz --feature exon

There output now shows the differences in the exons.

	156527 20662 177189

Let's check the coding sequences (CDS)

	latumapic --fasta ce270.fa.gz --file1 ce270.gff3.gz --file2 ce270.gtf.gz --feature CDS

They are identical. Phew.

	177189 0 177189

As it turns out, I removed a bunch of stuff from the GFF3 file because
it was so huge, and in doing so, I removed some of the non-coding
transcripts present in the GTF. How can you see that? Read on.

## 3. Examining genome annotation with `calfo` ##

Let's get an overview of the C. elegans genome annotation at release 270
with `calfo`. we'll create reports for both gff3 and gtf.

	calfo --fasta ce270.fa.gz --gff ce270.gff3.gz --title ce270 --html ce270gff3.html
	calfo --fasta ce270.fa.gz --gff ce270.gff3.gz --title ce270 --html ce270gtf.html

Open the html pages in your favorite browser. You'll notice that the
GFF3 and GTF files have different feature types. For example, GTF
includes start and stop codons. Also, the GFF3 has 'mRNA' features while
the GTF has 'transcript' features. When annotation files are read into
`grimoire`, the features are converted to GFF3 as much as possible, so
both 'mRNA' and 'transcript' become 'mRNA' internally.

## 4. Visualizing annotation files with `kandi` ##

One of the bioinformatician's most useful debugging tools is their
powers of observation. However, it is often difficult to observe genome
annotations when they are represented as thousands of lines of text.
That's why there are genome browsers like Jbrowse, Gbrowse, IGV, IGB,
Apollo, etc. Genome browsers are usually the best way to examine
annotation at arbitrary resolution. But setting them up can be a bit of
a pain because the files have to be indexed properly or uploaded to a
website. So we're going to use the really simple genome browser `kandi`,
which you can find in `grimoire/bin`. This `vi`-like browser is
terminal-based and driven with the keyboard, not mouse. `kandi` can be
useful in a debugging or tutorial setting but is not intended for use on
whole genomes. Let's open up the annotation files and compare them in
`kandi`.

	kandi ce270.fa.gz ce270.gff3.gz ce270.gtf.gz

The program starts out in 'gene' view. On the left-hand quarter of the
display, you should see a purple glyph that is present in the GTF file
that isn't in the GFF file. If you hit 't' to go to 'transcript' view,
you can see that the feature is glpyh is colored yellow. Coding segments
are green and this isn't coding. If you move the cursor over that yellow
glyph and hit 'return', you can zoom in on it. If you zoom in more, you
will eventually see its DNA.

## 5. Creating training and testing sets with `haman` ##

When building a gene-finder or other sequence decoder, we need some kind
of training set. One convenient source is a closely related genome that
has already been annotated, In our case, we're going to use the previous
WS270 annotations. If we want to know how well our decoder performs,
we also need a test set. We're going to do the very simple thing of
splitting our data into halves: one half for the training, one set for
the testing.

	haman --fasta ce270.fa.gz --gff ce270.gff3.gz --segment gene --split 2 --out set

The previous command splits our data set into 4 files. Those marked
set-0 will be used for training while set-1 is for testing. To make this
abundantly clear, let's alias them.

	ln -s set-0.fa train.fa
	ln -s set-0.gff3 train.gff3
	ln -s set-1.fa test.fa
	ln -s set-1.gff3 test.gff3

Take a quick look at the training set with `kandi`. You'll find that
each gene now has its own piece of DNA with no other genes on it. But
there may be more than one transcript per gene. The third one down has
two transcripts, for example. Also, by default, there are 100 bp
flanking each gene. You may want more or less than that, which you can
do with the `--padding` parameter in `haman`.

## 6. Building an HMM with `milwa` ##

For debugging purposes, you should run `dumapic` to see your HMM.

For fun, you can generate sequences with `morlis`.

## 7. Decoding sequences with `halito` ##

## 8. Comparing predictions with `latumapic` ##

Use `kandi` too if you want to look at some differences.

## 9. Tuning models to improve accuracy ##

+ Lexicalized emissions
+ Optimal lengths of features
+ Cross validation

## 10. Other genomes ##

+ `at10.fa.gz` from TAIR10
+ `at10.gff.gz` from TAIR10
+ `at11.bed.gz` from Araport11

+ `ce271.gtf.gz` from WS271 via Ensembl
