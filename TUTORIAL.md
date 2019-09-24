Tutorial: grimoire
==================

In this tutorial, you will learn how to train HMMs, decode sequences,
and compare annotations. All the programs you will need are in
`grimoire/bin` and all the data files are in `grimoire/data`. Here's the
road map for this tutorial.

1. Setting up
2. Summarizing genome annotation with `calfo`
3. Visualizing annotation files with `kandi`
4. Creating training and testing sets with `haman`
5. Building an HMM with `milwa`
6. Visualizing an HMM with `dumapic`
6. Decoding sequences with `halito`
7. Comparing predictions with `morlis`
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
correctly, run the unit tests in the grimoire root directory.

	python3 setup.py test

If any of the tests fail, it may be because one of the python libraries
is missing or the graphviz executables are not in your `PATH`. Fix those
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

	ln -s $GRIMOIRE/data/ce*.gz .

## 2. Summarizing genome annotation with `calfo` ##

One of the first problems any bioinformatician faces is how to deal with
"other peoples' data". In this case, we're talking about genome
annotation in some standard-ish format. WormBase is the primary source
for C. elegans files. But you can also pick them up from Ensemble, for
example. Even at WormBase you can find more than one format for the
annotation. The GFF3 file contains _everything_, which may be more than
you want while the GTF file contains just the genes.

Examine the annotation files with `zless` or whatever (but if that
'whatever' is you uncompressing the file and then opening it with Word,
you should consider a change in career). You'll notice that the GFF3 and
GTF files have different names for the same feature (e.g. mRNA is the
same as transcript) and also different ways of describing the same
feature (e.g. CDS features contain stop codons in GFF3 but not GTF).

Let's get an overview of the C. elegans genome annotation at release 270
with `calfo`. We'll create reports for both gff3 and gtf.

	calfo --fasta ce270.fa.gz --gff ce270.gff3.gz --title ce270gff --html ce270gff3.html
	calfo --fasta ce270.fa.gz --gff ce270.gtf.gz --title ce270gtf --html ce270gtf.html

Open the html pages in your favorite browser. In Figure 2, you can see
that the files specify different feature types. But the figures
afterward all use the same terminology (i.e. mRNA instead of
transcript). That's because grimoire considers GFF3 to be its native
format and automatically converts GTF (and some other formats) to GFF3.

## 3. Visualizing annotation files with `kandi` ##

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
display at about 23000, you should see a purple glyph that is present in
the GTF file that isn't in the GFF file. If you hit 't' to go to
'transcript' view, you can see that the is glpyh is colored yellow.
Coding segments are green and this isn't coding. If you move the cursor
over that yellow glyph and hit 'return', you can zoom in on it. If you
zoom in more (with the '+' key) you will eventually see the DNA.

So the GFF3 and GTF files are pretty similar but not identical. For our
purposes, it doesn't matter, so we'll use the GFF3 file from now on.

## 4. Creating training and testing sets with `haman` ##

When building a gene-finder or other sequence decoder, we need some kind
of training set. One convenient source is a closely related genome that
has already been annotated, In our case, we're going to use the previous
WS270 annotations. If we want to know how well our decoder performs, we
also need a test set. We're going to do the very simple thing of
splitting our data into halves: one half for the training, one half for
the testing.

	haman --fasta ce270.fa.gz --gff ce270.gff3.gz --segment gene --split 2 --out set

The previous command splits our data set into 4 files. Those marked
set-0 will be used for training while set-1 is for testing. To make this
abundantly clear, let's alias them.

	ln -s set-0.fa train.fa
	ln -s set-0.gff3 train.gff3
	ln -s set-1.fa test.fa
	ln -s set-1.gff3 test.gff3

Take a quick look at the training set with `kandi`.

	kandi set-0*

You'll find that each gene now has its own piece of DNA with about 100
bp upstream and downstream (negative strand genes have been flipped and
have 101 bp rather than 100 - for internal debugging reasons). You may
want more or less than that, which you can do with the `--padding`
parameter in `haman`. If you look at a few genes, you'll not that there
may be more than one transcript per gene. The third gene in the list has
two transcripts, for example. In the real, messy biological world, a
gene may produce several transcripts, some of which may be quite rare.
However, in the computer world, this complexity is usually distilled
down to the simple rule that a gene creates exactly one transcript. This
makes training, testing, and evaluation much simpler. In this tutorial,
we will continue on with that tradition. However, this is an
oversimplification of the underlying biology and while grimoire is
capable of doing more complex things, those are outside the scope of
this tutorial. 

## 5. Building an HMM with `milwa` ##

There are several simple models that can be built with `milwa`. We're
going to build a model of the splice donor site and some flanking
sequence.

	milwa --fasta train.fa --gff train.gff3 --model don --canonical --first --hmm donor.hmm

In the command line above, `--model don` indicates we want to build the
donor model, `--canonical` means we only want canonical sequences (e.g.
donor sites starting with 'GT', `--first` means we only want the first
transcript if there are more than one, and `--hmm donor.hmm` specifies
the name of the output file.

Examine the `donor.hmm` file with `less` or whatever and you'll see that
it is formatted as a JSON document.

## 6. Visualizing an HMM with `dumapic` ##

When you create your own HMMs, it's easy to mess up the state
connections. So it's a good idea to visualize the state connection
diagram with `dumapic`.

	dumapic --hmm donor.hmm --svg donor.svg

You can view the `donor.svg` file with a variety of web browsers and
graphics programs. ImageMagick works well for converting to png or pdf.

HMMs are generative models, so in that spirit, grimoire includes a
program, `mogref`, to generate random sequences consistent with a model.
Feel free to skip this next step as it's just included 'for fun'.

	mogref --hmm donor.hmm --fasta fake.fa --gff fake.gff --count 10 --length 200 --seed 1

Note that the names of the features in the GFF file are not actually
following the GFF3 specification.

## 7. Decoding sequences with `halito` ##

The HMM we built with `milwa` modeled splice donor sites, but we don't
have a collection of splice donor sites to decode. We'll create one now
with `milwa` but this time we will save the sequences and not the HMM.

	milwa --fasta test.fa --gff test.gff3 --model don --canonical --first --source donors

You will now have two new files in your working directory: `donors.fa`
and `donors.gff`. Inspect these with `less` to make sure they look as
expected.

Now it's finally to decode some sequences with the HMM we built. To do
that, we use `halito`. You can run this multi-threaded, but the HMM is
so simple and the sequences are so short that it isn't worth the
overhead.

	halito --fasta donors.fa --hmm donor.hmm > out.gff

## 8. Comparing predictions with `morlis` ##

To compare the predictions `out.gff` with the test set `donors.gff` we
use `morlis`. Right now, development of `morlis` is highly volatile as
there are many ways one can interpret correctness. For this reason, the
following code may or may not do anything useful.

	morlis --fasta donors.fa --file1 donors.gff --file2 out.gff

## 9. Tuning models to improve accuracy ##

The HMM we built in this tutorial was somewhat terrible. The reason for
that is that each state was specified with a naive emission model. HMMs
become much more accurate when state emissions have context. A clear
example of this is that CDS states need a context of at least 2 in order
to prevent in-frame stop codons (i.e. the probability of emitting an 'A'
from the 3rd state of a codon should be zero if the previous emissions
were 'TA').

To improve the model, we might also want to make the donor site longer
to capture more of the surrounding sequence context. How long the donor
site should be depends on the specific genome.

The data set here was very small. There were only ~65 sequences in each
of the training and testing sets. When using such small data sets,
splitting them 50/50 isn't a great idea. It's better to use
cross-validation or jack-knifing.

