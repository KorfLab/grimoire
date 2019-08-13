Tutorial: grimoire
==================

In this tutorial, you will learn how to train HMMs, decode sequences,
and compare annotations. All the programs you will need are in
`grimoire/bin` and all the data files are in `grimoire/data`. Here's the
road map for this tutorial.

1. Setting up
2. Examining genome annotation with `calfo`
3. Visualizing annotation files with `kalki`
4. Creating training and testing sets with `haman`
5. Building an HMM with `milwa`
6. Decoding sequences with `halito`
7. Comparing predictions with `latumapic`
8. Tuning models to improve accuracy

Disclaimers
-----------

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

---

OK, if you're still around, let's start the tutorial.

---

1. Setup
--------

It is assumed that you either `pip3` installed the grimoire package
(probably in a virtual environment) or you are working with a cloned git
repo and have set your `PYTHONPATH` to include `grimoire` and your
`PATH` to include `grimoire/bin`. To make sure the pacakge is installed
correctly, run the unit tests.

	python3 setup.py test

If any of the tests fail, it may be because some library is missing. Fix
those before proceeding. Don't continue the tutorial unless all the unit
tests pass cleanly.

In the `grimoire/data` directory you will find the files we will be
using in this tutorial. These sequence and annotation files represent 1%
of each genome.

+ `at10.fa.gz` from TAIR10
+ `at10.gff.gz` from TAIR10
+ `at11.bed.gz` from Araport11
+ `ce270.fa.gz` from WS270
+ `ce270.gff.gz` from WS270
+ `ce270.gtf.gz` from WS270
+ `ce271.gtf.gz` from WS271 via Ensembl

To keep things tidy, create a working diretory somewhere and make
symlinks from the files above to your working directory.

2. Examining genome annotation with `calfo`
-------------------------------------------

This is now the sanity check as well as the reporter


3. Visualizing annotation files with `kalki`
--------------------------------------------

Or is this a really bad idea

4. Creating training and testing sets with `haman`
--------------------------------------------------

5. Building an HMM with `milwa`
-------------------------------

For debugging purposes, you should run `dumapic`.
For fun, you can generate sequences with `morlis`.

6. Decoding sequences with `halito`
-----------------------------------

7. Comparing predictions with `latumapic`
-----------------------------------------

8. Tuning models to improve accuracy
--------------------------------

+ Lexicalized emissions
+ Optimal lengths of features


