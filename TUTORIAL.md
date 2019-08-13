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

+ You can modify file permissions and set environment variables
+ You can install python modules (probably in virtual environments)
+ You can clone Git repositories
+ You will not email the authors asking for Unix/Python/Git help

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
around these problems.

+ You have worked with FASTA files before
+ You have worked with GFF files before
+ You have used BLAST or similar programs on the command line
+ You can write scripts in Python/Perl to munge files

OK, if none of this scared you off, it's time to start the tutorial.

---

Overview
--------

1. Setup
2. Checking annotation files with `kalki`
3. Reporting genome summary with `calfo`

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

2. Checking annotation files with `kalki`
-----------------------------------------

While we may wish all genome annotation was provided in one standard
format, this isn't the case. Even where such standards exist, such as
GFF3, it's unwise to assume that all GFF3 is formatted identically or
properly. Before we can do anything useful with the grimoire tools, we
have to check the annotation file(s) to make sure they will work with
the various grimoire tools.

In the `grimoire/data` directory you will find the files we will be
using in this tutorial. These sequence and annotation files represent 1%
of the genome.

+ `at10.fa.gz` from TAIR10
+ `at10.gff.gz` from TAIR10
+ `at11.bed.gz` from Araport11
+ `ce270.fa.gz` from WS270
+ `ce270.gff.gz` from WS270
+ `ce270.gtf.gz` from WS270
+ `ce271.gtf.gz` from WS271 via Ensembl





3. Examining genome annotation with `calfo`
-------------------------------------------



4. Creating training and testing sets with `haman`
--------------------------------------------------

5. Building an HMM with `milwa`
-------------------------------

For debugging purposes, you should run `dumapic`.
For fun, you can generate sequences with `morlis`.

6. Decoding sequences with `halito`
-----------------------------------

7. Comparing annotations with `latumapic`
-----------------------------------------

8. Tuning models to improve accuracy
--------------------------------

+ Lexicalized emissions
+ Optimal lengths of features


