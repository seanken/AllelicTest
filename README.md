# AllelicTest

This code implements a privacy preserving GWAS method discussed. In particular, for a given number of SNPs (mret), a given privacy budget (epsilon) and a given a genotype file (.bed file, where the .bed file includes phenotype information), the algorithm returns m_ret SNPs with high values of the allelic test statistic, but does so in an epsilon differentially private way. This implementation uses the algorithm introduced in the manuscript "Realizing Privacy Preserving Genome-wide Association Studies" (published in Bioinformatics, 2016) to implement the neighbor based method for picking high scoring SNPs. Details are given in the manuscript

# How to run the code

Assume we have a .bed file, filname.bed, as well as the corresponding .bed file, etc. Assume that you have a privacy budget of epsilon=1 (see paper for details), and want a list of mret=3 high scoring SNPs. To run the program type

python FastNeigh.py filename 1 3

Our implementation requires numpy, scipy, and pysnptools to be installed.

# Files

PrivPick.py: Implements the internals of our algorithm, as well as the score and noise based algorithms for returning high scoring SNPs privately.

FastNeigh.py: A simple user interface for playing with PrivPick.py.

loadFile.py: Used to load the file of interest.

Note that this implementation is mostly meant as a proof of concept, so care should be taken when using it in practice.