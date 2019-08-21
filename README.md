# Deep_binning
A program for unsupervised deep learning of contigs features. The result of Deep_binning could be used for binning process for other binning tools.

# Files description 

Deep_binning.py is the Python script for feature learning.

kmer.py is the Python script to count k-mer features for contigs in fasta format.

*.CSV files are output files of Deep_binning of metagenomes used in our paper.

# Basic usage

Deep_binning is writen in Python and respective TensorFlow should be installed

Python Deep_binning.py -i inputfile -n num_hidden -e iteration_number -l learning_rate

# Inputfile preparation 

The inputfile could be the k-mer feature or coverage-feature written in .CSV format, in which each column represents a k-mer and each row is a contig. K-mer feature calculation could refer to kmer.py in this repository. Coverage feature could refer to CONCOCT packages in GitHub or you can contact us for help.

# Parameters recomendation

-i: inputfile must be in csv file and the first column is the contigs namelist

-n: number of hidden layer's neurons, we recomend to set this value with about the half value of raw feature dimension. Take the 136 dimensional 4-mer frequency as an example, values between 40-70 are all ok for -n   

-e: iteration number, 100 or 200

-l: learning rate, 0.01 or 0.001

Different combinations of -n, -e and -l should to tried to train the Deep_binning. The standard of a well-trained model is simple, when the reconstruction error convergence is observed, the outputfile of Deep_binning is ready for use of subsequent binning process.

# Output:
The output file is named 'inputfile'_encoder.csv file with the representations of each contig. This file has the same structure with the input file, but has less number of columns, which is equal to the number of hidden layer's neurons you have set.

# Support
If you are having trouble running Deep_binning or interpretting any results, please don't hesitate to contact us: liuyun313@jlu.edu.cn.
