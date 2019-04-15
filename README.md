# Deep_binning
A program for unsupervised deep learning of contigs features. The result of Deep_binning could be used for binning process for other binning tools.

# Basic usage

Deep_binning is writen in Python and respective TensorFlow should be installed

Python encoder.py -i inputfile -n num_hidden -e iteration_number -l learning_rate

# Parameters recomendation

-i: inputfile must be in csv file and the first column is the contigs namelist

-n: number of hidden layers, we recomend to set this value with about the half value of raw feature dimension. Take the 136 dimensional 4-mer frequency as an example, values between 40-70 are all ok for -n   

-e: iteration number, 100 or 200

-l: learning rate, 0.01 or 0.001

Different combinations of -n, -e and -l should to tried to train the Deep_binning. The standard of a well-trained model is simple, when the reconstruction error convergence is observed, the outputfile of Deep_binning is ready for use of subsequent binning process.

# Output:
inputfile_encoder.csv file with the representations of each contig.

# Support
If you are having trouble running Deep_binning or interpretting any results, please don't hesitate to contact us: liuyun313@jlu.edu.cn.
