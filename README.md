# Batch Adjustment by Reference Alignment (BARA): Improved prediction performance in biological test sets with batch effects

A development r-package implementing the BARA algorithm can be found at [robingradin/bara](https://github.com/robingradin/bara)

## About
Batch effects are a common problem associated with biological data acquisition techniques.
They are especially troublesome in predictive settings were fixed training sets are used
to classify subsequently acquired test sets. In these settings, the unknown biological 
conditions are completely correlated with batch.

BARA was developed as an analytical tool to help adjust for batch effects between
training sets and test sets. It uses a few reference samples to make the adjustment 
between training sets and test sets in a compressed data space spanned by the training set.

This repository contains the code used to perform the analyses described in
Gradin et. al. (submitted manuscript). Because the expression matrices are too large to 
be shared via GitHub, the processed datasets can be downloaded via the link below.
Alternatively, the data processing can be manually performed by downloading the
individual datasets specified in the manuscript. The processed datasets is contained
within a named list. Each entry in the list contains one expression matrix (exprs)
and one annotation matrix (pheno).

 * [Link to processed datasets at Dropbox](https://www.dropbox.com/s/h4exq3urg37cstm/processed_datasets.RDS?dl=0)

To generate the data for the comparison between the methods, first set the
working directory to this repository. Then, run the script
`src/assessnormalization/run_all_scripts.R`. To generate the results examining 
the effect of varying the number of reference samples used in BARA, run 
the script `src/bara_ref_samples/bara.R`.
