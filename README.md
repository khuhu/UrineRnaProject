Urine Rna Project
============================

The Rmd file for the analysis gives the basic format/workflow for how whole-urine (DRE) expression data from multiplexed targeted sequencing was analyzed. Study was based on the original MiPs (Michigan Prostate Score) but utilized IonTorrent panel instead of TME - note same type of normalization by KLK3 transcript.

All input for different analyses was based on KLK3 normalized counts. Cross-validation was performed based on a 2:1 split of the training and validation data- based on proportion of samples in classification groups: (1) benign/grade-group 1 (2) grade-group 3 or higher. Random forest model via R package VSURF was used for variable selection - not for model building as interdisciplinary team favored some type of regression model. Glmnet was used (caret), as it's a regression model, but also employs machine learning by optimizing both L1 and L2 (Ridge and Lasso) distances through randomization of hyper parameters. Model buidling used five-fold cross-validation of the training set.

All finalized models (used in PMID 33812851) are presented as Robjs.