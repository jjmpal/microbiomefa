# Gut microbiome and fibrillatio atriorum

This repository contains our R code for gut microbiome and FA analysis. The file *microbiomefa.Rmd* is the main analysis while different helper functions are defined in separate *R* files.

The code calculates 

- cross-sectional microbial diversity (alpha and beta) at species level
- draws cross-sectional PCoA plots at species level
- checks cross-sectional associations between top bacterial genera and FA
- Cox models for participants that don't have FA at baseline
  - associations between top genera and FA
  - associations between baseline PCA-axes and FA
  - associations between baseline alpha diversity and FA
