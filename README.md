# MODY GENES <img src="https://avatars1.githubusercontent.com/u/70200041?s=400&u=4a50dc780ca9f22de9dd34530836b48972965275&v=4" align="right" width="120" />
<!-- badges: start -->
[![DOI](https://zenodo.org/badge/290171567.svg)](https://zenodo.org/badge/latestdoi/290171567)
<!-- badges: end -->

This repository contains all code and scripts used to perform the analyses described in the paper "Pathogenic variants in actionable MODY genes associate with type 2 diabetes".

The repository includes both [bioinformatics](bioinformatics) and [biostatistics](biostatistics) analyses.

**Note: the code has not been tested outside of its original environment and depends on data not publicly available.**

## Bioinformatics

The bioinformatics analysis was performed using the scripts in the [bioinformatics](bioinformatics) directory.  
The different analysis steps are listed in the plain text file [bioinformatics/template_analysis](bioinformatics/template_analysis), which is provided to the [bioinformatics/start_multiserveur.pl](bioinformatics/start_multiserveur.pl) script.  
Then the [bioinformatics/start_sample.pl](bioinformatics/start_sample.pl) script is launched with all the different steps of the analysis.

All the different scripts to perform each steps can be found in the [bioinformatics/Include](bioinformatics/Include) directory.  
These scripts use local databases:  
- MySQL databases with informations about the different samples (*e.g.*, sample ID, sex, demultiplexing index and lane, raw data path) and the dbSNP annotations (version 135).  
- A MongoDB database with the dbNSFP annotations (version 3.0).  
- A Redis database with mutation counting informations.

Tools, captures and sequencing adaptors paths can be passed by the files in the [bioinformatics/Config](bioinformatics/Config) directory.

Perl version used is Perl 5, version 16, subversion 3 (v5.16.3) on CentOS Linux 7 (Core, x86_64).

## Biostatistics

The statistical analyses were performed using the scripts in the [biostatistics](biostatistics) directory.  
All scripts are numbered according to the order of use.  
The directory [biostatistics/utils](biostatistics/utils) contains additional scripts needed in the analyses.

The R version used is the 3.6.3 (2020-02-29) on Debian GNU/Linux 9 (stretch, x86_64) available as a Docker image ([umr1283/stat:R363](https://hub.docker.com/r/umr1283/stat)).

## Contact

For questions and other discussion, please contact the authors.
