# PRNI 2017 Tutorial on Site Effects

The sample correlation matrix between fMRI BOLD time-series serves as a an important starting point for many types of network models for functional connectivity (FC) analyses. In recent years, functional connectivity research increasingly combines data collected across multiple sites. However, the presence of systematic measurement errors due to scanner and session specific artifacts can bias the entries of the correlation matrix, thus creating non-biological differences across sites. In addition to session effects, different sites possess demographic variations that introduce additional changes to correlation matrices. Thus, both session specific as well as site-specific effects on fMRI FC are instances of batch effects in multi-site fMRI network analyses. This tutorial will introduce nonparametric metrics to visualize and quantify site effects in functional connectivity as well as methods to mitigate such effects. In particular, this tutorial will emphasize how reducing session level artifacts also reduce site effects.

### Citation

Narayan, Manjari and Maron-Katz, Adi, *Tutorial: Methods to Diagnose and Ameliorate Site Effects in BOLD Functional Connectivity*, 7th International Workshop on Pattern Recognition in Neuroimaging, 2017

[![DOI](https://zenodo.org/badge/94818742.svg)](https://zenodo.org/badge/latestdoi/94818742)


### Description

This repository contains iPython notebooks that run on the octave kernel with matlab code snippets. They accompany the slides for the two lectures below.

### Authors

- Manjari Narayan 
- Adi Maron-Katz

### Folder Organization

    |--prni2017-site-effects/
      |-- setup.m
      |-- README.md
      |-- data/
      |-- notebooks/
      |-- src/
          |-- external/ 

### Usage

The notebooks can be run on cocalc.com

- Create a free account on https://cocalc.com
- Go to bit.ly/prnisite2017
- Copy the `prni2017-site-effects.zip` into your own personal project. 
(Note: You probably cannot download the large zip file without a paid subscription)

### Installation

**Requirements**: 

MATLAB 2015b + or Octave Version 4.0 (Other versions not guaranteed, but should be compatible)

**Installaion**: 

- Download and unzip the latest release. 
- Set `USE_OCTAVE` variable to false if necessary and run `setup.m` before using notebooks. 

## Tutorials

#### Tutorial 1: Visualizing and Detecting site effects [[Slides]](https://www.swipe.to/8776dXkd8mc98shXPZnWg32S6d8g) 
    
The sample correlation matrix between fMRI BOLD time-series serves as a an important starting point for many types of network models for functional connectivity (FC) analyses. In recent years, functional connectivity research increasingly combines data collected across multiple sites. However, the presence of systematic measurement errors due to scanner and session specific artifacts can bias the entries of the correlation matrix, thus creating non-biological differences across sites. In addition to session effects, different sites possess demographic variations that introduce additional changes to correlation matrices. Thus, both session specific as well as site-specific effects on fMRI FC are instances of batch effects in multi-site fMRI network analyses. This tutorial will introduce nonparametric metrics to visualize and quantify site effects in functional connectivity and briefly review related literature.

#### Tutorial 2: Methods to Ameliorate site effects on a correlation matrix [[Slides]](https://www.swipe.to/8776dXkd8mc98shXPZnWg32S6d8g?p=5vD0QlJHr)

The need to address batch effects is a common problem in many areas including genomics. We will briefly review methods from the neuroimaging and genomics literature to address this problem. We then focus on two relevant statistical techniques to ameliorate inter-site differences in the BOLD correlation matrix. The first approach is successive normalization that is commonly used in genomics. Successive normalization is an iterative estimator that effectively z-scores both rows and columns of the data matrix. This procedure ensures that the BOLD observations or rows of the data matrix across subjects within the same site and between sites are standardized and thus on the same "playing field". We will discuss the relationship between this approach and global signal removal. The second approach uses conditional correlation estimator such that the usual correlation matrix can be decomposed into a true BOLD correlation matrix and a nuisance correlation matrix. Consequently, depending on the level of analysis, this approach can produce an estimate of either a subject or session-specific or a site-specific nuisance correlation matrix. As a by-product of this method, we also obtain a formal nuisance-to-signal ratio that quantifies the contribution of the nuisance correlation to the uncorrected correlation matrix. While Pearson correlation statistics between pairs of fMRI time-series yield the maximum likelihood estimator for the correlation matrix under a multivariate normal model, the two techniques introduced in this tutorial amount to employing alternative statistical estimators for the correlation matrix under a simple matrix-variate model or a factor model, respectively. This tutorial will be accompanied by preprocessed data from ABIDE and open source code for participants to gain an understanding of the benefits and limitations of the two techniques.

