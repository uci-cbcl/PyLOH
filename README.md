README for PyLOH 1.0
======================


Introduction
============

Estimation of tumor purity is a crucial step in understanding
the complex copy number landscapes of heterogeneous tumor 
samples through next-generation sequencing data. A prominent 
problem in tumor purity estimation is how to disambiguate 
multiple compatible solutions, as distinct combinations of 
tumor purity and ploidy can fit the sequencing data equivalently
well. Current methods "break the tie" relying on either prior 
information or numerical heuristics rather than comprehensively 
explaining the sequencing data itself. Such approaches will 
lead to fundamental bias when the ground truth significantly
deviates from the prior information. Here we propose a full
probabilistic model-based method PyLOH that leverage the 
cluster pattern of loss of heterozygosity observed in paired
cancer sequencing data to disambiguate multiple compatible
solutions. We also introduce a novel visualization method 
"BAF heat map" to to characterize the cluster pattern of LOH.

Install
=======

Prerequisites
-------------

Python 2.7 are required to run PyLOH 1.0, [Python 2.7.3](http://www.python.org/download/releases/2.7.3/) is recommended. 



