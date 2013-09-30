README for PyLOH 1.0
======================


INTRODUCTION
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



INSTALL
=======



Prerequisites
-------------
* Although not mandatory, Linux system is recommended.

* Python (2.7). [Python 2.7.3](http://www.python.org/download/releases/2.7.3/) is recommended.

* [Numpy](http://www.numpy.org/)(>=1.6.1). You can download the source of Numpy from [here](http://sourceforge.net/projects/numpy/files/).

* [Scipy](http://www.scipy.org/)(>=0.10). You can download the source of Scipy from [here](http://sourceforge.net/projects/scipy/files/).

* [Pysam](https://code.google.com/p/pysam/)(>=0.7). To install Pysam, you also need to install [Cython](http://cython.org/) first. 

* [matplotlib](http://matplotlib.org/)(>=1.2.0) is required to plot BAF heat map.


Altough not required by PyLOH, [samtools](http://samtools.sourceforge.net/) can be useful for creating bam, bam index and fasta index 
files which are required by the pysam module of PyLOH. 

Install from source
-------------------
Download the compressed source file PyLOH-*.tar.gz and do as follows:

```
$ tar -xzvf PyLOH-*.tar.gz
$ cd PyLOH-*
$ python setup.py install
```

If you prefer to install PyLOH other than the default directory, you can also use this command:
```
$ python setup.py install --prefix /home/yili/
```

There are also `config/` and `bin/` folders under PyLOH-*. The `config/` folder contains example priors and the `bin/` folder contains useful 
utilities, such as the R code to run [BICseq](http://compbio.med.harvard.edu/Supplements/PNAS11.html) and the python script to convert
BICseq results to BED file. You can copy these two folders somewhere easily accessible.



