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

There are also `config/` and `bin/` folders under PyLOH-*. The `config/` folder contains example priors and the `bin/` folder contains 
useful utilities, such as the R code to run [BICseq](http://compbio.med.harvard.edu/Supplements/PNAS11.html) and the python script to 
convert BICseq results to BED file. You can copy these two folders somewhere easily accessible.



USAGE
=====


Overview
--------
PyLOH is composed of three modules: 
* `preprocess`. Preprocess the reads aliments of paired normal-tumor samples in BAM format and produce the paired counts file, 
preprocessed segments file and preprocessed BAF heat map file as output.
 
* `run_model`. Take the paired counts file and preprocessed segments file as input, estimate tumor purity, the copy number and the
allele type of each segment.

* `postprocess`. Take the preprocessed BAF heat map file as input and plot the BAF heat map for each segment as output.

The general workflow of PyLOH is this
![alt tag](https://github.com/uci-cbcl/PyLOH/blob/gh-pages/images/workflow.png?raw=true)


Preprocess
----------
This part of README is based on [JoinSNVMix](https://code.google.com/p/joint-snv-mix/wiki/runningOld). To preprocess the paired 
cancer sequencing data, execute:
```
$ PyLOH.py preprocess REFERENCE_GENOME.fasta NORMAL.bam TUMOUR.bam BASENAME --segments_bed SEGMENTS.bed --min_depth 20 --min_base_qual 10 --min_map_qual 10 --process_num 10
```

**REFERENCE_GENOME.fasta** The path to the fasta file that the paired BAM files aligned to. Note that the index file should be generated 
for the reference genome. This can be done by running samtools as follows:

`$ samtools faidx REFERENCE_GENOME.fasta`

**NORMAL.bam** The BAM file for the normal sample. The BAM index file should be generated for this file and named NORMAL.bam.bai. This can
be done by running

`$ samtools index NORMAL.bam`

**TUMOUR.bam** The bam file for the tumour sample. As for the normal this file needs to be indexed.

**BASENAME** The base name of preprocessed files to be created.

**--segments_bed** Use the genome segmentation stored in SEGMENTS.bed. If not provided, use 22 autosomes as the segmentaion. 
But using automatic segmentation algorithm to generate SEGMENTS.bed is highly recommended, such as [BICseq](http://compbio.med.harvard.edu/Supplements/PNAS11.html).

**--min_depth** Minimum depth in both normal and tumor sample required to use a site in the analysis.

**--min_base_qual** Minimum base quality required for each base.

**--min_map_qual** Minimum mapping quality required for each base.

**--process_num** Number of processes to launch for preprocessing.


Run model
---------
After the paired cancer sequencing data is preprocessed, we can run the probabilistic model of PyLOH by execute:
```
$ PyLOH.py run_model BASENAME --allele_number_max 2 --max_iters 100 --stop_value 1e-7
```
**BASENAME** The base name of preprocessed files created in the preprocess step.

**--allele_number_max** The maximum copy number of each allele allows to take.

**--priors** Path to the file of the prior distribution. The prior file must be consistent with the --allele_number_max. If not provided,
use uniform prior, which is recommended.

**--max_iters** Maximum number of iterations for training.

**--stop_value** Stop value of the EM algorithm for training. If the change of log-likelihood is lower than this value, stop training.


Postprocess
-----------
Currently, the postprocess module is only for plotting the BAF heat map of each segment:
```
$ PyLOH.py BAF_heatmap BASENAME
```

**BASENAME** The base name of preprocessed files created in the preprocess step.


Output files
------------
**\*.PyLOH.counts** The preprocessed paired counts file. It contains the allelic counts information of sites, which are heterozygous 
loci in the normal genome. The definition of each column in a *.PyLOH.counts file is listed here:

| Column    | Definition                                         | 
| :-------- | :------------------------------------------------- | 
| seg_index | Index of each segment                              |      
| normal_A  | Count of bases match A allele in the normal sample |
| normal_B  | Count of bases match B allele in the normal sample |
| tumor_A   | Count of bases match A allele in the tumor sample  |
| tumor_B   | Count of bases match B allele in the tumor sample  |

**\*.PyLOH.segments** The preprocessed segments file. It contains the genomic information of each segment. The definition of each
column in a *.PyLOH.segments file is listed here:

| Column           | Definition                                                              | 
| :--------------- | :---------------------------------------------------------------------- | 
| seg_name         | Name of the segment                                                     |      
| chrom            | Chromosome of the segment                                               |  
| start            | Start position of the segment                                           |
| end              | End position of the segment                                             |
| normal_reads_num | Count of reads mapped to the segment in the normal sample               |
| tumor_reads_num  | Count of reads mapped to the segment in the normal sample               |
| LOH_frec         | Fraction of LOH sites in the segment                                    |
| LOH_status       | FALSE -> no LOH; TRUE -> significant LOH; UNCERTAIN -> medium level LOH |
| log2_ratio       | Log2 ratio between tumor_reads_num and normal_reads_num                 |

**\*.PyLOH.segments.extended** The extended segments file after run_model. There are two additional columns:

| Column           | Definition                            | 
| :--------------- | :-------------------------------------| 
| allele_type      | Estimated allele type of the segment  |      
| copy_number      | Estimated copy number of the segment  |  

**\*.PyLOH.purity** Estimated tumor purity.

**\*.PyLOH.heatmap.pkl** The preprocessed BAF heat map file in Python pickle format.

**\*.PyLOH.heatmap.plot** The folder of BAF heat maps plotted for each segment. A typical BAF heat map looks like this
![alt tag](https://github.com/uci-cbcl/PyLOH/blob/gh-pages/images/BAF_heamap_sample.png?raw=true)



OTHER
=====

BIC-seq related utilities
-------------------------
We highly recommend using automatic segmentation algorithm to partition the tumor genome, and thus prepare the segments file in BED format.
For exmaple, we used [BICseq](http://compbio.med.harvard.edu/Supplements/PNAS11.html) in the original paper. To run a BICseq analysis, you
can copy the commands in `bin/BICseq.R` and paste them in a R interative shell. Or you can also run the R script from the command line:
```
$ R CMD BATCH bin/BICseq.R
```
Note that,`normal.bam` and `tumor.bam` must be in the same directory where you run the command. The R script will output a segments file
`segments.BICseq`. Then you can use the other script `bin/BICseq2bed.py` to convert the segments file into BED format:
```
$ BICseq2bed.py segments.BICseq segments.bed --seg_length 1000000
```

**--seg_length** Only convert segments with length longer than the threshold.
