# GCTA Genome-wide Complex Trait Analysis

This is a rescued version of GCTA 0.93.9 as published in

By Yang J, Lee SH, Goddard ME, Visscher PM. Genome-wide complex trait analysis (GCTA): methods, data analyses, and interpretations. Methods Mol Biol. 2013;1019:215-36. doi: [10.1007/978-1-62703-447-0_9](https://doi.org/10.1007/978-1-62703-447-0_9).

See also [wikipedia](https://en.wikipedia.org/wiki/Genome-wide_complex_trait_analysis)

The "GCTA" software package covers the GREML estimation of SNP
heritability, but also includes other functionality:

* Estimate the genetic relationship from genome-wide SNPs;
* Estimate the inbreeding coefficient from genome-wide SNPs;
* Estimate the variance explained by all the autosomal SNPs;
* Partition the genetic variance onto individual chromosomes;
* Estimate the genetic variance associated with the X-chromosome;
* Test the effect of dosage compensation on genetic variance on the X-chromosome;
* Predict the genome-wide additive genetic effects for individual subjects and for individual SNPs;
* Estimate the LD structure encompassing a list of target SNPs;
* Simulate GWAS data based upon the observed genotype data;
* Convert Illumina raw genotype data into PLINK format;
* Conditional & joint analysis of GWAS summary statistics without individual level genotype data
* Estimating the genetic correlation between two traits (diseases) using SNP data
* Mixed linear model association analysis

AUTHOR: Jian Yang, Hong Lee, Mike Goddard and Peter Visscher

CONTACT: jian.yang@uq.edu.au

NOTE: A newer version of gcta is being developed at
http://cnsgenomics.com/software/gcta/#Download. It uses the non-free
Intel MKL.

# INSTALL

Dependencies:

1. C++ compiler and build tools
2. Eigen libraries and include files on path

You will need to specify the path of library EIGEN in the Makefile
(keyword EigenLib). Example

    make -f Makefile

# EXAMPLE

## Simulate phenotypes from existing genotype

The following command generates 10 phenotypes from existing plink genotypes

    ./gcta64 --bfile example/HLC --simu-qt --simu-hsq 0.5 --simu-causal-loci example/HLC.causal --simu-rep 10 --out HLC

where HLC.causal contains

    rs6688000 0.122
    rs10746502 0.50
    rs1542899 0.72
    rs12404358 -0.35

# LICENSE

Released under GNU General Public License, v2 (see LICENSE)

# DOCUMENTATION

[PDF](doc/GCTA_UserManual_v1.24.pdf) and (for the new software)
http://cnsgenomics.com/software/gcta/#Overview
