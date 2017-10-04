# GCTA Genome-wide Complex Trait Analysis

By Yang J, Lee SH, Goddard ME, Visscher PM. Genome-wide complex trait analysis (GCTA): methods, data analyses, and interpretations. Methods Mol Biol. 2013;1019:215-36. doi: [10.1007/978-1-62703-447-0_9](https://doi.org/10.1007/978-1-62703-447-0_9).

See also [wikipedia](https://en.wikipedia.org/wiki/Genome-wide_complex_trait_analysis)

The "GCTA" software package covers the GREML estimation of SNP
heritability, but also includes other functionality:

*        Estimate the genetic relationship from genome-wide SNPs;
*        Estimate the inbreeding coefficient from genome-wide SNPs;
*        Estimate the variance explained by all the autosomal SNPs;
*        Partition the genetic variance onto individual chromosomes;
*        Estimate the genetic variance associated with the X-chromosome;
*        Test the effect of dosage compensation on genetic variance on the X-chromosome;
*        Predict the genome-wide additive genetic effects for individual subjects and for individual SNPs;
*        Estimate the LD structure encompassing a list of target SNPs;
*        Simulate GWAS data based upon the observed genotype data;
*        Convert Illumina raw genotype data into PLINK format;
*        Conditional & joint analysis of GWAS summary statistics without individual level genotype data
*        Estimating the genetic correlation between two traits (diseases) using SNP data
*        Mixed linear model association analysis

AUTHOR: Jian Yang, Hong Lee, Mike Goddard and Peter Visscher

CONTACT: jian.yang@uq.edu.au

YEAR: 2009-2011

LICENSE: Released under GNU General Public License, v2 (see LICENSE)

DOCUMENTATION: ?

LIBRARY: require the library EIGEN (version 3.0.3) support (see http://eigen.tuxfamily.org/). You will need to specify the path of library EIGEN in the Makefile (keyword EigenLib).

COMPILATION: You will need a standard C/C++ compiler such as GNU gcc and type

      make -f Makefile
