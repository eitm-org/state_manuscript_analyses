# STATE_analyses
Tooling and notebooks for analyzing STATE data


Scripts to make manuscript plots are organized like this:

![STATE figs](https://github.com/user-attachments/assets/fa17e95d-353f-4f96-a7ea-ce51a9c78546)


Script                                     | Input                                       | Output
------------------------------------------ | ------------------------------------------- | --------------
scripts/state_manuscript_data_analysis.ipy |                                             |
scripts/all_paper_figs_data_cleaning.R     | output from state_manuscript_data_analysis, |csvs for plotting
scripts/all_paper_figs.Rmd                 | output from all_paper_figs_data_cleaning    | all figures for manuscript


# ClairS workflows

Nextflow modules and workflows that implement and adapt [ClairS](https://github.com/HKU-BAL/ClairS) (Zheng et. al., 2023) for extracting sSNVs for ONT long-read sequencing. 

# Post-processing pipeline
## Installation
There are three primary tools used as part of the post-processing pipeline - bcftoolss v1.17, GATK v4.4.0.0, VEP v110.1

### bcftools

bcftools can be downloaded from the [here](https://samtools.github.io/bcftools/)

There are a few dependencies for bcftools to function, they are listed below,

* [zlib](http://zlib.net)
* [libbz2](http://bzip.org/)
* [liblzma](http://tukaani.org/xz/)
* [libcurl](https://curl.haxx.se/)
* [libcrypto](https://www.openssl.org/)
* [gsl](https://www.gnu.org/software/gsl/)
* [libperl](http://www.perl.org/)

The following dependencies are needed to configure and build bcftools,

* GNU make
* C compiler (gcc or clang)
* autoheader
* autoconf

HTSlib is also required for bcftools but it is attached along with the bcftools package during download and can be used instead of a standalone installation.

Detailed installation  instructions can be found [here](https://raw.githubusercontent.com/samtools/bcftools/develop/INSTALL)


### Genome Analysis ToolKit 
The latest release of GATK v4 can be downloaded from [here](https://github.com/broadinstitute/gatk/releases)

Primary requirement is Java 8 / JDK 1.8 (Oracle and OpenJDK supported)

Some of the tools require R libraries and Python v3.6.2 or higher.

R packages: 
* gsalib
* ggplot2
* reshape
* gplots

There is a docker version of GATK with all required dependencies [here](https://hub.docker.com/r/broadinstitute/gatk/)

GATK does not require installation and the precompiled jar files can be invoked using the gatk wrapper script and the folder added to PATH. Detailed instructions can be found [here](https://gatk.broadinstitute.org/hc/en-us/articles/360036194592-Getting-started-with-GATK4)


### ENSEMBL Variant Effect Predictor 
VEP can be downloaded from [here](https://useast.ensembl.org/info/docs/tools/vep/script/vep_download.html)

The recommended installation method is via Docker or Singularity.

Primary requirements of VEP are as follows,
* gcc, g++ and make
* Perl version 5.10 or above recommended (tested on 5.10, 5.14, 5.18, 5.22, 5.26)
* Perl packages:
    * Archive::Zip
    * DBD::mysql (version <=4.050)
    * DBI

After installing the dependencies, run the installer PERL script. 
Detailed instructions can be found [here](https://useast.ensembl.org/info/docs/tools/vep/script/vep_download.html#installer)
