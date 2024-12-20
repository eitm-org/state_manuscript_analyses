# STATE_analyses
Tooling and notebooks for analyzing STATE data


Scripts to make manuscript plots are organized like this:

![STATE figs](https://github.com/user-attachments/assets/fa17e95d-353f-4f96-a7ea-ce51a9c78546)


Script                                     | Input                                       | Output
------------------------------------------ | ------------------------------------------- | --------------
scripts/state_manuscript_data_analysis.ipy |                                             |
scripts/all_paper_figs_data_cleaning.R     | output from state_manuscript_data_analysis, |csvs for plotting
scripts/all_paper_figs.Rmd                 | output from all_paper_figs_data_cleaning    | all figures for manuscript



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

## Environment setup
1. run `git clone https://github.com/eitm-org/state_manuscript_analyses.git && cd state_manuscript_analyses`/
2. run `make venv && make requirements`.
3. activate venv by running `source venv/bin/activate`.
4. run `make install`.


## Data Setup
1. Download and unzip the references data folder, HG002_fully_resolved folder, and samples.csv from zenodo [ineset link here]
2. Download funcotator_dataSources.v1.7.20200521s from [here](https://console.cloud.google.com/storage/browser/broad-public-datasets/funcotator/funcotator_dataSources.v1.7.20200521g?pageState=(%22StorageObjectListTable%22:(%22f%22:%22%255B%255D%22))) and place it inside of the refs folder previously downloaded from zenodo. 
3. Modify paths in `post_processing/scripts/constants.py` to point to the respecting data paths.
4. Download and unzip tensig.zip from [here][https://drive.google.com/file/d/1jZVpvFOP8lOLKY1pTt1m8AUK9kyD-3u3/view] and place *only* the constants.Rdata file in to `post_processing/scripts/tensig/` folder.
5. rerun `make install`.


## To replicate the post processing pipeline results:
1. With venv activated, make sure the current working directory is in `state_manuscript_analyses`.
2. Download and unzip STATE_HG002_vcfs.zip and STATE_vcfs_f3_region_filtered.zip and place them in your current working directory.
3. Start the `post_processing_wf/post_processing_wf.sh` from line 15 (`python post_processing_wf/scripts/funcotate.py`) and run until the end. The pipeline cannot be replicated from the begining because we cannot release the unfiltered vcfs due to them possibly containing germline variants. However the dataset did release can be processed for functional annotation and additional characterization.
