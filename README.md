# SRA RNASeq Pipeline

A standardized pipeline for downloading and processing RNA-Seq data. Basic run requires only an SRA study accession or a Synapse data folder ID. 

## Installation 

This is a pre-release of the pipeline designed for use on the UT Health San Antonio bioinformatics cluster (GNU/Linux x86_64 Ubuntu bionic 18.04) and has not been tested on any other systems or configurations. 

To install: 

1. Download repo -- `git clone https://github.com/millerh1/SRA_RNASeq_Pipeline.git`
2. Copy all contents to place on `PATH` while preserving directory structure:\

  Easy example:  
  `cp -r SRA_RNASeq_Pipeline/{SRA_RNASeq_Pipeline,SRA_RNASeq_Pipeline_Genome_Build,SRA_RNASeq_Pipeline_Files} ~/bin`  
  
  Alternative example:   
  `find SRA_RNASeq_Pipeline/ -maxdepth 1 -mindepth 1 -exec cp -r {} ~/bin \;`

### Dependencies 

To operate all capabilities of the pipeline, the following packages must be installed and visible to `PATH` (version recommended is the tested version). The easiest method to install them is through Anaconda3/Miniconda3:

#### General:

* `Anaconda3/miniconda3 >= 4.7.6`
* `csvkit >= 1.0.4`
* `synapseclient >= 1.9.3`
* `mashmap >= 2.0`
* `salmon >= 0.14.1`
* `STAR >= 2.7.0`
* `STAR-Fusion >= 1.6.0`
* `igv-reports >= 0.9.2`
* `tetoolkit >= 2.0.3 (requires python=2.7)`

#### R-Language:

* `R >= 3.6.0`
* `jsonlite >= 1.6`
* `httr >= 1.4`

#### Perl-Language (required only for STAR-Fusion):

* `perl >= 5.26.2`
* `JSON::XS >= 4.0.2`
* `CARP::Assert >= 0.21`
* `Types::Serialiser >= 1.0`
* `URI::escape >= 1.76`
* `common::sense >= 3.74`
* `Set::IntervalTree >= 0.11`
* `DB_file >= 1.852`

**NOTE**: `tetoolkit` requires a separate conda environment with `python v2.7` installed. Also, `STAR-Fusion` requires `STAR v2.7.0`. If `STAR v2.7.1` is installed (the default `conda install` behavior), then a separate environment should be created with both `STAR-Fusion` and `STAR v2.7.0` -- this can be passed to the pipeline with the `--SF_env=environmentName` argument. 

## Using SRA RNASeq Pipeline

### Building a genome directory

If the installation path is properly located on `PATH` -- then the only required line is:  

`SRA_RNASeq_Pipeline_Genome_Build`

This will automatically build the genome directory inside of the `SRA_RNASeq_Pipeline_Files` folder. This is the **preferred** genome building method because it automatically generates all the files and indices required by the pipeline in the location that is expected by running the pipeline with default settings. 

Full set of options:  

```
SRA_RNASeq_Pipeline_Genome_Build [ -s <both|human|mouse>] [-P num_threads] [ -g genome_directory]

  -s|--species                  str               Select 'human', 'mouse', or 'both' [default = 'both'].
  -g|--genome_directory         genome_dir        Choose genome build location [default = package_dir/genome_build]
  -P|--num_threads              N                 Specify number of threads to utilize for commands [default = 40]
  -h|--help                                       Display usage info
```

#### Output

The output is a folder containing the genome resources (fasta files, aligner indices, etc) required to run the full pipeline. `genome_build` folder structure:   

```
SRA_RNASeq_Pipeline_Files/genome_build
├── human
│   ├── Assembly_Files
│   │   └── decoyTranscripts
│   ├── ctat_genome_lib_build_dir
│   │   ├── __chkpts
│   │   └── ref_genome.fa.star.idx
│   ├── Fasta_Files
│   ├── Salmon_Transcripts_Index
│   └── STAR_Genome_Index
└── mouse
    ├── Assembly_Files
    │   └── decoyTranscripts
    ├── ctat_genome_lib_build_dir
    │   ├── __chkpts
    │   └── ref_genome.fa.star.idx
    ├── Fasta_Files
    ├── Salmon_Transcripts_Index
    └── STAR_Genome_Index
```

### Running the pipeline

To run the pipeline, users should first find an SRA or Synapse project they want to process data from, then input the accession in the pipeline. 

SRA Example:

`SRA_RNASeq_Pipeline -a SRP045832 -p adriamycin_project`

Synapse Example:

`SRA_RNASeq_Pipeline --synapseID=syn4228582 --username=userName -p iPSC_project`

Full usage info:

```
SRA_RNASeq_Pipeline [-a SRP_accession] [-p project_directory] [-t tests] [-g genome_directory] [-P num_threads]

  -a|--SRP_accession     SRP       SRA study accession to download data from. [e.g. SRP045672].
  -p|--project_dir       dir       Name of project directory. [default = 'RNASeq_SRA_Pipeline_Project']
  -t|--tests             string    's' [salmon] 't' [TECount] 'S' [Splicing-STAR] 'f' [STAR-Fusion] Default: 'stSf'
  -g|--genome_dir        dir       Name of genome directory. [default = package_dir/genome_build]
  -P|--num_threads       int       Specify number of threads. [default = 1]
  --python27_env         env       If 't' in tests, name of conda env with python 2.7 and TEtoolkit installed.
  --SF_env               env       Conda ENV with STAR-Fusion and STAR 2.7.0 [Due to conflict with 2.7.1].
  --fastp                          Use fastp to perform adapter trimming, filtering, and fastq QC.
  --force_spliceSE                 Not recommended: If 'S' in tests, single-end reads are run with splicing STAR.
  --force_fusionSE                 Not recommended: If 'f' in tests, single-end reads are run with STAR-fusion.
  --synapseID            synID     Parent of desired file directory.
  --synapseSpecies       string    If --synapseID is set, specify species. ['human' or 'mouse] Default: 'human'
  --synapseMergeFile     file      Mapping from downloaded files to sample names [required for synapse fastqs].
  --username             string    Synapse username for downloading data (if --synapse is specified).
  --help                           Display usage info
```

