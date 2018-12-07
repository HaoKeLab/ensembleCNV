## iPattern

### Installation

To request the iPattern package, please contact the corresponding author Dr. Stephen W. Scherer (stephen.scherer@sickkids.ca) of the [paper](https://www.ncbi.nlm.nih.gov/pubmed/20531469). For more information about iPattern, please refer to the [paper](https://www.ncbi.nlm.nih.gov/pubmed/?term=21552272).

After obtaining the package (e.g., ipn.0.581.tar.gz is the version we received), please follow the instructions in the iPattern tutorial enclosed in the package for installation and usage. Here we echo the installation instructions in their tutorial. 

#### Requirements
- R (2.7.1+)
- The R ppc package – it can be downloaded from http://www-stat.stanford.edu/~tibs/PPC/Rdist/index.html
- Python (2.5.5+)

#### Setup

- untar the file with “tar -zvxf ipn.0.581.tar.gz”
- setup environment:
  - set up environment variable IPNBASE, export IPNBASE='/your/path/to/ipn'
  - set up environment variable PYTHONPATH=$PYTHONPATH:'/your/path/to/ipn/ipnlib'

Note: the directory structure/name must be kept as it is. Changing the directory structure will break the iPattern pipeline, the pipeline finds all the necessary scripts based on IPNBASE and the directory structure. When PBS job submitting system is not available, you can use –noqsub option to run iPattern sequentially.

### Analysis workflow

The script `run_iPattern.R` contains a template to run iPattern in a cluster environment. The names of directories involved in analysis need to  specified by users. 

