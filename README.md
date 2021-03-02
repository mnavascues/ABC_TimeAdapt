## TimeAdapt

TimeAdapt makes a joint inference of demography and selection from longitudinal population genomic data. The inference is based on approximate Bayesian computation via random forest. It takes whole genome data and information about the ages of samples to perform simulations (using SLiM + pyslim/msprime).

### Requirements

The code has been tested with the following versions:

- Python 3.8.5
  - scikit-allel 1.3.2
  - argparse
  - msprime 0.7.4
  - numpy
  - pyslim 0.403
  - scipy 1.5.4
- R 3.6.1
  - argparser 0.6
  - extraDistr 1.8.11
  - rcarbon 1.3.1
- SLiM 3.4

### Configuration using conda

Creation of the environment from scratch:
```shell
$ conda create -n timeadaptenv python==3.8.5 r-base=3.6.1
$ conda activate timeadaptenv
$ pip install -r requirements.txt 
$ conda install -c r r-rcarbon=1.2.0
$ conda install -c r r-argparser=0.4
$ conda install -c r r-extraDistr=1.8.11
$ conda env export > timeadaptenv.yml
```

Creation of the environment via yml file:
```shell
$ conda env create -f timeadaptenv.yml
```

File requirements.R is also provided for an alternative way to install R packages.

### Run using snakemake

Run one batch of simulations with Snakefile (specify parameters within Snakefile, see below)
```shell
$ conda activate timeadaptenv
$ snakemake sim
```

Alternatively you can generate your own pipeline. For simulations, params.R (use "params.R -h" for help) generates files slim_\*.sh and pyslim_\*.sh with the SLiM and Python command lines that produce each simulation (SLiM must be run first).

## Input parameters (as used in the Snakefile)

- *project_name*: name of the analysis project, a folder with that name is created in results folder and all output written inside.
- *batch*: identifier of a batch of simulations for the project. A folder named with that identifier is created within the project folder and results from the simulations of that bacth written inside.
- *seed*: seed for the random number generator (integer)
- *sample_file*: file (with path) containing sample information (see below)
- *genome_file*: file (with path) containing genome information (see below)
- *generations_forward*: number of generations simulated in forward (by SLiM)
- *periods_forward*: number of periods in which the forward simulation is divided
- *gen_len_prior_sh1*: prior for generation length (shape1 of rescaled beta distribution)
- *gen_len_prior_sh2*: prior for generation length (shape2 of rescaled beta distribution)
- *gen_len_prior_min*: prior for generation length (minimum of rescaled beta distribution)
- *gen_len_prior_max*: prior for generation length (maximum of rescaled beta distribution)
- *num_of_sims*: number of simulation (in the batch)
- *pop_size_prior_min*: prior for population size (minimum of uniform distribution)
- *pop_size_prior_max*: prior for population size (maximum of uniform distribution)
- *mut_rate_prior_mean*: prior for mutation rate (mean of normal distribution)
- *mut_rate_prior_sd*: prior for mutation rate (sd of normal distribution)


### Input files

- *Sample information file*: Text file with information on the samples to be analysed. First line should be a header containing the following column names:
   - "sampleID": String identifying each individual. Do not use spaces or tabs.
   - "age14C": Point estimate of age from radiocarbon dating (for ancient samples, set NA for modern samples or samples with a known calendar year). It indicates time of death.
   - "age14Cerror": Error for age from radiocarbon dating.
   - "year": Calendar year for the time of sampling (for present and modern samples, such as museum/collection samples, set NA for ancient samples)
   - "coverage": Numeric value indicating sequencing depth (do NOT add "x", write "30" not "30x")
   - "damageRepair": Indicate whather the library for ancient samples used a damage repair step. Indicate TRUE or FALSE (indicate TRUE for modern samples)
   - "groups": Numerical code to indicate the groups of individuals for which summary statistics will be calculated. For each digit, individuals with the same number form one group, each digit defines a set of groups. Summary statistics comparing groups (such as Fst). Numbers indicating groups must be consecutive and start with 0.

- *Genome information file*:
