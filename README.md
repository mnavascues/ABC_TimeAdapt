## TimeAdapt

TimeAdapt makes (i.e. will eventually make) a joint inference of demography and selection from longitudinal population genomic data. The inference is based on approximate Bayesian computation via random forest. It takes whole genome data and information about the ages of samples to perform simulations (using SLiM + pyslim/msprime).

### Requirements

The code has been tested with the following versions (on Ubuntu 20.04):

- Python 3.8.5
  - scikit-allel 1.3.2
  - msprime 1.0.0b1
  - numpy
  - pyslim 0.403
  - scipy 1.5.4
  - pytest
- R 3.6.1
  - ini 0.3.1
  - extraDistr 1.8.11
  - rcarbon 1.3.1
  - testthat 2.1.1
- SLiM 3.4

### Usage

These instructions are written for myself and might need to be adapted to other configurations.

Creation of the environment from scratch:
```shell
$ conda create -n timeadapt python==3.8.5 r-base=3.6.1
$ conda activate timeadapt
$ conda install scikit-allel=1.3.2
$ conda install msprime=1.0.0b1
$ conda install numpy=1.20.2
$ conda install pyslim=0.600
$ conda install scipy=1.5.3
$ conda install pandas=1.2.3
$ conda install pytest=6.2.3
$ conda install -c r r-rcarbon=1.2.0
$ conda install -c r r-ini=0.3.1
$ conda install -c r r-extraDistr=1.8.11
$ conda install -c r r-testthat=2.1.1
$ conda install pip
$ python -m pip install msprime --pre --upgrade
$ conda env export > timeadapt.yml
```

Creation of the environment via yml file:
```shell
$ conda env create -f timeadapt.yml
```

Run (using snakemake) tests (testthat for R, pytest for Python)
```shell
$ conda activate timeadapt
$ snakemake test
```

Run one batch of simulations with Snakefile and get directed acyclic graph of pipeline
```shell
$ conda activate timeadapt
$ snakemake sim
$ snakemake --dag | dot -Tsvg > results/workflow_dag.svg
```

Alternatively you can create your own pipeline. For simulations, simulations.R generates files slim_\*.sh and pyslim_\*.sh with the SLiM and Python command lines that produce each simulation (SLiM must be run first, then Python).

Remove all results from project folder
```shell
$ snakemake clean
```


### Input parameters (in *.ini config file)

| Parameter name | type | description |
|---|---|---------------|
|**[Settings]**|||
| *project_name* | string | name of the analysis project, a folder with that name is created in results folder and all output written inside.|
| *batch* | integer | identifier of a batch of simulations for the project. A folder named with that identifier is created within the project folder and results from the simulations of that bacth written inside.|
| *sample_file* | string |  path + file name containing sample information (see below)|
| *genome_file* | string | path + file name containing genome information (see below)|
| *num_of_sims* | positive integer | number of simulation (in the batch)|
| *seed* | integer | seed for the random number generator |
|**[Model]**|||
| *generations_forward* | positive integer | number of generations simulated in forward (by SLiM)|
| *periods_forward* | positive integer | number of periods in which the forward simulation is divided|
|**[Priors]**|||
| *gen_len_prior_sh1* | | prior for generation length (shape1 of rescaled beta distribution)|
| *gen_len_prior_sh2*| | prior for generation length (shape2 of rescaled beta distribution)|
| *gen_len_prior_min*| | prior for generation length (minimum of rescaled beta distribution)|
| *gen_len_prior_max*| | prior for generation length (maximum of rescaled beta distribution)|
| *pop_size_prior_min* | positive integer | prior for population size (minimum of uniform distribution)|
| *pop_size_prior_max* | positive integer | prior for population size (maximum of uniform distribution)|
| *mut_rate_prior_mean* |  |prior for mutation rate (mean of normal distribution)|
| *mut_rate_prior_sd* |  | prior for mutation rate (sd of normal distribution)|


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

...


### How to cite TimeAdapt

TimeAdapt implements the method described by Pavinato et al. (2020). The code puts together several tools that should be acknowledged when using TimeAdapt: SLiM, pyslim, msprime (TODO: add full citation for these)
