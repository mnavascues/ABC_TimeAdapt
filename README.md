## TimeAdapt

TimeAdapt makes (*will eventually make*) a joint inference of demography and selection from longitudinal population genomic data. The inference is based on approximate Bayesian computation via random forest. It takes whole genome data and information about the ages of samples to perform simulations (using SLiM + pyslim/msprime).

### Requirements

TimeAdapt is a collection of scripts in Python and SLiM (3.6). They have been tested in an Ubuntu (20.04) machine using a Conda environment and using a Snakemake workflow to run them. The Conda environment was created with the following commands (on 12.Nov.2021):

```shell
conda create -n timeadapt python==3.8.10 r-base=3.6.3
conda activate timeadapt
pip install msprime==1.0.2
conda install pyslim=0.600
conda install scikit-allel=1.3.5
conda install dill=0.3.4
conda install rpy2=3.4.5
conda install -c r r-rcarbon=1.2.0
conda install pytest
conda install flake8
conda env export > timeadapt.yml
```

Create the environment via yml file:
```shell
conda env create -f timeadapt.yml
```

### Usage

Input files:

- *config file*: a ini file with options regarding the Setting, Model, Prior and Statistics to be used in your analysis

- *sample file*: a text file (in the form of a space separated table) with metadata of your samples (ID, age, sequencing coverage, etc)

- *genome file*: a text file with a description of the genome organization (chromosomes, recombination map)

- *data file*: a vcf file with the genetic data

To run different parts of the analysis with snakemake :

```shell
snakemake rule -C options_file='path/to/your/config_file.ini'
```
Where `rule` is one of the rules defines in the snakefile. For instance, running `snakemake reftable -C options_file='tests/config_project.ini'` will create small reference table using parameters in file `tests/config_project.ini`. Typically the user will use rule `reftable` to run simulations and create the reference table (parameters, summary statistics and latent variables), `afforestation` to grow random forests and ...

Before running your analysis is highly recommended to performs some tests. Unit tests (using pytest) can be run using `snakemake test`.

![Directed acyclic graph for the workflow using the test project configuration (`snakemake --dag | dot -Tsvg > workflow_dag.svg`)](workflow_dag.svg)

#### Cleaning old files

You can remove all files (from a specific project or from all projects) from your results folder:

```shell
snakemake clean_project -C options_file='path/to/your/config_file.ini'
snakemake clean
```


### Input parameters (in *.ini config file)

| Parameter name | type | description |
|---|---|---------------|
|**[Settings]**|||
| *verbose* | integer | Level of information to output on screen (0 = minimal; 1,2,3,... = increasing detailed output; -1 = none).|
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

Example:

```
[Settings]
seed = 1234567890
quiet = False
project = test
batch = 1
sample_file = tests/input/test_sample.txt
genome_file = tests/input/test_genome.txt
num_of_sims = 3

[Model]
generations_forward = 400
periods_forward = 8
periods_coalescence = 1

[Priors]
gen_len_prior_sh1 = 2
gen_len_prior_sh2 = 1.465967
gen_len_prior_min = 26
gen_len_prior_max = 30
pop_size_prior_min = 10
pop_size_prior_max = 200
mut_rate_prior_mean = 0.00000005
mut_rate_prior_sd = 0.5
```




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

TimeAdapt implements the method described by Pavinato *et al.* (2020) with some modifications. The code puts together several tools that should be acknowledged when using TimeAdapt: SLiM, pyslim, msprime (TODO: add full citation for these)


### Funding

This project has received funding from the European Union’s Horizon 2020 research and innovation programme under the Marie Skłodowska-Curie grant agreement No 791695 (TimeAdapt).

### Licence

TimeAdapt: joint inference of demography and selection

Copyright (C) 2021  Miguel de Navascués, Uppsala universitet, INRAE


This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.


This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.


You should have received a copy of the GNU General Public License along with this program. If not, see https://www.gnu.org/licenses/.





