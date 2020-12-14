## TimeAdapt

TimeAdapt makes a joint inference of demography and selection from longitudinal population genomic data. The inference is based on approximate Bayesian computation via random forest. It takes whole genome data and information about the ages of samples to perform simulations (using SLiM + pyslim/msprime).

### Requirements

- Python 3.8.5
  - scikit-allel 1.3.2
  - argparse
  - msprime 0.7.4
  - numpy
  - pyslim 0.403
  - scipy 1.5.4
- R 3.6.3
  - abcrf 1.8.1
  - argparser 0.6
  - extraDistr 1.8.11
  - rcarbon 1.3.1
- SLiM 3.4

### Configuration

```shell
$ conda create -n timeadaptenv python==3.8.5
$ conda activate timeadaptenv
$ pip install -r requirements.txt 
$ Rscript requirements.R
```

yml file created via:

```shell
$ conda env export > timeadaptenv.yml
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
