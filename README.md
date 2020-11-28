## TimeAdapt

TimeAdapt makes a joint inference of demography and selection from longitudinal population genomic data. The inference is based on approximate Bayesian computation via random forest. It takes whole genome data and information about the ages of samples to perform simulations (using SLiM + pyslim/msprime).

### Requirements

- Python 3.8.5
  - allel
  - argparse
  - msprime 0.7.4
  - numpy
  - pyslim 0.403
  - scipy
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

