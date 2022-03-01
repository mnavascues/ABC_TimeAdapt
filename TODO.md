## TO DO

### Model

- [ ] *SNP calling simulation* **high** Test new version of function. Remove variants only polymorphic in aDNA.
- [ ] *Verify calculation of Ne (latent variable) in SLiM* Is it using the right number of generations?
- [ ] *Selection.* **high**
  - [ ] *Selection on standing variation.* **high**
  - [ ] *Selection on background selection.* **high**
  - [ ] *Selection on new mutation.* **low**
- [ ] *Admixture* add a (ghost) population as a source of gene flow to focus population
- [ ] *Age of samples for generation time lower than a year* Modify scripts/input so that organisms with a generation time lower that a year can have multiple sampling times in the same year. As it is done now the minimum possible generation time is 1 year!
- [x] ~~*verify that number of periods in forward simulation is compatible with length of forward simulation*~~
- [x] ~~*Demography* **high** implement changes of population size in coalescent simulation~~


### Statistics

- [ ] *Script to perfomr ABCRF* **high**
- [ ] *Allow a wider choice of prior distributions*
- [ ] *Add 3 groups and 4 groups summary statistrics* **low** F3, F4, abba/baba
- [ ] *Add runs statistics for two populations.* **low** Calculate runs statistics by making random pairs of individuals from each population.
- [x] ~~*Store RoH distribution as individual values and not as a list within the sumstats dictionary* **high** So that when transforms to pandas dataframe it does not create multiple rows for a single simulation~~

### Software

- [ ] *Script to pool stats and parameters into a single reference table* **high**
  - [ ] *Pool parameters and latent variables*
  - [ ] *export reference table as a file*
  - [x] ~~*Pool summary statistics in a dataframe*~~
- [x] ~~*write genome map for slim in a single file for the project* to be used for all simulations~~
  - [ ] *verify the existence of the recombination map created by a previous batch* **low**
- [ ] **Add ascii plots** for high levels of verbosity. Package 'txtplot' 
- [x] ~~*Refactor getparams from R to Python.* **low** The only R package that will be kept is rcarbon for simulating age for subfossil samples: use rpy2.~~
- [x] ~~*Separate function to read project options from simulation options*~~
- [x] ~~*output latent variables for each simulation in a separate file* instead of writing all simulation on the same file. This should be more robust and make Snakemake workflow more linear~~

### Testing

- [ ] *Revise and move all test from timeadapt.py to test_timeadapt.py* **low**
- [ ] *coverage* learn how to calculate coverage and add it to test workflow
- [x] ~~*make pytest & rpy2 work on github actions* currently given an error on github but tests pass in local???: this was done by using a conda environment~~
- [x] ~~*Set a dataset test from bee example for integration testing* **high**~~
- [x] ~~*Change how pytest compare numpy arrays* **low** FROM assert (numpy_array == expected_result).all() TO assert list(numpy_array) == expected_result~~

### Documentation

- [ ] *References for software and packages used.* **high**
  - [ ] *SLiM.* 
  - [ ] *tskit/msprime/pyslim.* 
  - [ ] *scikit.allel.* 
  - [ ] *rcarbon.* 
  