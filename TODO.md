## TO DO

### Model

- [ ] *Selection.* **high**
  - [ ] *Selection on standing variation.* **high**
  - [ ] *Selection on background selection.* **high**
  - [ ] *Selection on new mutation.* **low**


### Statistics

- [ ] *Script to perfomr ABCRF* **high**
- [ ] *Allow a wider choice of prior distributions*
- [ ] *Add 3 groups and 4 groups summary statistrics* **low** F3, F4, abba/baba
- [ ] *Add runs statistics for two populations.* **low** Calculate runs statistics by making random pairs of individuals from each population.
- [x] ~~*Store RoH distribution as individual values and not as a list within the sumstats dictionary* **high** So that when transforms to pandas dataframe it does not create multiple rows for a single simulation~~

### Software

- [ ] *Script to pool stats and parameters into a single reference table* **high**
  - [x] ~~*Pool summary statistics in a dataframe*~~
  - [ ] *Pool parameters and latent variables*
  - [ ] *export reference table as a file*
- [ ] *write genome map for slim in a single file for the project*
- [ ] *Use json files instead of ini* **low** SLiM can import parameters from json
- [x] ~~*Refactor getparams from R to Python.* **low** The only R package that will be kept is rcarbon for simulating age for subfossil samples: use rpy2.~~
- [x] ~~*Separate function to read project options from simulation options*~~

### Testing

- [ ] *test that number of periods in forward simulation is compatible with length of forward simulation*
- [ ] *Set a dataset test from bee example for integration testing* **high**
- [ ] *Revise and move all test from timeadapt.py to test_timeadapt.py* **low**
- [ ] *Change how pytest compare numpy arrays* **low** FROM assert (numpy_array == expected_result).all() TO assert list(numpy_array) == expected_result
- [ ] *make pytest & rpy2 work on github workflow* currently send an error on github but tests pass in local

### Documentation

- [ ] *References for software and packages used.* **high**
  - [ ] *SLiM.* 
  - [ ] *tskit/msprime/pyslim.* 
  - [ ] *scikit.allel.* 
  - [ ] *rcarbon.* 
  