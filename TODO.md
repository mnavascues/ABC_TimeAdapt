## TODO

### Model

- [ ] *Selection.* **high**
  - [ ] *Selection on standing variation.* **high**
  - [ ] *Selection on background selection.* **high**
  - [ ] *Selection on new mutation.* **low**


### Statistics

- [ ] *Store RoH distribution as individual values and not as a list within the sumstats dictionary* **high** So that when transforms to pandas dataframe it does not create multiple rows for a single simulation
- [ ] *Script to perfomr ABCRF* **high**
- [ ] *Add runs statistics for two populations.* **low** Calculate runs statistics by making random pairs of individuals from each population.

### Software

- [ ] *Script to pool stats and parameters into a single reference table* **high**
- [ ] *Refactor getparams from R to Python.* **low** The only R package that will be kept is rcarbon for simulating age for subfossil samples: use rpy2.
- [ ] *Revise and move all test from timeadapt.py to test_timeadapt.py* **low**
- [x] ~~*Separate function to read project options from simulation options*~~

### Testing

- [ ] *Change how pytest compare numpy arrays* **low** FROM assert (numpy_array == expected_result).all() TO assert list(numpy_array) == expected_result

### Documentation

- [ ] *References for software and packages used.* **high**
  - [ ] *SLiM.* **high**
  - [ ] *tskit/msprime/pyslim.* **high**
  - [ ] *scikit.allel.* **high**
  - [ ] *rcarbon.* **low**
  