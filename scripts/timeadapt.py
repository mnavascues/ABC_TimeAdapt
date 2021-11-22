#    TimeAdapt: joint inference of demography and selection
#    Copyright (C) 2021  Miguel de Navascués, Uppsala universitet, INRAE
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
import configparser  # for reading ini files
import pandas as pd  # for reading "table" files 
import numpy as np
import scipy.stats as st
import random
import allel
import math
import tempfile # for creating temporal files on testing
import pytest

# Using R package rcarbon through rpy2
from rpy2.robjects.packages import importr
import rpy2.robjects as robjects
rcarbon = importr("rcarbon")

### PRINT INFO ######################################################################################
def print_info(script_name,verbose,project=None,batch=None,sim=None):
  if verbose >=1 :
    print("#########################################")
  if verbose >=0 :
    if batch is not None:
      project = None
      if sim is None :
        print("TimeAdapt - %s - batch %s" %(script_name,batch) )
      else : 
        print("TimeAdapt - %s - batch %s - simulation %s" %(script_name,batch,sim) )
    elif project is not None:
      print("TimeAdapt - %s - project %s" %(script_name,project) )
    else :
      print("TimeAdapt - %s" %script_name )
  if verbose >=1 :
    print("by Miguel de Navascués")
    print("INRAE & Uppsala universitet")
    print("miguel.navascues@inrae.fr")
    print("#########################################")

### GET PROJECT OPTIONS ###############################################################################
def get_project_options(proj_options_file):
  proj_options = configparser.ConfigParser()
  proj_options.read(proj_options_file)
  # Settings
  assert 'Settings' in proj_options,"Missing [Settings] section in file"
  assert 'project' in proj_options['Settings'],"Missing 'project' parameter in file"
  project        = proj_options.get('Settings','project')
  num_of_batches = proj_options.getint('Settings','num_of_batches')
  try:
    verbose      = proj_options.getint('Settings','verbose')
  except:
    verbose      = int(proj_options.getfloat('Settings','verbose'))
  genome_file    = proj_options.get('Settings','genome_file')
  sample_file    = proj_options.get('Settings','sample_file')
  num_of_sims    = proj_options.getint('Settings','num_of_sims')
  seed           = proj_options.getint('Settings','seed')
  # Model
  assert 'Model' in proj_options,"Missing [Model] section in project options file"
  generations_forward = proj_options.getint('Model','generations_forward')
  periods_forward = proj_options.getint('Model','periods_forward')
  times_of_change_forw = get_times_of_change(generations_forward,periods_forward)
  periods_coalescence = proj_options.getint('Model','periods_coalescence')
  times_of_change_back = get_times_of_change(100000,periods_coalescence,mode="exponential")
  # Priors
  assert 'Priors' in proj_options,"Missing [Priors] section in project options file"
  gen_len_sh1 = proj_options.getfloat('Priors','gen_len_sh1')
  gen_len_sh2 = proj_options.getfloat('Priors','gen_len_sh2')
  gen_len_min = proj_options.getfloat('Priors','gen_len_min')
  gen_len_max = proj_options.getfloat('Priors','gen_len_max')
  assert gen_len_max>=gen_len_min,"Maximum must be higher than minimum"
  pop_size_min = proj_options.getfloat('Priors','pop_size_min')
  pop_size_max = proj_options.getfloat('Priors','pop_size_max')
  assert pop_size_max>=pop_size_min,"Maximum must be higher than minimum"
  mut_rate_mean = proj_options.getfloat('Priors','mut_rate_mean')
  mut_rate_sd = proj_options.getfloat('Priors','mut_rate_sd')
  # Statistics
  # TODO
  return {"project":project,
          "num_of_batches":num_of_batches,
          "genome_file":genome_file, 
          "sample_file":sample_file, 
          "verbose":verbose,
          "num_of_sims":num_of_sims,
          "seed":seed,
          "generations_forward":generations_forward,
          "periods_forward":periods_forward,
          "periods_coalescence":periods_coalescence,
          "times_of_change_forw":times_of_change_forw,
          "times_of_change_back":times_of_change_back,
          "gen_len_sh1":gen_len_sh1,
          "gen_len_sh2":gen_len_sh2,
          "gen_len_min":gen_len_min,
          "gen_len_max":gen_len_max,
          "pop_size_min":pop_size_min,
          "pop_size_max":pop_size_max,
          "mut_rate_mean":mut_rate_mean,
          "mut_rate_sd":mut_rate_sd}

### GET TIMES OF CHANGE ######################################################################################
def get_times_of_change(total_length,number_of_periods,mode="regular",exponent=2):
  if number_of_periods==1: return None
  if total_length<=number_of_periods:
    raise ValueError('number of periods must be lower than length of simulation')
  if mode=="exponential":
    t=total_length/((pow(exponent,number_of_periods-1)-1)/(exponent-1))
    if t<1:
      raise ValueError('incompatible value, try less periods or longer simulation')
    times_of_change = [int(t)]
    for i in range(1,number_of_periods-1):
      times_of_change.append( times_of_change[i-1] + int(t*pow(exponent,i)) )
  else:
    if mode!="regular":print("Warning: wrong value for mode in get_times_of_change(); using mode='regular' instead")
    x=total_length/number_of_periods
    times_of_change = []
    while x<total_length-1:
      times_of_change.append(int(x))
      x = x + total_length/number_of_periods
  # the next two errors should not happen if total_length>number_of_periods
  if len(times_of_change)!=number_of_periods-1:
    raise ValueError('wrong number of times')
  if len(set(times_of_change))<len(times_of_change):
    raise ValueError('repeated values in times of change')
  return times_of_change

### GET SIM OPTIONS ######################################################################################
def get_sim_options(sim_options_file):
  sim_options = configparser.ConfigParser()
  sim_options.read(sim_options_file)
  sim          = sim_options.get('Simulation','sim')
  batch        = sim_options.get('Simulation','batch')
  ss           = [int(i) for i in sim_options.get("Sample","ss").split()]
  chrono_order = [int(i) for i in sim_options.get("Sample","chrono_order").split()]
  assert sum(ss)==len(chrono_order), "Verify number of samples, inconsistent sample size across simulation options file"
  N            = [int(i) for i in sim_options.get("Demography","N").split()]
  for i in N:
    assert i>0, "Verify population sizes, they can be only positive values (excluding zero)"
  mu           = sim_options.getfloat('Genome','mu')
  assert mu>=0, "Verify mutation rate, it can be only positive values and zero"
  ttratio      = sim_options.getfloat('Genome','ttratio')
  seq_error    = sim_options.getfloat('Genome','seq_error')
  seed_coal    = sim_options.getint('Seeds','seed_coal')
  seed_mut     = sim_options.getint('Seeds','seed_mut')

  return {"sim":sim,
          "batch":batch,
          "ss":ss, 
          "chrono_order":chrono_order, 
          "N":N, 
          "mu":mu, 
          "ttratio":ttratio, 
          "seq_error":seq_error, 
          "seed_coal":seed_coal, 
          "seed_mut":seed_mut}

### GET OPTIONS ######################################################################################
def get_options(proj_options_file,sim_options_file):
  options = get_project_options(proj_options_file)
  options.update(get_sim_options(sim_options_file))
  return options

### READ SAMPLE INFO ######################################################################################
def read_sample_info(sample_info_file):
  # TODO : make it work with only ancient data (i.e. no year column) or no ancient data (no age14C)
  info = pd.read_table(filepath_or_buffer=sample_info_file, sep=r'\s+',
                       converters={'groups': lambda x: str(x)})
  # check header of file
  header = set(info.columns.values)
  expected_header = set(["sampleID","age14C","age14Cerror","year","coverage","damageRepair","groups"])
  assert header==expected_header, "Sample file does not have the expected column names"
  
  sample_size = len(info)
  group_levels = len(info["groups"][1])
  groups = np.full((group_levels, sample_size), 0, "int")
  is_ancient = []
  is_modern = []
  total_ancient = 0
  for i, row in info.iterrows():
      if len(row["groups"]) != group_levels:
          print("Error: verify that all individuals are assigned to groups")
          # TODO : interrupt program here
      if np.isnan(row["year"]):
          is_ancient.append(True)
          is_modern.append(False)
          total_ancient = total_ancient + 1
      else:
          is_ancient.append(False)
          is_modern.append(True)
      for level in range(0, group_levels):
          groups[level, i] = row["groups"][level]
  t0 = info["year"].max()

  return {"id":info["sampleID"],
          "age14C":info["age14C"],
          "age14Cerror":info["age14Cerror"],
          "ageBCAD":info["year"],
          "t0":info["year"].max(skipna=True),
          "coverage":info["coverage"],
          "is_ancient":is_ancient, 
          "is_modern":is_modern,
          "is_dr":info["damageRepair"], 
          "total_ancient":total_ancient,
          "size":sample_size,
          "group_levels":group_levels, 
          "groups":groups}

### READ GENOME INFO ######################################################################################
def get_genome_info(genome_info_file):
  table = pd.read_table(filepath_or_buffer=genome_info_file, sep=r'\s+')
  # check header of file
  header = set(table.columns.values)
  expected_header = set(["Chromosome","Position","Recombination_rate"])
  assert header==expected_header, "Genome file does not have the expected column names"

  nchr = max(table["Chromosome"])
  chromosomes = list(table["Chromosome"])
  rates = list(table["Recombination_rate"])
  positions = list(table["Position"])
  chr_ends_index = list(np.append(np.where(np.diff(chromosomes)==1)[0],len(chromosomes)-1))
  chr_ends = list(map(positions.__getitem__,chr_ends_index))
  rescaling_values = np.append(0,np.cumsum(chr_ends))
  chromo = 0
  # sum lenght of previous chromosomes to positions
  for i in range(0,len(positions)):
    positions[i] = positions[i] + rescaling_values[chromo]
    if (i in chr_ends_index):
      chromo += 1
  chr_ends = list(map(positions.__getitem__,chr_ends_index))
  L=chr_ends[-1]-1
  # insert recombination rate between chromosomes
  slim_rates = rates[:]
  for chromo in range(0,nchr-1):
    positions.insert(chr_ends_index[chromo]+1+chromo,positions[chr_ends_index[chromo]+chromo]+1)
    rates.insert(chr_ends_index[chromo]+1+chromo,math.log(2))
    slim_rates.insert(chr_ends_index[chromo]+1+chromo,0.5)
  slim_positions = [int(x-1) for x in positions]
  positions = [int(x) for x in positions]
  chr_ends = [int(x) for x in chr_ends] 
  positions.insert(0,0) # insert first position

  return {"nchr":nchr,
          "chr_ends":chr_ends,
          "L":int(L),
          "msprime_r_map":{"rates":rates, "positions":positions},
          "slim_r_map":{"rates":slim_rates, "positions":slim_positions}}

### GET AGE PDF ######################################################################################
# wrapper of rcarbon.calibrate, /!\ outputs ages in BCAD not BP
def get_age_pdf(x,errors,calCurves):
  try:
    sample_size = len(x)
    x = robjects.IntVector(x)
    errors = robjects.IntVector(errors)
  except:
    sample_size = 1
  res = rcarbon.calibrate(x=x, errors=errors, calCurves=calCurves, verbose=False)
  age_pdf = []
  for i in range(0,sample_size):
    ageBCAD = list(rcarbon.BPtoBCAD(res[1][i][0]))
    age_pdf.append({"ageBCAD":ageBCAD,"PrDens":list(res[1][i][1])})
  return age_pdf

### GET SAMPLE AGES ######################################################################################
def get_sample_ages(sample,age_pdf,gen_len):
  sample_ages = [None]*sample["size"]
  ancient_counter=0
  for i in range(0,sample["size"]):
    if sample["is_modern"][i]: sample_ages[i] = sample["ageBCAD"][i]
    if sample["is_ancient"][i]:
      sample_ages[i] = random.choices(age_pdf[ancient_counter]["ageBCAD"],
                                      weights = age_pdf[ancient_counter]["PrDens"])
      
      ancient_counter += 1
    sample_ages[i] = int(abs(sample_ages[i]-sample["t0"])/gen_len)
  return sample_ages

### SAMPLE PARAMETER TRAJECTORY #######################################################################
def sample_param_trajectory(times,minimum,maximum,factor=10):
  trajectory = np.zeros(shape=times)
  trajectory[0] = st.loguniform.rvs(a=minimum, b=maximum, size=1)
  for i in range(1,times):
    trajectory[i] = trajectory[i-1] * st.loguniform.rvs(a=1/factor, b=1*factor, size=1)
    trajectory[i] = max( min(trajectory[i],maximum) , minimum )
  return trajectory  





























### MAKE EMPTY GENOTYPE ARRAY ····························································
def empty_genotype_array(n_loci, n_samples, ploidy=2, allele=-1):
  """
  Creates a genotype array with all values as missing (-1) for a given number
  of samples, loci and ploidy

  :return: empty_ga
  """
  empty_ga = allel.GenotypeArray(np.full((n_loci, n_samples, ploidy), allele, dtype='i1'), dtype='i1')
  return empty_ga
def test_empty_genotype_array():
  test_ga = empty_genotype_array(3, 4, ploidy=2, allele=-1)
  assert type(test_ga) is allel.model.ndarray.GenotypeArray
  assert len(test_ga) == 3
  assert all(test_ga[0,0] == [-1,-1])
  assert all(test_ga[2,3] == [-1,-1])
  test_ga = empty_genotype_array(3, 4, ploidy=2, allele=0)
  assert all(test_ga[1,1] == [0,0])
### end MAKE EMPTY GENOTYPE ARRAY ····························································

### SNP CALLING FROM SIMULATED READS (WITH SEQUENCING ERROR)  ····························
def snp_calling(true_genotype, f_num_reads, error_rate=0.005, reads_th=8, score_th=10, ratio_th=3, dr=True, transversion=True):
    """
    snp_calling function takes perfect simulated data from one locus of one 
    diploid individual and adds missing data and error according to the number 
    of reads of the site, error rate of the sequencing technology and, for 
    ancient DNA not sequenced from damage repair (dr) libraries, creates 
    missing data for transition SNPs (since they cannot be distinguished from
    aDNA damage)

    :param true_genotype:
    :param f_num_reads:
    :param error_rate:
    :param reads_th:
    :param score_th:
    :param ratio_th:
    :param dr:
    :param transversion:
    :return:
    """
    if dr is False and transversion is False:
        genotype_call = [-1, -1]
    elif f_num_reads >= reads_th:
        derived_count = sum(true_genotype)
        p_derived = derived_count / 2. * (1 - error_rate) + (1 - derived_count / 2.) * error_rate
        derived_reads = st.binom.rvs(f_num_reads, p_derived)
        ancestral_reads = f_num_reads - derived_reads
        if f_num_reads >= score_th:
            if derived_reads == 0:
                genotype_call = [0, 0]
            elif ancestral_reads == 0:
                genotype_call = [1, 1]
            else:
                genotype_call = [1, 1]
                ratio_of_scores = derived_reads / ancestral_reads
                if (ratio_of_scores >= 1 / ratio_th) & (ratio_of_scores <= ratio_th):
                    genotype_call = [0, 1]
                elif derived_reads > ancestral_reads:
                    genotype_call = [1, 1]
                else:
                    genotype_call = [0, 0]
        else:
            random_allele = st.binom.rvs(1, derived_reads / f_num_reads)
            genotype_call = np.full(2, random_allele)
    else:
        genotype_call = [-1, -1]
    return genotype_call
def test_snp_calling():
  np.random.seed(1234)
  genotype_call = snp_calling( [0, 1], 100, error_rate=0.005, reads_th=1,
                score_th=10, ratio_th=3, dr=True, transversion=True)
  assert genotype_call == [0,1]
  genotype_call = snp_calling( [0, 1], 1, error_rate=0.005, reads_th=10,
                score_th=10, ratio_th=3, dr=True, transversion=True)
  assert genotype_call == [-1,-1]
### end SNP CALLING FROM SIMULATED READS (WITH SEQUENCING ERROR)  ····························

### SIMULATE SEQUENCING  ····························
def sequencing(ts, ssize, ttr, seq_error, dr, cov):
  if len(cov) != ssize:
    msg = "Number of coverage values (length=" + str(len(cov)) + \
          ") and number of samples (ssize=" + str(ssize) + \
          ") do not match"
    raise ValueError(msg)

  geno_data = empty_genotype_array(n_loci=ts.num_sites,
                                   n_samples=ssize,
                                   ploidy=2)
  positions = []
  locus = 0
  for variant in ts.variants():
    positions.append(round(variant.position))
    #print(variant.position)
    var_genotypes = variant.genotypes
    # print("--------------------------------")
    # print(var_genotypes)
    num_reads = np.random.poisson(lam=cov, size=ssize)
    transversion_snp = True
    if np.random.random() < ttr / (ttr + 1):
      transversion_snp = False
    # print(num_reads)
    for i in range(0, 2 * ssize, 2):
      if len(variant.alleles)==2:
        genotype_call = snp_calling(true_genotype=var_genotypes[i:(i + 2)],
                                                  f_num_reads=num_reads[int(i / 2)],
                                                  error_rate=seq_error,
                                                  dr=dr[int(i / 2)],
                                                  transversion=transversion_snp)
      else:
        genotype_call = [-1, -1] # this removes all SNP with more than two alleles
      geno_data[locus, int(i / 2)] = genotype_call
    locus = locus + 1
    #print(locus)
  return geno_data, positions
#def test_sequencing():
  # np.random.seed(1234)
  # TODO create a ts and some test from it
### end SIMULATE SEQUENCING  ····························

### TEST DATA FOR SUMMARY STATISTISCS
#                                  ind 0   ind 1   ind 2   ind 3  
test_ga_A = allel.GenotypeArray([[[ 0, 0],[ 0, 0],[ 0, 0],[ 0, 0]], # locus 0
                                 [[ 0, 1],[ 0, 1],[ 0, 1],[ 1, 1]], # locus 1
                                 [[-1,-1],[-1,-1],[-1,-1],[-1,-1]], # locus 2
                                 [[ 1, 1],[ 1, 1],[ 0, 1],[ 1, 1]], # locus 3
                                 [[ 0, 1],[ 0, 0],[ 0, 1],[ 0,-1]]],# locus 4
                                 dtype='i1')
#              0   1   2   3   4
test_pos_A = (10,123,234,299,340)
#                                  ind 0   ind 1   ind 2   ind 3  
test_ga_B = allel.GenotypeArray([[[ 0, 1],[ 0, 1],[ 1, 1],[ 0, 0]], # locus 0
                                 [[ 0, 1],[ 0, 1],[ 0, 1],[ 1, 1]], # locus 1
                                 [[ 0, 1],[ 0, 1],[-1,-1],[ 1, 1]], # locus 2
                                 [[ 0, 1],[ 0, 0],[ 0, 1],[ 1, 1]], # locus 3
                                 [[ 0, 1],[ 0, 1],[ 0, 0],[ 1, 1]], # locus 4
                                 [[ 0, 1],[ 0, 1],[ 0, 0],[ 1, 1]], # locus 5
                                 [[ 1, 1],[ 0, 1],[ 0, 1],[ 1, 1]], # locus 6
                                 [[ 0, 1],[ 0, 1],[ 0, 0],[ 1, 1]], # locus 7
                                 [[ 1, 1],[ 0, 1],[ 0, 0],[ 0, 1]], # locus 8
                                 [[ 0, 1],[ 0, 1],[ 0, 0],[ 1, 1]], # locus 9
                                 [[ 1, 1],[ 0, 1],[ 0, 1],[ 1, 1]], # locus 10
                                 [[ 0, 1],[ 0, 1],[ 0, 0],[ 1, 1]], # locus 11
                                 [[ 0, 0],[ 0, 1],[ 0, 1],[ 0, 1]], # locus 12
                                 [[-1,-1],[-1,-1],[-1,-1],[ 1, 1]], # locus 13
                                 [[ 1, 1],[ 1, 1],[ 0, 0],[ 1, 1]], # locus 14
                                 [[ 1, 1],[ 1, 1],[ 0, 0],[ 1, 1]], # locus 15
                                 [[ 1, 1],[ 0, 1],[ 0, 0],[ 0, 1]], # locus 16
                                 [[ 1, 1],[ 1, 1],[ 1, 0],[ 1, 1]], # locus 17
                                 [[ 1, 1],[ 1, 1],[ 0, 0],[ 1, 1]], # locus 18
                                 [[ 0, 1],[ 0, 0],[ 0, 1],[-1,-1]]],# locus 19
                                 dtype='i1')
#             0  1  2  3  4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19
test_pos_B = (4,10,50,77,99,123,134,150,178,201,209,234,256,270,299,311,315,340,358,378)

#                                 ga,        pos, nchr, chr_ends, w_size, expected_S
testdata_S = [pytest.param(test_ga_A, test_pos_A,    1,    [400],     50,          3, id="A"),
              pytest.param(test_ga_B, test_pos_B,    1,    [400],     50,         19, id="B")]
#                                 ga,        pos, nchr, chr_ends, w_size, expected_Pi
testdata_Pi = [pytest.param(test_ga_A, test_pos_A,    1,    [400],     50, 0.00315476, id="A"),
               pytest.param(test_ga_B, test_pos_B,    1,    [400],     50, 0.02418452, id="B")]
#                                 ga,        pos, nchr, chr_ends, w_size,  expected_WT
testdata_WT = [pytest.param(test_ga_A, test_pos_A,    1,    [400],     50, 0.00289256, id="A"),
               pytest.param(test_ga_B, test_pos_B,    1,    [400],     50, 0.01831956, id="B")]







### SAVE RESULTS from st.describe() INTO DICT ····························
def save_moments_2_dict(moments,sumstats,sample_name,sep,stat_name):
  sumstats[sample_name+sep+"min"+stat_name] = float(np.ma.getdata(moments[1][0]))
  sumstats[sample_name+sep+"max"+stat_name] = float(np.ma.getdata(moments[1][1]))
  sumstats[sample_name+sep+"m"+stat_name] = float(np.ma.getdata(moments[2]))
  sumstats[sample_name+sep+"v"+stat_name] = float(np.ma.getdata(moments[3]))
  sumstats[sample_name+sep+"s"+stat_name] = float(np.ma.getdata(moments[4]))
  sumstats[sample_name+sep+"k"+stat_name] = float(np.ma.getdata(moments[5]))
### end SAVE RESULTS from st.describe() INTO DICT  ····························

### CALCULATE SUMMARY STATISTICS : SINGLE SAMPLE  ····························
def single_sample_sumstats(ga,pos,nchr,chr_ends,w_size,sumstats,name="",sep="_",quiet=True):
  ac = ga.count_alleles()
  sfs = allel.sfs_folded(ac)
  # total number of segregating sites (excluding monomorphic and all missing data)
  segsites = sum(sfs)-sfs[0]
  sumstats[name+sep+"S"]=segsites
  if quiet is False: print("Segregating sites: "+ str(segsites))
  # Observed heterozygosity and Inbreeding coefficient
  ho = st.describe(allel.heterozygosity_observed(ga), nan_policy='omit')
  save_moments_2_dict(ho,sumstats,name,sep,"Ho")
  fis = st.describe(allel.inbreeding_coefficient(ga), nan_policy='omit')
  save_moments_2_dict(fis,sumstats,name,sep,"Fis")
  if quiet is False: print("Ho: " + str(ho) + "; Fis: " + str(fis) )
  # pairwise genetic diversity
  total_pi = allel.sequence_diversity(pos, ac, start=1, stop=chr_ends[nchr-1])
  sumstats[name+sep+"Pi"]=total_pi
  w_pi, _, _, _ = allel.windowed_diversity(pos, ac, size=w_size, start= 1, stop=chr_ends[0])
  for chromo in range(1,nchr):
    temp, _, _, _ = allel.windowed_diversity(pos, ac, size=w_size, start= 1+chr_ends[chromo-1], stop=chr_ends[chromo])
    np.append(w_pi,temp)
  pi = st.describe(w_pi)
  save_moments_2_dict(pi,sumstats,name,sep,"Pi")
  if quiet is False: print("π: "+ str(total_pi) + " ; " + str(pi))
  # Watterson theta (from number of segregating sites)
  w_W_theta, _, _, _ = allel.windowed_watterson_theta(pos, ac, size=w_size, start=1, stop=chr_ends[0])
  for chromo in range(1,nchr):
    temp, _, _, _ = allel.windowed_watterson_theta(pos, ac, size=w_size, start=1+chr_ends[chromo-1], stop=chr_ends[chromo])
    np.append(w_W_theta,temp)
  W_theta = st.describe(w_W_theta)
  save_moments_2_dict(W_theta,sumstats,name,sep,"WT")
  if quiet is False: print("Watterson θ: "+ str(W_theta))
  # Tajima's D
  total_Taj_D = allel.tajima_d(ac, pos, start=1, stop=chr_ends[nchr-1])
  sumstats[name+sep+"TD"]=total_Taj_D
  w_Taj_D, _, _ = allel.windowed_tajima_d(pos, ac, size=w_size, start=1, stop=chr_ends[0])
  for chromo in range(1,nchr):
    temp, _, _ = allel.windowed_tajima_d(pos, ac, size=w_size, start=1+chr_ends[chromo-1], stop=chr_ends[chromo])
    np.append(w_Taj_D,temp)
  Taj_D = st.describe(w_Taj_D, nan_policy='omit')
  save_moments_2_dict(Taj_D,sumstats,name,sep,"TD")
  if quiet is False: print("Tajima's D: "+ str(total_Taj_D)+ " ; " + str(Taj_D))
  # distribution of sizes for naive runs of homozygosity (distance between heterozygous positions)
  roh_distribution = np.full(int(round(np.log10(w_size)))+1, 0)
  roh_distribution = roh_distribution + windowed_distribution_roh(ga, pos, w_size, start=1, stop=chr_ends[0])
  for chromo in range(1,nchr):
    roh_distribution + windowed_distribution_roh(ga, pos, w_size, start=1+chr_ends[chromo-1], stop=chr_ends[chromo])
  for i in range(0,len(roh_distribution)):
    sumstats[name+sep+"RoHD"+sep+str(i)]=roh_distribution[i]
  if quiet is False: print("Distribution of Runs of Homozygosity: "+ str(roh_distribution))
  return

@pytest.mark.parametrize("ga,pos,nchr,chr_ends,w_size,expected_S", testdata_S)
def test_single_sample_sumstats_S(ga,pos,nchr,chr_ends,w_size,expected_S):
  test_sumstats = {}
  single_sample_sumstats(ga, pos, nchr, chr_ends, w_size, test_sumstats)
  assert test_sumstats["_S"] == expected_S
  assert test_sumstats["_S"] >= 0

@pytest.mark.parametrize("ga,pos,nchr,chr_ends,w_size,expected_Pi", testdata_Pi)
def test_single_sample_sumstats_Pi(ga,pos,nchr,chr_ends,w_size,expected_Pi):
  test_sumstats = {}
  single_sample_sumstats(ga, pos, nchr, chr_ends, w_size, test_sumstats)
  assert test_sumstats["_Pi"] == pytest.approx(expected_Pi)
  assert test_sumstats["_Pi"] >= 0
  assert test_sumstats["_mPi"] == pytest.approx(expected_Pi)

@pytest.mark.parametrize("ga,pos,nchr,chr_ends,w_size,expected_WT", testdata_WT)
def test_single_sample_sumstats_WT(ga,pos,nchr,chr_ends,w_size,expected_WT):
  test_sumstats = {}
  single_sample_sumstats(ga, pos, nchr, chr_ends, w_size, test_sumstats)
  assert test_sumstats["_mWT"] >= 0
  assert test_sumstats["_mWT"] == pytest.approx(expected_WT)



def test_single_sample_sumstats():
  test_nchr = 1
  test_chr_end = [400]
  test_w_s = 50
  # 4 individuals (columns), 20 loci (rows)
  test_ga = allel.GenotypeArray([[[ 0, 1],[ 0, 1],[ 1, 1],[ 0, 0]], #  0
                                 [[ 0, 1],[ 0, 1],[ 0, 1],[ 1, 1]], #  1
                                 [[ 0, 1],[ 0, 1],[-1,-1],[ 1, 1]], #  2
                                 [[ 0, 1],[ 0, 0],[ 0, 1],[ 1, 1]], #  3
                                 [[ 0, 1],[ 0, 1],[ 0, 0],[ 1, 1]], #  4
                                 [[ 0, 1],[ 0, 1],[ 0, 0],[ 1, 1]], #  5
                                 [[ 1, 1],[ 0, 1],[ 0, 1],[ 1, 1]], #  6
                                 [[ 0, 1],[ 0, 1],[ 0, 0],[ 1, 1]], #  7
                                 [[ 1, 1],[ 0, 1],[ 0, 0],[ 0, 1]], #  8
                                 [[ 0, 1],[ 0, 1],[ 0, 0],[ 1, 1]], #  9
                                 [[ 1, 1],[ 0, 1],[ 0, 1],[ 1, 1]], # 10
                                 [[ 0, 1],[ 0, 1],[ 0, 0],[ 1, 1]], # 11
                                 [[ 0, 0],[ 0, 1],[ 0, 1],[ 0, 1]], # 12
                                 [[-1,-1],[-1,-1],[-1,-1],[ 1, 1]], # 13
                                 [[ 1, 1],[ 1, 1],[ 0, 0],[ 1, 1]], # 14
                                 [[ 1, 1],[ 1, 1],[ 0, 0],[ 1, 1]], # 15
                                 [[ 1, 1],[ 0, 1],[ 0, 0],[ 0, 1]], # 16
                                 [[ 1, 1],[ 1, 1],[ 1, 0],[ 1, 1]], # 17
                                 [[ 1, 1],[ 1, 1],[ 0, 0],[ 1, 1]], # 18
                                 [[ 0, 1],[ 0, 0],[ 0, 1],[-1,-1]]],# 19
                                 dtype='i1')
  #           0  1  2  3  4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19
  test_pos = (4,10,50,77,99,123,134,150,178,201,209,234,256,270,299,311,315,340,358,378)
  test_sumstats = {}
  test_name = "test"
  single_sample_sumstats(test_ga, test_pos, test_nchr, test_chr_end, test_w_s, test_sumstats, test_name)
  assert test_sumstats["test_maxTD"] == pytest.approx(1.71931236)
  #assert (test_sumstats["test_RoHD"] == [3,8,0]).all()
  assert test_sumstats["test_RoHD_0"] == 3
  assert test_sumstats["test_RoHD_1"] == 8
  assert test_sumstats["test_RoHD_2"] == 0
  test_sumstats = {}
  single_sample_sumstats(test_ga, test_pos, test_nchr, test_chr_end, 390, test_sumstats, test_name)
  #assert (test_sumstats["test_RoHD"] == [3,25,1,0]).all()
  assert test_sumstats["test_RoHD_0"] == 3
  assert test_sumstats["test_RoHD_1"] == 25
  assert test_sumstats["test_RoHD_2"] == 1
  assert test_sumstats["test_RoHD_3"] == 0



### end CALCULATE SUMMARY STATISTICS : SINGLE SAMPLE  ····························


### CALCULATE SUMMARY STATISTICS : ROH  ····························

def roh(gv,pos,missing_threshold=3):
  num_of_alt = gv.to_n_alt(fill=-1)
  missing_sites = np.where(num_of_alt==-1)[0]
  if missing_threshold > np.size(missing_sites) :
    het_sites = np.where(num_of_alt==1)[0]
    roh_limits = list(map(pos.__getitem__, list(het_sites)))
    roh = np.diff(roh_limits)
  else :
    roh = np.array([])
  return roh

def test_roh():
  # one individual, 11 loci
  test_ga = allel.GenotypeArray([[[0,1]],
                                 [[0,1]],
                                 [[0,0]],
                                 [[0,1]],
                                 [[0,1]],
                                 [[0,0]],
                                 [[0,1]],
                                 [[0,1]],
                                 [[1,1]],
                                 [[0,1]],
                                 [[0,1]]], dtype='i1')
  test_pos = (5,10,15,20,100,200,300,1300,3000,10000,20000)
  test_roh = roh(test_ga,test_pos)
  assert (test_roh == [5,10,80,200,1000,8700,10000]).all()
  test_ga = allel.GenotypeArray([[[ 0, 1]],
                                 [[ 0, 1]],
                                 [[ 0, 0]],
                                 [[-1,-1]],
                                 [[ 0, 1]],
                                 [[ 0, 0]],
                                 [[ 0, 1]],
                                 [[-1,-1]],
                                 [[-1,-1]],
                                 [[ 0, 1]],
                                 [[ 0, 1]]], dtype='i1')
  test_roh = roh(test_ga,test_pos)
  assert np.size(test_roh) == 0
  test_roh = roh(test_ga,test_pos,4)
  assert (test_roh == [5,90,200,9700,10000]).all()


def distribution_roh(ga, pos, w_start, w_stop, number_of_bins):
  # w_start and w_stop are the limits of the windows as positions in the genotype array, not positions in the genome
  if number_of_bins>1:
    roh_distribution = np.full(number_of_bins, 0)
  else:
    msg = "Negative value or zero. Number of bins has to be a positive integer"
    raise ValueError(msg)
  for ind in range(0,ga.n_samples):
    gv = ga[w_start:w_stop,ind]
    roh_lenghts = roh(gv, pos[w_start:w_stop])
    for roh_size in range(0,number_of_bins):
      roh_distribution[roh_size] += sum( roh_lenghts >= 10**(roh_size) ) - sum( roh_lenghts >= 10**(roh_size+1) )
  if (roh_distribution < 0).any():
    msg = "Negative value. Number of observations of RoH bin sizes has to be zero or higher"
    raise ValueError(msg)
  return roh_distribution

def test_distribution_roh():
  # 3 individuals, 11 loci
  test_ga = allel.GenotypeArray([[[ 0, 1],[ 0, 1],[ 0, 0]], #  0
                                 [[ 0, 1],[ 0, 1],[ 0, 1]], #  1
                                 [[ 0, 0],[-1,-1],[ 1, 1]], #  2
                                 [[ 0, 1],[ 1, 1],[ 0, 1]], #  3
                                 [[ 0, 1],[ 0, 1],[ 0, 1]], #  4
                                 [[ 0, 0],[ 0, 0],[-1,-1]], #  5
                                 [[ 0, 1],[ 0, 1],[ 0, 1]], #  6
                                 [[ 0, 1],[ 0, 1],[ 0, 1]], #  7
                                 [[ 1, 1],[ 0, 0],[ 0, 1]], #  8
                                 [[ 0, 1],[ 0, 1],[ 0, 1]], #  9
                                 [[ 0, 1],[ 0, 1],[ 0, 0]]],# 10
                                 dtype='i1')
  #            0   1   2   3    4    5    6     7     8      9     10
  test_pos = [ 5, 10, 15, 20, 100, 200, 300, 1300, 3000, 10000, 20000]
  test_start = 0
  test_end = 11
  number_of_bins = 5
  d_roh = distribution_roh(test_ga, test_pos, test_start, test_end, number_of_bins)
  assert (d_roh == [2,5,3,7,2]).all()


def windowed_distribution_roh(ga, pos, size, start, stop):
  # start and stop are the limts (in bp) of the genome whre the RoH are computed
  number_of_bins = int(round(np.log10(size)))+1
  roh_distribution = np.full(number_of_bins, 0)
  # verify that positions are monotonically increasing
  if np.all(pos[1:] > pos[:-1]):
    windows = allel.position_windows(pos, size=size, start=start, stop=stop, step=size)
    locs = allel.window_locations(pos, windows)
    for window_start, window_stop in locs :
      roh_distribution = roh_distribution + distribution_roh(ga, pos, window_start, window_stop, number_of_bins)
  else:
    msg = "Wrong order. Vector of positions has to be monotonically increasing"
    raise ValueError(msg)
  return roh_distribution

def test_windowed_distribution_roh():
  # 3 individuals, 11 loci
  test_ga = allel.GenotypeArray([[[ 0, 1],[ 0, 1],[ 0, 0]], #  0
                                 [[ 0, 1],[ 0, 1],[ 0, 1]], #  1
                                 [[ 0, 0],[-1,-1],[ 1, 1]], #  2
                                 [[ 0, 1],[ 1, 1],[ 0, 1]], #  3
                                 [[ 0, 1],[ 0, 1],[ 0, 1]], #  4
                                 [[ 0, 0],[ 0, 0],[-1,-1]], #  5
                                 [[ 0, 1],[ 0, 1],[ 0, 1]], #  6
                                 [[ 0, 1],[ 0, 1],[ 0, 1]], #  7
                                 [[ 1, 1],[ 0, 0],[ 0, 1]], #  8
                                 [[ 0, 1],[ 0, 1],[ 0, 1]], #  9
                                 [[ 0, 1],[ 0, 1],[ 0, 0]]],# 10
                                 dtype='i1')
  #            0   1   2   3    4    5    6     7     8      9     10
  test_pos = [ 5, 10, 15, 20, 100, 200, 300, 1300, 3000, 10000, 20000]
  test_size = 5000
  d_roh = windowed_distribution_roh(test_ga, test_pos, test_size, 1, 20000)
  assert (d_roh == [2,5,3,4,0]).all()

### CALCULATE SUMMARY STATISTICS : TWO SAMPLE  ····························
def two_samples_sumstats(ga,pair_of_groups,pos,nchr,chr_ends,w_size,sumstats,name="",sep="_"):
  a, b, c = allel.weir_cockerham_fst(g       = ga,
                                     subpops = pair_of_groups)
  fst = np.sum(a) / (np.sum(a) + np.sum(b) + np.sum(c))
  sumstats[name+sep+"Fst"]=fst
  fst_per_variant = (np.sum(a, axis=1) / (np.sum(a, axis=1) + np.sum(b, axis=1) + np.sum(c, axis=1)))
  moments_fst_per_variant = st.describe(fst_per_variant, nan_policy='omit')
  save_moments_2_dict(moments_fst_per_variant,sumstats,name,sep,"l_Fst")
  fst_per_window, _, _ = allel.windowed_weir_cockerham_fst(pos     = pos,
                                                           g       = ga,
                                                           subpops = pair_of_groups,
                                                           size    = w_size,
                                                           start   = 1,
                                                           stop    = chr_ends[0])
  for chromo in range(1,nchr):
    temp, _, _ = allel.windowed_weir_cockerham_fst(pos     = pos,
                                                   g       = ga,
                                                   subpops = pair_of_groups,
                                                   size    = w_size,
                                                   start   = 1+chr_ends[chromo-1],
                                                   stop    = chr_ends[chromo])
    np.append(fst_per_window,temp)
  moments_fst_per_window = st.describe(fst_per_window)
  save_moments_2_dict(moments_fst_per_window,sumstats,name,sep,"w_Fst")
  return

