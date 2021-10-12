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
import allel
import math
import tempfile # for creating temporal files on testing
import pytest

### GET OPTIONS ··············································································
def get_options(proj_options_file,sim_options_file):
  proj_options = configparser.ConfigParser()
  proj_options.read(proj_options_file)
  sim_options = configparser.ConfigParser()
  sim_options.read(sim_options_file)

  project      = proj_options.get('Settings','project')
  batch        = proj_options.get('Settings','batch')
  genome_file  = proj_options.get('Settings','genome_file')
  sample_file  = proj_options.get('Settings','sample_file')
  sim          = sim_options.get('Simulation','sim')
  ss           = [int(i) for i in sim_options.get("Sample","ss").split()]
  chrono_order = [int(i) for i in sim_options.get("Sample","chrono_order").split()]
  N            = [int(i) for i in sim_options.get("Demography","N").split()]   
  mu           = sim_options.getfloat('Genome','mu')
  ttratio      = sim_options.getfloat('Genome','ttratio')
  seq_error    = sim_options.getfloat('Genome','seq_error')
  seed_coal    = sim_options.getint('Seeds','seed_coal')
  seed_mut     = sim_options.getint('Seeds','seed_mut')

  return project, batch, sim, genome_file, sample_file, ss, chrono_order, N, mu, ttratio, seq_error, seed_coal, seed_mut
def test_get_options():
  project, batch, sim, genome_file, sample_file, ss, chrono_order, N, mu, ttratio, seq_error, seed_coal, seed_mut = \
           get_options(proj_options_file = "tests/input/config_project.ini", sim_options_file  = "tests/input/sim_1.ini")
  assert project == "test"
  assert batch == "1"
  assert sim == "1"
  assert genome_file == "tests/input/test_genome.txt"
  assert sample_file == "tests/input/test_sample.txt"
  assert ss == [4,1,1,1,1,1,1,1,1,1,1,1,1,1]
  assert chrono_order == [0,1,2,3,10,9,4,8,6,5,7,12,13,16,11,15,14]
  assert N == [11,85,200,200,200,37,10,30,71]
  assert N[0] == 11
  assert seed_coal == 1066941363
  assert seed_mut == 4083239485
### end GET OPTIONS ··········································································

### READ SAMPLE INFO FILE .........
def read_sample_info(sample_info_file):
    """
    Read text file with information on sample. Format of the file:
    --------------------------------------------------------------------------
    sampleID           age14C  age14Cerror  year  coverage  damageRepair  groups
    B_Ju_hoan_North-4  NA      NA           2010  40.57     TRUE          00
    S_Ju_hoan_North-1  NA      NA           2010  46.49     TRUE          00
    BallitoBayA        1980    20           NA    12.94     FALSE         11
    BallitoBayB        2110    30           NA    1.25      TRUE          12
    --------------------------------------------------------------------------

    :param sample_info_file: path of file
    :return:
    """
    # TODO : change damageRepair for ancientDamage (which should be more intuitive)
    # TODO : make it work with only ancient data (i.e. no year column)

    info = pd.read_table(filepath_or_buffer=sample_info_file, sep=r'\s+',
                         converters={'groups': lambda x: str(x)})

    sample_id = info["sampleID"]
    coverage = info["coverage"]
    is_dr = info["damageRepair"]

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

    return sample_id, coverage, is_ancient, is_modern, is_dr, total_ancient, \
           sample_size, group_levels, groups
def test_read_sample_info():
    _, temporary_file_name = tempfile.mkstemp()
    with open(temporary_file_name, 'w') as f:
        f.write("sampleID age14C age14Cerror year coverage damageRepair groups\n")
        f.write("modern   NA     NA          2010 30.03    TRUE         0\n")
        f.write("ancient  1980   20          NA   10.01    TRUE         1\n")
    sample_id, coverage, is_ancient, is_modern, is_dr, total_ancient, \
    sample_size, group_levels, \
    groups = read_sample_info(sample_info_file=temporary_file_name)
    assert sample_id[0] == "modern"
    assert sample_id[1] == "ancient"
    assert coverage[0] == 30.03
    assert is_ancient[0] is False
    assert is_ancient[1] is True
    assert is_modern[0] is True
    assert is_modern[1] is False
    assert is_dr[0]
    assert is_dr[1]
    assert total_ancient == 1
    assert sample_size == 2
    assert group_levels == 1
### end READ SAMPLE INFO FILE .........

### GET GENOME MAP ····································································
def get_genome_map(gf):
  table = pd.read_table(filepath_or_buffer=gf, sep=r'\s+')
  nchr = max(table["Chromosome"])
  rates = table["Recombination_rate"]
  lengths = table["Length"]
  ends = pd.Series(lengths).cumsum()
  return nchr, rates, ends
def get_recombination_map(gf):
  table = pd.read_table(filepath_or_buffer=gf, sep=r'\s+')
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
  # insert recombination rate between chromosomes
  for chromo in range(0,nchr-1):
    positions.insert(chr_ends_index[chromo]+1+chromo,positions[chr_ends_index[chromo]+chromo]+1)
    rates.insert(chr_ends_index[chromo]+1+chromo,math.log(2))
  positions.insert(0,0) # insert first position
  return nchr, chr_ends, rates, positions
def test_get_recombination_map():
    _, temporary_file_name = tempfile.mkstemp()
    with open(temporary_file_name, 'w') as f:
      f.write("Chromosome	Position	Recombination_rate\n")
      f.write("1	 1000000	1E-08\n")
      f.write("1	 2000000	1E-07\n")
      f.write("1	10000000	1E-08\n")
      f.write("2	 1000000	1E-08\n")
      f.write("2	 2000000	1E-10\n")
      f.write("2	10000000	1E-07\n")
      f.write("3	 2000000	1E-08\n")
      f.write("3	 4000000	1E-09\n")
      f.write("3	10000000	1E-08\n")
      f.write("4	 3000000	1E-08\n")
      f.write("4	 5000000	1E-10\n")
      f.write("4	10000000	1E-07\n")
    num_of_chr, chr_ends, rates, positions = get_recombination_map(temporary_file_name)
    assert num_of_chr==4
    assert (chr_ends == [10000000,20000000,30000000,40000000])
    assert (rates == [1E-08,1E-07,1E-08,math.log(2),
                      1E-08,1E-10,1E-07,math.log(2),
                      1E-08,1E-09,1E-08,math.log(2),
                      1E-08,1E-10,1E-07])
    assert (positions == [       0, 1000000, 2000000,10000000,
                          10000001,11000000,12000000,20000000,
                          20000001,22000000,24000000,30000000,
                          30000001,33000000,35000000,40000000])
def test_get_genome_map():
  num_of_chr, chrom_rates, chrom_ends = get_genome_map(gf="tests/input/human_genome.txt")
  assert num_of_chr == 22
  assert (chrom_rates == [1.14856E-08,1.10543E-08,1.12796E-08,1.12312E-08,1.12809E-08,
                   1.12229E-08,1.17646E-08,1.14785E-08,1.17807E-08,1.33651E-08,
                   1.17193E-08,1.30502E-08,1.09149E-08,1.11973E-08,1.38358E-08,
                   1.48346E-08,1.58249E-08,1.5076E-08,1.82201E-08,1.71783E-08,
                   1.30452E-08,1.4445E-08]).all()
  assert (chrom_ends == [248956422,491149951,689445510,879660065,1061198324,
                   1232004303,1391350276,1536488912,1674883629,1808681051,
                   1943767673,2077042982,2191407310,2298451028,2400442217,
                   2490780562,2574038003,2654411288,2713028904,2777473071,
                   2824183054,2875001522]).all()
### end GET GENOME MAP ································································

### MAKE RECOMBINATION MAP ····························································
def make_rec_map(nchr, c_rates, c_ends):
  pos = [0]
  rat = []
  for i in range(nchr-1):
    # print(i)
    pos.append(c_ends[i])
    pos.append(c_ends[i]+1)
    rat.append(c_rates[i])
    rat.append(math.log(2))
  pos.append(c_ends[nchr-1])
  rat.append(c_rates[nchr-1])
  return pos, rat
def test_make_rec_map():
  positions, rates = make_rec_map(nchr=3, c_rates=[1e-8,2e-8,3e-8], c_ends=[10,20,30])
  assert positions == [0,10,11,20,21,30]
  assert rates ==  pytest.approx([1e-8,math.log(2),2e-8,math.log(2),3e-8])
### end MAKE RECOMBINATION MAP ························································

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

def save_moments_2_dict(moments,sumstats,sample_name,sep,stat_name):
  sumstats[sample_name+sep+"min"+stat_name] = np.ma.getdata(moments[1][0])
  sumstats[sample_name+sep+"max"+stat_name] = np.ma.getdata(moments[1][1])
  sumstats[sample_name+sep+"m"+stat_name] = np.ma.getdata(moments[2])
  sumstats[sample_name+sep+"v"+stat_name] = np.ma.getdata(moments[3])
  sumstats[sample_name+sep+"s"+stat_name] = np.ma.getdata(moments[4])
  sumstats[sample_name+sep+"k"+stat_name] = np.ma.getdata(moments[5])

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
    temp, _, _, _ = allel.windowed_tajima_d(pos, ac, size=w_size, start=1+chr_ends[chromo-1], stop=chr_ends[chromo])
    np.append(w_Taj_D,temp)
  Taj_D = st.describe(w_Taj_D, nan_policy='omit')
  save_moments_2_dict(Taj_D,sumstats,name,sep,"TD")
  if quiet is False: print("Tajima's D: "+ str(total_Taj_D)+ " ; " + str(Taj_D))
  # distribution of sizes for naive runs of homozygosity (distance between heterozygous positions)
  roh_distribution = np.full(int(round(np.log10(w_size)))+1, 0)
  roh_distribution = roh_distribution + windowed_distribution_roh(ga, pos, w_size, start=1, stop=chr_ends[0])
  for chromo in range(1,nchr):
    roh_distribution + windowed_distribution_roh(ga, pos, w_size, start=1+chr_ends[chromo-1], stop=chr_ends[chromo])
  sumstats[name+sep+"RoHD"]=roh_distribution
  if quiet is False: print("Distribution of Runs of Homozygosity: "+ str(roh_distribution))
  return
  
def test_single_sample_sumstats():
  # 4 individuals (columns), 5 loci (rows)
  test_ga = allel.GenotypeArray([[[ 0, 0],[ 0, 0],[ 0, 0],[ 0, 0]],
                                 [[ 0, 1],[ 0, 1],[ 0, 1],[ 1, 1]],
                                 [[-1,-1],[-1,-1],[-1,-1],[-1,-1]],
                                 [[ 1, 1],[ 1, 1],[ 0, 1],[ 1, 1]],
                                 [[ 0, 1],[ 0, 0],[ 0, 1],[ 0,-1]]], dtype='i1')
  test_pos = (10,123,234,299,340)
  test_nchr = 1
  test_chr_end = [400]
  test_w_s = 50
  test_sumstats = {}
  single_sample_sumstats(test_ga, test_pos, test_nchr, test_chr_end, test_w_s, test_sumstats)
  assert test_sumstats["_S"]==3
  assert test_sumstats["_Pi"] == pytest.approx(0.003154762)
  assert test_sumstats["_minPi"] == pytest.approx(0.)
  assert test_sumstats["_maxPi"] == pytest.approx(0.01071429)
  assert test_sumstats["_mPi"] == pytest.approx(0.003154762)
  assert test_sumstats["_minWT"] == pytest.approx(0.)
  assert test_sumstats["_maxWT"] == pytest.approx(0.0077135)
  assert test_sumstats["_mWT"] == pytest.approx(0.002892563)
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
  assert (test_sumstats["test_RoHD"] == [3,8,0]).all()
  test_sumstats = {}
  single_sample_sumstats(test_ga, test_pos, test_nchr, test_chr_end, 390, test_sumstats, test_name)
  assert (test_sumstats["test_RoHD"] == [3,25,1,0]).all()



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


