import configparser  # for reading ini files
import pandas as pd  # for reading "table" files 
import math

### GET OPTIONS ··············································································
def get_options(proj_options_file,sim_options_file):
  proj_options = configparser.ConfigParser()
  proj_options.read(proj_options_file)
  sim_options = configparser.ConfigParser()
  sim_options.read(sim_options_file)

  project     = proj_options.get('Settings','project')
  batch       = proj_options.get('Settings','batch')
  genome_file = proj_options.get('Settings','genome_file')

  sim         = sim_options.get('Simulation','sim')

  N  = [int(i) for i in sim_options.get("Demography","N").split()]   

  return project, batch, sim, genome_file, N


def test_get_options():
  project, batch, sim, genome_file, N = get_options(proj_options_file = "tests/input/config_project.ini",
                                   sim_options_file  = "tests/input/sim_1.ini")
  assert project == "test"
  assert batch == "1"
  assert sim == "1"
  assert genome_file == "tests/input/human_genome.txt"
  assert N == [20,10]
  assert N[0] == 20
 

### end GET OPTIONS ··········································································

### GET GENOME MAP ····································································
def get_genome_map(gf):
  table = pd.read_table(filepath_or_buffer=gf, sep=r'\s+')
  nchr = max(table["Chromosome"])
  rates = table["Recombination_rate"]
  lengths = table["Length"]
  ends = pd.Series(lengths).cumsum()
  return nchr, rates, ends
  

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
  assert rates ==  [1e-8,math.log(2),2e-8,math.log(2),3e-8]


### end MAKE RECOMBINATION MAP ························································




