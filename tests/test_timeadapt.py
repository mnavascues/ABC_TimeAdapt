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
import importlib.util
spec = importlib.util.spec_from_file_location("timeadapt", "scripts/timeadapt.py")
timeadapt = importlib.util.module_from_spec(spec)
spec.loader.exec_module(timeadapt)
# import timeadapt
import tempfile # for creating temporal files on testing
import allel
import numpy as np
import pandas as pd
import math
import pytest

# TEST GET OPTIONS #############################################################################################

_, temp_config_file_1 = tempfile.mkstemp()
with open(temp_config_file_1, 'w') as f:
  f.write("[Settings]\n")
  f.write("project=test\n")
  f.write("verbose=0\n")
  f.write("project_dir=results/test\n")
  f.write("[Model]\n")
  f.write("periods_coalescence=10\n")
  f.write("times_of_change_back=196 587 1370 2935 6067 12329 24853 49902 100000\n")
  f.write("[Sample]\n")
  f.write("size=17\n")
  f.write("coverage=40.57 46.49 34.22 37.24 12.94 1.25 0.01 1.1 2.3 0.8 1 11.09647 0.227619 11.58128 4.998657 0.130293 22.19364\n")
  f.write("is_damaged=FALSE FALSE FALSE FALSE FALSE FALSE TRUE TRUE TRUE TRUE TRUE FALSE TRUE FALSE FALSE TRUE FALSE\n")
  f.write("group_levels=3\n")
  f.write("groups=aaa aaa aaa aaa bbb bbc bbc bbc bbc bbb bbb bcd bbc bcc bcd bcd bcd\n")
  f.write("[Genome]\n")
  f.write("nchr=4\n")
  f.write("chr_ends=10000000 20000000 30000000 40000000\n")
  f.write("msprime_r_map_positions=0 1000000 2000000 3000000 4000000 10000000 10000001 11000000 12000000 14000000 18000000 20000000 20000001 22000000 24000000 26000000 30000000 30000001 33000000 35000000 37000000 40000000\n")
  f.write("msprime_r_map_rates=1e-08 1e-07 1e-08 1e-09 1e-08 0.693147180559945 1e-08 1e-10 1e-09 1e-08 1e-07 0.693147180559945 1e-08 1e-09 1e-08 1e-09 0.693147180559945 1e-08 1e-10 1e-08 1e-07\n")

_, temp_sim_file_1 = tempfile.mkstemp()
with open(temp_sim_file_1, 'w') as f:
  f.write("[Simulation]\n")
  f.write("batch=1\n")
  f.write("sim=1\n")
  f.write("[Sample]\n")
  f.write("ss=4 1 1 2 1 1 1 1 1 1 1 1 1\n")
  f.write("msprime_ts=0 18 41 68 71 72 82 106 134 193 267 278 343\n")
  f.write("chrono_order=0 1 2 3 10 9 4 8 6 5 7 12 13 16 11 15 14\n")
  f.write("[Demography]\n")
  f.write("N=58 33 12 37 200 136 130 174 46 12\n")
  f.write("[Genome]\n")
  f.write("mu=1.53224247401956e-07\n")
  f.write("ttratio=2\n")
  f.write("seq_error=0.005\n")
  f.write("[Seeds]\n")
  f.write("seed_coal=1543026934\n")
  f.write("seed_mut=1939636257\n")

result_config_1 = {"project" : "test",
                   "verbose" : 0,
                   "periods_coalescence" : 10,
                   "times_of_change_back" : [196, 587, 1370, 2935, 6067, 12329, 24853, 49902, 100000],
                   "nchr" : 4,
                   "chr_ends" : [10000000, 20000000, 30000000, 40000000],
                   "msprime_r_map" : {"positions" : [0, 1000000, 2000000, 3000000, 4000000, 10000000, 10000001, 11000000, 12000000, 14000000, 18000000, 20000000, 20000001, 22000000, 24000000, 26000000, 30000000, 30000001, 33000000, 35000000, 37000000, 40000000],
                                      "rates" : [1e-08, 1e-07, 1e-08, 1e-09, 1e-08, 0.693147180559945, 1e-08, 1e-10, 1e-09, 1e-08, 1e-07, 0.693147180559945, 1e-08, 1e-09, 1e-08, 1e-09, 0.693147180559945, 1e-08, 1e-10, 1e-08, 1e-07]},
                   "total_sample_size" : 17,
                   "coverage" : [40.57, 46.49, 34.22, 37.24, 12.94, 1.25, 0.01, 1.1, 2.3, 0.8, 1, 11.09647, 0.227619, 11.58128, 4.998657, 0.130293, 22.19364],
                   "is_damaged" : [False, False, False, False, False, False, True, True, True, True, True, False, True, False, False, True, False],
                   "group_levels" : 3,
                   "groups" : ["aaa", "aaa", "aaa", "aaa", "bbb", "bbc", "bbc", "bbc", "bbc", "bbb", "bbb", "bcd", "bbc", "bcc", "bcd", "bcd", "bcd"]}

result_sim_1 = {"sim" : "1",
                "batch" : "1",
                "ss" : [4, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1], 
                "chrono_order" : [0, 1, 2, 3, 10, 9, 4, 8, 6, 5, 7, 12, 13, 16, 11, 15, 14], 
                "N" : [58, 33, 12, 37, 200, 136, 130, 174, 46, 12], 
                "mu" : 1.53224247401956e-07, 
                "ttratio" : 2, 
                "seq_error" : 0.005, 
                "seed_coal" : 1543026934, 
                "seed_mut" : 1939636257}

result_config_1_sim_1 = dict(result_config_1)
result_config_1_sim_1.update(result_sim_1)

test_config_files = [pytest.param(temp_config_file_1, temp_sim_file_1, result_config_1_sim_1, id = "1_1")]
#test_config_files = [pytest.param(temp_config_file_1, temp_sim_file_1, result_config_1_sim_1, id="1_1"),
#                     pytest.param(temp_config_file_2, temp_sim_file_1, result_config_2_sim_1, id="2_1")]

@pytest.mark.parametrize("temp_config_file,temp_sim_file,expected_result", test_config_files)
def test_get_options(temp_config_file,temp_sim_file,expected_result):
  options = timeadapt.get_options(proj_options_file = temp_config_file, sim_options_file = temp_sim_file)
  
  assert type(options["project"]) is str
  assert options["project"] == expected_result["project"]
  assert type(options["batch"]) is str
  assert options["batch"] == expected_result["batch"]
  assert type(options["sim"]) is str
  assert options["sim"] == expected_result["sim"]
  assert type(options["verbose"]) is int
  assert options["verbose"] ==  expected_result["verbose"]
  assert type(options["ss"]) is list
  for i in options["ss"]:
    assert type(i) is int
  assert options["ss"] == expected_result["ss"]
  assert type(options["chrono_order"]) is list
  for i in options["chrono_order"]:
    assert type(i) is int
  assert options["chrono_order"] == expected_result["chrono_order"]
  assert type(options["N"]) is list
  for i in options["N"]:
    assert type(i) is int
  assert options["N"] == expected_result["N"]
  assert type(options["mu"]) is float
  assert options["mu"] == pytest.approx(expected_result["mu"])
  assert type(options["seed_coal"]) is int
  assert options["seed_coal"] == expected_result["seed_coal"]
  assert type(options["seed_mut"]) is int
  assert options["seed_mut"] == expected_result["seed_mut"]
  assert options["times_of_change_back"] == expected_result["times_of_change_back"]
  assert options["periods_coalescence"] == expected_result["periods_coalescence"]

# TEST GET TIMES OF CHANGE  ########################################################################################

def test_get_times_of_change():
    with pytest.raises(ValueError) as e:
        timeadapt.get_times_of_change(10,10)
    assert str(e.value) == 'number of periods must be lower than length of simulation' 

# TEST READ SAMPLE INFO #############################################################################################

_, temp_sample_file_1 = tempfile.mkstemp()
with open(temp_sample_file_1, 'w') as f:
  f.write("sampleID age14C age14Cerror year coverage damageRepair groups\n")
  f.write("modern   NA     NA          2010 30.03    TRUE         0\n")
  f.write("ancient  1980   20          NA   10.01    TRUE         1\n")

result_sample_1 = {"id":["modern","ancient"],
                   "age14C":[np.nan,1980],
                   "age14Cerror":[np.nan,20],
                   "ageBCAD":[2010,np.nan],
                   "t0":2010,
                   "coverage":[30.03,10.01],
                   "is_ancient":[False,True],
                   "is_modern":[True,False],
                   "is_dr":[True,True],
                   "total_ancient":1,
                   "size":2,
                   "group_levels":1,
                   "groups":np.array([0,1])}

test_sample_files = [pytest.param(temp_sample_file_1, result_sample_1, id="1")]

@pytest.mark.parametrize("sample_file,expected_result", test_sample_files)
def test_read_sample_info(sample_file,expected_result):
  sample_info = timeadapt.read_sample_info(sample_info_file=sample_file)
  assert list(sample_info["id"]) == expected_result["id"]
  assert list(sample_info["age14C"]) == pytest.approx(expected_result["age14C"],nan_ok=True)
  assert list(sample_info["age14Cerror"]) == pytest.approx(expected_result["age14Cerror"],nan_ok=True)
  assert list(sample_info["ageBCAD"]) == pytest.approx(expected_result["ageBCAD"],nan_ok=True)
  assert sample_info["t0"] == expected_result["t0"]
  assert list(sample_info["coverage"]) == expected_result["coverage"]
  assert list(sample_info["is_ancient"]) == expected_result["is_ancient"]
  assert list(sample_info["is_modern"]) == expected_result["is_modern"]
  assert list(sample_info["is_dr"]) == expected_result["is_dr"]
  assert sample_info["total_ancient"] == expected_result["total_ancient"]
  assert sample_info["size"] == expected_result["size"]
  assert sample_info["group_levels"] == expected_result["group_levels"]
  assert (sample_info["groups"] == expected_result["groups"]).all()

# TEST READ GENOME INFO #############################################################################################

_, temp_genome_file_1 = tempfile.mkstemp()
with open(temp_genome_file_1, 'w') as f:
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

_, temp_genome_file_2 = tempfile.mkstemp()
with open(temp_genome_file_2, 'w') as f:
  f.write("Chromosome	Position	Recombination_rate\n")
  f.write("1	 100000	1E-08\n")
  f.write("1	 200000	1E-07\n")
  f.write("1	1000000	1E-08\n")

result_genome_1 = {"nchr":4,
                   "chr_ends":[10000000,20000000,30000000,40000000],
                   "msprime_r_map":{"rates":[1E-08,1E-07,1E-08,math.log(2),
                                             1E-08,1E-10,1E-07,math.log(2),
                                             1E-08,1E-09,1E-08,math.log(2),
                                             1E-08,1E-10,1E-07],
                                    "positions":[       0, 1000000, 2000000,10000000,
                                                 10000001,11000000,12000000,20000000,
                                                 20000001,22000000,24000000,30000000,
                                                 30000001,33000000,35000000,40000000]},
                   "slim_r_map":{"rates":[1E-08,1E-07,1E-08,0.5,
                                          1E-08,1E-10,1E-07,0.5,
                                          1E-08,1E-09,1E-08,0.5,
                                          1E-08,1E-10,1E-07],
                                 "positions":[           999999, 1999999, 9999999,
                                              10000000,10999999,11999999,19999999,
                                              20000000,21999999,23999999,29999999,
                                              30000000,32999999,34999999,39999999]}}

result_genome_2 = {"nchr":1,
                   "chr_ends":1000000,
                   "msprime_r_map":{"rates":[1E-08,1E-07,1E-08],
                                    "positions":[0, 1000000, 2000000,10000000]},
                   "slim_r_map":{"rates":[1E-08,1E-07,1E-08],
                                 "positions":[999999, 1999999,9999999]}}


test_genome_files = [pytest.param(temp_genome_file_1, result_genome_1, id="1"),
                     pytest.param(temp_genome_file_2, result_genome_2, id="2")]

@pytest.mark.parametrize("genome_file,expected_result", test_genome_files)
def test_get_genome_info(genome_file,expected_result):
    genome_info = timeadapt.get_genome_info(genome_file)
    assert genome_info["nchr"]==genome_info["nchr"]
    assert genome_info["chr_ends"]==genome_info["chr_ends"]
    assert genome_info["msprime_r_map"]["rates"]==genome_info["msprime_r_map"]["rates"]
    assert genome_info["msprime_r_map"]["positions"]==genome_info["msprime_r_map"]["positions"]
    assert genome_info["slim_r_map"]["rates"]==genome_info["slim_r_map"]["rates"]
    assert genome_info["slim_r_map"]["positions"]==genome_info["slim_r_map"]["positions"]

### TEST SAMPLE PARAMETER TRAJECTORY #######################################################################

test_sample_trajectory = [pytest.param(10,1,1000,id="1"),
                          pytest.param(100,50,1000000,id="2")]
@pytest.mark.parametrize("times,mininum,maximum", test_sample_trajectory)
def test_sample_param_trajectory(times,mininum,maximum):
  traj = timeadapt.sample_param_trajectory(times,mininum,maximum)
  assert all(traj>=mininum)
  assert all(traj<=maximum)
  assert np.size(traj)==times





