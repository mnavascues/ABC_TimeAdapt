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
import timeadapt
import allel
import pytest

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


@pytest.mark.parametrize("ga,pos,nchr,chr_ends,w_size,expected_S", testdata_S)
def test_single_sample_sumstats_S(ga,pos,nchr,chr_ends,w_size,expected_S):
  test_sumstats = {}
  timeadapt.single_sample_sumstats(ga, pos, nchr, chr_ends, w_size, test_sumstats)
  assert test_sumstats["_S"] == expected_S
  assert test_sumstats["_S"] >= 0

@pytest.mark.parametrize("ga,pos,nchr,chr_ends,w_size,expected_Pi", testdata_Pi)
def test_single_sample_sumstats_Pi(ga,pos,nchr,chr_ends,w_size,expected_Pi):
  test_sumstats = {}
  timeadapt.single_sample_sumstats(ga, pos, nchr, chr_ends, w_size, test_sumstats)
  assert test_sumstats["_Pi"] == pytest.approx(expected_Pi)
  assert test_sumstats["_Pi"] >= 0
  assert test_sumstats["_mPi"] == pytest.approx(expected_Pi)

@pytest.mark.parametrize("ga,pos,nchr,chr_ends,w_size,expected_WT", testdata_WT)
def test_single_sample_sumstats_WT(ga,pos,nchr,chr_ends,w_size,expected_WT):
  test_sumstats = {}
  timeadapt.single_sample_sumstats(ga, pos, nchr, chr_ends, w_size, test_sumstats)
  assert test_sumstats["_mWT"] >= 0
  assert test_sumstats["_mWT"] == pytest.approx(expected_WT)

