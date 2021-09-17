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
import sys
import msprime
import pyslim
import timeadapt

## how to run:
#  python scripts/coalsim.py tests/input/config_project.ini tests/input/sim_1.ini

def main():
  # get options for project and simulation:
  project, batch, sim, genome_file, _, _, _, N, _, _, _, seed_coal, _ = \
           timeadapt.get_options(proj_options_file = sys.argv[1], sim_options_file = sys.argv[2])
  
  # get recombination map:
  _, rates, positions = timeadapt.get_recombination_map(gf = genome_file)
  rate_map = msprime.RateMap(position = positions, rate = rates) 

  # simulate with msprime
  # https://tskit.dev/msprime/docs/latest/intro.html
  msp_ts = msprime.sim_ancestry(samples            = N[0],
                                population_size    = N[0],
                                model              = "dtwf",
                                recombination_rate = rate_map,
                                random_seed        = seed_coal)

  # make tree-sequence a SLiM-tree-sequence
  slim_ts = pyslim.annotate_defaults(msp_ts, model_type="WF", slim_generation=1)
  
  # save tree
  slim_ts.dump("results/"+project+"/"+batch+"/coalsim_"+sim+".trees")


############################################################################################################
############################################################################################################
if __name__ == "__main__":
    main()
