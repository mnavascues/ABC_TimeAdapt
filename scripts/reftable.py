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
import timeadapt

def main():
  # get options for project and simulation:
  project, batch, genome_file, _, verbose = \
           timeadapt.get_project_options(proj_options_file = sys.argv[1])

  # print program name
  if verbose >=1 :
    print("#########################################")
    print("#                                       #")
  if verbose >=0 :
    print("#      TimeAdapt - reftable.py           #")
  if verbose >=1 :
    print("#      by Miguel de Navascués           #")
    print("#      INRAE & Uppsala universitet      #")
    print("#      miguel.navascues@inrae.fr        #")
    print("#                                       #")
    print("#########################################")

############################################################################################################
if __name__ == "__main__":
    main()

