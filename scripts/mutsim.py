import sys
import msprime
import pyslim
import timeadapt

def main():
  # get options for project and simulation:
  project, batch, sim, genome_file, sample_file, N, seed_coal, seed_mut = \
           timeadapt.get_options(proj_options_file = sys.argv[1], sim_options_file = sys.argv[2])

  # read tree sequence from SLiM output file
  treesq = pyslim.load("results/"+project+"/"+batch+"/forwsim_"+sim+".trees")
  #print(treesq)
  tree = treesq.first()
  print(tree.draw(format="unicode"))




############################################################################################################
if __name__ == "__main__":
    main()
