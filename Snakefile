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

import configparser
import numpy as np

# check if user specified a config file in command line (normal behaviour except for tests)
try:
  config["options_file"]
except:
  config["options_file"] = 'tests/config_project.ini'
options_file = config["options_file"]

# check if user specified a python path in command line (necessary in -some?- clusters)
try:
  config["python_path"]
except:
  config["python_path"] = 'python'
python_path  = config["python_path"]

# READ OPTIONS FILE
options      = configparser.ConfigParser()
options.read(options_file)

# SETTINGS
project        = options.get('Settings','project')
num_of_batches = options.getint('Settings','num_of_batches')
batches        = range(1,num_of_batches+1)
num_of_sims    = options.getint('Settings','num_of_sims')
sims           = range(1,num_of_sims+1)
sample_file    = options.get('Settings','sample_file')
genome_file    = options.get('Settings','genome_file')
seed           = options.getint('Settings','seed')
np.random.seed(seed)
seeds = np.random.randint(1, 2**32-1, num_of_batches)

# run simulations
rule afforestation:
    input:
        ref_table_sumstats = expand('results/{p}/{b}/ref_table_sumstats.pkl',p=project,b=batches),
        ref_table_latent_variables = expand('results/{p}/{b}/ref_table_latent_variables.pkl',p=project,b=batches),
        ref_table_params = expand('results/{p}/{b}/ref_table_params.pkl',p=project,b=batches)
    resources:
        runtime_min=10

# pool stats and params in a single reference table
rule reftable:
    input:
        ref_table_sumstats = expand('results/{p}/{b}/ref_table_sumstats.pkl',p=project,b=batches),
        ref_table_latent_variables = expand('results/{p}/{b}/ref_table_latent_variables.pkl',p=project,b=batches)
    resources:
        runtime_min=10
rule reftable_batch:
    input:
        script = 'scripts/reftable.py',
        sumstats_files = expand('results/{{p}}/{{b}}/sumstats_{s}.pkl',s=sims),
        latent_variables_files = expand('results/{{p}}/{{b}}/latent_variables_{s}.txt',s=sims)
    output:
        ref_table_sumstats = 'results/{p}/{b}/ref_table_sumstats.pkl',
        ref_table_latent_variables = 'results/{p}/{b}/ref_table_latent_variables.pkl'
    resources:
        runtime_min=10
    shell:
        '{python_path} {input.script} {options_file} {wildcards.b}'



# simulation of mutations with msprime
rule mutsim:
    input:
        script = 'scripts/mutsim.py',
        sim_ini = 'results/{p}/{b}/sim_{s}.ini',
        forwsim_trees = 'results/{p}/{b}/forwsim_{s}.trees'
    output:
        sumstats_files = temp('results/{p}/{b}/sumstats_{s}.pkl')
    resources:
        runtime_min=30
    shell:
        '{python_path} {input.script} {options_file} {input.sim_ini}'


# simulation with SLiM
rule forwsim:
    input:
        script='scripts/forwsim.slim',
        slim_options='results/{p}/{b}/sim_{s}.eidos',
        coalsim_trees='results/{p}/{b}/coalsim_{s}.trees'
    output:
        forwsim_trees = temp('results/{p}/{b}/forwsim_{s}.trees'),
        latent_variables_file = temp('results/{p}/{b}/latent_variables_{s}.txt')
    resources:
        runtime_min=30
    shell:
        'slim -l 0 -d "option_file=\'{input.slim_options}\'" {input.script}'

# simulation of tree sequence with msprime
rule coalsim:
    input:
        script='scripts/coalsim.py',
        sim_ini='results/{p}/{b}/sim_{s}.ini'
    output:
        coalsim_trees=temp('results/{p}/{b}/coalsim_{s}.trees')
    resources:
        runtime_min=120
    shell:
        '{python_path} {input.script} {options_file} {input.sim_ini}'

# read parameters and sample from priors
rule getparams:
    input:
        ref_table_params = expand('results/{p}/{b}/ref_table_params.pkl',p=project,b=batches),
        latent_variables_file = expand('results/{p}/{b}/latent_variables.txt',p=project,b=batches),
        slim_options = expand('results/{p}/{b}/sim_{s}.eidos',p=project,b=batches,s=sims),
        sim_ini = expand('results/{p}/{b}/sim_{s}.ini',p=project,b=batches,s=sims)
    resources:
        runtime_min=10
rule getparams_batch:
    input:
        script = 'scripts/getparams.py'
    params:
        seed = lambda wildcards, seeds=seeds: seeds[int(wildcards.b)-1]
    output:
        slim_options = temp(expand('results/{{p}}/{{b}}/sim_{s}.eidos',s=sims)),
        sim_ini = temp(expand('results/{{p}}/{{b}}/sim_{s}.ini',s=sims)),
        ref_table_params = 'results/{p}/{b}/ref_table_params.pkl'
    resources:
        runtime_min=10
    shell:
        '{python_path} {input.script} {options_file} {wildcards.b} {params.seed}'
    
# delete all files in project directory
rule clean_project:
    shell: 'rm -rf results/'+project

# delete all files in results directory
rule clean:
    shell: 'rm -rf results/'

# run tests
rule test:
    shell:
        '''
        flake8 scripts/timeadapt.py --count --select=E9,F63,F7,F82 --show-source --statistics
        pytest -v scripts/timeadapt.py
        pytest -v scripts/test_timeadapt.py
        '''



