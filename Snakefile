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
# (specify config file by adding "-C options_file='PATH/TO/YOUR/CONFIG_FILE.ini'" to command line)
try:
  config["options_file"]
except:
  config["options_file"] = 'tests/config_project.ini'
options_file = config["options_file"]

# check if user specified a batch number
# (specify config file by adding "-C options_file='PATH/TO/YOUR/CONFIG_FILE.ini'" to command line)
try:
  config["batch"]
except:
  config["batch"] = 1
batch = config["batch"]

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
num_of_sims    = options.getint('Settings','num_of_sims')
sims           = range(1,num_of_sims+1)
sample_file    = options.get('Settings','sample_file')
genome_file    = options.get('Settings','genome_file')
seed           = options.getint('Settings','seed')
#np.random.seed(seed)
#seeds = np.random.randint(1, 2**32-1, num_of_batches)
keep_files     = True


# run simulations
rule abc:
    input:
        ref_table_sumstats         = expand('results/{p}/{b}/ref_table_sumstats.RDS',         p = project, b = batch),
        ref_table_latent_variables = expand('results/{p}/{b}/ref_table_latent_variables.RDS', p = project, b = batch),
        ref_table_params           = expand('results/{p}/{b}/ref_table_params.RDS',           p = project, b = batch)
    output:
        forest_file = touch('afforestation.done')
    resources:
        runtime_min = 10

# pool stats and latent variables in a single reference table (per batch)
rule reftable:
    input:
        script                 = 'scripts/reftable.R',
        options_RDS            = expand('results/{p}/project_options.RDS',          p = project),
        sumstats_files         = expand('results/{p}/{b}/sumstats_{s}.txt',         p = project, b = batch, s = sims),
        latent_variables_files = expand('results/{p}/{b}/latent_variables_{s}.txt', p = project, b = batch, s = sims)
    output:
        ref_table_sumstats         = expand('results/{p}/{b}/ref_table_sumstats.RDS',         p = project, b = batch),
        ref_table_latent_variables = expand('results/{p}/{b}/ref_table_latent_variables.RDS', p = project, b = batch)
    resources:
        runtime_min = 10
    shell:
        'Rscript {input.script} {input.options_RDS} {batch}'

# all simulations
rule simulate:
    input:
        sumstats_files        = expand('results/{p}/{b}/sumstats_{s}.txt',         p = project, b = batch, s = sims),
        latent_variables_file = expand('results/{p}/{b}/latent_variables_{s}.txt', p = project, b = batch, s = sims)

# simulation of mutations with msprime
rule mutsim:
    input:
        script        = 'scripts/mutsim.py',
        sim_ini       = 'results/{p}/{b}/sim_{s}.ini',
        project_ini   = 'results/{p}/project_options.ini',
        forwsim_trees = 'results/{p}/{b}/forwsim_{s}.trees'
    output:
        #sumstats_files = temp('results/{p}/{b}/sumstats_{s}.txt')
        sumstats_files = 'results/{p}/{b}/sumstats_{s}.txt'
    resources:
        runtime_min = 30
    shell:
        '{python_path} {input.script} {input.project_ini} {input.sim_ini}'


# simulation with SLiM
rule forwsim:
    input:
        script        = 'scripts/forwsim.slim',
        slim_options  = 'results/{p}/{b}/sim_{s}.eidos',
        coalsim_trees = 'results/{p}/{b}/coalsim_{s}.trees'
    output:
        #forwsim_trees         = temp('results/{p}/{b}/forwsim_{s}.trees'),
        #latent_variables_file = temp('results/{p}/{b}/latent_variables_{s}.txt')
        forwsim_trees         = 'results/{p}/{b}/forwsim_{s}.trees',
        latent_variables_file = 'results/{p}/{b}/latent_variables_{s}.txt'
    resources:
        runtime_min = 30
    shell:
        'slim -l 0 -d "option_file=\'{input.slim_options}\'" {input.script}'

# simulation of tree sequence with msprime
rule coalsim:
    input:
        script      = 'scripts/coalsim.py',
        sim_ini     = 'results/{p}/{b}/sim_{s}.ini',
        project_ini = 'results/{p}/project_options.ini'
    output:
        #coalsim_trees = temp('results/{p}/{b}/coalsim_{s}.trees')
        coalsim_trees = 'results/{p}/{b}/coalsim_{s}.trees'
    resources:
        runtime_min = 120
    shell:
        '{python_path} {input.script} {input.project_ini} {input.sim_ini}'

# read parameters and sample from priors
rule getparams:
    input:
        script      = 'scripts/getparams.R',
        options_RDS = expand('results/{p}/project_options.RDS', p = project)
    output:
        #slim_options     = temp(expand('results/{p}/{b}/sim_{s}.eidos', p = project, b = batch, s = sims)),
        #sim_ini          = temp(expand('results/{p}/{b}/sim_{s}.ini', p = project,   b = batch, s = sims)),
        slim_options     = expand('results/{p}/{b}/sim_{s}.eidos', p = project, b = batch, s = sims),
        sim_ini          = expand('results/{p}/{b}/sim_{s}.ini', p = project, b = batch, s = sims),
        ref_table_params = expand('results/{p}/{b}/ref_table_params.RDS', p = project, b = batch)
    resources:
        runtime_min = 10
    shell:
        'Rscript {input.script} {input.options_RDS} {batch}'

# read project parameters
rule setproject:
    input:
        script = 'scripts/setproject.R',
        options_file = options_file
    output:
        expand('results/{p}/project_options.RDS', p = project),
        expand('results/{p}/project_options.ini', p = project)
    resources:
        runtime_min = 10
    shell:
        'Rscript {input.script} {input.options_file}'

# delete all files in project directory
rule clean_project:
    shell: 'rm -rf results/' + project

# delete all files in results directory
rule clean:
    shell: 'rm -rf results/'

# run tests
rule test:
    shell:
        '''
        echo "\n----------- RUNNING TESTTHAT (R) ------------------\n"
        Rscript -e "library(testthat); test_file(\'tests/test_timeadapt.R\')"
        echo "\n----------- RUNNING FLAKE8 (PYTHON) ---------------\n"
        flake8 scripts/timeadapt.py --count --select=E9,F63,F7,F82 --show-source --statistics
        echo "\n----------- RUNNING PYTEST (PYTHON) ---------------\n"
        #pytest -v scripts/timeadapt.py
        pytest -v tests/test_timeadapt.py
        '''

