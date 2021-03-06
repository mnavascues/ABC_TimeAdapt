import configparser

# READ OPTIONS FILE
options_file = 'tests/input/config_project.ini'
options      = configparser.ConfigParser()
options.read(options_file)

# SETTINGS
project     = options.get('Settings','project')
batch       = options.get('Settings','batch')
seed        = options.get('Settings','seed')
sample_file = options.get('Settings','sample_file')
genome_file = options.get('Settings','genome_file')
num_of_sims = options.getint('Settings','num_of_sims')

sims = range(1,num_of_sims+1)

# localrules: sim, getparams, clean, clean_project, clean_batch, test

# run simulations
rule sim:
    input: forwsim_trees = expand('results/{p}/{b}/forwsim_{s}.trees',p=project,b=batch,s=sims)

# simulation with msprime
rule forwsim:
    input:
        script='scripts/forwsim.slim',
        slim_options='results/{p}/{b}/sim_{s}.eidos',
        coalsim_trees='results/{p}/{b}/coalsim_{s}.trees'
    output:
        forwsim_trees='results/{p}/{b}/forwsim_{s}.trees'
    shell:
        'slim -d "option_file=\'{input.slim_options}\'" {input.script}'

# simulation with msprime
rule coalsim:
    input:
        script='scripts/coalsim.py',
        sim_ini='results/{p}/{b}/sim_{s}.ini'
    output:
        coalsim_trees='results/{p}/{b}/coalsim_{s}.trees'
    shell:
        'python {input.script} {options_file} {input.sim_ini}'


# read parameters and sample from priors
rule getparams:
    input:
        script='scripts/getparams.R'
    output:
        slim_options=expand('results/{p}/{b}/sim_{s}.eidos',p=project,b=batch,s=sims),
        sim_ini=expand('results/{p}/{b}/sim_{s}.ini',p=project,b=batch,s=sims) 
    shell:
        'Rscript {input.script} {options_file}'

# delete all files in batch directory
rule clean_batch:
    shell: 'rm -rf results/'+project+'/'+str(batch)

# delete all files in project directory
rule clean_project:
    shell: 'rm -rf results/'+project

# delete all files in results directory
rule clean_all:
    shell: 'rm -rf results/'

# run tests
rule test:
    shell:
        '''
        echo "----------- RUNNING TESTTHAT (R) ---------------\n"
        Rscript -e "library(testthat); test_file(\'tests/testthat/test_timeadapt.R\')"
        echo "\n----------- RUNNING PYTEST (PYTHON) ---------------\n"
        pytest -v scripts/timeadapt.py
        '''



