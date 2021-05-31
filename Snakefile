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

localrules: sim, getparams, clean, clean_project, clean_batch, test

# run simulations
rule sim:
    #input: sim_stats = expand('results/{p}/{b}/stats_{s}.txt',p=project,b=batch,s=sims)
    input: coalsim_trees = expand('results/{p}/{b}/coalsim_{s}.trees',p=project,b=batch,s=sims)

# simulation with msprime
rule coalsim:
    input:
        script='scripts/coalsim.py',
        sim_ini='results/{p}/{b}/sim_{s}.ini'
    output:
        'results/{p}/{b}/coalsim_{s}.trees'
    shell:
        'python {input.script} {options_file} {input.sim_ini}'


# read parameters and sample from priors
rule getparams:
    input:
        script='scripts/getparams.R'
    output:
        slim_command=expand('results/{p}/{b}/slim_{s}.sh',p=project,b=batch,s=sims) ,
        sim_ini=expand('results/{p}/{b}/sim_{s}.ini',p=project,b=batch,s=sims) 
        #pyslim_command=expand('results/{p}/{b}/pyslim_{s}.sh',p=project,b=batch,s=sims)  
        #slim_command=temp(expand('results/{p}/{b}/slim_{s}.sh',p=project,b=batch,s=sims)) ,
        #pyslim_command=temp(expand('results/{p}/{b}/pyslim_{s}.sh',p=project,b=batch,s=sims))  
    shell:
        'Rscript {input.script} {options_file}'

# delete all files in batch directory
rule clean_batch:
    shell: 'rm -rf results/'+project+'/'+str(batch)

# delete all files in project directory
rule clean_project:
    shell: 'rm -rf results/'+project

# delete all files in project directory
rule clean:
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



