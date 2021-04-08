import configparser

# READ OPTIONS FILE
options_file = 'testproject.ini'
options = configparser.ConfigParser()
options.read(options_file)

# SETTINGS
project = options.get('Settings','project')
batch = options.get('Settings','batch')
seed = options.get('Settings','seed')
sample_file = options.get('Settings','sample_file')
genome_file = options.get('Settings','genome_file')
num_of_sims = options.getint('Settings','num_of_sims')

sims = range(1,num_of_sims+1)

localrules: sim, getparams, clean, clean_batch, test

# run simulations
rule sim:
    input: sim_stats = expand('results/{p}/{b}/stats_{s}.txt',p=project,b=batch,s=sims)

# read parameters and sample from priors
rule getparams:
    input:
        'R/simulations.R'
    output:
        slim_command=expand('results/{p}/{b}/slim_{s}.sh',p=project,b=batch,s=sims) ,
        sim_ini=expand('results/{p}/{b}/sim_{s}.ini',p=project,b=batch,s=sims) 
        #pyslim_command=expand('results/{p}/{b}/pyslim_{s}.sh',p=project,b=batch,s=sims)  
        #slim_command=temp(expand('results/{p}/{b}/slim_{s}.sh',p=project,b=batch,s=sims)) ,
        #pyslim_command=temp(expand('results/{p}/{b}/pyslim_{s}.sh',p=project,b=batch,s=sims))  
    shell:
        'Rscript {input} {options_file}'

# simulate with SLiM
rule slim:
    input:
        'results/{p}/{b}/slim_{s}.sh'
    output:
        'results/{p}/{b}/slim_{s}.tree'
        #temp('results/{p}/{b}/slim_{s}.tree')
    shell:
        'bash {input}'

# simulate with msprime (via pyslim) and calculate summary stats
rule pyslim:
    input:
        # pyslim_command='results/{p}/{b}/pyslim_{s}.sh',
        sim_ini='results/{p}/{b}/sim_{s}.ini',
        slim_trees='results/{p}/{b}/slim_{s}.tree'
    output:
        'results/{p}/{b}/stats_{s}.txt'
    shell:
        'python python/msprimeNstats.py {options_file} {input.sim_ini}'

# delete all files in batch directory
rule clean_batch:
    shell: 'rm -rf results/'+project+'/'+str(batch)

# delete all files in project directory
rule clean:
    shell: 'rm -rf results/'+project

# run tests
rule test:
    shell:
        '''
        Rscript -e "library(testthat); test_file(\'tests/testthat/test_myfunctions.R\')"
        echo "\n\n\n"
        pytest -v python/myfun.py
        '''
