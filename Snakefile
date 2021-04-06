import configparser

# READ OPTIONS FILE
options_file = 'testproject.ini'
options = configparser.ConfigParser()
options.read(options_file)

# SETTINGS
project_name = options.get('Settings','project_name')
batch = options.get('Settings','batch')
seed = options.get('Settings','seed')
sample_file = options.get('Settings','sample_file')
genome_file = options.get('Settings','genome_file')
num_of_sims = options.getint('Settings','num_of_sims')

sims = range(1,num_of_sims+1)

localrules: getparams

# run simulations
rule sim:
    input: sim_stats = expand('results/{p}/{b}/stats_{s}.txt',p=project_name,b=batch,s=sims)

# read parameters and sample from priors
rule getparams:
    input:
        'R/simulations.R'
    output:
        slim_command=expand('results/{p}/{b}/slim_{s}.sh',p=project_name,b=batch,s=sims) ,
        pyslim_command=expand('results/{p}/{b}/pyslim_{s}.sh',p=project_name,b=batch,s=sims)  
        #slim_command=temp(expand('results/{p}/{b}/slim_{s}.sh',p=project_name,b=batch,s=sims)) ,
        #pyslim_command=temp(expand('results/{p}/{b}/pyslim_{s}.sh',p=project_name,b=batch,s=sims))  
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

# simulate with msprime (via pyslim) aNd calculate summary stats
rule pyslim:
    input:
        pyslim_command='results/{p}/{b}/pyslim_{s}.sh',
        slim_trees='results/{p}/{b}/slim_{s}.tree'
    output:
        'results/{p}/{b}/stats_{s}.txt'
    shell:
        'bash {input.pyslim_command}'

# delete all files in project/batch directory
rule clean:
    shell: 'rm -rf results/'+project_name+'/'+str(batch)

# delete all files in project/batch directory
rule clean_project:
    shell: 'rm -rf results/'+project_name

# run tests
rule test:
    shell:
        '''
        Rscript -e "library(testthat); test_file(\'tests/testthat/test_myfunctions.R\')"
        echo "\n\n\n"
        pytest -v python/myfun.py
        '''
