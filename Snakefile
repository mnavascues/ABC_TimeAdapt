# set options (default options are for integration test)
project_name = "test"
batch = 1
seed = 1234567890
sample_file = "tests/sample_info_test.txt"
genome_file = "tests/genome_info_test.txt"
generations_forward = 400
periods_forward = 8
gen_len_prior_sh1 = 2
gen_len_prior_sh2 = 1.465967
gen_len_prior_min = 26
gen_len_prior_max = 30
num_of_sims = 3
pop_size_prior_min = 10
pop_size_prior_max = 200
mut_rate_prior_mean = "0.00000005"
mut_rate_prior_sd = 0.5

sims = range(1,num_of_sims+1)

# read parameters and sample from priors
rule params:
    input:
        simR='R/simulations.R'
    output:
        slim_command=expand('results/{p}/{b}/slim_{s}.sh',p=project_name,b=batch,s=sims) ,
        pyslim_command=expand('results/{p}/{b}/pyslim_{s}.sh',p=project_name,b=batch,s=sims)  
        #slim_command=temp(expand('results/{p}/{b}/slim_{s}.sh',p=project_name,b=batch,s=sims)) ,
        #pyslim_command=temp(expand('results/{p}/{b}/pyslim_{s}.sh',p=project_name,b=batch,s=sims))  
    shell: 'Rscript {input.simR}             \
                    -q TRUE                  \
                    -d {seed}                \
                    -p {project_name}        \
                    -b {batch}               \
                    -i {sample_file}         \
                    -g {genome_file}         \
                    -f {generations_forward} \
                    -w {periods_forward}     \
                    -l {gen_len_prior_sh1} {gen_len_prior_sh2} {gen_len_prior_min} {gen_len_prior_max}\
                    -s {num_of_sims}         \
                    -n {pop_size_prior_min} {pop_size_prior_max}\
                    -u {mut_rate_prior_mean} {mut_rate_prior_sd}'

# simulate with SLiM
rule slim:
    input:
        'results/{p}/{b}/slim_{s}.sh'
    output:
        'results/{p}/{b}/slim_{s}.tree'
        #temp('results/{p}/{b}/slim_{s}.tree')
    threads: 4
    shell:
        'bash {input}'

# simulate with msprime (via pyslim) aNd calculate summary stats
rule pyslim:
    input:
        pyslim_command='results/{p}/{b}/pyslim_{s}.sh',
        slim_trees='results/{p}/{b}/slim_{s}.tree'
    output:
        'results/{p}/{b}/stats_{s}.txt'
    threads: 4
    shell:
        'bash {input.pyslim_command}'


rule sim:
    input: sim_stats = expand('results/{p}/{b}/stats_{s}.txt',p=project_name,b=batch,s=sims)
   

# run tests
rule test:
    shell: 'Rscript -e "library(testthat);  test_file(\'tests/testthat/test_myfunctions.R\',reporter=\'fail\')"'
