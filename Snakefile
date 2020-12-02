project_name = "test"
batch = 1
seed = 1234567890
sample_file = "data/sample_info_test.txt"
genome_file = "data/genome_info_test.txt"
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
        simR='R/params.R'
    output:
        slim_command=expand('results/{project_name}/{batch}/slim_{sims}.sh',
                             project_name=project_name,
                             batch=batch,
                             sims=sims)
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
        'results/{p}/{b}/slim{s}.tree'
    threads: 4
    shell:
        'bash {input}'

rule all:
    input:
       slim_trees = expand('results/{p}/{b}/slim{s}.tree',
                           p=project_name,
                           b=batch,
                           s=sims)
   

# run tests
rule test:
    shell: 'Rscript -e "library(testthat);  test_file(\'tests/testthat/test_myfunctions.R\')"'
