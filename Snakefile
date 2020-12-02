project_name = "test"
batch = 1
seed = 1234567890
sample_file = "data/sample_info_test.txt"
genome_file = "data/genome_info_test.txt"
generations_forward = 400
periods_forward = 8

# simulate from prior using SLiM
rule simulation:
    input:
        simR='R/simulation.R'
    shell: 'Rscript {input.simR} -q F -d {seed} -p {project_name} -b {batch} -i {sample_file} -g {genome_file} -f {generations_forward} -w {periods_forward}'

# run tests
rule test:
    shell: 'Rscript -e "library(testthat);  test_file(\'tests/testthat/test_myfunctions.R\')"'
