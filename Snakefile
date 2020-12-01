project_name = "test"
batch = 1
seed = 1234567890
sample_file = "data/sample_info_test.txt"
genome_file = "data/genome_info_test.txt"

# simulate from prior using SLiM
rule simulation:
    input:
        simR='src/simulation.R'
    shell: 'Rscript {input.simR} -q F -d {seed} -p {project_name} -b {batch} -i {sample_file} -g {genome_file}'
