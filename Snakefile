project_name = "test"
batch = 1
seed = 1234567890

# simulate from prior using SLiM
rule simulation:
    input:
        simR='src/simulation.R'
    shell: 'Rscript {input.simR} -q F -d {seed} -p {project_name} -b {batch}'
