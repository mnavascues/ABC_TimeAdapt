# simulate from prior using SLiM
rule simulation:
    input:
        simR='src/simulation.R'
    shell: 'Rscript {input.simR} -q F'
