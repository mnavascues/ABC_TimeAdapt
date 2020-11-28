# simulate from prior using SLiM
rule simulation:
    input:
        simR='src/simulation.R'
    shell: 'Rscript {input.simR} -q F -d 1234567890 -p "test" -b 1'
