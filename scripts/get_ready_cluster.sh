# Run on UPPMAX
###############

cd timeadapt
module load conda
source conda_init.sh
conda activate timeadapt
module load bioinfo-tools
module load snakemake
module load SLiM
# nohup snakemake sim --profile slurm &
jobinfo -u mnavas
