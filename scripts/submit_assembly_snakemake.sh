#!/bin/bash

#-----------------------------------------------------
#   Submitter: Zexuan Zhao
#   Email: zzhao127@umd.edu
#-----------------------------------------------------

# Slurm sbatch parameters:
# 1. Number of nodes (usually 1)
#SBATCH --ntasks=1
# 2. Number of cores
#SBATCH --cpus-per-task=40
# 3. Running time:
# 	Format:
# 		M (M minutes)
# 		M:S (M minutes, S seconds)
# 		H:M:S (H hours, M minutes, S seconds)
# 		D-H (D days, H hours)
# 		D-H:M (D days, H hours, M minutes)
# 		D-H:M:S (D days, H hours, M minutes, S seconds)
#SBATCH -t 7-0:0:0
# 4. Memory in MB
#SBATCH --mem=102400
# 5. Running mode:
# 	Allow other jobs running on the node assigned to your job: --oversubscribe
#	Don't allow: --exclusive
#SBATCH --oversubscribe
#6. Select partition to run
#	Don't change the setting below. If you want to run on debug mode,
#	use --partition=debug when calling `sbatch this_script`.
#SBATCH --partition=standard

# Avoid inherit variables from the process called this script.
#SBATCH --export=NONE
#

# Section to ensure we have the "module" command defined
unalias tap >& /dev/null
if [ -f ~/.bash_profile ]; then
	source ~/.bash_profile
elif [ -f ~/.profile ]; then
	source ~/.profile
fi

# Set SLURM_EXPORT_ENV to ALL.  This prevents the --export=NONE flag
# from being passed to mpirun/srun/etc, which can cause issues.
# We want the environment of the job script to be passed to all
# tasks/processes of the job
export SLURM_EXPORT_ENV=ALL

# Module load section
# First clear our module list
module purge
# Then load modules. Check available modules here: https://hpcc.umd.edu/hpcc/help/software.html#module
module load samtools
module load fastqc
module load blast-plus
module load bwa
# Conda environment loading:
eval "$(conda shell.bash hook)"
conda activate LotmariaGenomeAssembly


# Section to output information identifying the job, etc.
echo "Slurm job ${SLURM_JOBID} running on"
hostname
echo "To run on ${SLURM_NTASKS} CPU cores across ${SLURM_JOB_NUM_NODES} nodes"
echo "All nodes: ${SLURM_JOB_NODELIST}"
date
pwd
echo "Loaded modules are:"
module list

############################################################
############################################################
# Main program                                             #
############################################################
############################################################

# Your code
snakemake --cores 40

# Finishing up
echo "Job finished with exit code $ECODE"
date

# Exit with the cached exit code
exit $ECODE
