#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=4gb
#SBATCH -t 8:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=simmo536@umn.edu
#SBATCH -p msismall
#SBATCH -o %A_%a.out
#SBATCH -e %A_%a.err

cd ~/2024-WECarrayPhasing
module load matlab
matlab -nodisplay -r \
"iVar = ${SLURM_ARRAY_TASK_ID}; \
display(['iVar = ',num2str(iVar)]); \
SS = ${SS}; \
display(['SS = ',num2str(SS)]); \
study_hydElecPTO_checkValves"

# Commands to use
# sbatch --export=SS=1 --array=1-50 ~/2024-WECarrayPhasing/study_hydElecPTO_checkValves.sh
# dos2unix  study_hydElecPTO_checkValves.sh

