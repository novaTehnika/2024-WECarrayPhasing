#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=16gb
#SBATCH -t 24:00:00
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
D_m_base = ${DISP}*1e-6/(2*pi); \
display(['D_m_base = ',num2str(D_m_base),' m^3/rad']); \
NumWECs = ${NUMWECS}; \
display(['NumWECs = ',num2str(NumWECs)]); \
study_hydElecPTO_arrayTpDesign"

# Commands to use
# sbatch --export=SS=1,DISP=316,NUMWECS=2 --array=1-2000 ~/2024-WECarrayPhasing/study_hydElecPTO_arrayTeDesign.sh
# dos2unix  study_hydElecPTO_arrayTeDesign.sh

