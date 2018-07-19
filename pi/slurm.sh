#!/bin/bash
#
#SBATCH --job-name=simulate_pi
#SBATCH --partition=batch
#SBATCH --output=%A-%a.out
#SBATCH --error=%A-%a.err
#SBATCH --ntasks=1
#SBATCH --time=00:05
#SBATCH --mem-per-cpu=4096MB
#SBATCH --array=1-100

echo "trabajo \"${SLURM_JOB_NAME}\""
echo "    id: ${SLURM_JOB_ID}"
echo "    array_id: ${SLURM_ARRAY_TASK_ID}"
echo "    partición: ${SLURM_JOB_PARTITION}"
echo "    nodos: ${SLURM_JOB_NODELIST}"
echo
date +"inicio %F - %T"

echo "
--------------------------------------------------------------------------------
"

# INICIO VARIABLES IMPORTANTES
#
# NO TOCAR. No se aceptan reclamos en caso de modificar estas líneas. Deberán
# incluirlas siempre, hasta próximo aviso.
#
[ -r /etc/profile.d/odin-users.sh ] && . /etc/profile.d/odin-users.sh
#
# FIN VARIABLES IMPORTANTES

srun Rscript simulate.R $SLURM_ARRAY_TASK_ID

echo "
--------------------------------------------------------------------------------
"

date +"fin %F - %T"


