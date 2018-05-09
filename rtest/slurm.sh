#!/bin/bash

#SBATCH --job-name="tarea_en_r"
#SBATCH --nodes=4
#SBATCH --ntasks=10
#SBATCH --time=00:00:02
#SBATCH --mem=200MB
#SBATCH --error="tarea_en_r-%j.err"
#SBATCH --output="tarea_en_r-%j.out"
#SBATCH --open-mode=append


# Testear si la configuración es correcta. Suele ser más útil y cómodo
# especificarlo en la línea de comandos, es decir:
#
#  $ sbatch --test-only script.sh
#
##SBATCH --test-only

echo "trabajo \"${SLURM_JOB_NAME}\""
echo "    id: ${SLURM_JOB_ID}"
echo "    partición: ${SLURM_JOB_PARTITION}"
echo "    nodos: ${SLURM_JOB_NODELIST}"
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

# PROGRAMA A CORRER
#
#srun /PATH/AL/PROGRAMA

srun Rscript --vanilla compute.R data.in salida.out

echo "
--------------------------------------------------------------------------------
"

date +"fin %F - %T"
