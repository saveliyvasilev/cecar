#!/bin/bash

# Ejecuta un programa con uso de MPI, corriendo 16 tareas, distribuidas en 4
# nodos. Esto no asegura que se corran 4 tareas por nodo. Para eso, leer la
# documentación sobre la opción `--ntasks-per-node'.

# Correr con `sbatch /ruta/a/este/script

#SBATCH --job-name="test_mpi"
#SBATCH --nodes=8
#SBATCH --ntasks=16
#SBATCH --workdir=/home/svassiliev/cecar/mpi_test
#SBATCH --error="ejemplo-mpi-%j.err"
#SBATCH --output="ejemplo-mpi-%j.out"
##SBATCH --partition=freeriders
#SBATCH --time=10:00

echo "trabajo \"${SLURM_JOB_NAME}\""
echo "    id: ${SLURM_JOB_ID}"
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

# El path al programa es `/home/USUARIO/ejemplos/programa-mpi`. Como más arriba
# seteamos el directorio de trabajo a `/home/USUARIO/ejemplos`, el programa se
# encuentra en el directorio de trabajo, por lo que basta poner
# `./programa-mpi`.
srun ./a.out

echo "
--------------------------------------------------------------------------------
"

date +"fin %F - %T"
