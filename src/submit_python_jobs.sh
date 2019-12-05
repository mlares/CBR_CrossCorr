#!/bin/bash

# SLURM script for: CLEMENTE
 
## Las líneas #SBATCH configuran los recursos de la tarea
## (aunque parezcan estar comentadas)

# More info:
# http://homeowmorphism.com/articles/17/Python-Slurm-Cluster-Five-Minutes


## Nombre de la tarea
#SBATCH --job-name=CMB_corr

## Cola de trabajos a la cual enviar.
#SBATCH --partition=small

## tasks requested
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20

## STDOUT
#SBATCH -o submit_python_jobs.out

## STDOUT
#SBATCH -e submit_python_jobs.err

## Tiempo de ejecucion. Formato dias-horas:minutos.
#SBATCH --time 0-1:00

## Script que se ejecuta al arrancar el trabajo

## Cargar el entorno del usuario incluyendo la funcionalidad de modules
## No tocar
. /etc/profile

source /home/${USER}/.bashrc
conda activate

## Launch program

srun python /home/mlares/projects/E.26/CBR_CrossCorr/src/run_correlation.py ../set/config_small.ini


## launch script
## $>sbatch submit_python_jobs.sh