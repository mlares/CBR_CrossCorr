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


## nodes requested
#SBATCH -N 2

## tasks requested
#SBATCH -n 6

## cores requested
#SBATCH -c 3

## Memory in Mb requested
#SBATCH --mem=10

## STDOUT
#SBATCH -o submit_python_jobs.out

## STDOUT
#SBATCH -e submit_python_jobs.err



## Procesos a largar.
## Es OpenMP, o sea que un proceso en un nodo y muchos hilos.
#SBATCH --ntasks=10

## Hilos por proceso
## Poner el mismo valor acá que en OMP_NUM_THREADS/MKL_NUM_THREADS
#SBATCH --cpus-per-task=56

## Tiempo de ejecucion. Formato dias-horas:minutos.
#SBATCH --time 3-0:00

## Script que se ejecuta al arrancar el trabajo

## Cargar el entorno del usuario incluyendo la funcionalidad de modules
## No tocar
. /etc/profile

## Configurar OpenMP/MKL/etc con la cantidad de cores detectada.
# export OMP_NUM_THREADS=56
# export MKL_NUM_THREADS=56

## Load modules 
# module load clemente etc...

conda activate

## Launch program

srun python ../codes/recon/eboss_clustering/python/zelrec_cluster_4.py



