#!/bin/bash

#SBATCH -N 1
#SBATCH -n 1
#SBATCH --ntasks-per-node=1
#SBATCH --ntasks-per-core=1
#SBATCH -J P4mm_pola
#SBATCH --time=00:10:00
#SBATCH --mail-user=julien.gosteau@cemes.fr
#SBATCH --mem=180G
#SBATCH --acctg-freq=task=1


#Positionnement de l'environnement
module purge
module load intel/18.2
module load intelmpi/18.2
#var environnements
export OMP_NUM_THREADS=1
export MALLOC_MMAP_MAX_=0
export MALLOC_TRIM_THRESHOLD_=-1
export FORT_BUFFERED=true
export decfort_dump_flag=true
export I_MPI_DAPL_SCALABLE_PROGRESS=1
export I_MPI_DAPL_TRANSLATION_CACHE=1
export MKL_CBWR=Auto

#répertoire d'installation des binaires VASP 5.4.4
#export INSTALL_DIR=/usr/local/vasp/5.4.4/cpu/intelmpi/bin
export path=`pwd`
#export INSTALL_DIR=~

#création du répertoire de travail
#cd /tmpdir/${USER}
#WORKDIR=${SLURM_JOBID}
#mkdir ${WORKDIR}
cores=$SLURM_NTASKS
#On se positionne dans le répertoire de travail
#cd ${WORKDIR}

#Copie des fichiers input VASP
## ATTENTION ATTENTION ATTENTION ##
## ATTENTION ATTENTION ATTENTION ##
## ATTENTION ATTENTION ATTENTION ##
# La ligne ci-dessous est à adapter à votre cas
#cp -r /chemin_input_vasp/ .

## Definir le nom du binaire vasp à utiliser vasp_std vasp_gam vasp_ncl
#export binaire=vasp_relax_z123
#export binaire=vasp_std
#
#cp  ${INSTALL_DIR}/${binaire} .
export INTALL_DIR="/gpfs/work/p1229/gosteau/PROGRAMME_python/Programme/create_bands/PROCAR_split/"
cp $INSTALL_DIR/split.exe $path/.
#commande de lancement du binaire MPI (MPI pur, non multi-threadé)
time srun $(placement 1 1 ) $path/split.exe $path > out

#Récupération des données
#cp * /chemin_output_vasp/.

jobinfo $SLURM_JOBID
infoincidentjob
