#!/bin/bash

WORKDIR=/work/sharma/afung/TornadoSim/ISUSim12
FOAMDIR=/work/sharma/afung/openfoam_v2312
MEM=128G

rm mesh.sbatch
rm prep.sbatch
rm run.sbatch
rm enableTurb.sbatch

cat << EOF > mesh.sbatch
#!/bin/bash
#SBATCH -A sharma
#SBATCH -J foamMesh
#SBATCH -D $WORKDIR
#SBATCH -N 1
#SBATCH --partition=nova
#SBATCH -n 8

#SBATCH --output=log_mesh
#SBATCH --time=8:00:00
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --mail-user=afung@iastate.edu
#SBATCH --mem=$MEM

cd $WORKDIR

module purge
module load openmpi
# module load openfoam-intel/6.0
source $FOAMDIR/etc/bashrc

rm -rf constant/polyMesh

blockMesh
surfaceFeatureExtract

rm system/decomposeParDict
cp system/decomposeParDict_mesh system/decomposeParDict

decomposePar

mpirun -np 8 snappyHexMesh -parallel

reconstructParMesh -latestTime -mergeTol 1E-6

rm -rf constant/polyMesh
mv 3/polyMesh constant
rmdir 3
rm -rf processor*

checkMesh

sbatch prep.sbatch

EOF

cat << EOF > prep.sbatch
#!/bin/bash
#SBATCH -A sharma
#SBATCH -J foamRun
#SBATCH -D $WORKDIR
#SBATCH -N 1
#SBATCH --partition=nova
#SBATCH -n 1
#SBATCH --mem=$MEM

#SBATCH --output=log_prep
#SBATCH --time=2:00:00
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --mail-user=afung@iastate.edu

cd $WORKDIR

module purge
module load openmpi

source $FOAMDIR/etc/bashrc

topoSet

rm 0/*
cp 0_MT/* 0

setExprFields

rm system/decomposeParDict
cp system/decomposeParDict_run system/decomposeParDict

sbatch run.sbatch

EOF

cat << EOF > run.sbatch
#!/bin/bash

#SBATCH -A sharma
#SBATCH -J foamRun
#SBATCH -D $WORKDIR
#SBATCH -N 1
#SBATCH --partition=nova
#SBATCH -n 64
#SBATCH --mem=$MEM

#SBATCH --output=log
#SBATCH --time=54:00:00
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --mail-user=afung@iastate.edu

cd $WORKDIR

module purge
module load openmpi

source $FOAMDIR/etc/bashrc

decomposePar

mpirun -np 64 porousSimpleFoam -parallel

EOF


echo "Enter r to run, p to prep and run, or m to mesh, prep, and run"
read opt

if [ "$opt" = "m" ]
then
	echo "Submitting FOAM Mesh Job"
	sbatch ./mesh.sbatch
fi

if [ "$opt" = "p" ]
then
	echo "Submitting FOAM Prep Job"
	sbatch ./prep.sbatch
fi

if [ "$opt" = "r" ]
then
	echo "Submitting FOAM Run Job"
	sbatch ./run.sbatch
fi

