#BSUB -P BIP240
#BSUB -J a0
#BSUB -W 2:00 
#BSUB -nnodes 105
#BSUB -alloc_flags "maximizegpfs smt1 gpumps gpudefault"
### End BSUB Options and begin shell commands

mkdir -p plots
mkdir -p data
mkdir -p dump
mkdir -p e2e_pos
mkdir -p com_pos
mkdir -p link_pos
mkdir -p info

module load gcc
module load fftw
module load cuda

module load miniforge3/24.3.0-0

python3 live.polymer.py

#number of resources
#nres=`echo "$LSB_DJOB_NUMPROC-1" | bc`
nres=1050

LMP_PATH=/gpfs/alpine2/proj-shared/bip240/lammps_build/lammps-2Aug2023/build_summit

echo "variable xx uloop $nres" > in.variables
echo >> in.variables
printf  "variable vseed universe " >> in.variables
for (( x=0; x<$nres; x++ ))
do
    printf "%i " $RANDOM >> in.variables

done

jsrun -n $nres -a4 -c4 -g0 -r10  $LMP_PATH/lmp_summit -in in_summit.lammps -p `echo $nres`x4
