nres=1

mkdir -p plots
mkdir -p data
mkdir -p dump
mkdir -p com_pos
mkdir -p link_pos
mkdir -p mon_pos
mkdir -p info

echo "variable xx uloop $nres" > in.variables
echo >> in.variables
printf  "variable vseed universe " >> in.variables
for (( x=0; x<$nres; x++ ))
do
    printf "%i " $RANDOM >> in.variables

done

time mpirun -n 16 lmp_mpi -in input.lammps > out.run

rm log.*
