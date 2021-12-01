#!/bin/bash

: "
File: run.sh
Author: J.A. Stevens
Email: j.a.stevens@rug.nl
Github: https://github.com/jan-stevens
Description:
            This script performs free energy of solvation calculation for
            a given set of bead types / pH values.
"

# Useful variables
basedir=$(pwd)
GMX="gmx"
version="gromacs-2021.3"
threads=6

# Load Gromacs if available
if ! command -v /usr/local/$version/bin/GMXRC &> /dev/null
then
    echo "$version could not be found"
    exit 0
else
    echo "Gromacs is found"
    source /usr/local/$version/bin/GMXRC
fi

# Load alchemical_analysis
# Ref: https://github.com/MobleyLab/alchemical-analysis
if ! command -v alchemical_analysis &> /dev/null
then
    echo "Alchemical_analysis could not be found"
    echo "The simulation will continue by using the Gromacs implementation of thermal integration"
else
    echo "Alchemical_analysis is found"
fi

# List of bead types / pH values for which to perform FE calculations
list=( C1-3.0 C1-7.5 C4-3.0 C4-7.5 C4h-3.0 C4h-7.5 N4a-3.0 N4a-7.5 SC2-3.0 SC2-7.5 SC3-3.0
   SC3-7.5 SC5-3.0 SC5-7.5 SC6-3.0 SC6-7.5 SN4a-3.0 SN4a-7.5 SX2-3.0 SX2-7.5 TC2-3.0 TC2-7.5
   TC4-3.0 TC4-7.5 TC5-3.0 TC5-7.5 TN3d-3.0 TN3d-7.5 X2-3.0 X2-7.5
)
# Note: Easily parallelizable by splitting this list and running on multiple machines.

# Create output directory
mkdir -p Output ; cd Output

# Iterate over all the bead types / pH values for which the FE is calculated.
for run in ${list[@]}
do
        mkdir -p $run ; cd $run
        start=`date +%s.%N`
        rundir=$(pwd)
        bead=${run::-4}

    {
        # Prepare the correct input .top and .gro from templates
        sed -e "s/MOLE/$bead/g" -e "s/VALUE/${run: -3}/" "$basedir/System/system.top" > "system.top"
        sed "s/MOLE/$(printf '%+4s' "$bead")/g" "../../System/start.gro" > "start.gro"

        # Loop over the \lambda values used for the FE calculation and prepare / perform MD
        # simulation that samples the respective \lambda value
        for state in {0..29}
        do
           mkdir -p $state
           cd $state
           mkdir -p Minimization
           cd Minimization

           sed -e "s/sedstate/$state/" -e "s/MOLE/${bead}B/" $basedir/mdp_files/min.mdp > min.mdp
           $GMX grompp -f min.mdp -p $rundir/system.top -c $rundir/start.gro -o min.tpr -maxwarn 10
           $GMX mdrun -nt $threads -rdd 1.4 -v -deffnm min

           rm -rf \#*
           cd ../
           mkdir -p NVT
           cd NVT

           sed -e "s/sedstate/$state/" -e "s/MOLE/${bead}B/" $basedir/mdp_files/nvt.mdp > nvt.mdp
           $GMX grompp -f nvt.mdp -p $rundir/system.top -c ../Minimization/min.gro -maxwarn 10 -o nvt.tpr
           $GMX mdrun -nt $threads -pin on -nb gpu -rdd 1.4 -v -deffnm nvt

           rm -rf \#*
           cd ../
           mkdir -p NPT
           cd NPT

           sed -e "s/sedstate/$state/" -e "s/MOLE/${bead}B/" $basedir/mdp_files/npt.mdp > npt.mdp
           $GMX grompp -f npt.mdp -p $rundir/system.top -c ../NVT/nvt.gro -maxwarn 10 -o npt.tpr
           $GMX mdrun -nt $threads -nb gpu -pin on -rdd 1.4 -v -deffnm npt

           rm -rf \#*
           cd ../
           mkdir -p Dynamics
           cd Dynamics

           sed -e "s/sedstate/$state/" -e "s/MOLE/${bead}B/" $basedir/mdp_files/md.mdp > md.mdp
           gmx grompp -f md.mdp -p $rundir/system.top -c ../NPT/npt.gro -maxwarn 10 -o md.tpr
           gmx mdrun -v -nt $threads -nb gpu -pin on -rdd 1.4 -dlb yes -cpi state.cpt -deffnm md
           cp md.xvg ../../dhdl.$state.xvg

           rm -rf \#*
           cd ../../
        done

        # Perform final FE calculation using Alchemical_analysis package and Gromacs implementation
        alchemical_analysis -p dhdl -s 500
        gmx bar -f dhdl.*.xvg -o -oi -oh | grep total > final.dat
        rm -rf \#*

    # Print the total runtime of one FE calculation
    end=`date +%s.%N`
    runtime=$( echo "$end - $start" | bc -l )
    echo "
          ####################################################################\n
          ###   The Calculations for $bead took $runtime to complete.   ###\n
          ####################################################################\n"
    } 2>&1 | tee run.dmp # These brackets log the output of the FE calculation to a dumpfile

    cd ..
done
