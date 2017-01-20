#!/bin/bash
#PBS -N generateZJJ
#PBS -o generateZJJ_out
#PBS -e generateZJJ_err
#PBS -l mem=400mb
#PBS -l cput=10000:00:00
#PBS -l walltime=10000:00:00
#PBS -m bea
#PBS -M edson.carquin@usm.cl
cd $PBS_O_WORKDIR

use root5
use gcc49

./bin/mg5_aMC generateZJJ.dat
