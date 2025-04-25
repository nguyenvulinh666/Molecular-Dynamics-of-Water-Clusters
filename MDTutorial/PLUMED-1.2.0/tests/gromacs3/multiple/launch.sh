#!/bin/bash

EXE=/share/home33/carlo/Templates/source/gromacs-3.3.3/src/kernel/mdrun

cd w1/
${EXE} > gromacs.out.1.save &
cd ../w2
${EXE} > gromacs.out.2.save 
cd ../
