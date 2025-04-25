Copyright-computer modelling lab., University of Modena and Reggio Emilia.
Author- Dr. Pedone Alfonso January 2005

The cluster program is available free of charge to academic
institutions. Copies should be obtained from the author or through the
DL_POLY package.

No claim is made that the program is free from errors, the user is
responsible for checking the validity of their results.

Details of the algorithm and its application can be found into the
paper: "A tool for the prediction of crystalline phases obtained
from controlled crystallization of glasses".
G.Lusvardi,G.Malavasi,L.Menabue,M.C.Menziani,A.Pedone and U.Segre
J.Phys.Chem.B, 109, 21586-21592, 2005 web Release Date: 22-oct-2005
DOI: 10.1021/jp0546857

Aim of the program:

This is a FORTRAN interactive program able to search for clusters of
ions (crystallization nucleus) with stoichiometry similar to a crystal
phase that you are looking for into the structure of glasses modelled
by MD simulation.

The source data is assumed to be formatted and compatible with the
CONFIG and REVCON file written by the subroutine REVIVE.f of the
DL_POLY package.

An explanation of the CONFIG file format is available on line to the
following link:
http://www.cse.clrc.ac.uk/msi/software/DL_POLY/MANUALS/USRMAN3/node113.html

However, an input file for the cluster program is attached with the program.

After you run the executable, the program asks you the following questions:

1.Name of the REVCON or CONFIG input file you want to analyze.
2.Name of the output file 
3.The stoichiometry of the crystal phase you are looking for.

For example, if you are searching for a crystallization nucleus of the
phase NaCaPO4, in a multi-component glass made of O,Si,Na,Ca and P
ions, the program will ask you:

>>Insert number of O ions
>>4.0           (this should be your answer)
>>Insert number of Si ions
>>0.0           (the crystal phase does not contain Si atoms)
>>Insert number of Na ions
>>1.0
... and so on.


OUTPUT file:

Record 1: Number of the atom that lies at the centre of the spherical cluster.
Record 2-3: Cartesian coordinates of the central atom
Record 4: Radius of the cluster
Record 5-n:In this records the composition and minimum stoichiometry of the cluster with smaller displacement found is listed.
Record n+1:displacement of the cluster

The subsequent records contain the same information for bigger clusters.  

