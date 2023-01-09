# PP_N_Body - Parallel Barnes Hut Solution to N Body Problem

MPI assignment for the Parallel Programming Practical (PPP)

Content:

- nbody		Contains sequential version of the N-body algorithm
			Create parallel version of the N-body algorithm there 
			and name it nbody-par.

- docs		Place there your report.

- bin		Contains nbody-sanity-check (comparison with expected output)
			Do not submit programs which do not go through sanity-checks successfully.

DAS-5 MPI Run:
prun -v -1 -np 4 -script ${PRUN_ETC}/prun-openmpi ./nbody-par

## Cost Zone
