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

## Work Flow

No Cost zone

1. Build Tree

2.Dispacth work to worker
work = 	compute_forces(world);
    	compute_velocities(world);
        compute_positions(world);
3.Worker compute
4.Collect work

1. Worker compute local force. Local astro bodies accumlates one force


prun -v -16 -np 8 -sge-script /cm/shared/package/reserve.slurm/etc/prun-openmpi nbody/nbody-par 10000 0 nbody.ppm 100

## Communicate overhead

## Load Imbalance
Root worker takes extra work to distribute work

