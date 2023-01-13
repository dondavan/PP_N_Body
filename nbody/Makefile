CFLAGS	+= -Wall
CFLAGS	+= -O3
CFLAGS	+= -g2
export TMPDIR=/tmp

all: nbody-seq

nbody-seq: nbody-seq.c
	gcc $(CFLAGS) -o nbody-seq nbody-seq.c -lm

nbody-par: nbody-par.c
	mpicc $(CFLAGS) -o nbody-par nbody-par.c -lm

par-run:
	mpirun -np 4 ./nbody-par 32 2 ../nbody.ppm 1000000

seq-run:
	./nbody-seq 32 2 ../nbody.ppm 1000000

clean:
	rm -f *.o nbody-seq nbody-par *~ *core
