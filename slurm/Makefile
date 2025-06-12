GCC = mpicc
CFLAGS = -O3 -I../vertexSet

VSFOLDER = ../vertexSet/

target = $(VSFOLDER)main
SOURCES = $(VSFOLDER)vertexSet.c $(VSFOLDER)vertexSetIterator.c
HEADERS = $(VSFOLDER)vertexSet.h $(VSFOLDER)bitmap_custom.h $(VSFOLDER)vertexSetIterator.h
 
build = ../build

test_it: 
	mpirun -np 1 ./$(build)/test.o 0

test_red_exactHalfing: 
	sbatch test_exactly.job

test_red_approxHalfing: 
	sbatch test_approx.job

perf_red_exactHalfing:
	sbatch performance_exactly.job

perf_red_approxHalfing: 
	sbatch performance_approx.job

copy_to_cluster:
	scp -r ../vertexSet/ hydra:~/.
	scp -r ../slurm hydra:~/.

copy_from_cluster:
	scp -r hydra:~/data/ ../
