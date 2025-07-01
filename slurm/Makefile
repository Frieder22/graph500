GCC = mpicc
CFLAGS = -O3 -I../vertexSet

VSFOLDER = ../vertexSet/

target = $(VSFOLDER)main
SOURCES = $(VSFOLDER)vertexSet.c $(VSFOLDER)vertexSetIterator.c
HEADERS = $(VSFOLDER)vertexSet.h $(VSFOLDER)bitmap_custom.h $(VSFOLDER)vertexSetIterator.h
 
build = ../build

test: $(VSFOLDER)test.c $(SOURCES) $(HEADERS)
	$(GCC) $(CFLAGS) -o $(build)/test.o $(VSFOLDER)test.c $(SOURCES)

perf: $(VSFOLDER)performance.c
	$(GCC) $(CFLAGS) -o $(build)/perf.o $(VSFOLDER)performance.c $(SOURCES)

load_MPI:
	spack load openmpi@4.1.6

test_red_exactHalfing: test
	sbatch test_exactly.job

test_red_approxHalfing: test
	sbatch test_approx.job

test_red_dense: perf
	sbatch test_dense.job

perf_red_exactHalfing: perf
	sbatch performance_exactly.job

perf_red_approxHalfing: perf
	sbatch performance_approx.job

perf_red_dense: perf
	sbatch performance_dense.job

copy_to_cluster:
	scp -r ../vertexSet/ hydra:~/.
	scp -r ../slurm hydra:~/.

copy_from_cluster:
	scp -r hydra:~/data/ ../

.PHONY: clean_Data

clean_Data:
	rm ../data/performanceData/*