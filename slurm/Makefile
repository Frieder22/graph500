CONFIG = default.config
include ${CONFIG}

GCC = mpicc
CFLAGS = -O3 -I../vertexSet

VSFOLDER = ../vertexSet/

target = $(VSFOLDER)main
SOURCES = $(VSFOLDER)vertexSet.c $(VSFOLDER)vertexSetIterator.c
HEADERS = $(VSFOLDER)vertexSet.h $(VSFOLDER)bitmap_custom.h $(VSFOLDER)vertexSetIterator.h $(VSFOLDER)testStructure.h
 
build = ../build

test: $(VSFOLDER)test.c $(SOURCES) $(HEADERS)
	$(GCC) $(CFLAGS) -o $(build)/test.o $(VSFOLDER)test.c $(SOURCES)

perf: $(VSFOLDER)performance.c $(SOURCES) $(HEADERS)
	$(GCC) $(CFLAGS) -o $(build)/perf.o $(VSFOLDER)performance.c $(SOURCES)

submit_slurm: test perf
	sbatch run.job

copy_to_cluster:
	scp -r ../vertexSet/ hydra:~/.
	scp -r ../slurm hydra:~/.

init_cluster:
	ssh -t hydra 'bash --login -i -c "mkdir -p build && mkdir data/performanceData/intelmpi && mkdir data/performanceData/openmpi && mkdir data/performanceData/mpich"'

run: copy_to_cluster
	ssh -t hydra 'bash --login -i -c "${MPI_LOAD} && cd slurm && make submit_slurm"'

get_Data:
	scp -r hydra:~/data/ ../

.PHONY: clean_Data

clean_Data:
	rm ../data/performanceData/*