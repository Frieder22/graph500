GCC = mpicc
CFLAGS = -O3 -I ../vertexSet -I cJSON

VSFOLDER = ../vertexSet/

target = $(VSFOLDER)main
SOURCES = $(VSFOLDER)vertexSet.c $(VSFOLDER)vertexSetIterator.c cJSON/cJSON.c
HEADERS = $(VSFOLDER)vertexSet.h $(VSFOLDER)bitmap_custom.h $(VSFOLDER)vertexSetIterator.h $(VSFOLDER)testStructure.h
 
build = ../build

perf: executeExperimentRuns.c $(SOURCES) $(HEADERS)
	$(GCC) $(CFLAGS) -o $(build)/executeExperimentRuns.o executeExperimentRuns.c $(SOURCES) 

run_local: perf
	mpirun -np 4 $(build)/executeExperimentRuns.o runs.json

init_cluster:
	ssh -t hydra 'bash --login -i -c "mkdir -p build && mkdir -p data && mkdir -p performance"'
	scp -r cJSON hydra:~/performance/.
	scp executeExperimentRuns.c hydra:~/performance/.
	scp Makefile hydra:~/performance/.

copy_to_cluster:
	scp -r ../vertexSet/ hydra:~/.
	scp run.job hydra:~/performance/
	scp runs.json hydra:~/performance/

run: 
	python setupExperiment.py
	make copy_to_cluster
	ssh -t hydra 'bash --login -i -c "cd performance && sbatch run.job"'


.PHONY: clean_Data

clean_Data:
	echo Not implemented