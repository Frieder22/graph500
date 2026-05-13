CFLAGS = -Drestrict=__restrict__ -O3 -DGRAPH_GENERATOR_MPI -DREUSE_CSR_FOR_VALIDATION -I../aml -fcommon -I../vertexSet
# added -fcommon to avoid dublicated variable conflicts
LDFLAGS = -lpthread 
MPICC = mpicc

all: graph500_reference_bfs_sssp graph500_reference_bfs 
#graph500_custom_bfs graph500_custom_bfs_sssp

GENERATOR_SOURCES = ../generator/graph_generator.c ../generator/make_graph.c ../generator/splittable_mrg.c ../generator/utils.c
VERTEXSET_SOURCES = ../vertexSet/vertexSet.c ../vertexSet/vertexSetIterator.c
SOURCES = main.c utils.c validate.c ../aml/aml.c csr_custom.c
HEADERS = common.h csr_reference.h bitmap_reference.h csr_custom.h

path:
	echo $(PWD)/../vertexSet

graph500_reference_bfs_sssp: bfs_reference.c $(SOURCES) $(HEADERS) $(GENERATOR_SOURCES) csr_reference.c sssp_reference.c
	$(MPICC) $(CFLAGS) -DSSSP $(LDFLAGS) -o graph500_reference_bfs_sssp bfs_reference.c sssp_reference.c csr_reference.c $(SOURCES) $(GENERATOR_SOURCES) -lm

graph500_reference_bfs: bfs_reference.c $(SOURCES) $(HEADERS) $(GENERATOR_SOURCES) csr_reference.c
	$(MPICC) $(CFLAGS) $(LDFLAGS) -o graph500_reference_bfs bfs_reference.c csr_reference.c $(SOURCES) $(GENERATOR_SOURCES) -lm

graph500_custom_bfs: bfs_custom.c $(SOURCES) $(HEADERS) $(GENERATOR_SOURCES) $(VERTEXSET_SOURCES) csr_reference.c 
	$(MPICC) $(CFLAGS) $(LDFLAGS) -o graph500_custom_bfs bfs_custom.c csr_reference.c $(SOURCES) $(GENERATOR_SOURCES) $(VERTEXSET_SOURCES) -lm

graph500_custom_bfs_sssp: bfs_custom.c sssp_custom.c $(SOURCES) $(HEADERS) $(GENERATOR_SOURCES) csr_reference.c
	$(MPICC) $(CFLAGS) $(LDFLAGS) -o graph500_custom_bfs_sssp bfs_custom.c sssp_custom.c csr_reference.c $(SOURCES) $(GENERATOR_SOURCES) -lm

run: graph500_custom_bfs
	mpirun -np 4 graph500_custom_bfs 12

run_ref: graph500_reference_bfs
	mpirun -np 4 graph500_reference_bfs 7

clean:
	-rm -f graph500_*
