#include "common.h"
#include "csr_custom.h"
#include <assert.h>


void divideTuplegraph_divisible(tuple_graph* const tg){
    // only rank 0 holds correct edge list
	// distribute edges equally to the processes
	int64_t nlocaledges = tg->nglobaledges / size; // every process has correct nglobaledges
												   // we assume number of ranks devide the edges as a whole
	packed_edge* local_edges;
	local_edges = (packed_edge*) malloc(sizeof(packed_edge) * nlocaledges);
	MPI_Scatter(tg->edgememory,
				nlocaledges,
				packed_edge_mpi_type, // mpi_Type is already provided by framework
				local_edges,
				nlocaledges,
				packed_edge_mpi_type,
				0,  // rank 0 has all edges
				MPI_COMM_WORLD);
	tg->edgememory = local_edges;
	tg->nlocaledeges = nlocaledges;				
}

size_t* getVertexSpacing(tuple_graph* tg){
	size_t* vertexCount;
	vertexCount = (size_t*) malloc(tg->nglobalverts * sizeof(size_t));
	
	// initialize with 0
	for (size_t i = 0; i < tg->nglobalverts; i++) {
		vertexCount[i] = 0;
	}

	// count occurence of each vertex
	packed_edge edge;	
	for (size_t i = 0; i < tg->nlocaledeges; i++) {
		edge = tg->edgememory[i];
		vertexCount[edge.v0_low]++;
		vertexCount[edge.v1_low]++;
	}

	// perform scan to get starting positions in CSR
	for (size_t i = 1; i < tg->nglobalverts; i++) {
		vertexCount[i] += vertexCount[i-1];
	}

	assert(vertexCount[tg->nglobalverts - 3] == 2*tg->nlocaledeges); // Error in counting the edges 
	return vertexCount;

}

uint64_t* getDataArray(){

}


void createDistributedGraph(const tuple_graph* const tg, distributedGraph_CSR* const graph){
    graph->nGlobalEdges = tg->nglobaledges;
    

    // devide edges equally on threads
    divideTuplegraph_divisible(tg);
	graph->nLocaledges = tg->nlocaledeges;

	// get indices for each vertex
	graph->indices = getVertexSpacing(tg);		

	// fill data array
	printf("Hello\n");
}
