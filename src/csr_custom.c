#include "common.h"
#include "csr_custom.h"
#include <assert.h>

void divideTuplegraph_divisible(tuple_graph* const tg){
    // Assumption to work properly: rank 0 holds correct edge list
	//								number of edges is devisible by size of COMM_WORLD

	// distribute edges equally to the processes
	uint64_t nlocaledges = tg->nglobaledges / size; // every process has correct nglobaledges
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

uint64_t* getVertexSpacing(const tuple_graph* const tg){
	uint64_t* vertexCount;
	vertexCount = (uint64_t*) malloc((tg->nglobalverts + 1) * sizeof(uint64_t));
	
	// initialize with 0
	for (uint64_t i = 0; i < tg->nglobalverts + 1; i++) {
		vertexCount[i] = 0;
	}

	// count occurence of each vertex
	packed_edge edge;	
	for (uint64_t i = 0; i < tg->nlocaledeges; i++) {
		edge = tg->edgememory[i];
		vertexCount[edge.v0_low + 1]++;
		vertexCount[edge.v1_low + 1]++;
	}

	// perform scan to get starting positions in CSR
	for (uint64_t i = 1; i < tg->nglobalverts + 1; i++) {
		vertexCount[i] += vertexCount[i-1];
	}

	assert(vertexCount[tg->nglobalverts] == 2*tg->nlocaledeges); // Error in counting the edges 
	return vertexCount;

}

void setDataArray(const tuple_graph* const tg, distributedGraph_CSR* const graph){
	uint32_t* vertexCount;
	uint32_t* data;
	uint32_t vertexOffset, vertex1, vertex2;
	packed_edge edge;
	
	vertexCount = (uint32_t*) malloc(tg->nglobalverts * sizeof(uint32_t));
	data = (uint32_t*) malloc(graph->indices[tg->nglobalverts-1] * sizeof(uint32_t));

	// initialize with 0
	for (uint32_t i = 0; i < tg->nglobalverts; i++) {
		vertexCount[i] = 0;
	}

	// fill in data
	for (uint64_t i = 0; i < tg->nlocaledeges; i++) {
		// find corresponding edge
		edge = tg->edgememory[i];

		// set data in for outgoing direction
		vertex1 = edge.v0_low;
		vertex2 = edge.v1_low;
		vertexOffset = graph->indices[vertex1];
		data[vertexCount[vertex1] + vertexOffset] = vertex2;
		vertexCount[vertex1]++;

		// set data in for ingoing direction
		vertexOffset = graph->indices[vertex2];
		data[vertexCount[vertex2] + vertexOffset] = vertex1;
		vertexCount[vertex2]++;
	}

	// put in struct
	graph->data = data;

	free(vertexCount);
}


void printNeighbors(distributedGraph_CSR* const graph, uint32_t vertex, int checkRank){
	packed_edge edge;
	if (rank == checkRank) {
		printf("Neighbours of vertex %ld (for rank %d):\n", vertex, rank);
		for (uint64_t i = graph->indices[vertex]; i < graph->indices[vertex+1]; i++){
			printf("%ld ", graph->data[i]);
		}
		printf("\n");
	}
}

void createDistributedGraph(const tuple_graph* const tg, distributedGraph_CSR* const graph){
    // set nGlobaledges
	graph->nGlobalEdges = tg->nglobaledges;
    
    // divide edges equally on threads
    divideTuplegraph_divisible(tg);
	graph->nLocaledges = tg->nlocaledeges;

	// get indices for each vertex
	graph->indices = getVertexSpacing(tg);		

	// put in data
	setDataArray(tg, graph);
}
/*
void getNeighbours(distributedGraph_CSR* const graph, uint32_t vertex,  uint32_t *start, uint32_t *end){
	*start = &(graph->data[graph->indices[vertex]]);
	*end = &(graph->data[graph->indices[vertex + 1]]);
}
*/

void freeDistributedGraph(distributedGraph_CSR* const graph){
	free(graph->data);
	free(graph->indices);
	free(graph);
}