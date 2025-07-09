#include "common.h"
#include "csr_custom.h"
#include <assert.h>

void divideTuplegraph_divisible(const tuple_graph* const tg, tuple_graph* const tg_split){
    // Assumption to work properly: rank 0 holds correct edge list
	//								every rank knows the number of global edges

	// Only works, if there are less than 2147483647 edges (2^32)
	assert(tg->nglobaledges <= INT_MAX);

	// find number of elements and displacements
	int nlocaledges[size];
	int offset[size+1];
	offset[0] = 0;
	for (size_t i = 1; i < size + 1; i++) {
		offset[i] = i * tg->nglobaledges / size;
		nlocaledges[i-1] = offset[i] - offset[i-1];
	}

	// Scatter edgelist amongst ranks
	packed_edge* local_edges;
	local_edges = (packed_edge*) malloc(sizeof(packed_edge) * nlocaledges[rank]);
	MPI_Scatterv(tg->edgememory, nlocaledges, offset,
				packed_edge_mpi_type, local_edges, nlocaledges[rank],
				packed_edge_mpi_type, 0, MPI_COMM_WORLD);
	
	// Set relevant info in new split variant
	tg_split->edgememory = local_edges;
	tg_split->nlocaledeges = nlocaledges[rank];
	tg_split->nglobaledges = tg->nglobaledges;

	// set to NULL for good measure ;)
	local_edges = NULL;		
}

uint64_t* getVertexSpacing(const tuple_graph* const tg, const uint64_t nglobalverts){
	uint64_t* vertexCount;
	vertexCount = (uint64_t*) malloc((nglobalverts + 1) * sizeof(uint64_t));
	
	// initialize with 0
	for (uint64_t i = 0; i < nglobalverts + 1; i++) {
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
	for (uint64_t i = 1; i < nglobalverts + 1; i++) {
		vertexCount[i] += vertexCount[i-1];
	}

	assert(vertexCount[nglobalverts] == 2*tg->nlocaledeges && "Error in counting the edges"); 
	return vertexCount;

}

void setDataArray(const tuple_graph* const tg, distributedGraph_CSR* const graph, const int64_t nglobalverts){
	uint32_t* vertexCount;
	uint32_t* data;
	uint32_t vertexOffset, vertex1, vertex2;
	packed_edge edge;

	vertexCount = (uint32_t*) malloc(nglobalverts * sizeof(uint32_t));
	data = (uint32_t*) malloc(graph->indices[nglobalverts] * sizeof(uint32_t));

	// initialize with 0
	for (uint32_t i = 0; i < nglobalverts; i++) {
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
		printf("Neighbours of vertex %d (for rank %d):\n", vertex, rank);
		for (uint64_t i = graph->indices[vertex]; i < graph->indices[vertex+1]; i++){
			printf("%d ", graph->data[i]);
		}
		printf("\n");
	}
}

void createDistributedGraph(const tuple_graph* const tg, distributedGraph_CSR* const graph, int64_t nglobalverts){
    // set nGlobaledges
	graph->nGlobalEdges = tg->nglobaledges;


	// set nGlobalVerts
	graph->nGlobalVerts = nglobalverts;
    
    // divide edges equally on threads
	tuple_graph tg_split;
    divideTuplegraph_divisible(tg, &tg_split);
	graph->nLocaledges = tg->nlocaledeges;

	// get indices for each vertex
	graph->indices = getVertexSpacing(&tg_split, nglobalverts);

	// put in data
	setDataArray(&tg_split, graph, nglobalverts);

	// free allocated memory from tg_split 
	free(tg_split.edgememory);
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