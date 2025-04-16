#include "common.h"
#include "csr_custom.h"


int64_t divideTuplegraph_divisible(tuple_graph* const tg){
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
    return nlocaledges;

}


void distributedGraph(tuple_graph* const tg, distributedGraph_CSR* const graph){
    graph->nGlobalEdges = tg->nglobaledges;
    
    // devide edges equally on threads
    uint64_t nlocaledges = divideTuplegraph_divisible(tg);


									


}
