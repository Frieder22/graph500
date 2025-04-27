#if !defined(CSR_CUSTOM_H)
#define CSR_CUSTOM_H

#include "common.h"


typedef struct{
    size_t nGlobalEdges;
    size_t nLocaledges;
    uint32_t* data;
    uint64_t* indices;
} distributedGraph_CSR;

/**
 * Divides the edges of an tuple graph evenly onto the ranks. Only works, if the number of
 * edges is divisble by the size of COMM_WORLD
 */
void divideTuplegraph_divisible(tuple_graph* const tg);

/**
 * Creates a CSR array to have fast access the edges of a node. The edges are distributed 
 * evenly on the ranks. (edge oriented distribution to ranks and vertex oriented organization
 * in each rank)
 */
void createDistributedGraph(const tuple_graph* const tg, distributedGraph_CSR* const graph);

/**
 * Find neighbours of a vertex.
 * @return start and pointer of array
 */
void getNeighbours(distributedGraph_CSR* const graph, uint32_t vertex, uint32_t* start,  uint32_t* end); 

/**
 * Correct cleanup of distributed graph struct
 * @todo implementation
 */
void freeDistributedGraph(distributedGraph_CSR* const graph);

#endif // CSR_CUSTOM_H
