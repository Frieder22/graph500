#if !defined(CSR_CUSTOM_H)
#define CSR_CUSTOM_H

#include "common.h"


typedef struct{
    size_t nGlobalEdges;
    size_t nLocaledges;
    uint64_t data;
    size_t* reference;
} distributedGraph_CSR;

void divideTuplegraph_divisible(tuple_graph* const tg);

void createDistributedGraph(const tuple_graph* const tg, distributedGraph_CSR* const graph);


#endif // CSR_CUSTOM_H
