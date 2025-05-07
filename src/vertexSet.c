
#include "vertexSet.h"
#include <stddef.h>
#include <stdint.h>
#include <mpi.h>


void Vertexset_Init(Vertexset* vs, uint32_t maxsize){
    vs->size = 0;
    vs->maxsize = maxsize;

    // init bitarray
    size_t size_bitarray = (maxsize + 1) / sizeof(uint64_t);
    vs->size_bitarray = size_bitarray;
    vs->bitarray = (uint64_t*) malloc(size_bitarray * sizeof(uint64_t));

    // init sparse array 
    // allocate to much memory to fill up whole array with all possible
    // vertices (no dublicates)
    vs->sparseArray = (uint32_t*) malloc(maxsize * sizeof(uint32_t));

    //calculate size, where dense variant is more efficient
    vs->sizeCrit = size_bitarray * sizeof(uint64_t) / sizeof(32);
};

void Vertexset_Add(uint32_t vertex, Vertexset* vs);

bool Vertexset_TransformToDense(Vertexset* vs);

bool Vertexset_TransformToSparse(Vertexset* vs);


void Vertexset_Free(Vertexset* vs);

