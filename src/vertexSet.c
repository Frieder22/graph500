
#include "vertexSet.h"
#include "bitmap_custom.h"
#include <stddef.h>
#include <stdint.h>
#include <mpi.h>


void Vertexset_Init(Vertexset* vs, uint32_t maxsize){
    vs->size = 0;
    vs->maxsize = maxsize;

    // init bitarray
    size_t size_bitarray = (maxsize + 1) / sizeof(unsigned long long);
    vs->size_bitarray = size_bitarray;
    vs->bitarray = (unsigned long long*) malloc(size_bitarray * sizeof(unsigned long long));

    // init sparse array 
    // allocate to much memory to fill up whole array with all possible
    // vertices (no dublicates)
    vs->sparseArray = (uint32_t*) malloc(maxsize * sizeof(uint32_t));
    vs->sparseArray[0] = NULL;

    //calculate size, where dense variant is more efficient
    vs->sizeCrit = size_bitarray * sizeof(unsigned long long) / sizeof(32);

    // set default as sparse
    vs->isdense = false;
};

void Vertexset_Add(uint32_t vertex, Vertexset* vs){
    if (vs->isdense) {
        SET_VISITED_CUSTOM((vs->bitarray), vertex);
    }
    
};

bool Vertexset_Contains(uint32_t vertexx, Vertexset* vs);

bool Vertexset_TransformToDense(Vertexset* vs);

bool Vertexset_TransformToSparse(Vertexset* vs);

void Vertexset_Allreduce(Vertexset* vs, int VERTEXSET_OP);


void Vertexset_Deinit(Vertexset* vs){
    free(vs->bitarray);
    free(vs->sparseArray);
};

