#if !defined(VERTEXSET)
#define VERTEXSET

#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>

typedef struct Vertexset{
    size_t size;
    size_t maxsize;
    size_t sizeCrit;
    bool isdense;

    uint64_t *bitarray;
    uint32_t size_bitarray;

    uint32_t *sparseArray;
} Vertexset;

void Vertexset_Init(Vertexset* vs, uint32_t maxsize);

void Vertexset_Add(uint32_t vertex, Vertexset* vs);

bool Vertexset_TransformToDense(Vertexset* vs);

bool Vertexset_TransformToSparse(Vertexset* vs);


void Vertexset_Free(Vertexset* vs);



#endif // VERTEXSET
