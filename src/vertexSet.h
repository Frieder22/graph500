#if !defined(VERTEXSET)
#define VERTEXSET

#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>
#include <bitmap_custom.h>

typedef struct Vertexset{
    size_t size;
    size_t maxsize;
    size_t sizeCrit;
    bool isdense;

    unsigned long long *bitarray;
    uint32_t size_bitarray;

    uint32_t *sparseArray;
} Vertexset;

void Vertexset_Init(Vertexset* vs, uint32_t maxsize);

void Vertexset_Add(uint32_t vertex, Vertexset* vs);

bool Vertexset_Contains(uint32_t vertexx, Vertexset* vs);

bool Vertexset_TransformToDense(Vertexset* vs);

bool Vertexset_TransformToSparse(Vertexset* vs);

void Vertexset_Allreduce(Vertexset* vs, int VERTEXSET_OPERATION);

void Vertexset_Deinit(Vertexset* vs);


#define VERTEXSET_OR 0
#define VERTEXSET_AND 1
#define VERTEXSET_XOR 2

#endif // VERTEXSET
