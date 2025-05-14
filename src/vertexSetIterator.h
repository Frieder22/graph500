#if !defined(VERTEXITERATOR)
#define VERTEXITERATOR

#include "vertexSet.h"

typedef struct vertexSetIterator {
    // general information
    Vertexset* vs;
    bool isdense;
    uint32_t currVertex;

    // for dense representation
    size_t wordIdx;
    int wordShift;
    uint64_t word;

    // for sparse representation
    size_t sparseIdx;
} vertexSetIterator;

void vertexSetIterator_Init(vertexSetIterator* it, Vertexset* vs);

bool vertexSetIterator_Has_next(vertexSetIterator* it);

uint32_t vertexSetIterator_Next(vertexSetIterator* it);

void vertexSetIterator_Reset(vertexSetIterator* it);

void vertexSetIterator_Deinit(vertexSetIterator* it);

#endif // VERTEXITERATOR
