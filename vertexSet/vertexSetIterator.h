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


/**
 * Initialized Iterator for a given vertexSet
 * @param vs vertexSet, the iterator is for
 * @param it iterator, that should be inititalized
 */
void vertexSetIterator_Init(vertexSetIterator* it, const Vertexset* vs);

/**
 * Checks, if the iterator can return another value (true)
 * or is at the end of its set (false)
 * @param it corresponding iterator
 * @return bool, that indicates, if another element is
 *         available
 */
bool vertexSetIterator_Has_next(vertexSetIterator* it);

/**
 * Returns Next item of iterator.
 * @param it corresponding iterator
 * @return returns next element
 */
uint32_t vertexSetIterator_Next(vertexSetIterator* it);

/**
 * Resets the iterator to point to first element.
 * @param it corresponding iterator
 */
void vertexSetIterator_Reset(vertexSetIterator* it);


#endif // VERTEXITERATOR
