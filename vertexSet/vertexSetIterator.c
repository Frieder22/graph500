#include "vertexSet.h"
#include "vertexSetIterator.h"
#include "bitmap_custom.h"

void vertexSetIterator_Init(vertexSetIterator* const it, Vertexset* const vs) {
    // set corresponding vertex set
    it->vs = vs;

    // set corresponding representation
    it->isdense = vs->isdense;

    // setting for dense representation
    if (it->isdense) {
        it->wordIdx = 0;
        it->wordShift = 0;
        if (vs->approx_halfing) {
            it->word = vs->bitArray[(vs->size_bitarray - vs->indexShift)%vs->size_bitarray];
        } else {
            it->word = vs->bitArray[0];
        }

    } else {
        // setting for sparse representation
        it->sparseIdx = 0;
    } 
};

bool vertexSetIterator_Has_next(vertexSetIterator* const it){
    if (it->isdense) {
        if (!it->vs->approx_halfing){
            for ( size_t i = it->wordIdx; i < it->vs->size_bitarray; i++){
                // find a word thats not 0
                if (it->vs->bitArray[i] != 0ULL){
                    // find the 1s in the word
                    for (int j = it->wordShift; j < ulong_bits; j++) {
                        if (it->word &1UL == 1) { //found a 1
                            // save current vertex
                            it->currVertex = i * ulong_bits + j;
                            
                            //shift one to right
                            it->word >>= 1; 
                            // save progress
                            it->wordIdx = i;
                            it->wordShift = j + 1;
                            return true;
                        }
                        it->word >>= 1; //shift one to right
                    }
                }
                it->wordIdx = i + 1;
                it->wordShift = 0;
                it->word = it->vs->bitArray[i+1];
            }
            return false;
        } else {
            for ( size_t wordIdx = it->wordIdx; wordIdx < it->vs->size_bitarray; wordIdx++){
                int i = (wordIdx - it->vs->indexShift + it->vs->size_bitarray) % it->vs->size_bitarray;
                // find a word thats not 0
                if (it->vs->bitArray[i] != 0ULL){
                    // find the 1s in the word
                    for (int j = it->wordShift; j < ulong_bits; j++) {
                        if (it->word &1UL == 1) { //found a 1
                            // save current vertex
                            it->currVertex = wordIdx * ulong_bits + j;
                            
                            //shift one to right
                            it->word >>= 1; 
                            // save progress
                            it->wordIdx = wordIdx;
                            it->wordShift = j + 1;
                            return true;
                        }
                        it->word >>= 1; //shift one to right
                    }
                }
                it->wordIdx = wordIdx + 1;
                it->wordShift = 0;
                it->word = it->vs->bitArray[(i+1)%it->vs->size_bitarray];
            }
            return false;
        }
        
        
    } else {
        return it->sparseIdx < it->vs->sizeSparse;
    }
};

uint32_t vertexSetIterator_Next(vertexSetIterator* const it){
    if (it->isdense) {
        return it->currVertex;
    } else {
        // increase of sparseIdxx and return of value in same line
        return it->vs->sparseArray[it->sparseIdx++];
    }
    
};

void vertexSetIterator_Reset(vertexSetIterator* const it){
    if (it->isdense) {
        Vertexset *vs = it->vs;
        it->wordIdx = 0;
        it->wordShift = 0;
        if (vs->approx_halfing) {
            it->word = vs->bitArray[(vs->size_bitarray - vs->indexShift)%vs->size_bitarray];
        } else {
            it->word = vs->bitArray[0];
        }
    } else {
        it->sparseIdx = 0;
    } 
};

