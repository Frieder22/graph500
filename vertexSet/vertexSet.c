#include "vertexSet.h"
#include "bitmap_custom.h"
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <mpi.h>
#include <assert.h>
#include <stdio.h>
#include <utils.h>


void Vertexset_Init(Vertexset* const vs, const uint32_t maxsize, const MPI_Comm MPI_COMM, bool shiftPattern){
    vs->maxsize = maxsize;
    
    // init bitArray
    size_t size_bitarray = (maxsize + ulong_bits - 1) / ulong_bits;
    vs->size_bitarray = size_bitarray;
    vs->bitArray = (unsigned long long*) malloc(size_bitarray * sizeof(unsigned long long));
    vs->bitBuffer = (unsigned long long*) malloc(size_bitarray * sizeof(unsigned long long));
    // init to beeing 0
    Bitmap_Clean(vs->bitArray, vs->size_bitarray);
    Bitmap_Clean(vs->bitBuffer, vs->size_bitarray);

    
    // init MPI stuff
    vs->MPI_COMM = MPI_COMM;
    MPI_Comm_rank(MPI_COMM, &vs->mpi_rank);
    MPI_Comm_size(MPI_COMM, &vs->mpi_size);

    // check if systems definition of ULL is compatible with this implementation
    assert(sizeof(unsigned long long) == 8); // 8 byte (=64 bit)

    // calculate size, where dense variant is seen as more efficient (mean strategy)
    vs->sizeCrit = ((maxsize + ulong_bits - 1)  * (size_t) (vs->mpi_size - 1)) / (32* vs->mpi_size * Vertexset_Log2_floor(vs->mpi_size));
    
    // init sparse array 
    // allocate to much memory to fill up whole array with all possible
    // vertices (no duplicates)
    vs->sizeSparse = 0;
    vs->sparseArray = (uint32_t*) malloc((vs->sizeCrit * 2 + 2) * sizeof(uint32_t));
    vs->sparseBuffer = (uint32_t*) malloc(vs->sizeCrit * 2 * sizeof(uint32_t));
    // buffer for informations of other ranks
    vs->sizesAll = (int*) malloc(vs->mpi_size * sizeof(int));
    vs->displ= (int*) malloc((vs->mpi_size + 1) * sizeof(int));

    // set default as sparse
    vs->isdense = false;

    // set usage of approximate halfing
    vs->shiftPattern = shiftPattern;

    if (vs->shiftPattern) {
        // calculate start and end of each block
        int blockIndices_input[vs->mpi_size + 1];
        int block_nElements_input[vs->mpi_size];
        
        int *blockIndices_buffer, *block_nElements_buffer;
        
        blockIndices_buffer = (int*) malloc(sizeof(int)*(vs->mpi_size + 1));
        block_nElements_buffer = (int*) malloc(sizeof(int)* (vs->mpi_size));

        // calculate indices for input (same for every process)
        blockIndices_input[0] = 0;
        for (size_t i = 1; i < vs->mpi_size + 1; i++) {
            blockIndices_input[i] = (vs->size_bitarray * i)  / vs->mpi_size;
            block_nElements_input[i-1] = blockIndices_input[i] - blockIndices_input[i-1];
            assert(blockIndices_input[i]>=0);
        }
        int indexShift = blockIndices_input[vs->mpi_rank];
    
        // calculate indices for R buffer
        blockIndices_buffer[0] = 0;
        for (size_t i = 1; i < vs->mpi_size + 1; i++) {
            blockIndices_buffer[i] = blockIndices_buffer[i - 1] + block_nElements_input[(i - 1 + vs->mpi_rank) % vs->mpi_size ];
            block_nElements_buffer[i-1] = blockIndices_buffer[i] - blockIndices_buffer[i-1];
            assert(blockIndices_buffer[i]>=0);

        }

        vs->block_Indices = blockIndices_buffer;
        vs->block_nElements = block_nElements_buffer;
        vs->indexShift = indexShift;
    }
};

void Vertexset_SetCritsizeStrategy(Vertexset* const vs, critSizeStrategies citSizeStrategy){
    //calculate critical Size depending on the stategies
    size_t bitsInDenseArray = vs->size_bitarray * ulong_bits;
    switch (citSizeStrategy){
    case CRITSIZE_MAX:
        vs->sizeCrit = bitsInDenseArray / 32;
        break;
    case CRITSIZE_MIN:
        vs->sizeCrit = bitsInDenseArray / (32 * vs->mpi_size);
        break;
    case CRITSIZE_MEAN:
        vs->sizeCrit = (bitsInDenseArray * (size_t) (vs->mpi_size - 1)) / (32* vs->mpi_size * Vertexset_Log2_floor(vs->mpi_size));
        break;
    default:
        return;
    }
    

    // adapt allocated memory for sparse buffer
    uint32_t* temp;
    temp = (uint32_t*) realloc(vs->sparseArray, (vs->sizeCrit * 2 + 2) * sizeof(uint32_t));
    if (!temp){
        free(temp);
        assert(false && "couldn't realloc for Crit size strategy");
    }
    vs->sparseArray = temp;

    temp = (uint32_t*) realloc(vs->sparseBuffer, vs->sizeCrit * 2 * sizeof(uint32_t));
    if (!temp){
        free(temp);
        assert(false && "couldn't realloc for Crit size strategy");
    }
    vs->sparseBuffer = temp;
}

void Vertexset_SetCritsize(Vertexset* const vs, int sizeCrit){
    vs->sizeCrit = sizeCrit;

    // adapt allocated memory for sparse buffer
    vs->sparseArray = (uint32_t*) realloc(vs->sparseArray, (vs->sizeCrit * 2 + 2) * sizeof(uint32_t));
    vs->sparseBuffer = (uint32_t*) realloc(vs->sparseBuffer, vs->sizeCrit * 2 * sizeof(uint32_t));
}

void Vertexset_Add(Vertexset* const vs, const uint32_t vertex) {
    assert(vertex < vs->maxsize);
    if (vs->isdense) {
        if (!vs->shiftPattern) {
            Bitmap_Set(vs->bitArray, vertex);
        } else {
            Bitmap_Set_Shifted(vs->bitArray, vertex, vs->indexShift, vs->size_bitarray);
        }
    } else {
        vs->sparseArray[vs->sizeSparse++] = vertex;
        if (vs->sizeSparse >= vs->sizeCrit) {
            Vertexset_TransformToDense(vs);
        }
    }
};

bool Vertexset_Contains(const Vertexset* const vs, const uint32_t vertex) {
    assert(vertex < vs->maxsize);
    if (vs->isdense) {
        if (!vs->shiftPattern) {
            return Bitmap_Test(vs->bitArray, vertex);
        } else {
            return Bitmap_Test_Shifted(vs->bitArray, vertex, vs->indexShift, vs->size_bitarray);
        }  
    } else {
        for (size_t i = 0; i < vs->sizeSparse; i++) {
            if (vs->sparseArray[i] == vertex) {
                return true;
            }
        }
        return false;
    }
};

void Vertexset_Clean(Vertexset* const vs) {
    Bitmap_Clean(vs->bitArray, vs->size_bitarray);
    vs->sizeSparse = 0;
    vs->isdense = false;
};

bool Vertexset_TransformToDense(Vertexset* const vs) {
    // only do, if it's sparse
    if (!vs->isdense){
        uint32_t vertex;
        Bitmap_Clean(vs->bitArray, vs->size_bitarray);
        for (size_t i = 0; i < vs->sizeSparse; i++){
            vertex = vs->sparseArray[i];
            if (!vs->shiftPattern) {
                Bitmap_Set(vs->bitArray, vertex);
            } else {
                Bitmap_Set_Shifted(vs->bitArray, vertex, vs->indexShift, vs->size_bitarray);
            }
        }
        vs->isdense = true;
    }
    
};

bool Vertexset_TransformToSparse(Vertexset* const vs) {
    // only do, if it is dense
    if (vs->isdense){
        unsigned long long word;
        size_t countAdded = 0;
        if (vs->shiftPattern) {            
            for (size_t i = 0; i < vs->size_bitarray; i++){
                // find a word thats not 0
                if (vs->bitArray[i] != 0ULL){
                    word = vs->bitArray[i]; 
                    // find the 1s in the word
                    for (int j = 0; j < ulong_bits; j++) {
                        if (word &1 == 1) { //found a 1
                            // add to sparse list
                            vs->sparseArray[countAdded] = i * ulong_bits + j;
                            countAdded ++;
                        }
                        word >>= 1; //shift one to right
                    }
                }
            }
        } else {
            uint32_t vert;   
            for (size_t i = 0; i < vs->size_bitarray; i++){
                // find a word thats not 0
                if (vs->bitArray[i] != 0ULL){
                    word = vs->bitArray[i]; 
                    // find the 1s in the word
                    for (int j = 0; j < ulong_bits; j++) {
                        if (word &1 == 1) { //found a 1
                            // add to sparse list
                            vert = ((i + vs->indexShift)% vs->size_bitarray)  * ulong_bits + j;
                            if (vert < vs->maxsize){
                                vs->sparseArray[countAdded] = vert;
                                countAdded ++;
                            }
                        }
                        word >>= 1; //shift one to right
                    }
                }
            }
        }
        
        vs->sizeSparse = countAdded;
        vs->isdense = false;    
    }
}

void Vertexset_Union_Butterfly(Vertexset* const vs, const int VERTEXSET_OPERATION) {
    assert(VERTEXSET_OPERATION == VERTEXSET_OR); // no other operator implemented
    assert((vs->mpi_size & (vs->mpi_size - 1)) == 0); // communicator must be size of 2^k
    assert(!vs->shiftPattern); //works only for non shifed entries
    
    // find ciritical size
    size_t criticalSize = vs->sizeCrit;

    // find block indices
    int blockIdx[vs->mpi_size + 1];
    blockIdx[0] = 0;
    for (size_t i = 1; i < vs->mpi_size + 1; i++) {
        blockIdx[i] = vs->size_bitarray * i / vs->mpi_size;
    }
    
    // do reduce_scatter
    int commNeighbor, commNeighborRev;
    int startBlock;
    int startIndex, startIndexRecv;
    int nElements;
    int tag;
    MPI_Status status;
    MPI_Request req;
    int recvCount;
    int nIterations = Vertexset_Log2_floor(vs->mpi_size);
    int reverseRank = Vertexset_reverseBits(vs->mpi_rank, nIterations);
    for (int shift = vs->mpi_size / 2 ; shift >= 1; shift/=2) {
        // do bitflip with LOR to find neighbor
        commNeighbor = reverseRank ^ shift;
        commNeighborRev = Vertexset_reverseBits(commNeighbor, nIterations);

        // decide, which strategy to choose from
        if (!vs->isdense) {
            // send sparse array
            MPI_Isend(vs->sparseArray, vs->sizeSparse, MPI_INT32_T, commNeighborRev, 100, vs->MPI_COMM, &req);
        } else {
            // find start indices of send
            startBlock = (commNeighbor/shift) * shift;
            startIndex = blockIdx[startBlock];
            nElements = blockIdx[startBlock + shift] - startIndex;

            assert(nElements>=0);
            // send dense array
            MPI_Isend(vs->bitArray + startIndex, nElements, MPI_UNSIGNED_LONG_LONG, commNeighborRev, 200, vs->MPI_COMM, &req);
        }
        
        
        // get info, if recieved message is in dense or sparse format
        MPI_Probe(commNeighborRev, MPI_ANY_TAG, vs->MPI_COMM, &status);
        tag = status.MPI_TAG;
        
        if (tag==100) {
            // sparse array is recieved
            // get length of message
            MPI_Get_count(&status, MPI_INT32_T, &recvCount);
            // check if transformation to dense is needed
            if (!vs->isdense) {
                // recieve sparse format
                MPI_Recv(vs->sparseArray + vs->sizeSparse, recvCount, MPI_INT32_T, commNeighborRev, 100, vs->MPI_COMM, MPI_STATUS_IGNORE);
                vs->sizeSparse += recvCount;
                MPI_Wait(&req, MPI_STATUS_IGNORE);

                // transform to dense, if too big
                if (vs->sizeSparse > vs->sizeCrit) {
                    uint32_t vertex;
                    Bitmap_Clean(vs->bitArray, vs->size_bitarray);
                    for (size_t i = 0; i < vs->sizeSparse; i++){
                        vertex = vs->sparseArray[i];
                        Bitmap_Set(vs->bitArray, vertex);
                    }
                    vs->isdense = true;
                }
                
            } else {
                // fill into dense format
                MPI_Recv(vs->sparseBuffer, recvCount, MPI_INT32_T, commNeighborRev, 100, vs->MPI_COMM, MPI_STATUS_IGNORE);
                MPI_Wait(&req, MPI_STATUS_IGNORE);
                uint32_t vert;
                for (size_t i = 0; i < recvCount; i++) {
                    vert = vs->sparseBuffer[i];
                    Bitmap_Set(vs->bitArray, vert);
                }            
            }


        } else if (tag==200) {
            // dense array is recieved          
            // get length of message
            MPI_Get_count(&status, MPI_UNSIGNED_LONG_LONG, &recvCount);
            
            startBlock = (reverseRank/shift) * shift;
            startIndex = blockIdx[startBlock];
            
            // choose reduction style depenmding on set state (dense/sparse)            
            if (vs->isdense) {
                // Recieve into buffer
                MPI_Recv(vs->bitBuffer, recvCount, MPI_UNSIGNED_LONG_LONG, commNeighborRev, 200, vs->MPI_COMM, MPI_STATUS_IGNORE);
                // own set is dense, so just logic OR is performed
                MPI_Wait(&req, MPI_STATUS_IGNORE);

                // do reduction
                assert(recvCount == blockIdx[startBlock + shift]-startIndex);
                for (int i = startIndex; i < blockIdx[startBlock + shift]; i++) {
                    vs->bitArray[i] |= vs->bitBuffer[i-startIndex];
                }

            } else {
                // own set is sparse
                // recieve directly into main bitarray (no buffering is needed)
                // also we can assume the memory of the bitarrays is 0 everywhere.
                MPI_Recv(vs->bitArray + startIndex, recvCount, MPI_UNSIGNED_LONG_LONG, commNeighborRev, 200, vs->MPI_COMM, MPI_STATUS_IGNORE);

                // no wait is needed, as we used sparse array in the send and now only dense arrays are needed
                
                // reduction
                uint32_t vert;
                for (size_t i = 0; i < vs->sizeSparse; i++) {
                    vert = vs->sparseArray[i];
                    Bitmap_Set(vs->bitArray, vert);
                }
                vs->isdense=true;
            }
            
            
        } else {
            assert("recieved wrong tag" && false);
        }      
    }
    
    // In case only sparse communication was performed, we don't need allgather
    if (!vs->isdense) {
        return;
    }
    

    // do allGather
    for (int shift = 1; shift < vs->mpi_size; shift*=2)  {
        // do bitflip with LOR to find neighbor
        commNeighbor = reverseRank ^ shift;
        commNeighborRev = Vertexset_reverseBits(commNeighbor, nIterations);

        // find start indices of send
        startBlock = (reverseRank/shift) * shift;
        startIndex = blockIdx[startBlock];
        nElements = blockIdx[startBlock + shift] - startIndex;
        
        // calculate start index and recv count
        startBlock = (commNeighbor/shift) * shift;
        startIndexRecv = blockIdx[startBlock];
        recvCount = blockIdx[startBlock + shift]-startIndexRecv;
        
        // Sending and recieval of blocks
        MPI_Sendrecv(vs->bitArray + startIndex, nElements, MPI_UNSIGNED_LONG_LONG, commNeighborRev, 101,
                     vs->bitArray + startIndexRecv, recvCount, MPI_UNSIGNED_LONG_LONG, commNeighborRev, 101,
                    vs->MPI_COMM, MPI_STATUS_IGNORE);
    }

    // must be dense, if allgather was needed
    vs->isdense = true;
};

void Vertexset_Union_Shift(Vertexset* const vs, const int VERTEXSET_OPERATION) {
    assert(VERTEXSET_OPERATION == VERTEXSET_OR); //no other version is implemented
    assert(vs->shiftPattern);  //works only for shifted entries

    int criticalSize = vs->sizeCrit;
    int to;
    int from;
    int nElements;
    int startIndex;
    int tag;
    MPI_Status status;
    MPI_Request reqDense, reqSparse;
    int recvCount;        
    int shift, shift_old, shift_next;
    bool inBuff;
    int sizeBuff, lastBlockSize;
    int sizeSparse, tail;

    // find number of iterations
    int iterations;
    if ((vs->mpi_size & (vs->mpi_size - 1)) == 0) {
        // just the log, if mpi_size == 2^k
        iterations = Vertexset_Log2_floor(vs->mpi_size);
    } else {
        // ceil(log), if mpi_size =/= 2^k
        iterations = Vertexset_Log2_floor(vs->mpi_size) + 1;
    }

    // allocate memory to save skip sequence
    int skipSequence[iterations + 1];
    skipSequence[iterations] = vs->mpi_size;


    // shifts for first iteration
    shift_old = vs->mpi_size;
    shift = 1 << Vertexset_Log2_floor(vs->mpi_size-1);

    // sparse array setup
    sizeSparse = vs->sizeSparse;
    vs->sparseArray[sizeSparse] = sizeSparse;
    tail = vs->mpi_size;
    inBuff = false;

    // do adapted Reducescatter
    for (int i = 0; i < iterations; i++) {
        skipSequence[iterations - i - 1] = shift;
        // find communication neighbors
        to = (vs->mpi_rank + shift) % vs->mpi_size;
        from = (vs->mpi_rank - shift + vs->mpi_size) % vs->mpi_size;
        

        // decide, which variant should be send
        if (!vs->isdense) {
            // Handle tail location. Either in buffer or ontop of pile
            if (tail>shift || tail == 0) {
                // tail is bigger than shift and must be sent
                if (inBuff) {
                    // buffered block must be copied into pile
                    MPI_Wait(&reqSparse, MPI_STATUS_IGNORE);
                    memcpy(vs->sparseArray +  sizeSparse, vs->sparseBuffer, (sizeBuff+1)*sizeof(uint32_t));
                    sizeSparse += sizeBuff;
                    inBuff = false;
                }
                tail -= shift;
            } else {
                // tail is smaller or equal to shift and must be buffered
                if (!inBuff) {
                    // last block is is copied from pile to buffer
                    sizeBuff = vs->sparseArray[sizeSparse];
                    memcpy(vs->sparseBuffer, vs->sparseArray + sizeSparse - sizeBuff, (sizeBuff+1)*sizeof(uint32_t));
                    sizeSparse -= sizeBuff;
                    inBuff=true;
                }
            }
            // send sparse array
            MPI_Isend(vs->sparseArray, sizeSparse + 1 - inBuff, MPI_UINT32_T, to, inBuff, vs->MPI_COMM, &reqSparse);
        } else {
            // find start indices of dense send
            startIndex = vs->block_Indices[shift];
            
            // find the size of dense send
            nElements = vs->block_Indices[shift_old] - startIndex;

            // send dense array
            MPI_Isend(vs->bitArray + startIndex, nElements, MPI_UNSIGNED_LONG_LONG, to, 200, vs->MPI_COMM, &reqDense);
        }

        // get info, if recieved message is in dense or sparse format
        MPI_Probe(from, MPI_ANY_TAG, vs->MPI_COMM, &status);
        tag = status.MPI_TAG;
        
        if (tag < 200) {
            // sparse array is recieved
            // get length of message
            MPI_Get_count(&status, MPI_UINT32_T, &recvCount);
            
            if (!vs->isdense) {
                // reCieve directly in place
                // send sparse array
                MPI_Recv(vs->sparseArray + sizeSparse + 1 - inBuff, recvCount, MPI_UINT32_T, from, MPI_ANY_TAG, vs->MPI_COMM, MPI_STATUS_IGNORE);

                if (!inBuff) {
                    MPI_Wait(&reqSparse, MPI_STATUS_IGNORE);
                    // find new last block size of reveived data
                    lastBlockSize = vs->sparseArray[sizeSparse + recvCount];
                    
                    // overwrite old lastBlockSize with valid element that is not from new last block
                    vs->sparseArray[sizeSparse] = vs->sparseArray[sizeSparse + recvCount - lastBlockSize - 1];
                    
                    // overwrite valid element that is not from new last block with element from new last block
                    vs->sparseArray[sizeSparse + recvCount - lastBlockSize - 1] = vs->sparseArray[sizeSparse + recvCount - 1];

                    // place new lastBlockSize in correct position
                    vs->sparseArray[sizeSparse + recvCount - 1] = vs->sparseArray[sizeSparse + recvCount]; 
                }
                sizeSparse += recvCount - 1 + inBuff;
                
                // transform into dense, if too long
                if (sizeSparse + (inBuff*sizeBuff) > criticalSize) { //also count buffered elements!
                    // add vertices from sparse array
                    uint32_t vert;
                    for (size_t i = 0; i < sizeSparse; i++) {
                        vert = vs->sparseArray[i];
                        Bitmap_Set_Shifted(vs->bitArray, vert, vs->indexShift, vs->size_bitarray);
                    }
                    // add vertives from buffer                
                    if (inBuff) {
                        for (size_t i = 0; i < sizeBuff; i++) {
                            vert = vs->sparseBuffer[i];
                            Bitmap_Set_Shifted(vs->bitArray, vert, vs->indexShift, vs->size_bitarray);
                        }
                    }
                    vs->isdense = true;
                }
            } else {
                // fill into dense format
                inBuff = tag;
                MPI_Recv(vs->sparseArray, recvCount, MPI_UINT32_T, from, tag, vs->MPI_COMM, MPI_STATUS_IGNORE);
                
                uint32_t vert;
                // add received vertices from sparse array
                sizeSparse = recvCount - 1 + inBuff;          
                MPI_Wait(&reqDense, MPI_STATUS_IGNORE);
                for (size_t i = 0; i < sizeSparse; i++) {
                    vert = vs->sparseArray[i];
                    Bitmap_Set_Shifted(vs->bitArray, vert, vs->indexShift, vs->size_bitarray);
                }
            }
        } else {
            // dense array is recieved          
            // get length of message (can be done implicitly)
            recvCount = vs->block_Indices[shift_old - shift];
            // Recieve into buffer.
            MPI_Recv(vs->bitBuffer, recvCount, MPI_UNSIGNED_LONG_LONG, from, 200, vs->MPI_COMM, MPI_STATUS_IGNORE);

            if(vs->isdense){
                // Wait until send buffer message is safe
                MPI_Wait(&reqDense, MPI_STATUS_IGNORE);

                // do reduction           
                for (int i = 0; i < recvCount; i++) {
                    vs->bitArray[i] |= vs->bitBuffer[i];
                }
            } else {
                // swap bitmap buffer and bit array
                unsigned long long* temp = vs->bitArray;
                vs->bitArray = vs->bitBuffer;
                vs->bitBuffer = temp;
                temp = NULL;

                uint32_t vert;
                // add vertices from sparse array
                for (size_t i = 0; i < sizeSparse; i++) {
                    vert = vs->sparseArray[i];
                    Bitmap_Set_Shifted(vs->bitArray, vert, vs->indexShift, vs->size_bitarray);
                }
                // add vertices from bufferReverse
                if (inBuff) {
                    for (size_t i = 0; i < sizeBuff; i++) {
                        vert = vs->sparseBuffer[i];
                        Bitmap_Set_Shifted(vs->bitArray, vert, vs->indexShift, vs->size_bitarray);
                    }
                }
                vs->isdense = true;
            }
        }
        // update shift
        shift_old = shift;
        shift >>=1; // only exactly halfing
    }
    
    // decide, if allgather is needed.
    // In case only sparse communication was performed, we don't need allgather
    if (!vs->isdense) {
        if (inBuff) {
            // buffered block must be copied into pile
            memcpy(vs->sparseArray +  sizeSparse, vs->sparseBuffer, sizeBuff*sizeof(uint32_t));
            sizeSparse += sizeBuff;
            inBuff = false;
        }
        vs->sizeSparse = sizeSparse;
        return;
    }
    

    // do allGather
    for (int i=0; i < iterations; i++)  {
        // find shifts
        shift = skipSequence[i];
        shift_next = skipSequence[i+1];        
        
        // find communication neighbors
        to = (vs->mpi_rank - shift + vs->mpi_size) % vs->mpi_size;
        from = (vs->mpi_rank + shift) % vs->mpi_size;
        
        recvCount = vs->block_Indices[shift_next] - vs->block_Indices[shift];
        MPI_Sendrecv(vs->bitArray, vs->block_Indices[shift_next - shift], MPI_UNSIGNED_LONG_LONG, to, 101,
                     vs->bitArray + vs->block_Indices[shift], recvCount, MPI_UNSIGNED_LONG_LONG, from,101,
                     vs->MPI_COMM, MPI_STATUS_IGNORE);
    }
};

void Union_Allreduce_Naive(Vertexset* const vs, const int VERTEXSET_OPERATION){
    assert(VERTEXSET_OPERATION == VERTEXSET_OR); // no other operator implemented
    int displ [vs->mpi_size];

    int criticalSize = vs->sizeCrit;

    if (vs->isdense) {
        vs->sizeSparse = criticalSize + 1;
    }
    

    // find how many entries are proposed.
    MPI_Allgather(&(vs->sizeSparse), 1, MPI_INT, vs->sizesAll, 1, MPI_INT, MPI_COMM_WORLD);
    
    // calculate displacements with exclusive scan
    displ[0] = 0;
    for (int j = 1; j < vs->mpi_size + 1; j++){
        displ[j] = displ[j-1] + vs->sizesAll[j-1];
    }

    if (displ[vs->mpi_size] < criticalSize) {
        MPI_Allgatherv(vs->sparseArray, vs->sizeSparse, MPI_INT32_T, vs->sparseBuffer, vs->sizesAll, displ, MPI_INT32_T, MPI_COMM_WORLD);

        // switch pointer
        uint32_t *temp;
        temp = vs->sparseBuffer;
        vs->sparseBuffer = vs->sparseArray;
        vs->sparseArray = temp;
        temp = NULL;

        // set new sparse size
        vs->sizeSparse = displ[vs->mpi_size];
    } else {
        Vertexset_TransformToDense(vs);
        MPI_Allreduce(MPI_IN_PLACE, vs->bitArray, vs->size_bitarray, MPI_UINT64_T, MPI_BOR, MPI_COMM_WORLD);
    }
};

void Vertexset_Allreduce_Butterfly(Vertexset* const vs, const int VERTEXSET_OPERATION){
    assert(VERTEXSET_OPERATION == VERTEXSET_OR); // no other operator implemented
    assert((vs->mpi_size & (vs->mpi_size - 1)) == 0); // communicator must be size of 2^k

    // find block indices
    int blockIdx[vs->mpi_size + 1];
    blockIdx[0] = 0;
    for (size_t i = 1; i < vs->mpi_size + 1; i++) {
        blockIdx[i] = vs->size_bitarray * i / vs->mpi_size;
    }
    
    // do reduce_scatter
    int commNeighbor;
    int startBlock;
    int startIndex;
    int startIndexRecv;
    int nElements;
    int tag;
    MPI_Status status;
    MPI_Request req;
    int recvCount;
    for (int shift = vs->mpi_size / 2 ; shift >= 1; shift/=2) {
        // do bitflip with LOR to find neighbor
        commNeighbor = vs->mpi_rank ^ shift;

        // find start indices of send message
        startBlock = (commNeighbor/shift) * shift;
        startIndex = blockIdx[startBlock];
        nElements = blockIdx[startBlock + shift] - startIndex;

        // get indices of recieve message
        startBlock = (vs->mpi_rank/shift) * shift;
        startIndexRecv = blockIdx[startBlock];
        recvCount = blockIdx[startBlock + shift]-startIndexRecv;

        MPI_Sendrecv(vs->bitArray + startIndex, nElements, MPI_UNSIGNED_LONG_LONG, commNeighbor, 200,
                    vs->bitBuffer, recvCount, MPI_UNSIGNED_LONG_LONG, commNeighbor, 200, vs->MPI_COMM, MPI_STATUS_IGNORE);

        for (int i = startIndexRecv; i < blockIdx[startBlock + shift]; i++) {
            vs->bitArray[i] |= vs->bitBuffer[i-startIndexRecv];
        }
        
    }

    // do allGather
    for (int shift = 1; shift < vs->mpi_size; shift*=2)  {
        // do bitflip with LOR to find neighbor
        commNeighbor = vs->mpi_rank ^ shift;

        // find start indices of send
        startBlock = (vs->mpi_rank/shift) * shift;
        startIndex = blockIdx[startBlock];
        nElements = blockIdx[startBlock + shift] - startIndex;
        
        // calculate start index and recv count
        startBlock = (commNeighbor/shift) * shift;
        startIndexRecv = blockIdx[startBlock];
        recvCount = blockIdx[startBlock + shift]-startIndexRecv;
        
        // Sending and recieval of blocks
        MPI_Sendrecv(vs->bitArray + startIndex, nElements, MPI_UNSIGNED_LONG_LONG, commNeighbor, 101,
                    vs->bitArray + startIndexRecv, recvCount, MPI_UNSIGNED_LONG_LONG, commNeighbor, 101,
                    vs->MPI_COMM, MPI_STATUS_IGNORE);
    } 
};

void Vertexset_Allreduce_Butterfly_ReverseBits(Vertexset* const vs, const int VERTEXSET_OPERATION){
    assert(VERTEXSET_OPERATION == VERTEXSET_OR); // no other operator implemented
    assert((vs->mpi_size & (vs->mpi_size - 1)) == 0); // communicator must be size of 2^k

    // find block indices
    int blockIdx[vs->mpi_size + 1];
    blockIdx[0] = 0;
    for (size_t i = 1; i < vs->mpi_size + 1; i++) {
        blockIdx[i] = vs->size_bitarray * i / vs->mpi_size;
    }
    
    // do reduce_scatter
    int commNeighbor;
    int commNeighborRev;
    int startBlock;
    int startIndex;
    int startIndexRecv;
    int nElements;
    MPI_Request req;
    int recvCount;
    int nIterations = Vertexset_Log2_floor(vs->mpi_size);
    int reverseRank = Vertexset_reverseBits(vs->mpi_rank, nIterations);
    for (int shift = vs->mpi_size / 2 ; shift >= 1; shift/=2) {
        // do bitflip with LOR to find neighbor
        commNeighbor = reverseRank ^ shift;
        commNeighborRev = Vertexset_reverseBits(commNeighbor, nIterations);

        // find start indices of send message
        startBlock = (commNeighbor/shift) * shift;
        startIndex = blockIdx[startBlock];
        nElements = blockIdx[startBlock + shift] - startIndex;

        // get indices of recieve message
        startBlock = (reverseRank/shift) * shift;
        startIndexRecv = blockIdx[startBlock];
        recvCount = blockIdx[startBlock + shift]-startIndexRecv;

        MPI_Sendrecv(vs->bitArray + startIndex, nElements, MPI_UNSIGNED_LONG_LONG, commNeighborRev, 200,
            vs->bitBuffer, recvCount, MPI_UNSIGNED_LONG_LONG, commNeighborRev, 200, vs->MPI_COMM, MPI_STATUS_IGNORE);
            
        for (int i = startIndexRecv; i < blockIdx[startBlock + shift]; i++) {
            vs->bitArray[i] |= vs->bitBuffer[i-startIndexRecv];
        }
        
    }

    // do allGather
    for (int shift = 1; shift < vs->mpi_size; shift*=2)  {
        // do bitflip with LOR to find neighbor
        commNeighbor = reverseRank ^ shift;
        commNeighborRev = Vertexset_reverseBits(commNeighbor, nIterations);

        // find start indices of send
        startBlock = (reverseRank/shift) * shift;
        startIndex = blockIdx[startBlock];
        nElements = blockIdx[startBlock + shift] - startIndex;
        
        // calculate start index and recv count
        startBlock = (commNeighbor/shift) * shift;
        startIndexRecv = blockIdx[startBlock];
        recvCount = blockIdx[startBlock + shift]-startIndexRecv;
        
        // Sending and recieval of blocks
        MPI_Sendrecv(vs->bitArray + startIndex, nElements, MPI_UNSIGNED_LONG_LONG, commNeighborRev, 101,
                    vs->bitArray + startIndexRecv, recvCount, MPI_UNSIGNED_LONG_LONG, commNeighborRev, 101,
                    vs->MPI_COMM, MPI_STATUS_IGNORE);
    }  

};

void Vertexset_Allreduce_Shift(Vertexset* const vs, const int VERTEXSETOPERATION){
    assert(VERTEXSETOPERATION ==VERTEXSET_OR);
    // find number of iterations
    int iterations;
    if ((vs->mpi_size & (vs->mpi_size - 1)) == 0) {
        // just the log, if mpi_size == 2^k
        iterations = Vertexset_Log2_floor(vs->mpi_size);
    } else {
        // ceil(log), if mpi_size =/= 2^k
        iterations = Vertexset_Log2_floor(vs->mpi_size) + 1;
    }

    // allocate memory to save skip sequence
    int skipSequence[iterations + 1];
    skipSequence[iterations] = vs->mpi_size;
    

    // do reduce_scatter
    int to;
    int from;
    int recvCount;
    int nElements;
    int startIndex;
    int shift = Vertexset_nearestLog2(vs->mpi_size);
    int shift_old = vs->mpi_size;
    int shift_next;
    for (int i = 0; i < iterations; i++) {
        // Save skip sequence
        skipSequence[iterations - i - 1] = shift;

        // find communication neighbors
        to = (vs->mpi_rank + shift) % vs->mpi_size;
        from = (vs->mpi_rank - shift + vs->mpi_size) % vs->mpi_size;
        
        // find start indices of dense send
        startIndex = vs->block_Indices[shift];
        
        // find the size of dense send
        nElements = vs->block_Indices[shift_old] - startIndex;

        // find the size that is recieved          
        recvCount = vs->block_Indices[shift_old - shift];

        MPI_Sendrecv(vs->bitArray + startIndex, nElements, MPI_UNSIGNED_LONG_LONG, to, 200,
                        vs->bitBuffer, recvCount, MPI_UNSIGNED_LONG_LONG, from, 200, vs->MPI_COMM, MPI_STATUS_IGNORE);
        // do reduction           
        for (int i = 0; i < recvCount; i++) {
            vs->bitArray[i] |= vs->bitBuffer[i];
        }

        // update shift
        shift_old = shift;
        shift >>= 1; // exactly halfing
        
    }
    

    // do allGather
    for (int i=0; i < iterations; i++)  {
        // find shifts
        shift = skipSequence[i];
        shift_next = skipSequence[i+1];        
        
        // find communication neighbors
        to = (vs->mpi_rank - shift + vs->mpi_size) % vs->mpi_size;
        from = (vs->mpi_rank + shift) % vs->mpi_size;
        
        recvCount = vs->block_Indices[shift_next] - vs->block_Indices[shift];
        MPI_Sendrecv(vs->bitArray, vs->block_Indices[shift_next - shift], MPI_UNSIGNED_LONG_LONG, to, 101,
                    vs->bitArray + vs->block_Indices[shift], recvCount, MPI_UNSIGNED_LONG_LONG, from,101,
                    vs->MPI_COMM, MPI_STATUS_IGNORE);
    }
    vs->isdense = true; 
}

void Vertexset_Allgather_Butterfly(Vertexset* vs, int VERTEXSET_OPERATION){
    assert(VERTEXSET_OPERATION == VERTEXSET_OR); // no other operator implemented
    assert((vs->mpi_size & (vs->mpi_size - 1)) == 0); // communicator must be size of 2^k

    
    // do Allgather
    int commNeighbor;
    MPI_Status status;
    MPI_Request req;
    int recvCount;
    for (int shift = vs->mpi_size / 2 ; shift >= 1; shift/=2) {
        // do bitflip with LOR to find neighbor
        commNeighbor = vs->mpi_rank ^ shift;
       
        // send sparse array
        MPI_Isend(vs->sparseArray, vs->sizeSparse, MPI_UINT32_T, commNeighbor, 100, vs->MPI_COMM, &req);
        
        // get info, if recieved message is in dense or sparse format
        MPI_Probe(commNeighbor, 100, vs->MPI_COMM, &status);
        // get length of message
        MPI_Get_count(&status, MPI_UINT32_T, &recvCount);

        // recieve sparse format
        MPI_Recv(vs->sparseArray + vs->sizeSparse, recvCount, MPI_UINT32_T, commNeighbor, 100, vs->MPI_COMM, MPI_STATUS_IGNORE);
        MPI_Wait(&req, MPI_STATUS_IGNORE);

        vs->sizeSparse += recvCount;
    }
}

void Vertexset_Allgather_Butterfly_ReverseBits(Vertexset* vs, int VERTEXSET_OPERATION){
    assert(VERTEXSET_OPERATION == VERTEXSET_OR); // no other operator implemented
    assert((vs->mpi_size & (vs->mpi_size - 1)) == 0); // communicator must be size of 2^k
    
    int nIterations = Vertexset_Log2_floor(vs->mpi_size);
    int reverseRank = Vertexset_reverseBits(vs->mpi_rank, nIterations);

    // do Allgather
    int commNeighbor, commNeighborRev;
    MPI_Status status;
    MPI_Request req;
    int recvCount;
    for (int shift = vs->mpi_size / 2 ; shift >= 1; shift/=2) {
        // do bitflip with LOR to find neighbor
        commNeighbor = reverseRank ^ shift;
        commNeighborRev = Vertexset_reverseBits(commNeighbor, nIterations);

       
        // send sparse array
        MPI_Isend(vs->sparseArray, vs->sizeSparse, MPI_INT32_T, commNeighborRev, 100, vs->MPI_COMM, &req);
        
        // get info, if recieved message is in dense or sparse format
        MPI_Probe(commNeighborRev, MPI_ANY_TAG, vs->MPI_COMM, &status);
        // get length of message
        MPI_Get_count(&status, MPI_INT32_T, &recvCount);
        // recieve sparse format
        MPI_Recv(vs->sparseArray + vs->sizeSparse, recvCount, MPI_INT32_T, commNeighborRev, 100, vs->MPI_COMM, MPI_STATUS_IGNORE);
        vs->sizeSparse += recvCount;
        MPI_Wait(&req, MPI_STATUS_IGNORE);

    }
}

void Vertexset_Allgather_Shift(Vertexset* vs, int  VERTEXSET_OPERATION){
    // safety measure
    if (vs->mpi_size < 2) {
        return;
    }
    
    int shift = 1 << Vertexset_Log2_floor(vs->mpi_size-1);
    int to, from;
    int nBlocks;
    MPI_Request req;
    MPI_Status status;
    int recvCount;
    int lastBlockSize;
    int sizeBuff;
    int tail = vs->mpi_size;
    bool inBuff = false;
    int sizeSparse = vs->sizeSparse;

    // set size of last block
    vs->sparseArray[sizeSparse] = sizeSparse;
    while (shift >=1 ) {
        to = (vs->mpi_rank - shift + vs->mpi_size) % vs->mpi_size;
        from = (vs->mpi_rank + shift) % vs->mpi_size;
        if (tail>shift || tail == 0) {
            // tail is bigger than shift and must be sent
            if (inBuff) {
                // buffered block must be copied into pile
                MPI_Wait(&req, MPI_STATUS_IGNORE);
                memcpy(vs->sparseArray +  sizeSparse, vs->sparseBuffer, (sizeBuff+1)*sizeof(uint32_t));
                
                sizeSparse += sizeBuff;
                inBuff = false;
            }
            tail -= shift;
        } else {
            // tail is smaller or equal to shift and must be buffered
            if (!inBuff) {
                // last block is is copied from pile to buffer
                sizeBuff = vs->sparseArray[sizeSparse];
                memcpy(vs->sparseBuffer, vs->sparseArray + sizeSparse - sizeBuff, (sizeBuff+1)*sizeof(uint32_t));
                sizeSparse -= sizeBuff;
                inBuff=true;
            }
            
        }

        
        MPI_Isend(vs->sparseArray, sizeSparse + 1 - inBuff, MPI_UINT32_T, to, 100, vs->MPI_COMM, &req);
        MPI_Probe(from, MPI_ANY_TAG, vs->MPI_COMM, &status);
        MPI_Get_count(&status, MPI_UINT32_T, &recvCount);
        MPI_Recv(vs->sparseArray + sizeSparse + 1 - inBuff, recvCount, MPI_UINT32_T, from, 100, vs->MPI_COMM, MPI_STATUS_IGNORE);
        
        if (!inBuff) {
            MPI_Wait(&req, MPI_STATUS_IGNORE);
            // find new last block size of reveived data
            lastBlockSize = vs->sparseArray[sizeSparse + recvCount];

            // overwrite old lastBlockSize with valid element that is not from new last block
            vs->sparseArray[sizeSparse] = vs->sparseArray[sizeSparse + recvCount - lastBlockSize - 1];
            
            // overwrite valid element that is not from new last block with element from new last block
            vs->sparseArray[sizeSparse + recvCount - lastBlockSize - 1] = vs->sparseArray[sizeSparse + recvCount - 1];

            // place new lastBlockSize in correct position
            vs->sparseArray[sizeSparse + recvCount - 1] = vs->sparseArray[sizeSparse + recvCount]; 
        }
        
        
        sizeSparse += recvCount - 1 + inBuff;
        shift >>= 1; //exactly halfing
    }
    
    // copy block on pile, if it is buffered after last iteration
    if (inBuff) {
        // buffered block must be copied into pile
        memcpy(vs->sparseArray +  sizeSparse, vs->sparseBuffer, (sizeBuff)*sizeof(uint32_t));
        sizeSparse += sizeBuff;
        inBuff = false;
    }    
    vs->sizeSparse = sizeSparse;


}

void Vertexset_PrintSet(Vertexset* const vs){
    if(vs->isdense){
        Vertexset_TransformToSparse(vs);
    }
    int max_elements = 30;
    printf("Elements in rank %d: ", vs->mpi_rank);
    for (size_t i = 0; i < vs->sizeSparse; i++) {
        printf("%d ", vs->sparseArray[i]);
        if (i > max_elements) {
            printf("...");
            break;
        }
        
    }
    printf("\n");    
};

void Vertexset_Deinit(Vertexset* const vs){
    free(vs->bitArray);
    free(vs->sparseArray);
    free(vs->bitBuffer);
    free(vs->sparseBuffer);
    free(vs->displ);
    free(vs->sizesAll);
};