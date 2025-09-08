#include "vertexSet.h"
#include "bitmap_custom.h"
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <mpi.h>
#include <assert.h>
#include <stdio.h>

static inline int Vertexset_Log2_floor(int x){
    int i;
    for (i = 0; i < 32; i++) {
        x = x >> 1;
        if (x <= 0){
            break;
        }
    }
    return i;
};


void Vertexset_Init(Vertexset* const vs, const uint32_t maxsize, const MPI_Comm MPI_COMM){
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
    // calculate size, where dense variant is more efficient
    vs->sizeCrit = size_bitarray * sizeof(unsigned long long) / sizeof(uint32_t);
    
    // init sparse array 
    // allocate to much memory to fill up whole array with all possible
    // vertices (no duplicates)
    vs->sizeSparse = 0;
    vs->sparseArray = (uint32_t*) malloc(vs->sizeCrit * 2 * sizeof(uint32_t));
    vs->sparseBuffer = (uint32_t*) malloc(vs->sizeCrit * 2 * sizeof(uint32_t));
    // buffer for informations of other ranks
    vs->sizesAll = (int*) malloc(vs->mpi_size * sizeof(int));
    vs->displ= (int*) malloc((vs->mpi_size + 1) * sizeof(int));

    // set default as sparse
    vs->isdense = false;

    // set usage of approximate halfing
    vs->approx_halfing = (vs->mpi_size & (vs->mpi_size - 1)) != 0;

    if (vs->approx_halfing) {
        // calculate start and end of each block
        int blockIndices_input[vs->mpi_size + 1];
        int block_nElements_input[vs->mpi_size];
        
        int *blockIndices_buffer, *block_nElements_buffer;
        
        blockIndices_buffer = (int*) malloc(sizeof(int)*(vs->mpi_size + 1));
        block_nElements_buffer = (int*) malloc(sizeof(int)* (vs->mpi_size));

        // calculate indices for input (same for every thread)
        blockIndices_input[0] = 0;
        for (size_t i = 1; i < vs->mpi_size + 1; i++) {
            blockIndices_input[i] = vs->size_bitarray * i / vs->mpi_size;
            block_nElements_input[i-1] = blockIndices_input[i] - blockIndices_input[i-1];
        }
        int indexShift = blockIndices_input[vs->mpi_rank];
    
        // calculate indices for R buffer
        blockIndices_buffer[0] = 0;
        for (size_t i = 1; i < vs->mpi_size + 1; i++) {
            blockIndices_buffer[i] = blockIndices_buffer[i - 1] + block_nElements_input[(i - 1 + vs->mpi_rank) % vs->mpi_size ];
            block_nElements_buffer[i-1] = blockIndices_buffer[i] - blockIndices_buffer[i-1];
        }

        vs->block_Indices = blockIndices_buffer;
        vs->block_nElements = block_nElements_buffer;
        vs->indexShift = indexShift;
    }
    

};

void Vertexset_Add(Vertexset* const vs, const uint32_t vertex) {
    assert(vertex < vs->maxsize);
    if (vs->isdense) {
        if (!vs->approx_halfing) {
            Bitmap_Set(vs->bitArray, vertex);
        } else {
            Bitmap_Set_Shifted(vs->bitArray, vertex, vs->indexShift, vs->size_bitarray);
        }
    } else {
        assert(vs->sizeSparse <= vs->maxsize);
        vs->sparseArray[vs->sizeSparse++] = vertex;
        if (vs->sizeSparse >= vs->sizeCrit) {
            Vertexset_TransformToDense(vs);
        }
    }
};

bool Vertexset_Contains(const Vertexset* const vs, const uint32_t vertex) {
    assert(vertex < vs->maxsize);
    if (vs->isdense) {
        if (!vs->approx_halfing) {
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
};

bool Vertexset_TransformToDense(Vertexset* const vs) {
    // only do, if it's sparse
    if (!vs->isdense){
        uint32_t vertex;
        Bitmap_Clean(vs->bitArray, vs->size_bitarray);
        for (size_t i = 0; i < vs->sizeSparse; i++){
            vertex = vs->sparseArray[i];
            if (!vs->approx_halfing) {
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
        if (vs->approx_halfing) {            
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

void Vertexset_Allreduce(Vertexset* const vs, const int VERTEXSET_OP){
    // no other variant is implemented yet
    // use only sparse communication
    assert(false && "Automatic algorithm selection is not implemented");
}


void Vertexset_Allreduce_Pure(Vertexset* const vs, const int VERTEXSET_OPERATION){
    assert(VERTEXSET_OPERATION == VERTEXSET_OR); //no other version is implemented
    assert(!vs->approx_halfing);
    // get size information from other ranks
    int size_int = (int) vs->sizeSparse;
    MPI_Allgather(&size_int, 1, MPI_INT, vs->sizesAll, 1, MPI_INT, vs->MPI_COMM);

    // calculate displacements with exxclusive scan
    vs->displ[0] = 0;
    for (int i = 1; i < vs->mpi_size + 1; i++){
        vs->displ[i] = vs->displ[i-1] + vs->sizesAll[i-1];
    }

    if (vs->isdense) {
        MPI_Allreduce(MPI_IN_PLACE, vs->bitArray, vs->size_bitarray, MPI_UINT64_T, MPI_BOR, MPI_COMM_WORLD);
    } else {
        // Ceck, if resulting array would be too big
        if(vs->displ[vs->mpi_size] < vs->maxsize){
            // Do dense communication
            Vertexset_TransformToDense(vs);
            Vertexset_Allreduce_Pure(vs, VERTEXSET_OPERATION);
            vs->sizeSparse = 1;
            return;
        }


        // gather of all subarrays
        MPI_Allgatherv(vs->sparseArray, vs->sizeSparse, MPI_UINT32_T, vs->sparseBuffer, vs->sizesAll, (vs->displ), MPI_UINT32_T, vs->MPI_COMM);

        // switch pointer
        vs->sparseArray = vs->sparseBuffer;
        vs->sparseBuffer = NULL;
    }
    // set size
    vs->sizeSparse = vs->displ[vs->mpi_size];
    assert(vs->sizeSparse >= 0);
};

void Vertexset_Allreduce_Exact_Halfing(Vertexset* const vs, const int VERTEXSET_OPERATION) {
    assert(VERTEXSET_OPERATION == VERTEXSET_OR); // no other operator implemented
    assert((vs->mpi_size & (vs->mpi_size - 1)) == 0); // communicator must be size of 2^k
    assert(!vs->approx_halfing); //works only for non shifed entries
    
    // find ciritical size
    size_t criticalSize = vs->sizeCrit;

    // find block indices
    int blockIdx[vs->mpi_size + 1];
    blockIdx[0] = 0;
    for (int i = 1; i < vs->mpi_size + 1; i++) {
        blockIdx[i] = vs->size_bitarray * i / vs->mpi_size;
    }
    
    // do reduce_scatter
    int commNeighbor;
    int startBlock;
    int startIndex, startIndexRecv;
    int nElements;
    int tag;
    MPI_Status status;
    MPI_Request req;
    int recvCount;
    for (int shift = vs->mpi_size / 2 ; shift >= 1; shift/=2) {
        // do bitflip with LOR to find neighbor
        commNeighbor = vs->mpi_rank ^ shift;

        // decide, which strategy to choose from
        if (!vs->isdense) {
            // send sparse array
            MPI_Isend(vs->sparseArray, vs->sizeSparse, MPI_INT32_T, commNeighbor, 100, vs->MPI_COMM, &req);
        } else {
            // find start indices of send
            startBlock = (commNeighbor/shift) * shift;
            startIndex = blockIdx[startBlock];
            nElements = blockIdx[startBlock + shift] - startIndex;

            // send dense array
            MPI_Isend(vs->bitArray + startIndex, nElements, MPI_UNSIGNED_LONG_LONG, commNeighbor, 200, vs->MPI_COMM, &req);
        }
        
        
        // get info, if recieved message is in dense or sparse format
        MPI_Probe(commNeighbor, MPI_ANY_TAG, vs->MPI_COMM, &status);
        tag = status.MPI_TAG;
        
        if (tag==100) {
            // sparse array is recieved
            // get length of message
            MPI_Get_count(&status, MPI_INT32_T, &recvCount);
            // check if transformation to dense is needed
            if (!vs->isdense) {
                // recieve sparse format
                MPI_Recv(vs->sparseArray + vs->sizeSparse, recvCount, MPI_INT32_T, commNeighbor, 100, vs->MPI_COMM, MPI_STATUS_IGNORE);
                vs->sizeSparse += recvCount;

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
                MPI_Recv(vs->sparseBuffer, recvCount, MPI_INT32_T, commNeighbor, 100, vs->MPI_COMM, MPI_STATUS_IGNORE);
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
            
            
            // choose reduction style depenmding on set state (dense/sparse)            
            if (vs->isdense) {
                // Recieve into buffer
                MPI_Recv(vs->bitBuffer, recvCount, MPI_UNSIGNED_LONG_LONG, commNeighbor, 200, vs->MPI_COMM, MPI_STATUS_IGNORE);
                // own set is dense, so just logic OR is performed
                MPI_Wait(&req, MPI_STATUS_IGNORE);

                // do reduction
                startBlock = (vs->mpi_rank/shift) * shift;
                startIndex = blockIdx[startBlock];
                assert(recvCount == blockIdx[startBlock + shift]-startIndex);
                for (int i = startIndex; i < blockIdx[startBlock + shift]; i++) {
                    vs->bitArray[i] |= vs->bitBuffer[i-startIndex];
                }

            } else {
                // own set is sparse
                // recieve directly into main bitarray (no buffering is needed)
                // also we can assume the memory of the bitarrays is 0 everywhere.
                startBlock = (vs->mpi_rank/shift) * shift;
                startIndex = blockIdx[startBlock];
                Bitmap_Clean(vs->bitArray, vs->size_bitarray);
                MPI_Recv(vs->bitArray + startIndex, recvCount, MPI_UNSIGNED_LONG_LONG, commNeighbor, 200, vs->MPI_COMM, MPI_STATUS_IGNORE);

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

    // must be dense, if allgather was needed
    vs->isdense = true;
};


void Vertexset_Allreduce_Approximate_Halfing(Vertexset* const vs, const int VERTEXSET_OPERATION) {
    assert(VERTEXSET_OPERATION == VERTEXSET_OR); //no other version is implemented
    assert(vs->approx_halfing);  //works only for shifted entries

    int criticalSize = vs->sizeCrit;

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
    int sendNeighbor;
    int recvNeighbor;
    int nElements;
    int startIndex;
    int tag;
    MPI_Status status;
    MPI_Request req;
    int recvCount;
    int shift = vs->mpi_size;
    int shift_old, shift_next;
    for (int i = 0; i < iterations; i++) {
        // update shift
        shift_old = shift;
        shift = shift - shift/2; // ceil(shift/2)
        skipSequence[iterations - i - 1] = shift;
        
        // find communication neighbors
        sendNeighbor = (vs->mpi_rank + shift) % vs->mpi_size;
        recvNeighbor = (vs->mpi_rank - shift + vs->mpi_size) % vs->mpi_size;
        

        // decide, which variant should be send
        if (!vs->isdense) {
            // send sparse array
            MPI_Isend(vs->sparseArray, vs->sizeSparse, MPI_INT32_T, sendNeighbor, 100, vs->MPI_COMM, &req);
        } else {
            // find start indices of dense send
            startIndex = vs->block_Indices[shift];
            
            // find the size of dense send
            nElements = vs->block_Indices[shift_old] - startIndex;

            // send dense array
            MPI_Isend(vs->bitArray + startIndex, nElements, MPI_UNSIGNED_LONG_LONG, sendNeighbor, 200, vs->MPI_COMM, &req);
        }

        
        // get info, if recieved message is in dense or sparse format
        MPI_Probe(recvNeighbor, MPI_ANY_TAG, vs->MPI_COMM, &status);
        tag = status.MPI_TAG;
        
        if (tag==100) {
            // sparse array is recieved
            // get length of message
            MPI_Get_count(&status, MPI_INT32_T, &recvCount);
            
            if (!vs->isdense) {
                MPI_Recv(vs->sparseArray + vs->sizeSparse, recvCount, MPI_INT32_T, recvNeighbor, 100, vs->MPI_COMM, MPI_STATUS_IGNORE);
                vs->sizeSparse += recvCount;

                // transform to dense, if too big
                if (vs->sizeSparse > vs->sizeCrit) {
                    uint32_t vertex;
                    Bitmap_Clean(vs->bitArray, vs->size_bitarray);
                    for (size_t i = 0; i < vs->sizeSparse; i++){
                        vertex = vs->sparseArray[i];
                        Bitmap_Set_Shifted(vs->bitArray, vertex, vs->indexShift, vs->size_bitarray);
                    }
                    vs->isdense = true;
                }
            
            } else {
                // fill into dense format
                MPI_Recv(vs->sparseBuffer, recvCount, MPI_INT32_T, recvNeighbor, 100, vs->MPI_COMM, MPI_STATUS_IGNORE);
                MPI_Wait(&req, MPI_STATUS_IGNORE);
                uint32_t vert;
                for (size_t i = 0; i < recvCount; i++) {
                    vert = vs->sparseBuffer[i];
                    Bitmap_Set_Shifted(vs->bitArray, vert, vs->indexShift, vs->size_bitarray);
                }  
            }
            
        } else if (tag==200) {
            // dense array is recieved          
            // get length of message (can be done explicitly)
            recvCount = vs->block_Indices[shift_old - shift];
            
            // choose reduction style depenmding on set state (dense/sparse)            
            if (vs->isdense) {
                // Recieve into buffer.
                MPI_Recv(vs->bitBuffer, recvCount, MPI_UNSIGNED_LONG_LONG, recvNeighbor, 200, vs->MPI_COMM, MPI_STATUS_IGNORE);

                // Wait until send buffer message is safe
                MPI_Wait(&req, MPI_STATUS_IGNORE);

                // do reduction           
                for (int i = 0; i < recvCount; i++) {
                    vs->bitArray[i] |= vs->bitBuffer[i];
                }
            } else {
                // own set is sparse
                // recieve directly into main bitarray (no buffering is needed)
                // also we can assume the memory of the bitarrays is 0 everywhere.
                Bitmap_Clean(vs->bitArray, vs->size_bitarray);
                MPI_Recv(vs->bitArray, recvCount, MPI_UNSIGNED_LONG_LONG, recvNeighbor, 200, vs->MPI_COMM, MPI_STATUS_IGNORE);

                // no wait is needed, as we used sparse array in the send and now only dense arrays are needed
                
                // reduction
                uint32_t vert;
                for (size_t i = 0; i < vs->sizeSparse; i++) {
                    vert = vs->sparseArray[i];
                    Bitmap_Set_Shifted(vs->bitArray, vert, vs->indexShift, vs->size_bitarray);
                }
                vs->isdense=true;
            }
            
        } else {
            assert("recieved wrong tag" && false);
        }
    }
    
    // decide, if allgather is needed.
    // In case only sparse communication was performed, we don't need allgather
    if (!vs->isdense) {
        return;
    }


    // do allGather
    for (int i=0; i < iterations; i++)  {
        // find shifts
        shift = skipSequence[i];
        shift_next = skipSequence[i+1];        
        
        // find communication neighbors
        sendNeighbor = (vs->mpi_rank - shift + vs->mpi_size) % vs->mpi_size;
        recvNeighbor = (vs->mpi_rank + shift) % vs->mpi_size;
        
        recvCount = vs->block_Indices[shift_next] - vs->block_Indices[shift];
        MPI_Sendrecv(vs->bitArray, vs->block_Indices[shift_next - shift], MPI_UNSIGNED_LONG_LONG, sendNeighbor, 101,
                     vs->bitArray + vs->block_Indices[shift], recvCount, MPI_UNSIGNED_LONG_LONG, recvNeighbor,101,
                     vs->MPI_COMM, MPI_STATUS_IGNORE);
    }
    vs->isdense = true;
};

void Vertexset_Allreduce_Ring_Comm(Vertexset* const vs, const int VERTEXSET_OPERATION){
    assert(VERTEXSET_OPERATION == VERTEXSET_OR); // no other operator implemented
    
    // find block indices
    int blockIdx[vs->mpi_size + 1];
    blockIdx[0] = 0;
    for (int i = 1; i < vs->mpi_size + 1; i++) {
        blockIdx[i] = vs->size_bitarray * i / vs->mpi_size;
    }

    // find ciritical size
    size_t criticalSize = vs->sizeCrit;

    // prepare, if vertexset is sparse
    if (!(vs->isdense)) {
        Bitmap_Clean(vs->bitArray, vs->size_bitarray);
        int count=0;
        uint32_t vert;
        for (size_t i = 0; i < vs->sizeSparse; i++) {
            vert = vs->sparseArray[i];
            if (count <= (criticalSize/vs->mpi_size)){
                // add to sparse list and update bitarray
                if (!(Bitmap_Test(vs->bitArray, vert))) {
                    // add to correct sparse list, while avoiding dublicates
                    vs->sparseBuffer[count++] = vert;
                    Bitmap_Set(vs->bitArray, vert);
                }
            } else {
                // update bitarray (without caring about sparse)
                Bitmap_Set(vs->bitArray, vert);
            }
        }
        // set correctly counted size
        vs->sizeSparse = count;
        
        // switching pointer to original
        if (count <= (criticalSize/vs->mpi_size)) {
            // only neccary, when sparse variant is needed   
            uint32_t *temp;
            temp = vs->sparseArray;
            vs->sparseArray = vs->sparseBuffer;
            vs->sparseBuffer = temp;
            temp = NULL;        
        }
    }

    if (vs->isdense) {
        vs->sizeSparse = criticalSize + 1;
    }

    // do ring communication
    int sendBlock, nElementsDense, startIdx;
    int recvNeigh, sendNeigh;
    int recvCount, tag;
    MPI_Request req;
    MPI_Status status;
    for (size_t iteration = 0; iteration < 2*vs->mpi_size; iteration++) {
        // find, which block to send
        sendBlock = (vs->mpi_rank + iteration) % vs->mpi_size;

        // find start and size of block
        startIdx = blockIdx[sendBlock];
        nElementsDense = blockIdx[sendBlock + 1] - startIdx;

        // find communication neighbors
        recvNeigh = (vs->mpi_rank - 1 + vs->mpi_size) % vs->mpi_size;
        sendNeigh = (vs->mpi_rank + 1) % vs->mpi_size;


        // devide, if send sparse or dense array
        if (vs->sizeSparse <= criticalSize/vs->mpi_size) {
            // do sparse send
            MPI_Isend(vs->sparseArray, vs->sizeSparse, MPI_INT32_T, sendNeigh, 100, vs->MPI_COMM, &req);            
        } else {
            // do dense send
            MPI_Isend(vs->bitArray + startIdx, nElementsDense, MPI_UNSIGNED_LONG_LONG, sendNeigh, 200, vs->MPI_COMM, &req);
        }

        // find which kind of message is recieved
        MPI_Probe(recvNeigh, MPI_ANY_TAG, vs->MPI_COMM, &status);
        tag = status.MPI_TAG;
        if (tag==100) {
            // sparse recieve
            MPI_Get_count(&status, MPI_INT32_T, &recvCount);

            // recieve into buffer
            MPI_Recv(vs->sparseBuffer, recvCount, MPI_INT32_T, recvNeigh, 100, vs->MPI_COMM, MPI_STATUS_IGNORE);

            // wait for save usage of sparse array
            MPI_Wait(&req, MPI_STATUS_IGNORE);

            int count=vs->sizeSparse;
            uint32_t vert;
            // do reduction (append non dublicates and also update bitarray)
            for (size_t j = 0; j < recvCount; j++)  {
                vert = vs->sparseBuffer[j];
                if (count <= criticalSize/vs->mpi_size || iteration == vs->mpi_size - 1) {
                    // set bitmap, and add to sparse,if needed
                    if (!Bitmap_Test(vs->bitArray, vert)) {
                        // save in sparse array until its needed in next round
                        Bitmap_Set(vs->bitArray, vert);
                        vs->sparseArray[count++] = vert;                        
                    }
                } else {
                    Bitmap_Set(vs->bitArray, vert);
                }
            }

            // set count
            vs->sizeSparse = count;

        } else if (tag == 200) {
            // do dense recieve
            MPI_Get_count(&status, MPI_UNSIGNED_LONG_LONG, &recvCount);

            // recieve int dense buffer
            MPI_Recv(vs->bitBuffer, recvCount, MPI_UNSIGNED_LONG_LONG, recvNeigh, 200, vs->MPI_COMM, MPI_STATUS_IGNORE);

            // wait for save usage of dense array
            MPI_Wait(&req, MPI_STATUS_IGNORE);

            // do reduction
            startIdx = blockIdx[(recvNeigh + iteration) % vs->mpi_size];
            for (size_t j = 0; j < recvCount; j++) {
                vs->bitArray[startIdx + j] |= vs->bitBuffer[j];
            }

            vs->sizeSparse = criticalSize;
        }
        // last iteration defines format of output
        if(iteration == vs->mpi_size - 1){
            if (tag==100) {
                vs->isdense = false;
            } else if(tag==200) {
                vs->isdense = true;
            }    
        }
    }
};

void Vertexset_Allreduce_Dense(Vertexset* const vs, const int VERTEXSET_OPERATION){
    assert(VERTEXSET_OPERATION == VERTEXSET_OR); // no other operator implemented
    assert((vs->mpi_size & (vs->mpi_size - 1)) == 0); // communicator must be size of 2^k
    assert(vs->isdense);

    // find block indices
    int blockIdx[vs->mpi_size + 1];
    blockIdx[0] = 0;
    for (int i = 1; i < vs->mpi_size + 1; i++) {
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

        // send dense array
        //MPI_Isend(vs->bitArray + startIndex, nElements, MPI_UNSIGNED_LONG_LONG, commNeighbor, 200, vs->MPI_COMM, &req);

        // get indices of recieve message
        startBlock = (vs->mpi_rank/shift) * shift;
        startIndexRecv = blockIdx[startBlock];
        recvCount = blockIdx[startBlock + shift]-startIndexRecv;
        

        //MPI_Probe(commNeighbor, MPI_ANY_TAG, vs->MPI_COMM, &status);
        //MPI_Get_count(&status, MPI_LONG_LONG, &recvCount);
        MPI_Sendrecv(vs->bitArray + startIndex, nElements, MPI_UNSIGNED_LONG_LONG, commNeighbor, 200,
                    vs->bitBuffer, recvCount, MPI_UNSIGNED_LONG_LONG, commNeighbor, 200, vs->MPI_COMM, MPI_STATUS_IGNORE);

        

        // Recieve into buffer
        // MPI_Recv(vs->bitBuffer, recvCount, MPI_UNSIGNED_LONG_LONG, commNeighbor, 200, vs->MPI_COMM, MPI_STATUS_IGNORE);
        // MPI_Wait(&req, MPI_STATUS_IGNORE);
        // do reduction
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



void Vertexset_PrintSet(Vertexset* const vs){
    if(vs->isdense){
        Vertexset_TransformToSparse(vs);
    }
    int max_elements = 30;
    printf("Elements in rank %d: ", vs->mpi_rank);
    for (size_t i = 0; i < vs->sizeSparse; i++) {
        printf("%d ", vs->sparseArray[i]);
        if (i > max_elements) {
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



