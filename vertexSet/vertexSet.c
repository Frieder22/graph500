#include "vertexSet.h"
#include "bitmap_custom.h"
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <mpi.h>
#include <assert.h>
#include <stdio.h>

int Vertexset_Log2_floor(int x){
    int i;
    for (i = 0; i < 32; i++) {
        x = x >> 1;
        if (x <= 0){
            break;
        }
    }
    return i;
};

void Vertexset_Update_Sizes(int** sizesAll, int** buffer, int rank_diff, int mpi_size) {
    int *temp;
    for (int i = 0; i < mpi_size; i++) {
        (*buffer)[i] = (*sizesAll)[i] + (*sizesAll)[ (i + mpi_size - 1) % mpi_size];
    }
    temp = *sizesAll;
    *sizesAll = *buffer;
    *buffer = temp;
};

// from https://www.geeksforgeeks.org/qsort-function-in-c/
// Changed to sort descending
int comp_desc(const void *a, const void *b) {
    return (*(int *)b - *(int *)a);
};
int Vertexset_Find_critical_iteration(const Vertexset* vs){
    const int mpi_size = vs->mpi_size;
    const int *sizes = vs->sizesAll;
    
    // copy sizes array
    int sizes_copy[mpi_size];
    memcpy(sizes_copy, sizes, mpi_size);

    // sort sizes_copy in descending order
    qsort(sizes_copy, mpi_size, sizeof(int), comp_desc);

    // perform inclusive scan until critical size is reached
    int count;
    for (count = 0; count < mpi_size - 1; count++) {
        sizes_copy[count + 1] += sizes_copy[count];
        if (sizes_copy[count + 1] >  vs->sizeCrit) {
            break;
        }
    }
    
    // #count ranks are needed in worst case to overshoot
    // we need to switch in worst case floor(log_2(count)) iteration rounds
    return Vertexset_Log2_floor(count);
}


void Vertexset_Init(Vertexset* vs, uint32_t maxsize, MPI_Comm MPI_COMM){
    vs->maxsize = maxsize;
    
    // init bitArray
    size_t size_bitarray = (maxsize + (sizeof(unsigned long long)*8)) / (sizeof(unsigned long long)*8);
    vs->size_bitarray = size_bitarray;
    vs->bitArray = (unsigned long long*) malloc(size_bitarray * sizeof(unsigned long long));
    vs->bitBuffer = (unsigned long long*) malloc(size_bitarray * sizeof(unsigned long long));
    
    // init MPI stuff
    vs->MPI_COMM = MPI_COMM;
    MPI_Comm_rank(MPI_COMM, &vs->mpi_rank);
    MPI_Comm_size(MPI_COMM, &vs->mpi_size);
    
    // init sparse array 
    // allocate to much memory to fill up whole array with all possible
    // vertices (no duplicates)
    vs->sizeSparse = 0;
    vs->sparseArray = (uint32_t*) malloc(maxsize * sizeof(uint32_t));
    vs->sparseBuffer = (uint32_t*) malloc(maxsize * sizeof(uint32_t));
    // buffer for informations of other ranks
    vs->sizesAll = (int*) malloc(vs->mpi_size * sizeof(int));
    vs->displ= (int*) malloc((vs->mpi_size + 1) * sizeof(int));

    // check if systems definition of ULL is compatible with this implementation
    assert(sizeof(unsigned long long) == 8); // 8 byte (=64 bit)
    // calculate size, where dense variant is more efficient
    vs->sizeCrit = size_bitarray * sizeof(unsigned long long) / sizeof(uint32_t);

    // set default as sparse
    vs->isdense = false;
};

void Vertexset_Add(Vertexset* vs, uint32_t vertex) {
    assert(vertex < vs->maxsize);
    if (vs->isdense) {
        Bitmap_Set(vs->bitArray, vertex);
    } else {
        assert(vs->sizeSparse + 1 < vs->maxsize);
        vs->sparseArray[vs->sizeSparse] = vertex;
    }
    vs->sizeSparse++;
};

bool Vertexset_Contains(Vertexset* vs, uint32_t vertex) {
    assert(vertex < vs->maxsize);
    if (vs->isdense) {
        return Bitmap_Test(vs->bitArray, vertex);
    } else {
        for (size_t i = 0; i < vs->sizeSparse; i++) {
            if (vs->sparseArray[i] == vertex) {
                return true;
            }
        }
        return false;
    }
};

void Vertexset_Clean(Vertexset* vs) {
    if (vs->isdense) {
        Bitmap_Clean(vs->bitArray, vs->size_bitarray);
        vs->sizeSparse = 0;
    } else {
        vs->isdense = false;            
        vs->sizeSparse = 0;
    }
};

bool Vertexset_TransformToDense(Vertexset* vs) {
    // only do, if it's sparse
    if (!vs->isdense){
        uint32_t vertex;
        Bitmap_Clean(vs->bitArray, vs->size_bitarray);
        for (size_t i = 0; i < vs->sizeSparse; i++){
            vertex = vs->sparseArray[i];
            Bitmap_Set(vs->bitArray, vertex);
        }
        vs->isdense = true;
    }
    
};

bool Vertexset_TransformToSparse(Vertexset* vs) {
    // only do, if it is dense
    if (vs->isdense){
        unsigned long long word;
        size_t countAdded = 0;
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
        vs->sizeSparse = countAdded;
        vs->isdense = false;    
    }
}

void Vertexset_Allreduce(Vertexset* vs, int VERTEXSET_OP){
    // no other variant is implemented yet
    // use only sparse communication
    Vertexset_TransformToSparse(vs);
    Vertexset_Allreduce_Pure(vs, VERTEXSET_OP);
}


void Vertexset_Allreduce_Pure(Vertexset* vs, int VERTEXSET_OPERATION){
    assert(VERTEXSET_OPERATION == VERTEXSET_OR); //no other version is implemented
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

void Vertexset_Allreduce_Exact_Halfing(Vertexset* vs, int VERTEXSET_OPERATION) {
    assert(VERTEXSET_OPERATION == VERTEXSET_OR); // no other operator implemented
    assert((vs->mpi_size & (vs->mpi_size - 1)) == 0); // communicator must be size of 2^k
    
    // find ciritical size
    size_t criticalSize = vs->sizeCrit;

    // prepare, if vertexset is sparse
    if (!(vs->isdense)) {
        Bitmap_Clean(vs->bitArray, vs->size_bitarray);
        int count=0;
        int32_t vert;
        for (size_t i = 0; i < vs->sizeSparse; i++) {
            vert = vs->sparseArray[i];
            if (count <= criticalSize/2){
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
        if (count< criticalSize/2) {
            // only neccary, when sparse variant is needed   
            uint32_t *temp;
            temp = vs->sparseArray;
            vs->sparseArray = vs->sparseBuffer;
            vs->sparseBuffer = temp;
            temp = NULL;        
        }
    }

    // prepare vertexset, if it is dense
    if (vs->isdense){
        vs->sizeSparse = criticalSize + 1;
    }



    // find block indices
    int blockIdx[vs->maxsize + 1];
    blockIdx[0] = 0;
    for (int i = 1; i < vs->mpi_size + 1; i++) {
        blockIdx[i] = vs->size_bitarray * i / vs->mpi_size;
    }
    
    // do reduce_scatter
    int commNeighbor;
    int startBlock;
    int startIndex;
    int nElements;
    int tag;
    MPI_Status status;
    MPI_Request req;
    int recvCount;
    for (int shift = vs->mpi_size / 2 ; shift >= 1; shift/=2) {
        // do bitflip with LOR to find neighbor
        commNeighbor = vs->mpi_rank ^ shift;

        // decide, which strategy to choose from
        criticalSize /= 2;
        if (vs->sizeSparse < criticalSize) {
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

            MPI_Recv(vs->sparseBuffer, recvCount, MPI_INT32_T, commNeighbor, 100, vs->MPI_COMM, MPI_STATUS_IGNORE);

            MPI_Wait(&req, MPI_STATUS_IGNORE);
            int count=vs->sizeSparse;
            int32_t vert;
            // do reduction (append non dublicates and also update bitarray)
            for (size_t i = 0; i < recvCount; i++)  {
                vert = vs->sparseBuffer[i];
                if (count <= criticalSize/2 || shift==1) {
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

        } else if (tag==200) {
            // dense array is recieved          
            // get length of message
            MPI_Get_count(&status, MPI_UNSIGNED_LONG_LONG, &recvCount);
            
            // Recieve into buffer
            MPI_Recv(vs->bitBuffer, recvCount, MPI_UNSIGNED_LONG_LONG, commNeighbor, 200, vs->MPI_COMM, MPI_STATUS_IGNORE);
            MPI_Wait(&req, MPI_STATUS_IGNORE);

            // do reduction
            startBlock = (vs->mpi_rank/shift) * shift;
            startIndex = blockIdx[startBlock];
            assert(recvCount == blockIdx[startBlock + shift]-startIndex);
            for (int i = startIndex; i < blockIdx[startBlock + shift]; i++) {
                vs->bitArray[i] |= vs->bitBuffer[i-startIndex];
            }
            
            // set count, that is over critical size
            vs->sizeSparse = criticalSize + 1;
            
        } else {
            assert("recieved wrong tag" && false);
        }      
    }
    
    // In case only sparse communication was performed, we don't need allgather
    if (vs->sizeSparse <= criticalSize) {
        vs->isdense = false;
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
        
        //if (shift==2)
        //printf("nElements: %d\n", nElements);
        MPI_Isend(vs->bitArray + startIndex, nElements, MPI_UNSIGNED_LONG_LONG, commNeighbor, 101, vs->MPI_COMM, &req);
        
        // Check, if recieved message is in dense format
        MPI_Probe(commNeighbor, 101, vs->MPI_COMM, &status);
        MPI_Get_count(&status, MPI_LONG_LONG, &recvCount);
        
        // Recieve into buffer
        MPI_Recv(vs->bitBuffer, recvCount, MPI_LONG_LONG, commNeighbor, 101, vs->MPI_COMM, MPI_STATUS_IGNORE);
        
        // Send buffer can be used again
        MPI_Wait(&req, MPI_STATUS_IGNORE);
        
        // do fill in
        startBlock = (commNeighbor/shift) * shift;
        startIndex = blockIdx[startBlock];
        assert(recvCount == blockIdx[startBlock + shift]-startIndex);
        for (int i = startIndex; i < blockIdx[startBlock + shift]; i++) {
            vs->bitArray[i] = vs->bitBuffer[i-startIndex];
        }  
    }

    // must be dense, if allgather was needed
    vs->isdense = true;
    

};


void Vertexset_Allreduce_Approximate_Halfing(Vertexset* vs, int VERTEXSET_OPERATION) {
    assert(VERTEXSET_OPERATION == VERTEXSET_OR); //no other version is implemented

    int criticalSize = vs->sizeCrit;
    // transform into dense, if critical size is hit
    if (vs->sizeSparse >= criticalSize) {
        Vertexset_TransformToDense(vs);
    }

    // set size to maximum, if dense
    if (vs->isdense) {
        vs->sizeSparse = criticalSize + 1;
    }
    
    // calculate start and end of each block
    int blockIndices_input[vs->mpi_size + 1];
    int blockIndices_buffer[vs->mpi_size + 1];
    int block_nElements_input[vs->mpi_size];
    int block_nElements_buffer[vs->mpi_size];
    
    // calculate indices for input (same for every thread)
    blockIndices_input[0] = 0;
    for (size_t i = 1; i < vs->mpi_size + 1; i++) {
        blockIndices_input[i] = vs->size_bitarray * i / vs->mpi_size;
        block_nElements_input[i-1] = blockIndices_input[i] - blockIndices_input[i-1];
    }

    // calculate indices for R buffer
    blockIndices_buffer[0] = 0;
    for (size_t i = 1; i < vs->mpi_size + 1; i++) {
        blockIndices_buffer[i] = blockIndices_buffer[i - 1] + block_nElements_input[(i - 1 + vs->mpi_rank) % vs->mpi_size ];
        block_nElements_buffer[i-1] = blockIndices_buffer[i] - blockIndices_buffer[i-1];
    }
    
    // prepare, if vertexset is sparse
    if (!(vs->isdense)) {
        Bitmap_Clean(vs->bitArray, vs->size_bitarray);
        int count=0;
        int32_t vert;
        for (size_t i = 0; i < vs->sizeSparse; i++) {
            vert = vs->sparseArray[i];
            if (count <= criticalSize/2){
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
        if (count < criticalSize/2) {
            // only neccary, when sparse variant is needed   
            uint32_t *temp;
            temp = vs->sparseArray;
            vs->sparseArray = vs->sparseBuffer;
            vs->sparseBuffer = temp;
            temp = NULL;        
        }
    }

    // filling shifted dense buffer
    int32_t indexShift = blockIndices_input[vs->mpi_rank];
    for (size_t i = 0; i < vs->size_bitarray; i++) {
        vs->bitBuffer[i] = vs->bitArray[(i+indexShift) % vs->size_bitarray];
    }
    // Attention: from now on bitBuffer holds data and bitArray is used as buffer



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
    int nElementsDense;
    int startIndex;
    int tag;
    MPI_Status status;
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
        
        // find start indices of dense send
        startIndex = blockIndices_buffer[shift];

        // find the size of dense send
        nElementsDense = blockIndices_buffer[shift_old] - startIndex;

        // find critical size
        criticalSize = nElementsDense * sizeof(unsigned long long) / sizeof(uint32_t);
        
        // decide, which variant should be send
        // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (false && vs->sizeSparse < criticalSize) {
            // send sparse array
            MPI_Send(vs->sparseArray, vs->sizeSparse, MPI_INT32_T, sendNeighbor, 100, vs->MPI_COMM);
        } else {
            // send dense array
            MPI_Send(vs->bitBuffer + startIndex, nElementsDense, MPI_UNSIGNED_LONG_LONG, sendNeighbor, 200, vs->MPI_COMM);
        }

        
        // get info, if recieved message is in dense or sparse format
        MPI_Probe(recvNeighbor, MPI_ANY_TAG, vs->MPI_COMM, &status);
        tag = status.MPI_TAG;
        
        if (tag==100) {
            // sparse array is recieved
            // get length of message
            MPI_Get_count(&status, MPI_INT32_T, &recvCount);
            
            MPI_Recv(vs->sparseBuffer, recvCount, MPI_INT32_T, recvNeighbor, 100, vs->MPI_COMM, MPI_STATUS_IGNORE);
            int count=vs->sizeSparse;
            int32_t vert;
            // do reduction (append non dublicates and also update bitarray)
            for (size_t i = 0; i < recvCount; i++)  {
                vert = vs->sparseBuffer[i];
                
                // estimate critical size for next iteration and only append, when necessary
                if (count <= criticalSize/2 || shift==1) {
                    // set bitmap, and add to sparse,if needed
                    if (!Bitmap_Test_Shifted(vs->bitBuffer, vert, indexShift, vs->size_bitarray)) {
                        // save in sparse array until its needed in next round
                        Bitmap_Set_Shifted(vs->bitBuffer, vert, indexShift, vs->size_bitarray);
                        vs->sparseArray[count++] = vert;                        
                    }
                } else {
                    Bitmap_Set_Shifted(vs->bitBuffer, vert, indexShift, vs->size_bitarray);
                }
            }

            // set count
            vs->sizeSparse = count;
            
        } else if (tag==200) {
            // dense array is recieved          
            // get length of message (can be done exxplicitly)
            recvCount = blockIndices_buffer[shift_old - shift];

            // Recieve into buffer. Attention: Here is bitArray used as the buffer!!!
            MPI_Recv(vs->bitArray, recvCount, MPI_UNSIGNED_LONG_LONG, recvNeighbor, 200, vs->MPI_COMM, MPI_STATUS_IGNORE);

            // do reduction           
            for (int i = 0; i < recvCount; i++) {
                vs->bitBuffer[i] |= vs->bitArray[i];
            }
            
            // set count, that is over critical size
            vs->sizeSparse = criticalSize + 1;
            
        } else {
            assert("recieved wrong tag" && false);
        }
    }
    
    // decide, if allgather is needed.
    // common ground: originial common ground divided by (2^iterations)
    if (vs->sizeSparse <= vs->sizeCrit >> iterations) {
        vs->isdense = false;
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
        
        
        MPI_Send(vs->bitBuffer, blockIndices_buffer[shift_next - shift], MPI_UNSIGNED_LONG_LONG, sendNeighbor, 101, vs->MPI_COMM);
        

        recvCount = blockIndices_buffer[shift_next] - blockIndices_buffer[shift];
        
        // Recieve into buffer
        MPI_Recv(vs->bitBuffer + blockIndices_buffer[shift], recvCount, MPI_LONG_LONG, recvNeighbor, 101, vs->MPI_COMM, MPI_STATUS_IGNORE);
    }
    
    // transform shifted buffer back
    for (size_t i = 0; i < vs->size_bitarray; i++) {
        vs->bitArray[(i+indexShift) % vs->size_bitarray] = vs->bitBuffer[i];
    }

    vs->isdense = true;
};

void Vertexset_Allreduce_Dynamic(Vertexset* vs, int VERTEXSET_OPERATION){
    // Check if OR operation is used (other is not implemented yet)
    assert(VERTEXSET_OPERATION == VERTEXSET_OR);

    // Check if size of Communicatior is power of 2
    assert((vs->mpi_size & (vs->mpi_size - 1)) == 0);
    
    // all processes should have a sparse representation
    if (vs->isdense) {
        Vertexset_TransformToSparse(vs);
    }

    // calculate numver of communication rounds
    int max_iterations = Vertexset_Log2_floor(vs->mpi_size);
    
    // find the sizes of all the other ranks
    int size_int = (int) vs->sizeSparse;
    MPI_Allgather(&size_int, 1, MPI_INT, vs->sizesAll, 1, MPI_INT, vs->MPI_COMM);
    
    // find the iteration round, when commucication should switch to dense variant
    int crit_iteration = max_iterations; // = Vertexset_Find_critical_iteration(vs);
    
    // inits for actual communication
    int rank_diff, recv_rank, send_rank, round;
    int *buffer;
    buffer = (int*) malloc(sizeof(int)*vs->mpi_size);
    rank_diff = 1;
    
    // communication rounds sparse
    for (round = 0; round < crit_iteration; round++) {
        send_rank = (vs->mpi_rank + rank_diff) % vs->mpi_size;
        recv_rank = (vs->mpi_rank - rank_diff + vs->mpi_size) % vs->mpi_size;
        MPI_Sendrecv(vs->sparseArray,                           vs->sizesAll[vs->mpi_rank], MPI_UINT32_T, send_rank, round,
                     vs->sparseArray + vs->sizesAll[recv_rank], vs->sizesAll[recv_rank],    MPI_UINT32_T, recv_rank, round,
                    vs->MPI_COMM, MPI_STATUS_IGNORE);
        Vertexset_Update_Sizes(&vs->sizesAll, &buffer, rank_diff, vs->mpi_size);        
        
        // double rank shift
        rank_diff <<= 1;
    }
    vs->sizeSparse = vs->sizesAll[vs->mpi_rank];

    // transform to dense
    if (round < max_iterations) {
        Vertexset_TransformToDense(vs);
    }
    

    // communication rounds dense
    for (round = crit_iteration; round < max_iterations; round++) {
        send_rank = (vs->mpi_rank + rank_diff) % vs->mpi_size;
        recv_rank = (vs->mpi_rank - rank_diff + vs->mpi_size) % vs->mpi_size;
        MPI_Sendrecv(vs->bitArray, vs->size_bitarray, MPI_UINT64_T, send_rank, round,
                     vs->bitBuffer, vs->size_bitarray, MPI_UINT64_T, recv_rank, round,
                    vs->MPI_COMM, MPI_STATUS_IGNORE);

        // perform logical OR
        for (size_t i = 0; i < vs->size_bitarray; i++) {
            vs->bitArray[i] |= vs->bitBuffer[i];
        }
        
        // double the rank difference
        rank_diff <<= 1;
    }

    free(buffer);
};


void Vertexset_PrintSet(Vertexset* vs){
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


void Vertexset_Deinit(Vertexset* vs){
    free(vs->bitArray);
    free(vs->sparseArray);
};



