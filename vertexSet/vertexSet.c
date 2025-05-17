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
    size_t size_bitarray = (maxsize + 1) / sizeof(unsigned long long);
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
        vs->sparseArray[vs->sizeSparse] = vertex;
        vs->sizeSparse++;
    }
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
    vs->isdense = false;
    vs->sizeSparse = 0;
};


bool Vertexset_TransformToDense(Vertexset* vs) {
    // only do, if it's sparse
    if (!vs->isdense){
        Bitmap_Clean(vs->bitArray, vs->size_bitarray);
        for (size_t i = 0; i < vs->sizeSparse; i++){
            uint32_t vertex = vs->sparseArray[i];
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
    if (vs->isdense) {
        MPI_Allreduce(MPI_IN_PLACE, vs->bitArray, vs->size_bitarray, MPI_UINT64_T, MPI_BOR, MPI_COMM_WORLD);
    } else {
        // get size information from other ranks
        int size_int = (int) vs->sizeSparse;
        MPI_Allgather(&size_int, 1, MPI_INT, vs->sizesAll, 1, MPI_INT, vs->MPI_COMM);

        // calculate displacements with exxclusive scan
        vs->displ[0] = 0;
        for (int i = 1; i < vs->mpi_size + 1; i++){
            vs->displ[i] = vs->displ[i-1] + vs->sizesAll[i-1];
        }

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

        // set size
        vs->sizeSparse = vs->displ[vs->mpi_size];
    }
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



