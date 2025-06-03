#include "vertexSet.h"
#include "vertexSetIterator.h"
#include "mpi.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>


// for colored text output
// from https://stackoverflow.com/questions/3219393/stdlib-and-colored-output-in-c
#define ANSI_RED     "\x1b[31m"
#define ANSI_GREEN   "\x1b[32m"
#define ANSI_YELLOW  "\x1b[33m"
#define ANSI_BLUE    "\x1b[34m"
#define ANSI_MAGENTA "\x1b[35m"
#define ANSI_CYAN    "\x1b[36m"
#define ANSI_RESET   "\x1b[0m"

int rank, size;


void test_Iterator_sparse(){
    printf("Testing Iterator sparse...\n");
    Vertexset vs;
    Vertexset_Init(&vs, 500, MPI_COMM_WORLD);

    uint32_t addedVerts[] = {10, 65, 6,0, 9, 7, 499};
    int n = sizeof(addedVerts) / sizeof(addedVerts[0]);


    for (int i = 0; i < n; i++) {
        Vertexset_Add(&vs, addedVerts[i]);
    }

    vertexSetIterator it;
    vertexSetIterator_Init(&it, &vs);

    int count = 0;
    uint32_t element;
    printf("Testing Has_next() and Next()...\n");
    while (vertexSetIterator_Has_next(&it)) {
        element = vertexSetIterator_Next(&it);
        assert(element == addedVerts[count]);
        count++;
    }
    assert(count == n);

    vertexSetIterator_Reset(&it);
    count = 0;
    printf("Testing Reset()...\n");
    while (vertexSetIterator_Has_next(&it)) {
        element = vertexSetIterator_Next(&it);
        assert(element == addedVerts[count]);
        count++;
    }
    assert(count == n);


    Vertexset_Deinit(&vs);

    printf(ANSI_GREEN "Test Iterator Sparse SUCCESSFUL\n" ANSI_RESET);
    printf("--------------------------------------------\n\n");
};

// from https://www.geeksforgeeks.org/qsort-function-in-c/
int comp(const void *a, const void *b) {
    return (*(int *)a - *(int *)b);
};

void test_Iterator_dense(){
    printf("Testing Iterator dense...\n");
    Vertexset vs;
    Vertexset_Init(&vs, 500, MPI_COMM_WORLD);

    uint32_t addedVerts[] = {10, 7, 12, 32, 0, 333, 499};
    int n = sizeof(addedVerts) / sizeof(addedVerts[0]);

    for (int i = 0; i < n; i++) {
        Vertexset_Add(&vs, addedVerts[i]);
    }

    Vertexset_TransformToDense(&vs);

    vertexSetIterator it;
    vertexSetIterator_Init(&it, &vs);


    qsort(addedVerts, n, sizeof(uint32_t), comp);


    int count = 0;
    uint32_t element;
    printf("Testing Has_next() and Next()...\n");
    while (vertexSetIterator_Has_next(&it)) {
        element = vertexSetIterator_Next(&it);
        assert(element == addedVerts[count]);
        count++;
    }
    assert(count == n);

    vertexSetIterator_Reset(&it);
    count = 0;
    printf("Testing Reset()...\n");
    while (vertexSetIterator_Has_next(&it)) {
        element = vertexSetIterator_Next(&it);
        assert(element == addedVerts[count]);
        count++;
    }
    assert(count == n);

    Vertexset_Deinit(&vs);

    printf(ANSI_GREEN "Test Iterator Dense SUCCESSFUL\n" ANSI_RESET);
    printf("--------------------------------------------\n\n");
};

void test_Allreduce_Exactly_Halfing(){
    assert(size == 4); // only works correctly for 4 ranks
    bool isCorrect = true;

    Vertexset vs;
    Vertexset_Init(&vs, 10, MPI_COMM_WORLD);

    // every rank adds different amount of vertices
    for (size_t i = 0; i < rank; i++) {
        Vertexset_Add(&vs, rank + i);
    }
    // rank 0:  ()
    // rank 1:  (1)
    // rank 2:  (2, 3)
    // rank 3:  (3, 4, 5)
    // Reduced: (1, 2, 3, 4, 5)

    uint32_t reducedReference[] = {1,2,3,4,5};
    const int referenceSize = sizeof(reducedReference)/sizeof(uint32_t);

    Vertexset_Allreduce_Exact_Halfing(&vs, VERTEXSET_OR);

    vertexSetIterator it;
    vertexSetIterator_Init(&it, &vs);

    // Check if vertex set is subset of reference set
    uint32_t val;
    bool isInSet;
    while (vertexSetIterator_Has_next(&it)) {
        val = vertexSetIterator_Next(&it);
        isInSet = false;
        for (size_t i = 0; i < referenceSize; i++) {
            if (val == reducedReference[i]) {
                isInSet = true;
                break;
            }            
        }
        if (!isInSet) { // value can not be found in reference
            isCorrect = false;
            break;
        } 
    }

    // check if reference set is subset of vertexset
    uint32_t valRef;
    vertexSetIterator_Reset(&it);
    for (size_t i = 0; i < referenceSize; i++) {
        valRef = reducedReference[i];
        isInSet = false;
        while (vertexSetIterator_Has_next(&it)) {
            val = vertexSetIterator_Next(&it);
            if(val == valRef){
                isInSet = true;
                break;
            }
        }
        if (!isInSet) { // value can not be found in reference
            isCorrect = false;
            break;
        } 
    }
    

    // Check if other ranks had problems
    MPI_Allreduce(MPI_IN_PLACE, &isCorrect, 1, MPI_C_BOOL, MPI_LAND, MPI_COMM_WORLD);

    if (rank==0) {
        if (isCorrect) {
            printf(ANSI_GREEN "Test Allreduce_Exactly SUCCESSFUL\n" ANSI_RESET);
            printf("--------------------------------------------\n\n");
        } else {
            printf(ANSI_RED "Test Allreduce_Exactly FAILED\n" ANSI_RESET );
            printf("--------------------------------------------\n\n");
        }        
    }
}


int main(int argc, char *argv[]){
    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    assert(argc == 2);
    int testNumber = atoi(argv[1]);
    
    if (testNumber == 0){
        test_Iterator_dense();
        test_Iterator_sparse();
    }

    if (testNumber == 1){
        test_Allreduce_Exactly_Halfing();
    }

    MPI_Finalize();
    return 0;
}

