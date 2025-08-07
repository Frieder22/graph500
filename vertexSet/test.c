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

// from https://www.geeksforgeeks.org/qsort-function-in-c/
// comparator for qsort to sort in ascending order
int comp(const void *a, const void *b) {
    return (*(int *)a - *(int *)b);
};

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

bool test_Allreduce_single(void (*reduceFunc) (Vertexset*, int), int setSize, float filling){
    bool isCorrect = true;
    bool verbose = false;

    int insertions = (int) setSize * filling;

    Vertexset vs, control;
    Vertexset_Init(&vs, setSize, MPI_COMM_WORLD);
    Vertexset_Init(&control, setSize, MPI_COMM_WORLD);
    
    Vertexset_TransformToDense(&control);
    
    // all ranks have same seed
    srand(2);
    int num;
    for (size_t i = 0; i < insertions; i++) {
        num = rand()%setSize;
        Vertexset_Add(&control, num);
        if (i%size==rank) {
            Vertexset_Add(&vs, num);
        }
    }
    
    //Allreduce_Dense only works with dense arrays
    if(reduceFunc == Vertexset_Allreduce_Dense){
        Vertexset_TransformToDense(&vs);
    }

    // perform allreduce
    (*reduceFunc) (&vs, VERTEXSET_OR);

    vertexSetIterator it, it_control;
    vertexSetIterator_Init(&it, &vs);
    vertexSetIterator_Init(&it_control, &control);

    // Check if vertex set is subset of reference set
    int32_t val;
    while (vertexSetIterator_Has_next(&it)) {
        val = vertexSetIterator_Next(&it);

        if (!Vertexset_Contains(&control, val)) { // value can not be found in reference
            isCorrect = false;
            break;
        } 
    }

    if (verbose)
    if (rank==1){
        vertexSetIterator_Reset(&it);
        printf("rank %d: ", rank);
        while (vertexSetIterator_Has_next(&it)) {
            val = vertexSetIterator_Next(&it);
            printf("%d ", val);
        }
        printf("\n");
    }


    // check if reference set is subset of vertexset
    while (vertexSetIterator_Has_next(&it_control)) {
        val = vertexSetIterator_Next(&it_control);
        if (!Vertexset_Contains(&vs, val)) { // value can not be found in reference
            isCorrect = false;
            break;
        } 
    }

    if(verbose)
    if(rank==0){
        printf("\ncontrol: ");
        vertexSetIterator_Reset(&it_control);
        while (vertexSetIterator_Has_next(&it_control)) {
            val = vertexSetIterator_Next(&it_control);
            printf("%d ", val);
        }
        printf("\n");
    }
    
    if(verbose)
    if (isCorrect) {
        printf("Rank %d is correct!\n", rank);
    }
    

    // Check if other ranks had problems
    MPI_Allreduce(MPI_IN_PLACE, &isCorrect, 1, MPI_C_BOOL, MPI_LAND, MPI_COMM_WORLD);

    if (rank==0) {
        if (isCorrect) {
            printf(ANSI_GREEN "Test Allreduce SUCCESSFUL\n" ANSI_RESET);
            printf("--------------------------------------------\n\n");
        } else {
            printf(ANSI_RED "Test Allreduce FAILED\n" ANSI_RESET );
            printf("--------------------------------------------\n\n");
        }        
    }

    return isCorrect;
}

bool test_Allreduce_batch(void (*reduceFunc) (Vertexset*, int)){
    bool isCorrect = true;
    isCorrect &= test_Allreduce_single(reduceFunc, 1000000, 0.01);
    isCorrect &= test_Allreduce_single(reduceFunc, 100, 1.0);
    isCorrect &= test_Allreduce_single(reduceFunc, 1, 0.5);
    isCorrect &= test_Allreduce_single(reduceFunc, 63, 0.3);
    isCorrect &= test_Allreduce_single(reduceFunc, 64, 0.3);
    isCorrect &= test_Allreduce_single(reduceFunc, 65, 0.3);
    //isCorrect &= test_Allreduce_single(reduceFunc, 1<<25, 0.3);  
    if (rank==0) {
        if (isCorrect) {
            printf(ANSI_GREEN"\n--------------------------------------------\n");
            printf("|          All tests SUCCESSFUL            |\n");
            printf("--------------------------------------------\n\n"ANSI_RESET);
        } else {
            printf(ANSI_RED"\n--------------------------------------------\n");
            printf("|           Some test(s) FAILED            |\n");
            printf("--------------------------------------------\n\n"ANSI_RESET);
        }        
    } 
    return isCorrect;
}



int main(int argc, char *argv[]){
    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    assert(argc == 2);
    int testNumber = atoi(argv[1]);
    
    // test_it: test iterators
    switch (testNumber) {
        case 0: // Test iterator
            test_Iterator_dense();
            test_Iterator_sparse();
            break;
        
        case 1: // Test Exact halfing
            assert((size & (size - 1)) == 0); // only works for ranks = 2^k
            test_Allreduce_batch(Vertexset_Allreduce_Exact_Halfing);
            break;

        case 2: // Test Approximate halfing
            test_Allreduce_batch(Vertexset_Allreduce_Approximate_Halfing);
            break;

        case 3: // Test allreduce using ring topology
            test_Allreduce_batch(Vertexset_Allreduce_Ring_Comm);
            break;

        case 4: // Test dense allreduce
            test_Allreduce_batch(Vertexset_Allreduce_Dense);
            break;

        default:
            break;
    }

    MPI_Finalize();
    return 0;
}

