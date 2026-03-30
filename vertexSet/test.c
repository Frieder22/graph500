#include "vertexSet.h"
#include "vertexSetIterator.h"
#include "mpi.h"
#include "testStructure.h"
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
char *N_Nodes, *N_TasksPerNode;

extern enum Testcase Testcase;

// from https://www.geeksforgeeks.org/qsort-function-in-c/
// comparator for qsort to sort in ascending order
int comp(const void *a, const void *b) {
    return (*(int *)a - *(int *)b);
};

void printTestInfo(int testnumber){
    if(rank==0){
        printf("\n");
        switch (testnumber) {
        case ITERATOR:
            printf("Testing functionality of Iterator functions.\n");
            break;
        case EXACT_HALFING:
            printf("Testing functionality of Exact Halfing Reduce.\n");
            break;
        case APROXIMATE_HALFING:
            printf("Testing functionality of Aproximate Halfing Reduce.\n");
            break;
        case NAIVE:
            printf("Testing functionality of Naive Reduce.\n");
            break;
        case DENSE:
            printf("Testing functionality of Dense Reduce.\n");
            break;
        case VERTEXSET_TEST:
            printf("Testing functionality of Vertexset functions.\n");
            break;
        case ALLGATHER:
            printf("Testing functionality of Allgather.\n");
            break;     

        default:
            printf(ANSI_RED "Print function doesn't know this testcase!\n" ANSI_RESET);
            break;
        }
        if (N_Nodes) {
            printf("Using %s nodes with %s processes on each node (%d ranks).\n", N_Nodes, N_TasksPerNode, size);
        } else {
            printf("Using %d processes on one machine.\n", size);
        }
        printf("\n");
    }
}

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
    Vertexset vs, control;
    Vertexset_Init(&vs, 500, MPI_COMM_WORLD);
    Vertexset_Init(&control, 500, MPI_COMM_WORLD);

    uint32_t addedVerts[] = {10, 7, 12, 32, 0, 333, 499};
    int n = sizeof(addedVerts) / sizeof(addedVerts[0]);

    for (int i = 0; i < n; i++) {
        Vertexset_Add(&vs, addedVerts[i]);
        Vertexset_Add(&control, addedVerts[i]);
    }

    Vertexset_TransformToDense(&vs);
    vertexSetIterator it;
    vertexSetIterator_Init(&it, &vs);


    int count = 0;
    uint32_t element;
    printf("Testing Has_next() and Next()...\n");
    while (vertexSetIterator_Has_next(&it)) {
        element = vertexSetIterator_Next(&it);
        assert(Vertexset_Contains(&control, element));
        count++;
    }
    assert(count == n);

    vertexSetIterator_Reset(&it);
    count = 0;
    printf("Testing Reset()...\n");
    while (vertexSetIterator_Has_next(&it)) {
        element = vertexSetIterator_Next(&it);
        assert(Vertexset_Contains(&control, element));
        count++;
    }
    assert(count == n);

    Vertexset_Deinit(&vs);
    Vertexset_Deinit(&control);

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
    isCorrect &= test_Allreduce_single(reduceFunc, 100, 0.02);
    isCorrect &= test_Allreduce_single(reduceFunc, 1, 0.5);
    isCorrect &= test_Allreduce_single(reduceFunc, 63, 0.3);
    isCorrect &= test_Allreduce_single(reduceFunc, 64, 0.3);
    isCorrect &= test_Allreduce_single(reduceFunc, 65, 0.3);
    isCorrect &= test_Allreduce_single(reduceFunc, 1<<25, 0.3);  
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

bool test_Vertexset_sparse(){
    bool isCorrect = true;
    int maxsize = 500;
    const int n1 = 15;
    const int n2 = 7;
    srand(42);

    // Init of VS
    Vertexset vs;
    Vertexset_Init(&vs, maxsize, MPI_COMM_WORLD);
    
    // reference buffer
    uint32_t added[n1 + n2];
    uint32_t vert;

    // PART 1
    // Add random verts
    for (size_t i = 0; i < n1; i++) {
        vert = rand() % maxsize;

        // add in reference
        added[i] = vert;
        // add in VS
        Vertexset_Add(&vs, vert);
    }
    
    // Check for correct insertions
    for (size_t i = 0; i < n1; i++) {
        isCorrect &= (added[i] == vs.sparseArray[i]);
    }
    isCorrect &= (n1 == vs.sizeSparse);  // correct size

    // PART 2
    // Add more random verts
    for (size_t i = 0; i < n2; i++) {
        vert = rand() % maxsize;

        // add in reference
        added[n1 + i] = vert;
        // add in VS
        Vertexset_Add(&vs, vert);
    }

    // Check for correct insertions
    for (size_t i = 0; i < n1 + n2; i++) {
        isCorrect &= (added[i] == vs.sparseArray[i]);
    }
    isCorrect &= (n1 + n2 == vs.sizeSparse);  // correct size

    // Part 3
    // Check, if contains finds all inserted values
    for (size_t i = 0; i < n1 + n2; i++) {
        isCorrect &= Vertexset_Contains(&vs, added[i]);
    }
    
    // Check, if contains finds no missing values
    uint32_t complement[maxsize];
    for (size_t i = 0; i < maxsize; i++){
        complement[i] = i;
    }
    for (size_t i = 0; i < n1+n2; i++) {
        complement[added[i]] = maxsize; // mark already added vals
    }
    for (size_t i = 0; i < maxsize; i++) {
        vert = complement[i];
        if (vert != maxsize) {
            isCorrect &= (!Vertexset_Contains(&vs, vert));
        }
    }
    
    // Part 4
    // Check Cleanup
    Vertexset_Clean(&vs);
    isCorrect &= (vs.sizeSparse == 0);
    
    // Free Vertexset
    Vertexset_Deinit(&vs);

    // result output
    if (isCorrect) {
        printf(ANSI_GREEN "Test Vertexset_sparse SUCCESSFUL\n" ANSI_RESET);
        printf("--------------------------------------------\n\n");
    } else {
        printf(ANSI_RED "Test Vertexset_sparse FAILED\n" ANSI_RESET );
        printf("--------------------------------------------\n\n");
    }
}

bool test_Vertexset_dense(){
    bool isCorrect = true;
    int maxsize = 257;

    // Init of VS
    Vertexset vs;
    Vertexset_Init(&vs, maxsize, MPI_COMM_WORLD);

    // set to dense
    vs.isdense = 1;

    // clean VS
    Vertexset_Clean(&vs);
    
    // PART 1
    // Check Bitarray to be empty
    int nwords = vs.size_bitarray;
    for (size_t i = 0; i < nwords; i++) {
        isCorrect &= (vs.bitArray[i] == 0ULL);
    }

    // PART 2
    // Add a few numbers in VS, that are in the first word
    uint32_t verts[] = {1, 20, 3, 7, 63};

    unsigned long long word = 0ULL;
    for (size_t i = 0; i < sizeof(verts)/sizeof(uint32_t); i++) {
        // find correct word  
        word += 1ULL << verts[i];
        //add to VS
        Vertexset_Add(&vs, verts[i]);
    }
    isCorrect &= (word == vs.bitArray[0]);

    // PART 3
    // Add a few numbers in VS, that are in the second word
    uint32_t verts2[] = {65, 80, 100, 77, 111};

    unsigned long long word2 = 0ULL;
    for (size_t i = 0; i < sizeof(verts2)/sizeof(uint32_t); i++) {
        // find correct word
        word2 += 1ULL << (verts2[i] - ulong_bits);
        //add to VS
        Vertexset_Add(&vs, verts2[i]);
    }
    isCorrect &= (word2 == vs.bitArray[1]);
    // Check, if word 1 is unaltered
    isCorrect &= (word == vs.bitArray[0]);

    // PART 4:
    // Check correct allocation size for bitarray
    if (maxsize % ulong_bits == 0) {
        isCorrect &= (vs.size_bitarray == maxsize/ulong_bits);
    } else {
        isCorrect &= (vs.size_bitarray == maxsize/ulong_bits + 1);
    }

    // PART 5
    // check for correct contains functionality
    Vertexset_Add(&vs, maxsize - 1); // add edge cases
    Vertexset_Add(&vs, 0); // add edge cases
    bool isInSet;
    for (uint32_t vert = 0; vert < maxsize; vert++) {
        isInSet = false;
        for (size_t i = 0; i < sizeof(verts)/sizeof(uint32_t); i++) {
            isInSet |= (verts[i] == vert);
        }
        for (size_t i = 0; i < sizeof(verts2)/sizeof(uint32_t); i++) {
            isInSet |= (verts2[i] == vert);
        }
        isInSet |= (vert == 0);
        isInSet |= (vert == maxsize - 1);
        isCorrect &= (isInSet == Vertexset_Contains(&vs, vert));
    }
    
    // PART 6
    // Check again Bitarray to be empty
    Vertexset_Clean(&vs);
    for (size_t i = 0; i < nwords; i++) {
        isCorrect &= (vs.bitArray[i] == 0ULL);
    }

    // result output
    if (isCorrect) {
        printf(ANSI_GREEN "Test Vertexset_dense SUCCESSFUL\n" ANSI_RESET);
        printf("--------------------------------------------\n\n");
    } else {
        printf(ANSI_RED "Test Vertexset_dense FAILED\n" ANSI_RESET );
        printf("--------------------------------------------\n\n");
    }
}

bool test_Allgather_single(void (*allgatherFunc) (Vertexset*), int setSize, int filling){
    bool isCorrect = true;
    bool verbose = false;
    bool debug = false;

    if (debug && rank==0) printf(ANSI_RED"DEBUG mode!!\n"ANSI_RESET);

    int insertionsTotal = filling;

    Vertexset vs;
    Vertexset_Init(&vs, setSize, MPI_COMM_WORLD);
    uint32_t* control;
    control = (uint32_t*) malloc(insertionsTotal*sizeof(uint32_t));

    assert(insertionsTotal <= vs.sizeCrit);

    int num;
    if (debug) {
        // Fill arrays with only own rank ID
        insertionsTotal = size;
        for (size_t i = 0; i < size; i++) {
            num = i%size;
            control[i] = num;
            if (rank == num) {
                Vertexset_Add(&vs, num);
            }
        }
    } else {
        // Fill arrays with random numbers
        for (size_t i = 0; i < insertionsTotal; i++) {
            num = rand()%setSize;
            control[i] = num;
            if (i%size==rank) {
                Vertexset_Add(&vs, num);
            }
        }
    }
    

    // perform Allgather
    (*allgatherFunc) (&vs);

    // sort array and controll 
    qsort(vs.sparseArray, insertionsTotal, sizeof(uint32_t), comp);
    qsort(control, insertionsTotal, sizeof(uint32_t), comp);

    // compare elements
    for (size_t i = 0; i < insertionsTotal; i++){
        if (verbose && rank==0) {
            printf("arr, control: %d, %d\n", vs.sparseArray[i], control[i]);
        }
        
        if (vs.sparseArray[i]!= control[i]) {
            isCorrect=false;
        }
    }

    //compare size of lists
    if (vs.sizeSparse != insertionsTotal) {
        isCorrect = false;
        if (rank==0) {
            printf("Size doesn't match up!!!\n");
        }        
    }
    
    
    // Check if other ranks had problems
    MPI_Allreduce(MPI_IN_PLACE, &isCorrect, 1, MPI_C_BOOL, MPI_LAND, MPI_COMM_WORLD);
    
    // Priont resulrt
    if (rank==0) {
        if (isCorrect) {
            printf(ANSI_GREEN "Test Allgather SUCCESSFUL\n" ANSI_RESET);
            printf("--------------------------------------------\n\n");
        } else {
            printf(ANSI_RED "Test Allgather FAILED\n" ANSI_RESET );
            printf("--------------------------------------------\n\n");
        }
    }

    return isCorrect;
    
}


bool test_Allgather_batch(void (*allgatherFunc) (Vertexset*)){
    bool isCorrect = true;
    isCorrect &= test_Allgather_single(allgatherFunc, 500000, size*2+4);
    isCorrect &= test_Allgather_single(allgatherFunc, 500000, size + 7);
    isCorrect &= test_Allgather_single(allgatherFunc, 500000, size + 1);
    isCorrect &= test_Allgather_single(allgatherFunc, 1000000, 10*size + 3);

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
    // get the number of nodes via the enviroment viarable (provided by slurm)
    N_Nodes = (char*) getenv("SLURM_NNODES");
    // same for tasks per node
    N_TasksPerNode = (char*) getenv("SLURM_NTASKS_PER_NODE");

    printTestInfo(testNumber);

    // test_it: test iterators
    switch (testNumber) {
        case ITERATOR: // Test iterator
            test_Iterator_dense();
            test_Iterator_sparse();
            break;
        
        case EXACT_HALFING: // Test Exact halfing
            assert((size & (size - 1)) == 0); // only works for ranks = 2^k
            test_Allreduce_batch(Vertexset_Allreduce_Exact_Halfing);
            break;

        case APROXIMATE_HALFING: // Test Approximate halfing
            test_Allreduce_batch(Vertexset_Allreduce_Approximate_Halfing);
            break;

        case NAIVE: // Test allreduce using ring topology
            test_Allreduce_batch(Vertexset_Allreduce_Naive);
            break;

        case DENSE: // Test dense allreduce
            test_Allreduce_batch(Vertexset_Allreduce_Dense);
            break;

        case VERTEXSET_TEST: //Test functionality of vertexset
            test_Vertexset_sparse();
            test_Vertexset_dense();
            break;

        case ALLGATHER:
            test_Allgather_batch(Vertexset_Allgather);
        default:
            break;
    }

    MPI_Finalize();
    return 0;
}

