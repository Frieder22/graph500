#include "vertexSet.h"
#include "vertexSetIterator.h"
#include "mpi.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


// for colored text output
// from https://stackoverflow.com/questions/3219393/stdlib-and-colored-output-in-c
#define ANSI_RED     "\x1b[31m"
#define ANSI_GREEN   "\x1b[32m"
#define ANSI_YELLOW  "\x1b[33m"
#define ANSI_BLUE    "\x1b[34m"
#define ANSI_MAGENTA "\x1b[35m"
#define ANSI_CYAN    "\x1b[36m"
#define ANSI_RESET   "\x1b[0m"

#define OUTPUTPATH "../data/performanceData/"
#define N_DATAPOINTS 10

int rank, size;

struct performance {
    double timeSparse;
    double timeDense;
    double timeVertexSet;
} performance;


void performance_Allreduce(void (*reduceFunc) (Vertexset*, int)){
    bool isCorrect = true;
    bool verbose = false;


 
    // perform allreduce
    //(*reduceFunc) (&vs, VERTEXSET_OR);

    char filepath[100] = OUTPUTPATH;
    
    if (reduceFunc == Vertexset_Allreduce_Exact_Halfing) {
        char filename[] = "exact_Halfing";
        strcat(filepath, filename);
        strcat(filepath, ".txt");
    }

    if (reduceFunc == Vertexset_Allreduce_Approximate_Halfing) {
        char filename[] = "approx_Halfing";
        strcat(filepath, filename);
        strcat(filepath, ".txt");
    }
    
    FILE *f;
    if (rank==0){
        f = fopen(filepath, "w");
        if (f == NULL) {
            printf(ANSI_RED "File could not be opened!!\n" ANSI_RESET);
            return;
        }
        fprintf(f, "SetSize, VertexSet Time [s], Dense Time [s], Sparse Time [s]\n");
        printf("Results are written in %s\n", filepath);
    }
    int sizes[] = {20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000, 500000, 1000000, 2000000};
    int N = sizeof(sizes) / sizeof(int);

    int nDatapoints = 1;
    float filling = 0.1;
    int insertions;
    int32_t num;
    int setSize;

    Vertexset vs;
    
    size_t size_bitarray;
    unsigned long long *bitArray;
    
    int sparseSize, sparseSizeMax;
    int32_t *sparseArray, *sparseBuffer;
    int sizesAll[size];
    int displ[size + 1];

    double startTime;
    double endTime;

    // init Bitarray
    size_bitarray = (sizes[N-1] + (sizeof(unsigned long long)*8)) / (sizeof(unsigned long long)*8);
    bitArray = (unsigned long long*) malloc(size_bitarray * sizeof(unsigned long long));
    Bitmap_Clean(bitArray, size_bitarray);

    // init sparse array
    sparseSizeMax = (int)(sizes[N-1] * filling);
    sparseArray = (int32_t*) malloc(sparseSizeMax * sizeof(int32_t));
    sparseBuffer = (int32_t*) malloc(sparseSizeMax * sizeof(int32_t));

    assert(sparseArray != NULL);
    assert(sparseBuffer != NULL);
    
    srand(3);
    for (size_t i = 0; i < N; i++) {
        setSize = sizes[i];
        
        insertions = (int) (setSize * filling);
        for (size_t j = 0; j < N_DATAPOINTS; j++) {
            // init vertex set
            Vertexset_Init(&vs, setSize, MPI_COMM_WORLD);
            
            // fill arrays
            sparseSize = 0;
            Bitmap_Clean(bitArray, size_bitarray);
            for (size_t i = 0; i < insertions; i++) {
                num = rand()%setSize;
                if (i%size==rank) {
                    sparseArray[sparseSize++] = num;
                    Bitmap_Set(bitArray, num);
                    Vertexset_Add(&vs, num);
                }
            }

            // reduce Vertexset
            MPI_Barrier(MPI_COMM_WORLD);
            startTime = MPI_Wtime();
            (*reduceFunc) (&vs, VERTEXSET_OR);
            endTime = MPI_Wtime();
            performance.timeVertexSet = endTime - startTime;
            
            // reduce Bitmap
            MPI_Barrier(MPI_COMM_WORLD);
            startTime = MPI_Wtime();
            MPI_Allreduce(MPI_IN_PLACE, bitArray, size_bitarray, MPI_UINT64_T, MPI_BOR, MPI_COMM_WORLD);
            endTime = MPI_Wtime();
            performance.timeDense = endTime - startTime;
            
            // reduce Sparse Array
            MPI_Barrier(MPI_COMM_WORLD);
            startTime = MPI_Wtime();
            // get size information from other ranks
            MPI_Allgather(&sparseSize, 1, MPI_INT, sizesAll, 1, MPI_INT, MPI_COMM_WORLD);
            // calculate displacements with exclusive scan
            displ[0] = 0;
            for (int i = 1; i < size + 1; i++){
                displ[i] = displ[i-1] + sizesAll[i-1];
            }
            assert(displ[size] <= sparseSizeMax);
            MPI_Allgatherv(sparseArray, sparseSize, MPI_INT32_T, sparseBuffer, sizesAll, displ, MPI_INT32_T, MPI_COMM_WORLD);
            // put all numbers in (same) to perform logic OR (removes dublicates)
            for (size_t i = 0; i < displ[size]; i++) {
                Bitmap_Set(bitArray,sparseBuffer[i]);
            }
            endTime = MPI_Wtime();
            performance.timeSparse = endTime - startTime;

            // save results
            if (rank==0) {
                printf("true fillage: %d\n", displ[size]);
                fprintf(f, "%d,%lf,%lf,%lf\n", setSize, performance.timeVertexSet, performance.timeDense, performance.timeSparse);
            }

            // Deinit Vertexxset
            Vertexset_Deinit(&vs);
        }   
    }
    // Free fixxed arrays
    free(bitArray);
    free(sparseArray);
    
    if (rank==0) {
        printf("Performance Measurement done!\n");
        fclose(f);
    }
    
}


int main(int argc, char *argv[]){
    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    assert(argc == 2);
    int testNumber = atoi(argv[1]);
    

    // perf_red_exactHalfing: measure performance of allreduce 
    if (testNumber == 0){
        performance_Allreduce(Vertexset_Allreduce_Exact_Halfing);
    }

    // perf_red_approxHalfing: measure performance of allreduce 
    if (testNumber == 1){
        performance_Allreduce(Vertexset_Allreduce_Approximate_Halfing);
    }

    MPI_Finalize();
    return 0;
}

