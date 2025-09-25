#include "vertexSet.h"
#include "vertexSetIterator.h"
#include "mpi.h"
#include "testStructure.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>


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
#define N_DATAPOINTS 20
#define SETSIZE 500000
#define SETFILLING 0.1

int rank, size;
extern enum Testcase Testcase;

void fillInFunctionName(char* filepath, void (*reduceFunc) (Vertexset*, int)){
    if (reduceFunc == Vertexset_Allreduce_Exact_Halfing) {
        char filename[] = "exact_Halfing";
        strcat(filepath, filename);
    } else if (reduceFunc == Vertexset_Allreduce_Approximate_Halfing) {
        char filename[] = "approx_Halfing";
        strcat(filepath, filename);
    } else if (reduceFunc == Vertexset_Allreduce_Naive) {
        char filename[] = "naive";
        strcat(filepath, filename);
    } else if (reduceFunc == Vertexset_Allreduce_Dense) {
        char filename[] = "dense";
        strcat(filepath, filename);
    } else {
        printf(ANSI_RED "Performance testing for this function is not implemented!\n" ANSI_RESET);
        return;
    }
}

void fillInClusterConfig(char* filepath){
    // get the number of nodes via the enviroment viarable (provided by slurm)
    char *N_Nodes = getenv("SLURM_NNODES");
    // same for tasks per node
    char *N_TasksPerNode = getenv("SLURM_NTASKS_PER_NODE");
    
    // assume enviroment variables are only set by slurm
    if (N_Nodes){
        strcat(filepath, "_");
        strcat(filepath, N_Nodes);
        strcat(filepath, "x");
        strcat(filepath, N_TasksPerNode);
    } else {
        // identifier for local produced data
        strcat(filepath, "_local");
    }
}

void performance_sizeSeries(void (*reduceFunc) (Vertexset*, int), float filling){
    double performance[3];
    // Set up file IO 
    FILE *f;
    if (rank==0){
        char filepath[1000] = OUTPUTPATH;
        char fillingStr[20]; 
        sprintf(fillingStr, "_%f", filling);

        fillInFunctionName(filepath, reduceFunc);
        fillInClusterConfig(filepath);

        strcat(filepath, "_sizeSeries");
        strcat(filepath, fillingStr);
        strcat(filepath, ".txt");
    
        f = fopen(filepath, "w");
        if (f == NULL) {
            printf(ANSI_RED "File could not be opened!!\n" ANSI_RESET);
            return;
        }
        fprintf(f, "SetSize, VertexSet Time [s], Dense Time [s], Sparse Time [s]\n");
        printf("Size series performance measurement started...\n");
        printf("Results will be written in: \n%s\n", filepath);
    }

    // Definition of testing range
    int sizes[] = {20, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000, 500000, 1000000, 2000000, 5000000, 10000000};
    int N = sizeof(sizes) / sizeof(int);


    // general inits
    int insertions, setSize;
    int32_t num;
    float trueFillage;
    double startTime, endTime;
    srand(3);

    // allocate vertexset
    Vertexset vs;    
    
    // init Bitarray
    size_t size_bitarray;
    unsigned long long *bitArray;
    size_bitarray = (sizes[N-1] + (sizeof(unsigned long long)*8)) / (sizeof(unsigned long long)*8);
    bitArray = (unsigned long long*) malloc(size_bitarray * sizeof(unsigned long long));
    Bitmap_Clean(bitArray, size_bitarray);
    
    // init sparse array
    int sparseSize, sparseSizeMax;
    int32_t *sparseArray, *sparseBuffer;
    int sizesAll[size];
    int displ[size + 1];
    sparseSizeMax = (int)(sizes[N-1] * filling);
    sparseArray = (int32_t*) malloc(sparseSizeMax * sizeof(int32_t));
    sparseBuffer = (int32_t*) malloc(sparseSizeMax * sizeof(int32_t));

    // test performance of scale
    for (size_t i = 0; i < N; i++) {
        setSize = sizes[i];
        insertions = (int) (setSize * filling);

        
        // collect more than one data point for each configuration
        for (size_t j = 0; j < N_DATAPOINTS; j++) {
            // init vertex set
            Vertexset_Init(&vs, setSize, MPI_COMM_WORLD);

            size_bitarray = (setSize + (sizeof(unsigned long long)*8)) / (sizeof(unsigned long long)*8);

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
            performance[0] = endTime - startTime;
            
            // reduce Bitmap
            MPI_Barrier(MPI_COMM_WORLD);
            startTime = MPI_Wtime();
            MPI_Allreduce(MPI_IN_PLACE, bitArray, size_bitarray, MPI_UINT64_T, MPI_BOR, MPI_COMM_WORLD);
            endTime = MPI_Wtime();
            performance[1] = endTime - startTime;
            
            // reduce Sparse Array
            MPI_Barrier(MPI_COMM_WORLD);
            startTime = MPI_Wtime();
            // get size information from other ranks
            MPI_Allgather(&sparseSize, 1, MPI_INT, sizesAll, 1, MPI_INT, MPI_COMM_WORLD);
            // calculate displacements with exclusive scan
            displ[0] = 0;
            for (int j = 1; j < size + 1; j++){
                displ[j] = displ[j-1] + sizesAll[j-1];
            }
            assert(displ[size] <= sparseSizeMax);
            MPI_Allgatherv(sparseArray, sparseSize, MPI_INT32_T, sparseBuffer, sizesAll, displ, MPI_INT32_T, MPI_COMM_WORLD);
            // put all numbers in (same) to perform logic OR (removes dublicates)
            for (size_t j = 0; j < displ[size]; j++) {
                Bitmap_Set(bitArray,sparseBuffer[j]);
            }
            endTime = MPI_Wtime();
            performance[2] = endTime - startTime;

            // Find true fillage 
            sparseSize = 0;
            for (size_t j = 0; j < setSize; j++) {
                if (Bitmap_Test(bitArray, j)) {
                    sparseSize++;
                }
            }            

            // Find max runtime
            MPI_Allreduce(MPI_IN_PLACE, performance, 3, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

            // save results
            if (rank==0) {
                trueFillage = sparseSize/(1.0 * setSize);
                fprintf(f, "%d,%lf,%lf,%lf, %f\n", setSize, performance[0], performance[1], performance[2], trueFillage);
            }

            // Deinit Vertexxset
            Vertexset_Deinit(&vs);
        }   
    }
    // Free fixed arrays
    free(bitArray);
    free(sparseArray);
    free(sparseBuffer);
    
    if (rank==0) {
        printf(ANSI_GREEN "Size series performance measurement done!\n" ANSI_RESET);
        printf("-------------------------------------------------------------\n");
        fclose(f);
    }
}

void performance_fillingSeries(void (*reduceFunc) (Vertexset*, int), int setSize){
    // Set up file IO
    double performance[3];
    FILE *f;
    if (rank==0){
        char filepath[1000] = OUTPUTPATH;
        char setSizeStr[20]; 
        sprintf(setSizeStr, "_%d", setSize);

        fillInFunctionName(filepath, reduceFunc);
        fillInClusterConfig(filepath);

        strcat(filepath, "_fillingSeries");
        strcat(filepath, setSizeStr);
        strcat(filepath, ".txt");
    
        f = fopen(filepath, "w");
        if (f == NULL) {
            printf(ANSI_RED "File could not be opened!!\n" ANSI_RESET);
            return;
        }
        fprintf(f, "Fillage, VertexSet Time [s], Dense Time [s], Sparse Time [s]\n");
        printf("Filling series performance measurement started...\n");
        printf("Results will be written in: \n%s\n", filepath);
    }

    // Definition of testing range
    float fillings[] = {0.005, 0.01, 0.02, 0.03, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95};
    int N = sizeof(fillings) / sizeof(float);


    // general inits
    int insertions;
    int32_t num;
    float trueFillage, filling;
    double startTime, endTime;
    srand(3);

    // allocate vertexset
    Vertexset vs;    
    
    // init Bitarray
    size_t size_bitarray;
    unsigned long long *bitArray;
    size_bitarray = (setSize + (sizeof(unsigned long long)*8)) / (sizeof(unsigned long long)*8);
    bitArray = (unsigned long long*) malloc(size_bitarray * sizeof(unsigned long long));
    Bitmap_Clean(bitArray, size_bitarray);
    
    // init sparse array
    int sparseSize, sparseSizeMax;
    int32_t *sparseArray, *sparseBuffer;
    int sizesAll[size];
    int displ[size + 1];
    sparseSizeMax = (int)(setSize);
    sparseArray = (int32_t*) malloc(sparseSizeMax * sizeof(int32_t));
    sparseBuffer = (int32_t*) malloc(sparseSizeMax * sizeof(int32_t));

    // test performance of scale
    for (size_t i = 0; i < N; i++) {
        filling = fillings[i];
        insertions = (int) (setSize * filling);

        // collect more than one data point for each configuration
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
            performance[0] = endTime - startTime;
            
            // reduce Bitmap
            MPI_Barrier(MPI_COMM_WORLD);
            startTime = MPI_Wtime();
            MPI_Allreduce(MPI_IN_PLACE, bitArray, size_bitarray, MPI_UINT64_T, MPI_BOR, MPI_COMM_WORLD);
            endTime = MPI_Wtime();
            performance[1] = endTime - startTime;
            
            // reduce Sparse Array
            MPI_Barrier(MPI_COMM_WORLD);
            startTime = MPI_Wtime();
            // get size information from other ranks
            MPI_Allgather(&sparseSize, 1, MPI_INT, sizesAll, 1, MPI_INT, MPI_COMM_WORLD);
            // calculate displacements with exclusive scan
            displ[0] = 0;
            for (int j = 1; j < size + 1; j++){
                displ[j] = displ[j-1] + sizesAll[j-1];
            }
            assert(displ[size] <= sparseSizeMax);
            MPI_Allgatherv(sparseArray, sparseSize, MPI_INT32_T, sparseBuffer, sizesAll, displ, MPI_INT32_T, MPI_COMM_WORLD);
            // put all numbers in (same) to perform logic OR (removes dublicates)
            for (size_t j = 0; j < displ[size]; j++) {
                Bitmap_Set(bitArray,sparseBuffer[j]);
            }
            endTime = MPI_Wtime();
            performance[2] = endTime - startTime;

            // Find true fillage 
            sparseSize = 0;
            for (size_t j = 0; j < setSize; j++) {
                if (Bitmap_Test(bitArray, j)) {
                    sparseSize++;
                }
            }

            // Find max runtime
            MPI_Allreduce(MPI_IN_PLACE, performance, 3, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

            // save results
            if (rank==0) {
                trueFillage = sparseSize/(1.0 * setSize);
                fprintf(f, "%f,%lf,%lf,%lf, %f\n", filling, performance[0], performance[1], performance[2], trueFillage);
            }

            // Deinit Vertexxset
            Vertexset_Deinit(&vs);
        }   
    }
    // Free fixed arrays
    free(bitArray);
    free(sparseArray);
    free(sparseBuffer);
    
    if (rank==0) {
        printf(ANSI_GREEN "Fillage series performance measurement done!\n" ANSI_RESET);
            printf("-------------------------------------------------------------\n");
        fclose(f);
    }
}


int main(int argc, char *argv[]){
    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    assert(argc == 2);
    int testNumber = atoi(argv[1]);
    

    switch (testNumber) {
    // perf_red_exactHalfing: measure performance of allreduce 
    case EXACT_HALFING:
        if (rank==0) {
            printf("-------------------------------------------------------------\n");
            printf("Testing performance of Vertexset_Allreduce_Exact_Halfing...\n");
            printf("-------------------------------------------------------------\n");
        }
        performance_sizeSeries(Vertexset_Allreduce_Exact_Halfing, SETFILLING);
        performance_fillingSeries(Vertexset_Allreduce_Exact_Halfing, SETSIZE);
        break;

    // perf_red_approxHalfing: measure performance of allreduce 
    case APROXIMATE_HALFING:
        if (rank==0) {
            printf("-------------------------------------------------------------\n");
            printf("Testing performance of Vertexset_Allreduce_Approximate_Halfing...\n");
            printf("-------------------------------------------------------------\n");
        }        
        performance_sizeSeries(Vertexset_Allreduce_Approximate_Halfing, SETFILLING);
        break;

    // measure peprformance of ring comm reduce
    case NAIVE:
        if (rank==0) {
            printf("-------------------------------------------------------------\n");
            printf("Testing performance of Vertexset_Allreduce_Naive...\n");
            printf("-------------------------------------------------------------\n");
        }
        performance_sizeSeries(Vertexset_Allreduce_Naive, SETFILLING);
        performance_fillingSeries(Vertexset_Allreduce_Naive, SETSIZE);
        break;

    // measure peprformance of dense reduce
    case DENSE:
        if (rank==0) {
            printf("-------------------------------------------------------------\n");
            printf("Testing performance of Vertexset_Allreduce_Dense...\n");
            printf("-------------------------------------------------------------\n");
        }
        performance_sizeSeries(Vertexset_Allreduce_Dense, SETFILLING);
        performance_fillingSeries(Vertexset_Allreduce_Dense, SETSIZE);
        break;

    default:
        break;
    }

    MPI_Finalize();
    return 0;
}

