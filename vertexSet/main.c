#include "vertexSet.h"
#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"

int rank, size;

void printBinary(int num) {
    for (int i = sizeof(int) * 8 - 1; i >= 0; i--) {
        printf("%d", (num >> i) & 1);
    }
    printf("\n");
};

void printArr(int *arr, int len) {
    printf("Array from rank %d: ", rank);
    for (size_t i = 0; i < len; i++){
        printf("%d ", arr[i]);
    }
    printf("\n");
};

void predReduce(int *invec, int *inoutvec, int *len, MPI_Datatype *datatype){
    *inoutvec=20;
    /*
    if (*inoutvec == 0) {
        *inoutvec = *invec;
    }
    */

}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank==0) {    
    int32_t *a;
    int32_t *b;
    
    a = (int32_t*) malloc(15 * sizeof(int32_t));
    b = (int32_t*) malloc(15 * sizeof(int32_t));

    for (size_t i = 0; i < 15; i++) {
        a[i] = i;
        b[i] = i+15;
    }

    int32_t *temp;

    temp = a;
    a = b;
    b = temp;
    temp = NULL;

    free(a);
    free(b);
    free(temp);
    }
        
    /*
    int *buff;
    if (rank == 0) {
        int n = 10;
        buff = (int*) malloc(n * sizeof(int));
        
        MPI_Send(buff, n, MPI_INT, 1, 100, MPI_COMM_WORLD);
    }
    
    if (rank==1) {    
        MPI_Status status;
        int count;
        MPI_Probe(0, 100, MPI_COMM_WORLD, &status);
        MPI_Get_count(&status, MPI_INT, &count);
        
        printf("Count: %d", count);
    }
    */
    
    
    MPI_Finalize();
    return 0;
}



