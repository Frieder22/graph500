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

    int n = 3; // number of elements each rank should hold
    int len = n*size;
    int arr[n*size];
    int arrOUT[n*size];

    memset(arr, 0, sizeof(int)*len);
    memset(arrOUT, 0, sizeof(int)*len);

    for (int i = rank*n; i < n*(rank+1); i++) {
        arr[i] = i;
    }

    arr[1] = rank;

    printArr(arr, len);

    MPI_Op someOR;
    MPI_Op_create(predReduce, true, &someOR);

    MPI_Reduce(arr, arrOUT, len, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        printf("\n");
        printArr(arrOUT, len);
    }
    



    MPI_Finalize();
    return 0;
}



