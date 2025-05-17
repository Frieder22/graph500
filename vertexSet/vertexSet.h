#if !defined(VERTEXSET)
#define VERTEXSET

#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>
#include <bitmap_custom.h>
#include <mpi.h>

typedef struct Vertexset{
    // General information
    size_t maxsize;
    size_t sizeCrit;
    bool isdense;
    
    // Data for Bitarray
    unsigned long long *bitArray;
    unsigned long long *bitBuffer;
    uint32_t size_bitarray;
    
    // Data for sparse array
    size_t sizeSparse;
    uint32_t *sparseArray;
    uint32_t *sparseBuffer;
    int *sizesAll;
    int *displ;

    // mpi information
    MPI_Comm MPI_COMM;
    int mpi_rank;
    int mpi_size;
} Vertexset;

/**
 * Initializes Vertex Set object.
 * Runs in O(n)
 * n: maximum number of vertices
 * @param vs Vertex Set object, that needs to be initialized
 * @param maxsize Maximum number of vertices, that can be added
 * @param MPI_COMM MPI communicator, where reduce operation are possible
 */
void Vertexset_Init(Vertexset* vs, uint32_t maxsize, MPI_Comm MPI_COMM);

/**
 * Adds a vertex to vertexset.
 * Time complexities:
 * Dense: O(1)
 * Sparse: O(1)
 * @param vs corresponding Vertexset object 
 * @param vertex added vertex
 */
void Vertexset_Add(Vertexset* vs, uint32_t vertex);

/**
 * Checks, if a vertex is in vertexset.
 * Time complexities:
 * Dense: O(1)
 * Sparse: O(k)     ,k is number of elements in vs
 * @param vs corresponding Vertexset object 
 * @param vertex checked vertex
 */
bool Vertexset_Contains(Vertexset* vs, uint32_t vertex);

/**
 * Resets the memory of vertexset. After reset the VertexSet is
 * in sparse representation.
 * Sparse/Dense O(1)
 */
void Vertexset_Clean(Vertexset* vs);

/**
 * Transforms Vertexset to dense representation.
 * Time complexities:
 * Dense: O(1)
 * Sparse: O(k)     ,k is number of elements in vs
 * @param vs corresponding Vertexset object 
 */
bool Vertexset_TransformToDense(Vertexset* vs);

/**
 * Transforms Vertexset to sparse representation.
 * Time complexities:
 * Dense: O(n + k)     ,k is number of elements in vs
 * Sparse: O(1)
 * @param vs corresponding Vertexset object 
 */
bool Vertexset_TransformToSparse(Vertexset* vs);

/**
 * Performs an Allreduce operation with vertex sets from
 * other ranks.
 * Time complexities:
 * no guarantee
 * @param vs corresponding Vertexset object 
 * @param VERTEXSET_OPERATION which logical operation should be applied
 */
void Vertexset_Allreduce(Vertexset* vs, int VERTEXSET_OPERATION);

/**
 * Performs an Allreduce operation with vertex sets from
 * other ranks. Only sparse or dense representation are used during
 * the communciation. All ranks must have vs in the same kind of
 * representation.
 * Time complexities:
 * Dense: O(Allreduce(n))
 * Sparse: O(Allgather(k))
 * @param vs corresponding Vertexset object 
 * @param VERTEXSET_OPERATION which logical operation should be applied
 */
void Vertexset_Allreduce_Pure(Vertexset* vs, int VERTEXSET_OPERATION);

/**
 * Performs an Allreduce operation with vertex sets from
 * other ranks. At the beginning sparse representation is used and
 * for later communication rounds, switch to dense communication.
 * Time complexities, if vs is already in sparse representation:
 * best case (only sparse): O(log(p)*k*2^(log_2(p+1))-1)
 * worst case (only dense): O(log(p)*(n + n + k))
 * @param vs corresponding Vertexset object 
 * @param VERTEXSET_OPERATION which logical operation should be applied
 */
void Vertexset_Allreduce_Dynamic(Vertexset* vs, int VERTEXXSET_OPERATION);

/**
 * Prints the elements saved in Vertex set. At most 30 elements
 * are printed.
 * @param vs corresponding Vertexset object 
 */
void Vertexset_PrintSet(Vertexset* vs);

/**
 * Cleans up vs object.
 * @param vs corresponding Vertexset object 
 */
void Vertexset_Deinit(Vertexset* vs);


#define VERTEXSET_OR 0
#define VERTEXSET_AND 1
#define VERTEXSET_XOR 2

#endif // VERTEXSET
