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
    bool approx_halfing;
    
    // Data for Bitarray
    unsigned long long *bitArray;
    unsigned long long *bitBuffer;
    int size_bitarray;
    
    // Data for sparse array
    int sizeSparse;
    uint32_t *sparseArray;
    uint32_t *sparseBuffer;
    int *sizesAll;
    int *displ;

    // mpi information
    MPI_Comm MPI_COMM;
    int mpi_rank;
    int mpi_size;

    // block informations
    int *block_Indices;
    int *block_nElements;

    // index shift for approx_halfing
    int indexShift;
} Vertexset;

/**
 * Initializes Vertex Set object.
 * Runs in O(n)
 * n: maximum number of vertices
 * @param vs Vertex Set object, that needs to be initialized
 * @param maxsize Maximum number of vertices, that can be added
 * @param MPI_COMM MPI communicator, where reduce operation are possible
 */
void Vertexset_Init(Vertexset* const vs, const uint32_t maxsize, const MPI_Comm MPI_COMM);

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
bool Vertexset_Contains(const Vertexset* vs, uint32_t vertex);

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
 * other ranks. Each rank can have a different representation of 
 * VertexSet. The size of the communicator musst be some 2^k.
 * Time complexities:
 * Dense: O(Allreduce(n))
 * Sparse: O(Allgather(k))
 * @param vs corresponding Vertexset object 
 * @param VERTEXSET_OPERATION which logical operation should be applied
 */
void Vertexset_Allreduce_Exact_Halfing(Vertexset* vs, int VERTEXSET_OPERATION);

/**
 * Performs an Allreduce operation with vertex sets from
 * other ranks. Each rank can have a different representation of 
 * VertexSet.
 * Time complexities:
 * Dense: O(Allreduce(n))
 * Sparse: O(Allgather(k))
 * @param vs corresponding Vertexset object 
 * @param VERTEXSET_OPERATION which logical operation should be applied
 */
void Vertexset_Allreduce_Approximate_Halfing(Vertexset* vs, int VERTEXSET_OPERATION);

void Vertexset_Allreduce_Ring_Comm(Vertexset* vs, int VERTEXSET_OPERATION);

void Vertexset_Allreduce_Dense(Vertexset* vs, int VERTEXSET_OPERATION);


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
