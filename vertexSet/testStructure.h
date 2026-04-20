#if !defined(TESTSTRUCTURE)
#define TESTSTRUCTURE

#include "vertexSet.h"

typedef void (*vsFunc) (Vertexset*, int);

enum Testcase{
    ITERATOR,
    EXACT_HALFING,
    APROXIMATE_HALFING,
    NAIVE,
    DENSE,
    VERTEXSET_TEST,
    ALLGATHER
};

vsFunc getVertexsetFunc(int funcname){
    switch (funcname) {
    case EXACT_HALFING:
        return Vertexset_Allreduce_Exact_Halfing;
    case APROXIMATE_HALFING:
        return Vertexset_Allreduce_Approximate_Halfing;
    case NAIVE:
        return Vertexset_Allreduce_Naive;
    case DENSE:
        return Vertexset_Allreduce_Dense;
    case ALLGATHER:
        return Vertexset_Allgather;
    default:
        exit(1);
    }
}

#endif // TESTSTRUCTURE

