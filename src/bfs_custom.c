//Stub for custom BFS implementations

#include "common.h"
#include "aml.h"
#include "csr_reference.h"
#include "csr_custom.h"
#include "bitmap_reference.h"
#include <stdint.h>
#include <inttypes.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <limits.h>
#include <assert.h>
#include <stdint.h>

// added by me
#include <stdio.h>

//VISITED bitmap parameters
unsigned long *visited;
int64_t visited_size;

int64_t *pred_glob,*column;
int *rowstarts;
oned_csr_graph g_old;
distributedGraph_CSR graph;


// new and old local frontiers
size_t *front_new, *front_old;
// size of corresponding local frontiers
size_t front_new_size, front_old_size;

//user should provide this function which would be called once to do kernel 1: graph convert
void make_graph_data_structure(const tuple_graph* const tg) {
	//graph conversion, can be changed by user by replacing oned_csr.{c,h} with new graph format
	
	createDistributedGraph(tg, &graph);
	//convert_graph_to_oned_csr(tg, &g_old);
	
	if(rank == 0){
		uint32_t *start = START(52);
		uint32_t *end = END(52);
		printNeighbors(&graph, 52, 0);
		for (uint32_t *index = start; index != end; ++index) {
			printf("%ld ", *index);
		}
		printf("\n");
		
	}

	// create bitmap, where visited vertices are stored
	visited_size = (g_old.nlocalverts + ulong_bits - 1) / ulong_bits;
	visited = xmalloc(visited_size*sizeof(unsigned long));
	//user code to allocate other buffers for bfs

}

//user should provide this function which would be called several times to do kernel 2: breadth first search
//pred[] should be root for root, -1 for unrechable vertices
//prior to calling run_bfs pred is set to -1 by calling clean_pred
void run_bfs(int64_t root, int64_t* pred) {
	// allocate memory for frontiers
	front_new = (size_t*) malloc(sizeof(size_t) * g_old.nglobalverts);
	front_old = (size_t*) malloc(sizeof(size_t) * g_old.nglobalverts);
	front_new_size = 0;
	front_old_size = 0;

	if (VERTEX_OWNER(root)==rank) {
		pred[VERTEX_LOCAL(root)] = root;
		SET_VISITED(root);

		// add root to local front
		front_new[0] = root;
		front_new_size = 1;
	}
	
	pred_glob=pred;
	printf("Predecessor local= %ld", pred[0]);
	//user code to do bfs
}

//we need edge count to calculate teps. Validation will check if this count is correct
//user should change this function if another format (not standart CRS) used
void get_edge_count_for_teps(int64_t* edge_visit_count) {
	long i,j;
	long edge_count=0;
	for(i=0;i<g_old.nlocalverts;i++)
		if(pred_glob[i]!=-1) {
			for(j=g_old.rowstarts[i];j<g_old.rowstarts[i+1];j++)
				if(COLUMN(j)<=VERTEX_TO_GLOBAL(my_pe(),i))
					edge_count++;
		}
	aml_long_allsum(&edge_count);
	*edge_visit_count=edge_count;
}

//user provided function to initialize predecessor array to whatevere value user needs
void clean_pred(int64_t* pred) {
	int i;
	for(i=0;i<g_old.nlocalverts;i++) pred[i]=-1;
}

//user provided function to be called once graph is no longer needed
void free_graph_data_structure(void) {
	free_oned_csr_graph(&g_old);
	free(visited);
}

//user should change is function if distribution(and counts) of vertices is changed
size_t get_nlocalverts_for_pred(void) {
	return g_old.nlocalverts;
}  
