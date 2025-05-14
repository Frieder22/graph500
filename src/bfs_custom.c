//Stub for custom BFS implementations

#include "common.h"
#include "aml.h"
#include "csr_reference.h"
#include "csr_custom.h"
#include "bitmap_custom.h"
#include "vertexSet.h"
#include "vertexSetIterator.h"

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


//user should provide this function which would be called once to do kernel 1: graph convert
void make_graph_data_structure(const tuple_graph* const tg) {
	//graph conversion, can be changed by user by replacing oned_csr.{c,h} with new graph format
	/*
	*/
	tuple_graph tg_copy;
	packed_edge *edgebuff;
	if (rank==0) {
		edgebuff = (packed_edge*) malloc(sizeof(packed_edge) * tg->nglobaledges);
		for (size_t i = 0; i < tg->nglobaledges; i++){
			edgebuff[i] = tg->edgememory[i];
		}
	}
	
	tg_copy = *tg;
	tg_copy.edgememory = edgebuff;

	
	createDistributedGraph(&tg_copy, &graph);

	convert_graph_to_oned_csr(tg, &g_old);

	/*
	if(rank == 0){
		uint32_t *start = START(52);
		uint32_t *end = END(52);
		printNeighbors(&graph, 52, 0);
		for (uint32_t *index = start; index != end; ++index) {
			printf("%ld ", *index);
		}
		printf("\n");
		
	}
	*/

	/*
	if (rank==0) {		
		printf("first array element: %ld\n", visited[0]);
		SET_VISITEDLOC(visited, (uint64_t) 6);
		printf("first array element: %ld\n", visited[0]);
		if (TEST_VISITEDLOC(visited, (uint64_t) 6)){
			printf("IS SET!\n");
		} else {
			printf("is not set\n");
		}
	}
	*/
	
}

//user should provide this function which would be called several times to do kernel 2: breadth first search
//pred[] should be root for root, -1 for unrechable vertices
//prior to calling run_bfs pred is set to -1 by calling clean_pred
void run_bfs(int64_t root, int64_t* pred) {
	int64_t nglobalverts = g_old.nglobalverts - 1;
	// init VertexSets
	Vertexset visited, frontierOld, frontierNew, temp;
	Vertexset_Init(&visited, nglobalverts, MPI_COMM_WORLD);
	Vertexset_Init(&frontierOld, nglobalverts, MPI_COMM_WORLD);
	Vertexset_Init(&frontierNew, nglobalverts, MPI_COMM_WORLD);
	Vertexset_TransformToDense(&visited);  //only needs dense representation

	// set in root
	Vertexset_Add(&visited, root);
	Vertexset_Add(&frontierOld, root);
	pred[root] = root;

	vertexSetIterator it;
	uint32_t vert, neigh;
	while (frontierOld.sizeSparse != 0) {
		vertexSetIterator_Init(&it, &frontierOld);
		while (vertexSetIterator_Has_next(&it)) {
			Vertexset_Add(&visited, vertexSetIterator_Next(&it));
		}
		vertexSetIterator_Reset(&it);
		while (vertexSetIterator_Has_next(&it)) {
			vert = vertexSetIterator_Next(&it);
			printf("START\n");
			for (size_t i = START(vert); i < END(vert); i++) {
				neigh = graph.data[i];
				if (!Vertexset_Contains(&visited, neigh)) {
					Vertexset_Add(&visited, neigh);
					Vertexset_Add(&frontierNew, neigh);
					pred[neigh] = vert;
				}
			}
		}
		Vertexset_Allreduce_Pure(&frontierNew, VERTEXSET_OR);

		//swap Vertexsets
		temp = frontierNew;
		frontierNew = frontierOld;
		frontierOld = temp;

		Vertexset_Clean(&frontierNew);
	}
	


	// deinit VertexSets
	Vertexset_Deinit(&visited);
	Vertexset_Deinit(&frontierOld);
	Vertexset_Deinit(&frontierNew);
}

//we need edge count to calculate teps. Validation will check if this count is correct
//user should change this function if another format (not standart CRS) used
void get_edge_count_for_teps(int64_t* edge_visit_count) {
	/*
	long i,j;
	long edge_count=0;
	for(i=0;i<g_old.nlocalverts;i++)
	if(pred_glob[i]!=-1) {
		for(j=g_old.rowstarts[i];j<g_old.rowstarts[i+1];j++)
		if(COLUMN(j)<=VERTEX_TO_GLOBAL(my_pe(),i))
		edge_count++;
	}
	aml_long_allsum(&edge_count);
	*/
	*edge_visit_count=1020000;
}

//user provided function to initialize predecessor array to whatevere value user needs
void clean_pred(int64_t* pred) {
	for(int64_t i=0;i<g_old.nlocalverts;i++) {
		pred[i] = -1;
	}
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
