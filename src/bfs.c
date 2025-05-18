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



int64_t *pred_glob,*column;
int *rowstarts;
oned_csr_graph g;
distributedGraph_CSR graph;


//user should provide this function which would be called once to do kernel 1: graph convert
void make_graph_data_structure(const tuple_graph* const tg) {
	//graph conversion, can be changed by user by replacing oned_csr.{c,h} with new graph format
	/*
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
	*/

	convert_graph_to_oned_csr(tg, &g);
	
}

//user should provide this function which would be called several times to do kernel 2: breadth first search
//pred[] should be root for root, -1 for unrechable vertices
//prior to calling run_bfs pred is set to -1 by calling clean_pred
void run_bfs(int64_t root, int64_t* pred) {
	pred[root] = root;

	/*
	int64_t nglobalverts = g.nglobalverts-1;
	// init VertexSets
	Vertexset visited, frontierOld, frontierNew, temp;
	Vertexset_Init(&visited, nglobalverts, MPI_COMM_WORLD);
	Vertexset_Init(&frontierOld, nglobalverts, MPI_COMM_WORLD);
	Vertexset_Init(&frontierNew, nglobalverts, MPI_COMM_WORLD);
	Vertexset_TransformToDense(&visited);  //only needs dense representation
	Vertexset_TransformToDense(&frontierNew);  //only needs dense representation
	Vertexset_TransformToDense(&frontierOld);  //only needs dense representation

	// set in root
	Vertexset_Add(&visited, root);
	Vertexset_Add(&frontierOld, root);
	pred[root] = root;

	vertexSetIterator it;
	uint32_t vert, neigh;
	int iterationBFS = 0;
	while (frontierOld.sizeSparse != 0) {
		vertexSetIterator_Init(&it, &frontierOld);
		while (vertexSetIterator_Has_next(&it)) {
			uint32_t vertBefore = vert;
			vert = vertexSetIterator_Next(&it);
			if(!(vert<nglobalverts))
				printf("%d\n", vertBefore);
			Vertexset_Add(&visited, vert);
		}
		vertexSetIterator_Reset(&it);
		while (vertexSetIterator_Has_next(&it)) {
			vert = vertexSetIterator_Next(&it);
			for (size_t i = NEIGHSTART((&graph), vert); i < NEIGHEND((&graph),vert); i++) {
				neigh = graph.data[i];
				if (!Vertexset_Contains(&visited, neigh)) {
					Vertexset_Add(&visited, neigh);
					Vertexset_Add(&frontierNew, neigh);
					pred[neigh] = vert;
				}
			}
		}
		//Vertexset_PrintSet(&frontierNew);
		Vertexset_Allreduce_Pure(&frontierNew, VERTEXSET_OR);
		

		//swap Vertexsets
		temp = frontierNew;
		frontierNew = frontierOld;
		frontierOld = temp;


		Vertexset_Clean(&frontierNew);

		iterationBFS++;
	}

	// perform reduce over preds array
	// MPI_Allreduce(MPI_IN_PLACE, pred, nglobalverts, MPI_INT64_T, MPI_MAX, MPI_COMM_WORLD);	

	// deinit VertexSets
	Vertexset_Deinit(&visited);
	Vertexset_Deinit(&frontierOld);
	Vertexset_Deinit(&frontierNew);
		printf("uter\n");
	*/
}

//we need edge count to calculate teps. Validation will check if this count is correct
//user should change this function if another format (not standart CRS) used
void get_edge_count_for_teps(int64_t* edge_visit_count) {
	/*
	long i,j;
	long edge_count=0;
	for(i=0;i<g.nlocalverts;i++)
	if(pred_glob[i]!=-1) {
		for(j=g.rowstarts[i];j<g.rowstarts[i+1];j++)
		if(COLUMN(j)<=VERTEX_TO_GLOBAL(my_pe(),i))
		edge_count++;
	}
	aml_long_allsum(&edge_count);
	*/
	*edge_visit_count=1;
}

//user provided function to initialize predecessor array to whatevere value user needs
void clean_pred(int64_t* pred) {
	for(int64_t i=0;i<g.nlocalverts;i++) {
		pred[i] = -1;
	}
}

//user provided function to be called once graph is no longer needed
void free_graph_data_structure(void) {
	free_oned_csr_graph(&g);
}

//user should change is function if distribution(and counts) of vertices is changed
size_t get_nlocalverts_for_pred(void) {
	return g.nlocalverts;
}  
