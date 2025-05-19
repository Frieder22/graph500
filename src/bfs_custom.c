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

int64_t *column;
int64_t *pred_glob;
unsigned int * rowstarts;



size_t edgeCount;
oned_csr_graph g;
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

	
	convert_graph_to_oned_csr(tg, &g);
	createDistributedGraph(&tg_copy, &graph);

	free(edgebuff);	
}

//user should provide this function which would be called several times to do kernel 2: breadth first search
//pred[] should be root for root, -1 for unrechable vertices
//prior to calling run_bfs pred is set to -1 by calling clean_pred
void run_bfs(int64_t root, int64_t* pred) {
	int64_t  nglobalverts = g.nglobalverts - 1;
	int64_t* predslocal;
	predslocal  = (int64_t*) malloc(nglobalverts * sizeof(int64_t));
	for (size_t i = 0; i < nglobalverts; i++) {
		predslocal[i] = -1;
	}

	edgeCount = 0;

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

	predslocal[root] = root;

	vertexSetIterator it;
	uint32_t vert, neigh;
	int iterationBFS = 0;
	while (frontierOld.sizeSparse != 0) {
		vertexSetIterator_Init(&it, &frontierOld);
		while (vertexSetIterator_Has_next(&it)) {
			vert = vertexSetIterator_Next(&it);
			Vertexset_Add(&visited, vert);
		}
		vertexSetIterator_Reset(&it);
		while (vertexSetIterator_Has_next(&it)) {
			vert = vertexSetIterator_Next(&it);			
			for (size_t i = NEIGHSTART((&graph), vert); i < NEIGHEND((&graph),vert); i++) {
				neigh = graph.data[i];
				if (!Vertexset_Contains(&visited, neigh)) {
					edgeCount ++;
					Vertexset_Add(&visited, neigh);
					Vertexset_Add(&frontierNew, neigh);
					predslocal[neigh] = vert;
				}
			}
		}
		Vertexset_Allreduce_Pure(&frontierNew, VERTEXSET_OR);
		

		//swap Vertexsets
		temp = frontierNew;
		frontierNew = frontierOld;
		frontierOld = temp;


		Vertexset_Clean(&frontierNew);
		iterationBFS++;
	}


	// perform reduce over preds array
	MPI_Allreduce(MPI_IN_PLACE, predslocal, nglobalverts, MPI_INT64_T, MPI_MAX, MPI_COMM_WORLD);	
	// deinit VertexSets
	Vertexset_Deinit(&visited);
	Vertexset_Deinit(&frontierOld);
	Vertexset_Deinit(&frontierNew);

	/*
    // print preds list
	if(rank==0){
		for (size_t i = 0; i < nglobalverts; i++) {
			printf("%ld ", predslocal[i]);
		}
		printf("\n");
	}
	*/
	
	// put in correctly in output pred array
	for (size_t i = 0; i < nglobalverts; i++) {
		if (VERTEX_OWNER(i) == rank) {
			pred[VERTEX_LOCAL(i)] = predslocal[i];
		}		
	}

	pred_glob = pred;
	
}

//we need edge count to calculate teps. Validation will check if this count is correct
//user should change this function if another format (not standard CRS) used
void get_edge_count_for_teps(int64_t* edge_visit_count) {
	MPI_Allreduce(MPI_IN_PLACE, &edgeCount, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
	*edge_visit_count=edgeCount;
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
