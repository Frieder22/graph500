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
#include <stdio.h>

#define SHIFTPATTERN true
#define SHIFT 

int64_t *column;
int64_t *pred_glob;
unsigned int * rowstarts;

int64_t visited_size;
unsigned long *visited;
Vertexset *frontierOld, *frontierNew, *temp;
int64_t* predsAllRanks;

int64_t  nglobalverts;
size_t edgeCount;
oned_csr_graph g;
distributedGraph_CSR graph;



//user should provide this function which would be called once to do kernel 1: graph convert
void make_graph_data_structure(const tuple_graph* const tg) {
	//graph conversion, can be changed by user by replacing oned_csr.{c,h} with new graph format

	convert_graph_to_oned_csr(tg, &g);
	column=g.column;
	rowstarts=g.rowstarts;
 	nglobalverts = g.nglobalverts - 1;
	createDistributedGraph(tg, &graph, nglobalverts);

	// init bitmap
	visited_size = (nglobalverts + ulong_bits - 1) / ulong_bits;
	visited = malloc(visited_size*sizeof(unsigned long));

	// init VertexSets
	frontierNew = (Vertexset*) malloc(sizeof(Vertexset));
	frontierOld = (Vertexset*) malloc(sizeof(Vertexset));
	Vertexset_Init(frontierOld, nglobalverts, MPI_COMM_WORLD, SHIFTPATTERN);
	Vertexset_Init(frontierNew, nglobalverts, MPI_COMM_WORLD, SHIFTPATTERN);

	// init pred array
	predsAllRanks  = (int64_t*) malloc(nglobalverts * sizeof(int64_t));
}

void printPreds(int64_t *preds){
	if (rank==0) {
		printf("preds:\n");
		for (size_t i = 0; i < g.nglobalverts; i++) {
			printf("%d ,", preds[i]);
		}
		printf("\n");
	}
	
}

bool Vertexset_isFilled(Vertexset* vs){
	return vs->isdense + vs->sizeSparse;
}

void Bitmap_UnionWithVertexset(unsigned long* bitmap, Vertexset* vs){
	if (vs->isdense) {
		// dense variant
		for (size_t i = 0; i < vs->size_bitarray; i++) {
			bitmap[i] |= vs->bitArray[i];
		}
	} else {
		// sparse variant
		uint32_t vert;
		for (size_t i = 0; i < vs->sizeSparse; i++) {
			vert = vs->sparseArray[i];
			Bitmap_Set(bitmap, vert);
		}
	}
	
}

//user should provide this function which would be called several times to do kernel 2: breadth first search
//pred[] should be root for root, -1 for unrechable vertices
//prior to calling run_bfs pred is set to -1 by calling clean_pred
void run_bfs(int64_t root, int64_t* pred) {

	// set in root
	Bitmap_Clean(visited, visited_size);
	Bitmap_Set(visited, root);
	Vertexset_Add(frontierOld, root);

	predsAllRanks[root] = root;

	vertexSetIterator it;
	uint32_t vert, neigh;
	while (Vertexset_isFilled(frontierOld)) {
		vertexSetIterator_Init(&it, frontierOld);
		while (vertexSetIterator_Has_next(&it)) {
			vert = vertexSetIterator_Next(&it);
			Bitmap_Set(visited, vert);
		}
		vertexSetIterator_Init(&it, frontierOld);
		while (vertexSetIterator_Has_next(&it)) {
			vert = vertexSetIterator_Next(&it);			
			for (uint64_t i = NEIGHSTART((&graph), vert); i < NEIGHEND((&graph),vert); i++) {
				neigh = graph.data[i];
				if (!Bitmap_Test(visited, neigh)) {
					edgeCount ++;
					Bitmap_Set(visited, neigh);
					Vertexset_Add(frontierNew, neigh);
					predsAllRanks[neigh] = vert;
				}
			}
		}
		Vertexset_Union_Shift(frontierNew, VERTEXSET_OR);
		

		//swap Vertexsets
		temp = frontierNew;
		frontierNew = frontierOld;
		frontierOld = temp;
		
		Vertexset_Clean(frontierNew);
	}
	
	
	// perform reduce over preds array
	MPI_Allreduce(MPI_IN_PLACE, predsAllRanks, nglobalverts, MPI_INT64_T, MPI_MAX, MPI_COMM_WORLD);	

	
	// put in correctly in output pred array
	for (size_t i = 0; i < nglobalverts; i++) {
		if (VERTEX_OWNER(i) == rank) {
			pred[VERTEX_LOCAL(i)] = predsAllRanks[i];
		}		
	}
	pred_glob = pred;

}

//we need edge count to calculate teps. Validation will check if this count is correct
//user should change this function if another format (not standard CRS) used
void get_edge_count_for_teps(int64_t* edge_visit_count) {
	long i,j;
	long edge_count=0;
	for(i=0;i<g.nlocalverts;i++)
		if(pred_glob[i]!=-1) {
			for(j=rowstarts[i];j<rowstarts[i+1];j++)
				if(COLUMN(j)<=VERTEX_TO_GLOBAL(my_pe(),i))
					edge_count++;

		}

	aml_long_allsum(&edge_count);
	*edge_visit_count=edge_count;
}

//user provided function to initialize predecessor array to whatevere value user needs
void clean_pred(int64_t* pred) {
	for(int64_t i=0;i<g.nlocalverts;i++) {
		pred[i] = -1;
	}
	for (size_t i = 0; i < nglobalverts; i++) {
		predsAllRanks[i] = -1;
	}

}

//user provided function to be called once graph is no longer needed
void free_graph_data_structure(void) {
	free_oned_csr_graph(&g);
	free(visited);

	// deinit VertexSets
	free(frontierNew);
	free(frontierOld);
	Vertexset_Deinit(frontierOld);
	Vertexset_Deinit(frontierNew);

	// deinit preds for all rnks
	free(predsAllRanks);
}

//user should change is function if distribution(and counts) of vertices is changed
size_t get_nlocalverts_for_pred(void) {
	return g.nlocalverts;
}  
