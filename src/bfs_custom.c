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

int64_t *column;
int64_t *pred_glob;
unsigned int * rowstarts;



size_t edgeCount;
oned_csr_graph g;
distributedGraph_CSR graph;


//user should provide this function which would be called once to do kernel 1: graph convert
void make_graph_data_structure(const tuple_graph* const tg) {
	//graph conversion, can be changed by user by replacing oned_csr.{c,h} with new graph format

	convert_graph_to_oned_csr(tg, &g);
	createDistributedGraph(tg, &graph, g.nglobalverts -1);	
}

bool validate(int64_t root, int64_t *preds){
	// Indicator, that everything is correct
	bool isCorrect  = true;
	
	// Split Vertices onto ranks
	int32_t start_vert = graph.nGlobalVerts * (rank) / size;
	int32_t end_vert = graph.nGlobalVerts * (rank + 1) / size;
	
	// Phase 1: Check for cycles in preds and count level
	size_t level[graph.nGlobalVerts];
	{
		size_t i;
		int64_t pred;
		for (size_t vert = start_vert; vert < end_vert; vert++) {
			// find predecessor
			pred = preds[vert];

			// ignore verts, that have no predecessor
			if (pred == -1) {
				level[vert] = 0;
				continue;
			}

			// set root
			if (vert == root) {
				level[vert] = 0;
				continue;
			}
			
			i=1;
			while (i<graph.nGlobalVerts) {
				// check if pred is root
				if (pred == root) {
					level[vert] = i;
					break;
				}
				
				// Check if there is no pred
				if (pred == -1) {
					level[vert] = 0;
					break;
				}

				pred = preds[pred];
				i++;
			}
			if (!(i<graph.nGlobalVerts)) {
				printf("VALIDATION ERROR: Rank %d found a circular path in BFS-tree!\n", rank);
				isCorrect = false;
				break;
			}		
		}
	}

	// get displacements and counts from other ranks for MPI Collective
	int displ_mpi[size];
	int count_mpi[size];
	for (size_t r = 0; r < size; r++) {
		displ_mpi[r] = graph.nGlobalVerts * r / size;
		count_mpi[r] = graph.nGlobalVerts * (r+1) / size - displ_mpi[r];
	}
	// Communcicate calculated level to other ranks
	MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, level, count_mpi, displ_mpi, MPI_LONG, MPI_COMM_WORLD);
	

	// Phase 2: check, if counts of each pred-tree edge differs by exxactly one
	int64_t in;
	for (size_t out = start_vert; out < end_vert; out++) {
		in = preds[out];

		// only for connected edges and ones, that are not root
		if (in != -1 && out!= root) {
			int difference = level[out] - level[in];
			if (difference != 1) {
				printf("VALIDATION ERROR: Rank %d reported a predecessor with level difference not to 1!\n", rank);
				printf("difference: %d\n", difference);
				isCorrect = false;
				break;
			}
		}	
	}

	// Phase 3: Check, if every graph edge differs at most by one
	uint32_t neigh;
	uint64_t startIdx, endIdx;
	size_t level_vert, level_neigh;
	for (uint32_t vert = 0; vert < graph.nGlobalVerts; vert++) {
		for (uint64_t i = NEIGHSTART((&graph), vert); i < NEIGHEND((&graph), vert); i++) {
			neigh = graph.data[i];
			level_vert = level[vert];
			level_neigh = level[neigh];
			// if level == 0, vertex is either root or not reachable by BFS
			int difference = level_vert - level_neigh;
			if (difference > 1 || difference < -1) {
				printf("VALIDATION ERROR: Rank %d reported an edge with level differing more than 1!\n", rank);
				printf("difference: %d\n", difference);
				printf("root %d, pred1 %d, pred2 %d, levelRoot %d\n", root, preds[vert], preds[neigh], level[root]);
				isCorrect = false;
				break;
			}
		}
	}
	
	// Phase 4: Check, if all nodes were reached in BFS
	// <=> No included vertex should have a neighbor with pred[neigh] == -1
	for (uint32_t vert = 0; vert < graph.nGlobalVerts; vert++) {
		if (preds[vert] != -1) {
			for (uint64_t i = NEIGHSTART((&graph), vert); i < NEIGHEND((&graph), vert); i++) {
				neigh = graph.data[i];
				assert(neigh < graph.nGlobalVerts);
				if (preds[neigh] == -1) {
					printf("VALIDATION ERROR: Rank %d found a vertex, that was missed by BFS!\n", rank);
					printf("root %d, vertex %d, neigh %d\n", root, vert, neigh);
					printf("pred[vert] %d, pred[neigh] %d\n", preds[vert], preds[neigh]);
					isCorrect = false;
					break;
				}				
			}
			if (!isCorrect) {
				break;
			}
		}		
	}
	
	// Phase 5: Check, if all edges from tree appear in original graph
	int64_t pred;
	Vertexset checked;
	Vertexset_Init(&checked, graph.nGlobalVerts, MPI_COMM_WORLD);
	Vertexset_TransformToDense(&checked);
	for (size_t vert = 0; vert < graph.nGlobalVerts; vert++) {
		pred = preds[vert];
		
		// ignore not connected or root node
		if (pred == -1 || vert == root) {
			Vertexset_Add(&checked, vert);
		} else {
			// Check if any of the neighbors is pred
			for (uint64_t i = NEIGHSTART((&graph), vert); i < NEIGHEND((&graph), vert); i++) {
				neigh = graph.data[i];
				if (neigh == pred) {
					Vertexset_Add(&checked, vert);
				}
			}
		}
	}
	// Reduce with results from other ranks
	Vertexset_Allreduce_Pure(&checked, VERTEXSET_OR);
	// Check for missing entries
	for (size_t vert = start_vert; vert < end_vert; vert++) {
		if (!Vertexset_Contains(&checked, vert)) {
			printf("VALIDATION ERROR: Rank %d found a tree edge, that is not in original graph!\n", rank);
			isCorrect = false;
		}
	}
	Vertexset_Deinit(&checked);
	
	


	// Communicate status with other ranks
	MPI_Allreduce(MPI_IN_PLACE, &isCorrect, 1, MPI_C_BOOL, MPI_LAND, MPI_COMM_WORLD);
	assert(isCorrect);
	return isCorrect;
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
			for (uint64_t i = NEIGHSTART((&graph), vert); i < NEIGHEND((&graph),vert); i++) {
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

	// put in correctly in output pred array
	for (size_t i = 0; i < nglobalverts; i++) {
		if (VERTEX_OWNER(i) == rank) {
			pred[VERTEX_LOCAL(i)] = predslocal[i];
		}		
	}
	validate(root, predslocal);
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
