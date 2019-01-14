/* graph.c
 * 
 * Copyright (C) 2001 Havard Rue
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 * The author's contact information:
 *
 *       H{\aa}vard Rue
 *       Department of Mathematical Sciences
 *       The Norwegian University of Science and Technology
 *       N-7491 Trondheim, Norway
 *       Voice: +47-7359-3533    URL  : http://www.math.ntnu.no/~hrue  
 *       Fax  : +47-7359-3524    Email: havard.rue@math.ntnu.no
 *
 */
/*!
  \file graph.c
  \brief Functions for creating and handling general and regular graphs
  
  There are two ways of specifying a general graph.

  - <em> Reading from a file.</em>
  The simplest way to specify a general
  graph is to read the graph specifications from a file, using the
  function \c GMRFLib_read_graph(). This function allocates
  and initializes the members based on the information given on the
  input file, and customizes the graph for use in other library
  functions.
  - <em> Explicitly creating a \c GMRFLib_graph_tp - object.</em> 
  The user can create a general graph by allocating and
  initializing a variable of type \c GMRFLib_graph_tp, and then 
  (IMPORTANT!) calling the function \c GMRFLib_prepare_graph() to customize 
  the graph and check that the graph is consistently defined.  The user 
  should specify the the members \em n, \em nnbs and \em nbs in 
  \c GMRFLib_graph_tp. The member \em mothergraph_idx is initialized, 
  and needed, only if the graph is generated as a subgraph of another graph, 
  using \c GMRFLib_compute_subgraph().

  In addition to the functions for specifying general graphs, the library
  provides two separate functions for the specification of
  two-dimensional lattice graphs and one-dimensional linear graphs, and
  a separate setup for the generation of graphs for a weighted average
  model (see wa.c).

  - <em> A lattice graph:</em> A two-dimensional graph on a lattice is
  most easily specified using the function \c GMRFLib_make_lattice_graph(), 
  specifying the grid sizes \f$ n_{row} \f$  and
  \f$ n_{col} \f$ and the parameters \f$ m_{row} \f$ and \f$ m_{col} \f$, 
  defining the neighbourhood  \f$ (2 m_{row} + 1)\times (2 m_{col} + 1) \f$.
  - <em> A linear graph:</em> To specify a one-dimensional linear graph,
  i.e. an autoregressive AR(\em p)-model, use the function 
  \c GMRFLib_make_linear_graph(), specifying the number of nodes
  and the size of the neighbourhood.
*/

#include <assert.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#if !defined(__FreeBSD__)
#include <malloc.h>
#endif

#include "GMRFLib.h"
#include "GMRFLibP.h"

static const char RCSId[] = "$Id: graph.c 1 2013-03-28 13:54:24Z hanne $";

/*!
  \brief Creates an empty graph
  \param[in,out] graph A pointer to a \c GMRFLib_graph_tp pointer. 
  At output, \em graph points to an empty graph.
  \par Example:
  <tt> main(int argc, char* argv[])\n
  {\n
  ::::::::\n
  GMRFLib_graph_tp *graph;\n
  GMRFLib_make_empty_graph(&graph);\n
  ::::::::\n
  }
  </tt>
*/
int GMRFLib_make_empty_graph(GMRFLib_graph_tp **graph)
{
    /* 
       this function creates an empty graph
    */
    
    *graph = Calloc(1, GMRFLib_graph_tp); MEMCHK(*graph);

    /* 
       user variables
    */
    (*graph)->n               = 0;
    (*graph)->nbs             = NULL;
    (*graph)->nnbs            = NULL;

    /* 
       private variables
    */
    (*graph)->mothergraph_idx = NULL;
    
    return GMRFLib_SUCCESS;
}

/*!
  \brief Reads a graph from a file

  \param[in,out] graph At output, <em>(*graph)</em> has been initialized 
  by using the information on the file \em filename. 

  \param[in] filename The name of the file, formatted as described 
  below, containing the specification of the graph.

  \remarks The file is assumed to be of the format
  \verbinclude example_doxygen_file_format.txt
  Here, \em n, \em nbs and \em nnbs refer to the members of 
  \c GMRFLib_graph_tp, and <em>nn[i]</em> is the node number
  of the node with index \em i, running from \em 0 to <em>n-1</em>.  
  In words, the first row contains the number of nodes, \em n, and the
  successive rows, rows <em>i = 2,...,n+1</em>, list the node number,
  <em>nn[i]</em>, the number of neighbours, <em>nnbs[i]</em>, and
  the node numbers, <em>nbs[i][j]</em>, for the neighbours of each
  node \em i of the graph. The number of neighbours <em>nnbs[i]</em>
  might be zero.\n\n
  The function calls \c GMRFLib_prepare_graph() to
  customize the graph for further computations.

  \par Example:
  See \ref ex_graph

  \sa GMRFLib_prepare_graph, GMRFLib_print_graph
*/
int GMRFLib_read_graph(GMRFLib_graph_tp **graph, char *filename)
{
    /*
      read a graph in the following format

       N
       node[0] nnbs[0] nbs[node[0]][0] nbs[node[0]][1] ... nbs[node[0]][nnbs[0]-1] 
       node[1] nnbs[1] nbs[node[1]][0] nbs[node[1]][1] ... nbs[node[1]][nnbs[1]-1] 
       :
       node[N-1] nnbs[N-1] nbs[node[N-1]][0] nbs[node[N-1]][1] ... nbs[node[N-1]][nnbs[N-1]-1] 
     */
#define READ_ERROR() do { \
    if (1){\
        fprintf(stderr,"\n\n\t%s: error: file [%s]:\n",__GMRFLib_FuncName,filename);\
	fprintf(stderr,"\t\tfail to read tag[%1d] on line [%1d]\n", tag, lineno);\
    }    \
    if (fp) fclose(fp);\
    (*graph) = NULL;\
    GMRFLib_ERROR(GMRFLib_EREADGRAPH);} while(0)
 
    int *storage, n_neig_tot = 0, storage_indx;
    int i, j, lineno, tag, tnode;
    FILE *fp;

    if (!filename) return GMRFLib_SUCCESS;
    if (!(fp = fopen(filename, "r"))) GMRFLib_ERROR(GMRFLib_EOPENFILE);

    GMRFLib_make_empty_graph(graph);

    lineno = tag = 0;
    if (fscanf(fp, "%d", &((*graph)->n)) != 1) READ_ERROR();

    (*graph)->nnbs  = Calloc((*graph)->n, int);   MEMCHK((*graph)->nnbs);
    (*graph)->nbs   = Calloc((*graph)->n, int *); MEMCHK((*graph)->nbs);

    for(i=0;i< (*graph)->n;i++)
    {
	lineno = i+1;
	tag = 0;
	if (fscanf(fp, "%d", &tnode) != 1) READ_ERROR(); /* target node */

	if (tnode < 0 || tnode >= (*graph)->n)
	{
	    fprintf(stderr,"\n\n\t%s: error: file [%s]: line[%1d]: tag[%1d]\n",
		    __GMRFLib_FuncName, filename, lineno, tag); 
	    fprintf(stderr,"\t\tnode-number[%1d] is not in the range[0:%1d]\n", tnode, (*graph)->n);
	    if (fp) fclose(fp);
	    GMRFLib_ERROR(GMRFLib_EPARAMETER);
	}
	    
	tag++;
	if (fscanf(fp, "%d", &((*graph)->nnbs[tnode])) != 1) READ_ERROR();
	n_neig_tot += (*graph)->nnbs[tnode];
	if ((*graph)->nnbs[tnode])
	{
	    (*graph)->nbs[tnode] = Calloc((*graph)->nnbs[tnode], int); MEMCHK((*graph)->nbs[tnode]);
	    for(j=0;j<(*graph)->nnbs[tnode];j++)
	    {
		tag++;
		if (fscanf(fp, "%d", &((*graph)->nbs[tnode][j])) != 1) READ_ERROR();
	    }
	}
	else
	    (*graph)->nbs[tnode] = NULL;
    }
    fclose(fp);

    /*
      map the graph to a more computational convenient memory layout!
      use just one long vector to store all the neighbors.
     */
    storage      = Calloc(n_neig_tot, int);MEMCHK(storage);
    storage_indx = 0;
    for(i=0;i< (*graph)->n;i++)
	if ((*graph)->nnbs[i])
	{
	    memcpy(&storage[storage_indx], (*graph)->nbs[i], sizeof(int)* (*graph)->nnbs[i]);
	    free((void *)(*graph)->nbs[i]);
	    (*graph)->nbs[i] = &storage[storage_indx];
	    storage_indx += (*graph)->nnbs[i];
	}
	else
	    (*graph)->nbs[i] = NULL;

    if (GMRFLib_verify_graph_read_from_disc)
	if (GMRFLib_validate_graph(stderr, *graph) != GMRFLib_SUCCESS)
	    GMRFLib_ERROR(GMRFLib_EGRAPH);
    
    GMRFLib_prepare_graph(*graph);		  /* prepare the graph for computations */
    return GMRFLib_SUCCESS;
#undef READ_ERROR
}

/*!
  \brief Reads a graph from a file, binary format

  \param[in,out] graph At output, <em>(*graph)</em> has been initialized 
  by using the information on the file \em filename.

  \param[in] filename The name of the file, binary formatted similar to
  \c GMRFLib_read_graph() (without newlines), containing the specification
  of the graph.

  \sa GMRFLib_read_graph, GMRFLib_write_graph_binary
 */
int GMRFLib_read_graph_binary(GMRFLib_graph_tp **graph, char *filename)
{
    
    /*
      read a graph in the following format

       N
       node[0] nnbs[0] nbs[node[0]][0] nbs[node[0]][1] ... nbs[node[0]][nnbs[0]-1] 
       node[1] nnbs[1] nbs[node[1]][0] nbs[node[1]][1] ... nbs[node[1]][nnbs[1]-1] 
       :
       node[N-1] nnbs[N-1] nbs[node[N-1]][0] nbs[node[N-1]][1] ... nbs[node[N-1]][nnbs[N-1]-1] 
     */

#define READ_ERROR() do { \
    if (1){\
        fprintf(stderr,"\n\n\t%s: error: file [%s]:\n",__GMRFLib_FuncName,filename);\
	fprintf(stderr,"\t\tfail to read [%1d] bytes from byte [%1d]\n", nelm*sizeof(int), byte);\
    }    \
    if (fp) fclose(fp);\
    (*graph) = NULL;\
    GMRFLib_ERROR(GMRFLib_EREADGRAPH);} while(0)
 
    int *storage, n_neig_tot = 0, storage_indx;
    int i, byte, nelm, tnode;
    FILE *fp;

    if (!filename) return GMRFLib_SUCCESS;
    if (!(fp = fopen(filename, "r"))) GMRFLib_ERROR(GMRFLib_EOPENFILE);

    GMRFLib_make_empty_graph(graph);

    byte = 0;
    nelm = 1;
    //    if (fread(&((*graph)->n), sizeof(int), nelm, fp) != 1) READ_ERROR();
    byte += sizeof(int);
    
    (*graph)->nnbs  = Calloc((*graph)->n, int);   MEMCHK((*graph)->nnbs);
    (*graph)->nbs   = Calloc((*graph)->n, int *); MEMCHK((*graph)->nbs);

    for(i=0;i< (*graph)->n;i++)
    {
	nelm = 1;
	//	if (fread(&tnode, sizeof(int), nelm, fp) != 1) READ_ERROR(); /* target node */
	byte += sizeof(int);
	if (tnode < 0 || tnode >= (*graph)->n)
	{
	    fprintf(stderr,"\n\n\t%s: error: file [%s]: byte[%1d]\n",
		    __GMRFLib_FuncName, filename, byte); 
	    fprintf(stderr,"\t\tnode-number[%1d] is not in the range[0:%1d]\n", tnode, (*graph)->n);
	    if (fp) fclose(fp);
	    GMRFLib_ERROR(GMRFLib_EPARAMETER);
	}
	    
	nelm = 1;
	//	if (fread(&((*graph)->nnbs[tnode]), sizeof(int), nelm, fp) != 1) READ_ERROR();
	byte += sizeof(int);
	n_neig_tot += (*graph)->nnbs[tnode];
	if ((*graph)->nnbs[tnode])
	{
	    (*graph)->nbs[tnode] = Calloc((*graph)->nnbs[tnode], int); MEMCHK((*graph)->nbs[tnode]);
	    nelm = (*graph)->nnbs[tnode];
	    //   if (fread((*graph)->nbs[tnode], sizeof(int), nelm, fp) != nelm) READ_ERROR();
	    byte += sizeof(int)*((*graph)->nnbs[tnode]);
	}
	else
	    (*graph)->nbs[tnode] = NULL;
    }
    fclose(fp);

    /*
      map the graph to a more computational convenient memory layout!
      use just one long vector to store all the neighbors.
     */
    storage      = Calloc(n_neig_tot, int); MEMCHK(storage);
    storage_indx = 0;
    for(i=0;i< (*graph)->n;i++)
	if ((*graph)->nnbs[i])
	{
	    memcpy(&storage[storage_indx], (*graph)->nbs[i], sizeof(int)* (*graph)->nnbs[i]);
	    free((void *)(*graph)->nbs[i]);
	    (*graph)->nbs[i] = &storage[storage_indx];
	    storage_indx += (*graph)->nnbs[i];
	}
	else
	    (*graph)->nbs[i] = NULL;

    if (GMRFLib_verify_graph_read_from_disc)
	if (GMRFLib_validate_graph(stderr, *graph) != GMRFLib_SUCCESS)
	    GMRFLib_ERROR(GMRFLib_EGRAPH);
    
    GMRFLib_prepare_graph(*graph);		  /* prepare the graph for computations */
    return GMRFLib_SUCCESS;
#undef READ_ERROR
}

/*!
  \brief Prints the specification of a graph to standard
    output or a file
  \param[out] fp The \em FILE* on which to print the graph.
  \param[in] graph The graph to be printed.
  \sa GMRFLib_read_graph
 */
int GMRFLib_print_graph(FILE *fp, GMRFLib_graph_tp *graph)
{

    int i, j;

    if (!fp) return GMRFLib_SUCCESS;

    fprintf(fp, "graph has %1d nodes\n", graph->n);
    for(i=0;i<graph->n;i++)
    {
	fprintf(fp, "node %1d has %1d neighbors:", i, graph->nnbs[i]);
	for(j=0;j<graph->nnbs[i];j++) fprintf(fp, " %1d", graph->nbs[i][j]);
	fprintf(fp, "\n");
    }
    if (graph->mothergraph_idx)
	for(i=0;i<graph->n;i++)
	    fprintf(fp, "node %1d corresponds to mother-node %1d\n", i, graph->mothergraph_idx[i]);
    return GMRFLib_SUCCESS;
}

/*!
  \brief Write a graph to file in ascii format
  \param[in] filename The name of the file to store the graph in
    ascii format.
  \param[in] graph The graph to be written.
  \sa GMRFLib_read_graph.
 */
int GMRFLib_write_graph(char *filename, GMRFLib_graph_tp *graph)
{
    /* 
       write graph to file filename in the format so it can be read by 'read_graph'
    */
    
    int i, j;
    FILE *fp;

    if (!filename) return GMRFLib_SUCCESS;
    if (!graph) return GMRFLib_SUCCESS;
    
    fp = fopen(filename, "w");
    if (!fp) GMRFLib_ERROR(GMRFLib_EOPENFILE);

    fprintf(fp, "%1d\n", graph->n);
    for(i=0;i<graph->n;i++)
    {
	fprintf(fp, "%1d %1d", i, graph->nnbs[i]);
	for(j=0;j<graph->nnbs[i];j++) fprintf(fp, " %1d", graph->nbs[i][j]);
	fprintf(fp, "\n");
    }
    fclose(fp);
    
    return GMRFLib_SUCCESS;
}

/*!
  \brief Write a graph to file, binary format
  \param[in] filename The name of the file to store the graph in
    binary format. Can be read with \c GMRFLib_read_graph_binary().
  \param[in] graph The graph
  \sa GMRFLib_read_graph_binary, GMRFLib_write_graph.
 */
int GMRFLib_write_graph_binary(char *filename, GMRFLib_graph_tp *graph)
{
    /* 
       write graph to file filename in binary format so it can be read by 'read_graph_binary'
    */
    
    int i;
    FILE *fp;

    if (!filename) return GMRFLib_SUCCESS;
    if (!graph) return GMRFLib_SUCCESS;
    
    fp = fopen(filename, "w");
    if (!fp) GMRFLib_ERROR(GMRFLib_EOPENFILE);

    fwrite(&(graph->n), sizeof(int), 1, fp);
    for(i=0;i<graph->n;i++)
    {
	fwrite(&i, sizeof(int), 1, fp);
	fwrite(&(graph->nnbs[i]), sizeof(int), 1, fp);
	if (graph->nnbs[i]) fwrite(graph->nbs[i], sizeof(int), graph->nnbs[i], fp);
    }
    fclose(fp);
    
    return GMRFLib_SUCCESS;
}


/*!
  \brief Free a graph build with \c GMRFLib_read_graph()

  Frees the memory held by a \c GMRFLib_graph_tp -variable. 
  NOTE: To ensure safe memory
  handling, this function should only be applied to graphs generated
  by the graph-generating routines of the library, and NOT to
  user-specified graphs.

  \param[in,out] graph  A graph generated by a graph-generating
  library routine. At output, the graph and it's array members are
  all deallocated.
  \sa GMRFLib_read_graph, GMRFLib_make_linear_graph, 
  GMRFLib_make_lattice_graph, GMRFLib_prepare_graph
 */
int GMRFLib_free_graph(GMRFLib_graph_tp *graph)
{
    /* 
       free a graph build with ``GMRFLib_read_graph''
     */
    int i;

    if (!graph) return GMRFLib_SUCCESS;

    for(i=0;i<graph->n;i++)
	if (graph->nnbs[i])
	{
	    FREE(graph->nbs[i]);		  
	    break;				  /* new memory layout, only `free' the first!!! */
	}
    FREE(graph->nbs);
    FREE(graph->nnbs);
    FREE(graph->mothergraph_idx);
    FREE(graph);

    return GMRFLib_SUCCESS;
}

/*!
  \brief Return the number of non-zero elements in Q.
  \param[in] graph The graph.
 */
int GMRFLib_nQelm(GMRFLib_graph_tp *graph)
{
    /* 
       return the number of non-zero elements in Q
    */
    int nelm, i;

    for(i=0, nelm = graph->n; i < graph->n; i++) nelm += graph->nnbs[i];
    return nelm;
}

/*!
  \brief Return bit-number BITNO, bitno = 0, 1, 2, ..., 7.
 */
int GMRFLib_getbit(GMRFLib_uchar c, unsigned int bitno)
{
    /*
       return bit-number BITNO, bitno = 0, 1, 2, ..., 7
     */
    unsigned int zero = 0;
    return (int)((c >> bitno) & ~(~zero << 1));
}

/*!
  \brief Set bitno 0, 1, 2, ..., or 7, to TRUE.
 */
int GMRFLib_setbit(GMRFLib_uchar *c, unsigned int bitno)
{
    /*
      set bitno 0, 1, 2, ..., 7, to TRUE
    */
    unsigned int zero = 0;

    *c = *c | ((~(~zero << (bitno+1)) & (~zero << bitno)));
    return GMRFLib_SUCCESS;
}

/*!
  \brief Print all the bits in \c c to \c fp
 */
int GMRFLib_printbits(FILE *fp, GMRFLib_uchar c)
{
    /* 
       just print all the bits in C to FP
     */
    
    int j, nn = 8*sizeof(GMRFLib_uchar);
    fprintf(fp,"int=[%u] : ", c);

    for(j=0;j<nn;j++) fprintf(fp,"%1d", GMRFLib_getbit(c, (unsigned int)(nn-j-1)));
    fprintf(fp,"\n");

    return GMRFLib_SUCCESS;
}

/*!
  \brief Checks whether two nodes \em i and \em j are neighbours.
  \param[in] node, nnode The nodes \em i and \em j to be checked.
  \param[in] graph The graph to be checked.
  \return Returns \c GMRFLib_TRUE if the nodes \em node and \em nnode are neighbours, 
          otherwise it returns \c GMRFLib_FALSE.
  \par Example:
  \verbinclude example_doxygen_is_neighbour.txt
 */
int GMRFLib_is_neighb(int node, int nnode, GMRFLib_graph_tp *graph)
{
    /* 
       plain version.  return 1 if nnode is a neighbour of node, otherwise 0. assume that the
       nodes are sorted. (if node == nnode, then they are not neighbours.)

       Sat Nov 22 12:43:09 CET 2003 : I compared this version agains a hash-table variant, both
       using spmatrix and map_ii for each i. the last version was slower by a factor 2, the other
       even slower. it seems like this plain version is quite fast though...
    */
    
    int j, m, k;
    m = graph->nnbs[node];

    /* 
       fast return?
    */
    if (!m)                            return GMRFLib_FALSE;
    if (nnode < graph->nbs[node][0])   return GMRFLib_FALSE;
    if (nnode > graph->nbs[node][m-1]) return GMRFLib_FALSE;

    for(j=0;j<m;j++) 
    {
	k = graph->nbs[node][j];
	if (k > nnode)  return GMRFLib_FALSE;
	if (k == nnode) return GMRFLib_TRUE;
    }
    return GMRFLib_FALSE;
}

/*!
  \brief Prepare the graph by sort the vertices in increasing orders.

  Post-processes a user-specified \c GMRFLib_graph_tp object, such that 
  the resulting graph is of the format required by the library functions 
  operating on graphs.

  \param[in,out] graph  A graph explicitly specified by the user
  by creating and initializing a \c GMRFLib_graph_tp -object.  
  At output, the graph is customized.

  \note If the user spesify its own graph, then this function MUST used to prepare internal
  structures and ensure that they are correct, otherwise, strange errors can occure.

  \sa GMRFLib_read_graph
 */
int GMRFLib_prepare_graph(GMRFLib_graph_tp *graph)
{
    /* 
       prepare the graph by sort the vertices in increasing orders
     */
    GMRFLib_sort_nodes(graph);
    GMRFLib_make_nodes_unique(graph);
    
    return GMRFLib_SUCCESS;
}
int GMRFLib_make_nodes_unique(GMRFLib_graph_tp *graph)
{	
    /* 
       ensure the neigbours are unique. the neigbours must be sorted.
    */

    int i, j, k;
    if (!graph) return GMRFLib_SUCCESS;

    for(i=0;i<graph->n;i++)
	if (graph->nnbs[i])
	{
	    k = 0;
	    for(j=1;j<graph->nnbs[i];j++)
		if (graph->nbs[i][k] != graph->nbs[i][j]) graph->nbs[i][ ++k ] = graph->nbs[i][j];
	    graph->nnbs[i] = k+1;
	}

    return GMRFLib_SUCCESS;
}
int GMRFLib_intcmp(const void *a, const void *b)
{
    int ia, ib;
    ia = *((int *)a);
    ib = *((int *)b);
    if (ia > ib)  return  1;
    if (ia < ib)  return -1;
    return  0;
}
int GMRFLib_sort_nodes(GMRFLib_graph_tp *graph)
{	
    /* 
       sort the vertices in increasing order
    */

    int i;
    if (!graph) return GMRFLib_SUCCESS;
    for(i=0;i<graph->n;i++)
	if (graph->nnbs[i])
	    qsort(graph->nbs[i], (size_t)graph->nnbs[i], sizeof(int), GMRFLib_intcmp);

    return GMRFLib_SUCCESS;
}
int GMRFLib_compute_bandwidth(GMRFLib_graph_tp *graph, int *remap)
{
    int bandwidth=0, i, j, node;

    if (!graph) return GMRFLib_SUCCESS;
    for(i=0;i<graph->n;i++)
    {
	node = remap[i];
	for(j=0;j<graph->nnbs[i];j++)
	    bandwidth = MAX(bandwidth, node - remap[ graph->nbs[i][j] ]);
    }
    return bandwidth;
}
int GMRFLib_find_idx(int *idx, int n, int *iarray, int value)
{
    int i;
    for(i=0;i<n;i++)
	if (iarray[i] == value)
	{
	    if (idx) *idx = i;
	    return GMRFLib_SUCCESS;
	}
    return GMRFLib_EINDEX;
}

/*!
  \brief Validates the graph by checking that the members of a \c GMRFLib_graph_tp -object are defined consistently.

  \param[out] fp If \c !NULL, the function will write error messages stating the type of validation
  error on this file, in addition to the default standard output error message issued by
  GMRFLib_error_handler().  If \c NULL, only the default error message is issued.

  \param[in] graph The graph to be validated.

  \sa GMRFLib_error_handler
 */
int GMRFLib_validate_graph(FILE *fp, GMRFLib_graph_tp *graph)
{
    
    int i, j, jj, error=0;
    
    if (graph->n == 0) return GMRFLib_SUCCESS;
    if (graph->n < 0)
    {
	if (fp) fprintf(fp, "%s: error: graph->n = %1d < 0\n", __GMRFLib_FuncName, graph->n);
	GMRFLib_ERROR(GMRFLib_EPARAMETER);
    }
    
    /* 
       check nbs
     */
    for(i=0;i<graph->n;i++)
	for(j=0;j<graph->nnbs[i];j++)
	{
	    if (graph->nbs[i][j] < 0 || graph->nbs[i][j] >= graph->n)
	    {
		error++;
		if (fp)
		    fprintf(fp, "\n\n%s: error: nbs[%1d][%1d]=[%1d] is out of range\n", __GMRFLib_FuncName, i, j,
			    graph->nbs[i][j]); 
	    }
	    if (graph->nbs[i][j] == i)
	    {
		error++;
		if (fp)
		    fprintf(fp, "\n\n%s: error: nbs[%1d][%1d] = node[%1d]. not allowed!\n", __GMRFLib_FuncName, i, j, i);
	    }
	}
    if (error) return GMRFLib_EGRAPH;

    /* 
       check if i~j then j~i
     */

    for(i=0;i<graph->n;i++)
	for(j=0;j<graph->nnbs[i];j++)
	{
	    jj = graph->nbs[i][j];
	    if (GMRFLib_find_idx(NULL, graph->nnbs[jj], graph->nbs[jj], i) != GMRFLib_SUCCESS)
	    {
		if (fp)
		    fprintf(fp, "\n\n%s: error: node[%1d] has neighbor node[%1d], but not oposite\n",
			    __GMRFLib_FuncName, i, jj);  
		error++;
	    }
	}
    if (error) return GMRFLib_EGRAPH;

    /* 
       ok
     */
    if (0 && fp) fprintf(fp, "%s: graph is OK.\n", __GMRFLib_FuncName);
    return GMRFLib_SUCCESS;
}
int GMRFLib_remap_graph(GMRFLib_graph_tp **ngraph, GMRFLib_graph_tp *graph, int *remap)
{
    /* 
       return the remapped graph based on 'graph'. the returned graph has the identity mapping.
     */
    
    int i, j, k, nnb, indx, *hold = NULL;
    
    if (!graph) 
    {
	*ngraph = (GMRFLib_graph_tp *)NULL;
	return GMRFLib_SUCCESS;
    }
    
    GMRFLib_make_empty_graph(ngraph);
    (*ngraph)->n    = graph->n;
    (*ngraph)->nnbs = Calloc((*ngraph)->n, int);   MEMCHK((*ngraph)->nnbs);
    (*ngraph)->nbs  = Calloc((*ngraph)->n, int *); MEMCHK((*ngraph)->nbs);

    for(i=0;i<graph->n;i++)
    {
	k = remap[i];
	(*ngraph)->nnbs[k] = graph->nnbs[i];
	if ((*ngraph)->nnbs[k])
	{
	    (*ngraph)->nbs[k]  = Calloc((*ngraph)->nnbs[k], int);  MEMCHK((*ngraph)->nbs[k]);
	    for(j=0;j<(*ngraph)->nnbs[k];j++) (*ngraph)->nbs[k][j] = remap[graph->nbs[i][j]];
	}
	else
	    (*ngraph)->nbs[k] = NULL;
    }

    /* 
       rearrange into linear storage and free temporary storage
     */
    for(i=0, nnb=0;i<(*ngraph)->n;i++) nnb += (*ngraph)->nnbs[i];
    if (nnb)
    {
	hold = Calloc(nnb, int); MEMCHK(hold);
    }
    else
	hold = NULL;
    
    for(i=0, indx=0; i<(*ngraph)->n;i++) 
    {
	if ((*ngraph)->nnbs[i])
	{
	    memcpy(&hold[indx], (*ngraph)->nbs[i], (*ngraph)->nnbs[i]*sizeof(int));
	    FREE((*ngraph)->nbs[i]);
	    (*ngraph)->nbs[i] = &hold[indx]; 
	}
	else
	{
	    FREE((*ngraph)->nbs[i]);
	}
	indx += (*ngraph)->nnbs[i];
    }
    GMRFLib_prepare_graph(*ngraph);

    return GMRFLib_SUCCESS;
}

/*!
  \brief Returns a copy of a graph

  \param[out] graph_new A \c GMRFLib_graph_tp -object. At output, it contains a copy of \em graph.
  \param[in] graph_old A \c GMRFLib_graph_tp -object.
 */
int GMRFLib_copy_graph(GMRFLib_graph_tp **graph_new, GMRFLib_graph_tp *graph_old)
{
    if (!graph_old)
    {
	*graph_new = NULL;
	return GMRFLib_SUCCESS;
    }

    /* 
       a bit slow, but...
    */
    GMRFLib_compute_subgraph(graph_new, graph_old, NULL);

    /* 
       this mapping has to be the same if it is used in graph_old.
    */
    if (graph_old->mothergraph_idx)
    {
	assert((*graph_new)->mothergraph_idx);
	memcpy((*graph_new)->mothergraph_idx, graph_old->mothergraph_idx, graph_old->n*sizeof(int));
    }
    else
    {
	/* 
	   should not be used
	*/
	FREE((*graph_new)->mothergraph_idx);
    }
    return GMRFLib_SUCCESS;
}

/*!
  \brief Computes a subgraph of a graph by removing the nodes indicated by the \em remove_flag argument.

  \param[out] subgraph The subgraph generated by removing the nodes for which \em remove_flag is 1
    from the graph \em graph.
  \param[in] graph The graph from which to compute the subgraph.
  \param[in] remove_flag An array of length \em n, the number of nodes in \em graph, specifying
    which nodes to be removed creating the subgraphs. The nodes for which \em remove_flags is \em 0,
    IS included in the subgraph.

  \par Example:
  See \ref ex_graph
 */
int GMRFLib_compute_subgraph(GMRFLib_graph_tp **subgraph, GMRFLib_graph_tp *graph, char *remove_flag)
{

    /* 
       return a subgraph of graph, by ruling out those nodes for which remove_flag[i] is true,
       keeping those which remove_flag[i] false.
    */
    int i, j, nneig, nn, k, n_neig_tot, storage_indx, *sg_iidx, *storage, free_remove_flag=0;

    if (!graph) return GMRFLib_SUCCESS;
    GMRFLib_ASSERT(subgraph, GMRFLib_EINVARG);
    
    /* 
       to ease the interface: remove_flag = NULL is ok. then this just do a (slow) copy of the graph.
     */
    if (!remove_flag)				  
    {
	remove_flag = Calloc(graph->n, char); MEMCHK(remove_flag);
	free_remove_flag = 1;
    }

    GMRFLib_make_empty_graph(subgraph);
    
    for(i=0, nn=0;i<graph->n;i++) nn += (!remove_flag[i]);
    (*subgraph)->n = nn;
    if (!((*subgraph)->n)) return GMRFLib_SUCCESS;

    /* 
       create space
    */
    (*subgraph)->nnbs = Calloc((*subgraph)->n, int); MEMCHK((*subgraph)->nnbs);
    (*subgraph)->nbs  = Calloc((*subgraph)->n, int*); MEMCHK((*subgraph)->nbs);

    (*subgraph)->mothergraph_idx = Calloc((*subgraph)->n, int); MEMCHK((*subgraph)->mothergraph_idx); 
    sg_iidx                      = Calloc(graph->n, int); MEMCHK(sg_iidx);

    /* 
       make the mapping of nodes
    */
    for(i=0, k=0;i<graph->n;i++)
	if (!remove_flag[i])
	{
	    (*subgraph)->mothergraph_idx[k] = i;
	    sg_iidx[i]              = k;
	    k++;
	}
	else
	    sg_iidx[i] = -1;			  /* to force a failure if used wrong */

    /* 
       parse the graph and collect nodes not to be removed.
    */
    for(i=0, k=0, n_neig_tot = 0;i<graph->n;i++)
	if (!remove_flag[i])
	{
	    for(j=0, nneig=0;j<graph->nnbs[i];j++) nneig += (!remove_flag[graph->nbs[i][j]]);
	    n_neig_tot += nneig;
	    (*subgraph)->nnbs[k] = nneig;
	    
	    if (nneig > 0)
	    {
		(*subgraph)->nbs[k] = Calloc(nneig, int); MEMCHK((*subgraph)->nbs[k]);
		for(j=0, nneig=0;j<graph->nnbs[i];j++)
		    if (!remove_flag[graph->nbs[i][j]])
		    {
			(*subgraph)->nbs[k][nneig] = sg_iidx[graph->nbs[i][j]];
			nneig++;
		    }
	    }
	    else
		(*subgraph)->nbs[k] = NULL;

	    k++;
	}

    FREE(sg_iidx);
    
    /*
      map the graph to a more computational convenient memory layout!  use just one long vector to
      store all the neighbors.
    */
    if (n_neig_tot)
    {
	storage = Calloc(n_neig_tot, int);MEMCHK(storage);
	storage_indx = 0;
	for(i=0;i<(*subgraph)->n;i++)
	    if ((*subgraph)->nnbs[i])
	    {
		memcpy(&storage[storage_indx], (*subgraph)->nbs[i], sizeof(int)*(*subgraph)->nnbs[i]);
		FREE((*subgraph)->nbs[i]);
		(*subgraph)->nbs[i] = &storage[storage_indx];
		storage_indx += (*subgraph)->nnbs[i];
	    }
	    else
		(*subgraph)->nbs[i] = NULL;
    }
    
    GMRFLib_prepare_graph(*subgraph);		 

    if (free_remove_flag) FREE(remove_flag);	  /* if we have used our own */
    
    return GMRFLib_SUCCESS;
}
int GMRFLib_convert_to_mapped(double *destination,  double *source, GMRFLib_graph_tp *graph, int *remap)
{
    /* 
       convert from the real-world to the mapped world. source might be NULL.
     */
    
    int i;
    static int len      = 0;			  /* keep storage, never free! */
    static double *work = NULL;
    
    GMRFLib_ASSERT(destination, GMRFLib_EINVARG);
    GMRFLib_ASSERT(graph, GMRFLib_EINVARG);
    if (!(graph->n)) return GMRFLib_SUCCESS;

    if (destination && source && destination != source)
	for(i=0;i<graph->n;i++) destination[remap[i]] = source[i];
    else
    {
	if (graph->n > len)			  /* need more space? */
	{
	    /* 
	       some memchecker libs complain (default) about realloc a NULL pointer, allthough the
	       manual page says its ok.
	    */
	    
	    if (work)
		work = (double *)realloc(work, graph->n*sizeof(double));
	    else
		work = (double *)malloc(graph->n*sizeof(double));
	    MEMCHK(work);
	    len  = graph->n;
	}

	memcpy(work, destination, graph->n*sizeof(double));
	for(i=0;i<graph->n;i++) destination[remap[i]] = work[i];
    }
    return GMRFLib_SUCCESS;
}
int GMRFLib_convert_from_mapped(double *destination,  double *source, GMRFLib_graph_tp *graph, int *remap)
{
    /* 
       convert from the mapped-world to the real world. source might be NULL.
     */

    int i;
    static int len      = 0;			  /* keep storage, never free! */
    static double *work = NULL;
    
    GMRFLib_ASSERT(destination, GMRFLib_EINVARG);
    GMRFLib_ASSERT(graph, GMRFLib_EINVARG);
    if (!(graph->n)) return GMRFLib_SUCCESS;

    if (destination && source && destination != source)
	for(i=0;i<graph->n;i++) destination[i] = source[remap[i]];
    else
    {
	if (graph->n > len)			  /* need more space? */
	{
	    if (work)
		work = (double *)realloc(work, graph->n*sizeof(double));
	    else
		work = (double *)malloc(graph->n*sizeof(double));
	    MEMCHK(work);
	    len  = graph->n;
	}

	memcpy(work, destination, graph->n*sizeof(double));
	for(i=0;i<graph->n;i++) destination[i] = work[remap[i]];
    }
    return GMRFLib_SUCCESS;
}

/*!
  \brief  Compute \f$ \mbox{\boldmath $Q$}\mbox{\boldmath $x$} \f$.
  
  Computes \f$ \mbox{\boldmath $Q$} \mbox{\boldmath $x$} \f$, the matrix-vector product of the
  precision matrix <em>\b Q</em> and an array <em>\b x</em>. The elements of <em>\b x</em> should be
  in the original ordering of the nodes of the graph.

  \param[out] result At output, the length \em n array \em result contains 
  the value of the product  \f$ \mbox{\boldmath $Q$}\times \mbox{\boldmath $x$} \f$.
  \param[in] x The array <em>\b x</em>, as a length \em n array.
  \param[in] graph The graph on which <em>\b x</em> is defined.
  \param[in] Qfunc The function defining the precision matrix <em>\b Q</em>.
  \param[in] Qfunc_arg A \em character -pointer defining the
    address of a variable or data structure holding the arguments to
    the function \em Qfunc.

  \sa GMRFLib_Qfunc_tp, GMRFLib_xQx
 */
int GMRFLib_Qx(double *result, double *x,
	       GMRFLib_graph_tp *graph, GMRFLib_Qfunc_tp *Qfunc, char *Qfunc_arg)
{
    /* 
       compute result = Q*x.
     */
    
    int i, j, jj;

    GMRFLib_ASSERT(result, GMRFLib_EINVARG);
    GMRFLib_ASSERT(x, GMRFLib_EINVARG);
    GMRFLib_ASSERT(graph, GMRFLib_EINVARG);
    GMRFLib_ASSERT(Qfunc, GMRFLib_EINVARG);
    
    memset(result, 0, graph->n*sizeof(double));
    for(i=0;i<graph->n;i++)
    {
	result[i] += (*Qfunc)(i, i, Qfunc_arg)*x[i];
	for(j=0;j<graph->nnbs[i];j++)
	{
	    jj = graph->nbs[i][j];
	    result[i] += (*Qfunc)(i, jj, Qfunc_arg)*x[jj];
	}
    }
    return GMRFLib_SUCCESS;
}
/*!
  \brief  Compute \f$ \mbox{\boldmath $x^TQx$} \f$.
  
  Returns \f$ \mbox{\boldmath $x^TQx$} \f$: the scalar product of <em>\b x</em>, and the
  matrix-vector product of the precision matrix <em>\b Q</em> and an array <em>\b x</em>. The
  elements of <em>\b x</em> should be in the original ordering of the nodes of the graph.

  \param[in] x The array <em>\b x</em>, as a length \em n array.
  \param[in] graph The graph on which <em>\b x</em> is defined.
  \param[in] Qfunc The function defining the precision matrix <em>\b Q</em>.
  \param[in] Qfunc_arg A \em character -pointer defining the
    address of a variable or data structure holding the arguments to
    the function \em Qfunc.

  \sa GMRFLib_xQx
 */
double GMRFLib_xQx(double *x, GMRFLib_graph_tp *graph, GMRFLib_Qfunc_tp *Qfunc, char *Qfunc_arg)
{
    int i;
    double *y, res;

    if (!graph || graph->n <= 0 || !x) return 0.0;
    y = Calloc(graph->n, double); MEMCHK_RETVAL(y, -1.0);

    GMRFLib_Qx(y, x, graph, Qfunc, Qfunc_arg);
    for(i=0, res=0.0;i<graph->n;i++) res += y[i]*x[i];

    return res;
}
    
/*!  \brief Creates a graph on a lattice, customizing the graph for use in other library functions.

  \param[out] graph At output, <em>(*graph)</em> contains
  the specification of the lattice graph, customized for use in
  other library routines.
  \param[in] nrow, ncol The number of rows ( \f$ n_{row} \f$ ) and columns
  ( \f$ n_{col} \f$ ) in the \f$ n_{row}\times n_{col} \f$ lattice.
  \param[in] nb_row, nb_col Specification of the parameters \f$ nb_{row} \f$ 
  and \f$ nb_{col} \f$ of the \f$ (2 nb_{row}+1)\times(2 nb_{col}+1) \f$
  -neighbourhood of the lattice.
  \param[in] cyclic_flag If this argument is <em>>0</em>, the graph is made
  cyclic.

  \par Example:
  See \ref ex_graph

  \sa GMRFLib_lattice2node, GMRFLib_node2lattice
 */
int GMRFLib_make_lattice_graph(GMRFLib_graph_tp **graph,
			       int nrow, int ncol, int nb_row, int nb_col, int cyclic_flag)
{
    /* 
       make an lattice graph, with nodes from 0...n-1, and a (2 x nb_row +1) x (2 x nb_col + 1)
       neighborhood. if `cyclic_flag`, then make the graph cyclic. the prune_graph function might be
       useful to derive non-square neighborhoods.
    */
    int *hold = NULL, n, nnb, ir, ic, node, nnode, irow, icol;

    GMRFLib_ASSERT(nrow > 0, GMRFLib_EINVARG);
    GMRFLib_ASSERT(ncol > 0, GMRFLib_EINVARG);
    GMRFLib_ASSERT(nb_row >= 0, GMRFLib_EINVARG);
    GMRFLib_ASSERT(nb_col >= 0, GMRFLib_EINVARG);
    
    nb_row = MIN(nrow -1, nb_row);		  /* otherwise (i+nrow)%nrow will fail if cyclic */
    nb_col = MIN(ncol -1, nb_col);		  /* same here */
    n      = ncol*nrow;
    nnb    = (2*nb_row + 1)*(2*nb_col + 1);

    GMRFLib_make_empty_graph(graph);
    (*graph)->n    = n;
    (*graph)->nnbs = Calloc(n, int);  MEMCHK((*graph)->nnbs);
    (*graph)->nbs  = Calloc(n, int *);MEMCHK((*graph)->nbs);
    
    /* 
       if nh is large compared to n, the this graph may contain double (or more) occurences of one
       neigbor, but the `prepare_graph' function fix this.
    */
    if (nnb)
    {
	hold=Calloc(n*nnb, int); MEMCHK(hold);	  /* use a linear storage */
	for(node=0;node<n;node++) (*graph)->nbs[node] = &hold[node*nnb]; /* set pointers to it */

	for(node=0;node<n;node++)
	{
	    GMRFLib_node2lattice(node, &irow, &icol, nrow, ncol);
	    if (cyclic_flag)
	    {
		for(ir = irow - nb_row; ir <= irow + nb_row; ir++)
		    for(ic = icol - nb_col; ic <= icol + nb_col; ic++)
		    {
			GMRFLib_lattice2node(&nnode, (ir+nrow)%nrow, (ic+ncol)%ncol, nrow, ncol);
			if (nnode != node) (*graph)->nbs[node][ (*graph)->nnbs[node]++ ] = nnode;
		    }
	    }
	    else
	    {
		for(ir = MAX(0, irow - nb_row); ir <= MIN(nrow-1, irow + nb_row); ir++)
		    for(ic = MAX(0, icol - nb_col); ic <= MIN(ncol-1, icol + nb_col); ic++)
		    {
			GMRFLib_lattice2node(&nnode, ir, ic, nrow, ncol);
			if (nnode != node) (*graph)->nbs[node][ (*graph)->nnbs[node]++ ] = nnode;
		    }
	    }
	}
    }
    else
	for(node=0;node<n;node++) (*graph)->nbs[node] = NULL; /* if zero neighbors, the pointer should be zero */

    GMRFLib_prepare_graph(*graph);

    return GMRFLib_SUCCESS;
}
/*!
  \brief  Return the node number of a point in a lattice.

  In the case of a rectangular neighbourhood, the ordering of the
  nodes will be the one minimizing the bandwidth of the precision
  matrix of a GMRF on the graph.
  
  \param[out] node The node index \f$ (\in [0,n_{row}\times n_{col}-1]) \f$
    corresponding to the row and column indices of the graph.
  \param[in] irow Row index \f$ (\in [0,n_{row}-1]) \f$.
  \param[in] icol Column index \f$ (\in [0,n_{col}-1]) \f$.
  \param[in] nrow The number of rows of the lattice.
  \param[in] ncol The number of columns of the lattice.

  \remarks The node index \em node is computed so as to minimize the 
  bandwidth of the precision matrix <em>\b Q</em> of a GMRF on the graph 
  in the case of a rectangular neighbourhood, and will depend on the 
  relative sizes of \em nrow and \em ncol: \n
  If <em>(nrow > ncol)</em>, then <em>node = irow + icol*nrow </em>\n
  If <em>(nrow <= ncol)</em>, then <em>node = irow + icol*nrow </em>

  \par Example:
  See \ref ex_graph
  \sa GMRFLib_node2lattice, GMRFLib_make_lattice_graph
 */
int GMRFLib_lattice2node(int *node, int irow, int icol, int nrow, int ncol)
{
    /* 
       this defines the ordering of the nodes from the lattice. it does not depend on the nb_col and
       nb_row, but only on the size.
     */
    *node = (nrow > ncol ? irow + icol*nrow : icol + irow*ncol);
    return GMRFLib_SUCCESS;
}

/*!
  \brief  Compute the lattice indices of a node
  
  \param[in] node The node index \f$ (\in [0,n_{row}\times n_{col}-1]) \f$
    corresponding to the row and column indices of the graph.
  \param[out] irow Row index \f$ (\in [0,n_{row}-1]) \f$ corresponding
  to node number \em node in the graph.
  \param[out] icol Column index \f$ (\in [0,n_{col}-1]) \f$ corresponding
  to node number \em node in the graph.
  \param[in] nrow The number of rows of the lattice.
  \param[in] ncol The number of columns of the lattice.

  \remarks The lattice indices \em irow and \em icol are computed so as 
  to minimize the bandwidth of the precision matrix <em>\b Q</em> of a 
  GMRF on the graph in the case of a rectangular neighbourhood. They will 
  depend on the relative sizes of \em nrow and \em ncol: \n

  If <em>(nrow > ncol)</em>, then <em>irow = node - icol*nrow, 
  icol = node/nrow </em>\n
  If <em>(nrow <= ncol)</em>, then <em>irow = node/ncol, 
  icol = node - irow*ncol </em>

  \par Example:
  See \ref ex_graph
  \sa GMRFLib_lattice2node, GMRFLib_make_lattice_graph
 */
int GMRFLib_node2lattice(int node, int *irow, int *icol, int nrow, int ncol)
{
    /* 
       this is the inverse function of GMRFLib_lattice2node.
     */
    if (nrow > ncol)
    {
	*icol = node/nrow;
	*irow = node - (*icol)*nrow;
    }
    else
    {
	*irow = node/ncol;
	*icol = node - (*irow)*ncol;
    }
    return GMRFLib_SUCCESS;
}

/*!

  \brief Returns a new graph which is a copy of graph but where entries where the elements of the
  precision matrix <em>Q(i,j)=0, i~j,</em> are removed.

  \param[out] new_graph A \c GMRFLib_graph_tp -object. At
  output, it contains a copy of \em graph but all entries
  where <em>Q(i,j)=0</em> are removed.
  \param[in] graph At output, <em>(*graph)</em> has been initialized 
  by using the information on the file \em filename. 
  \param[in] Qfunc A function defining the elements of the \em \b Q -matrix.
  \param[in] Qfunc_arg A \em character -pointer defining the
  address of a variable or data structure holding the arguments to
  the function \em Qfunc.

  \sa GMRFLib_init_wa_problem, GMRFLib_init_nwa_problem.
 */
int GMRFLib_prune_graph(GMRFLib_graph_tp **new_graph, GMRFLib_graph_tp *graph,
			GMRFLib_Qfunc_tp *Qfunc, char *Qfunc_arg)
{

    /* 
       return a modified graph by removing entries where Q(i,j)=0, i~j.
     */
    int i, j, k, ii, *free_ptr = NULL, found;
    
    GMRFLib_ASSERT(graph, GMRFLib_EINVARG);
    GMRFLib_ASSERT(Qfunc, GMRFLib_EINVARG);
    
    GMRFLib_copy_graph(new_graph, graph);

    /* 
       this is a bit tricky. as GMRFLib_free_graph free's the first ptr where nnbs[i]>0, then we
       must make sure that this is the same ptr after pruning, as we modify `new_graph'. so, we need
       to store that ptr and make sure its ok after pruning.
     */
    for(i=0;i<(*new_graph)->n;i++)
	if ((*new_graph)->nnbs[i])
	{
	    free_ptr = (*new_graph)->nbs[i];
	    break;
	}

    for(i=0;i<graph->n;i++)
	if (graph->nnbs[i])
	{
	    for(j=0, k=0;j<graph->nnbs[i];j++)
	    {
		ii = graph->nbs[i][j];
		if ((*Qfunc)(i, ii, Qfunc_arg) != 0.0)
		{
		    (*new_graph)->nbs[i][k] = (*new_graph)->nbs[i][j];
		    k++;
		}
	    }
	    (*new_graph)->nnbs[i] = k;
	    if ((*new_graph)->nnbs[i] == 0) (*new_graph)->nbs[i] = NULL;
	}
    GMRFLib_prepare_graph(*new_graph);

    for(i=0, found = 0;i<(*new_graph)->n;i++)
	if ((*new_graph)->nnbs[i])
	{
	    if (graph->nbs[i] != free_ptr)
	    {
		int j; /* *pos = &((*new_graph)->nnbs[i]); */

		if (0) GMRFLib_msg(stdout, "\n\nNEW CODE HERE, HOPE ITS OK\n\n");

		for(j=0;j< (*new_graph)->nnbs[i]; j++) free_ptr[j] = (*new_graph)->nbs[i][j];
		(*new_graph)->nbs[i] = free_ptr;
	    }
	    found = 1;
	    break;				  
	}
    if (!found)
    {
	/* 
	   then all neighbors are gone, and we have to free free_ptr ourself.
	 */
	FREE(free_ptr);
    }

    return GMRFLib_SUCCESS;
}

/*!
  \brief Creates a linear graph, corresponding to an autoregressive AR(\em p)-model

  Make an linear graph, with nodes from <em>0...n-1</em>, and \em bw neighbors in each
  direction. hence, <em>bw = p</em>, makes the graph for an auto-regressive process or order \em
  p. if `cyclic_flag`, then make the graph cyclic.

  \param graph At output, (*graph) contains the specifications of the 
  linear graph, customized for use in other library routines.
  \param n The number of nodes in the graph.
  \param bw Specifies the bandwidth of the graph, such that each node 
  has \em bw neighbours in each direction. A bandwidth of \em bw 
  corresponds to an AR(\em bw)-model.
  \param cyclic_flag If this argument is > 0, the graph is made cyclic.

  \par Example:
  See \ref ex_graph
 */
int GMRFLib_make_linear_graph(GMRFLib_graph_tp **graph, int n, int bw, int cyclic_flag)
{
    int i, j, k, *hold = NULL;

    GMRFLib_ASSERT(n > 0, GMRFLib_EINVARG);
    GMRFLib_ASSERT(bw >= 0, GMRFLib_EINVARG);

    bw = MIN(n-1, MAX(0, bw));
    GMRFLib_make_empty_graph(graph);
    (*graph)->n    = n;
    (*graph)->nnbs = Calloc(n, int);  MEMCHK((*graph)->nnbs)
    (*graph)->nbs  = Calloc(n, int *);MEMCHK((*graph)->nbs);

    if (bw)
    {
	hold = Calloc(n*2*bw, int); MEMCHK(hold); /* use a linear storage */
	for(i=0;i<n;i++) (*graph)->nbs[i] = &hold[i*2*bw];	/* set pointers to it */

	if (cyclic_flag)
	{
	    for(i=0;i<n;i++)
	    {
		(*graph)->nnbs[i] = 2*bw;
		for(j=i-bw, k=0;j<=i+bw;j++) if (j != i) (*graph)->nbs[i][k++] = MOD(j, n);
	    }
	}
	else
	{
	    for(i=0;i<n;i++)
	    {
		(*graph)->nnbs[i] = (i - MAX(i-bw, 0)) + (MIN(n-1, i+bw) - i);
		for(j=i-bw, k=0;j<=i+bw;j++) if (j != i && LEGAL(j, n)) (*graph)->nbs[i][k++] = j;
		if (!k) (*graph)->nbs[i] = NULL;	  /* if zero neighbors, the pointer should be zero */
	    }
	}
    }
    else
	for(i=0;i<n;i++) (*graph)->nbs[i] = NULL;	  /* if zero neighbors, the pointer should be zero */

    GMRFLib_prepare_graph(*graph);
    return GMRFLib_SUCCESS;
}

/*!
  \brief Make the graph suitable to the RW1 or RW2 model with regular locations defined in def.

  \param[out] graph  The graph for the RW1 or RW2 model with regular locations

  \param[in] def The definition of the RW1 or RW2 model with regular locations

  \remark \c GMRFLib_make_rw_graph() is to be used in connection with \c GMRFLib_rw() both using the
  RW1 or RW2 model defined using \c GMRFLib_rwdef_tp.

  \remark There is an alternative set of tools for the case where the locations are irregular, see
  \c GMRFLib_make_crw_graph(), \c GMRFLib_crw() and \c GMRFLib_crwdef_tp.

  \sa GMRFLib_rw, GMRFLib_rwdef_tp
*/
int GMRFLib_make_rw_graph(GMRFLib_graph_tp **graph, GMRFLib_rwdef_tp *def)
{
    GMRFLib_ASSERT(def, GMRFLib_EINVARG);
    GMRFLib_ASSERT(def->order == 1 || def->order == 2,  GMRFLib_EINVARG);
    GMRFLib_ASSERT(def->n > def->order,  GMRFLib_EINVARG);
    GMRFLib_ASSERT(graph, GMRFLib_EINVARG);
 
    GMRFLib_make_linear_graph(graph, def->n, def->order, def->cyclic);
    return GMRFLib_SUCCESS;
}

/*!
  \brief Make the graph suitable to the (C)RW1 or CRW2 model with irregular locations defined in def.
  
  \param[out] graph  The graph for the (C)RW1 or CRW2 model with irregular locations
  
  \param[in] def The definition of the (C)RW1 or CRW2 model with irregular locations
  
  \remark \c GMRFLib_make_crw_graph() is to be used in connection with \c GMRFLib_crw() both using
  the (C)RW1 or CRW2 model defined using \c GMRFLib_crwdef_tp.
    
  \remark There is an alternative set of tools for the case where the locations are regular, see \c
  GMRFLib_make_rw_graph(), \c GMRFLib_rw() and \c GMRFLib_rwdef_tp.

  \sa GMRFLib_crw(), GMRFLib_crwdef_tp
*/
int GMRFLib_make_crw_graph(GMRFLib_graph_tp **graph, GMRFLib_crwdef_tp *def)
{
    int i, *hold, n;
    GMRFLib_graph_tp *gg;
    
    GMRFLib_ASSERT(def, GMRFLib_EINVARG);
    GMRFLib_ASSERT(def->order == 1 || def->order == 2,  GMRFLib_EINVARG);
    GMRFLib_ASSERT(def->n > def->order,  GMRFLib_EINVARG);
    GMRFLib_ASSERT(def->layout == GMRFLib_CRW_LAYOUT_PAIRS
		   || def->layout == GMRFLib_CRW_LAYOUT_BLOCK
		   || def->layout == GMRFLib_CRW_LAYOUT_SIMPLE,  GMRFLib_EINVARG);
    GMRFLib_ASSERT(graph, GMRFLib_EINVARG);

    n = def->n;
    
    if (def->order == 1 || (def->order == 2 && def->layout == GMRFLib_CRW_LAYOUT_SIMPLE))
    {
	GMRFLib_make_linear_graph(graph, n, def->order, 0);
	return GMRFLib_SUCCESS;
    }
    if (def->order == 2 && def->layout == GMRFLib_CRW_LAYOUT_PAIRS)
    {
	GMRFLib_make_linear_graph(graph, 2*n, 3, 0);
	return GMRFLib_SUCCESS;
    }

    /* 
       the rest is for the block'ed case for order 2. this is the tricky bit
    */
    GMRFLib_make_empty_graph(&gg);
    gg->n    = 2*n;
    gg->nnbs = Calloc(2*n, int);   MEMCHK(gg->nnbs);
    gg->nbs  = Calloc(2*n, int *); MEMCHK(gg->nnbs);
    hold     = Calloc(2*n*5, int); MEMCHK(hold);
    for(i=0;i<2*n;i++) gg->nbs[i] = &hold[i*5];	  /* this is ok */

    for(i=0;i<n;i++)
    {
	if (i == 0)
	{
	    gg->nnbs[i]     = 3;
	    gg->nbs[i][0]   = i+n;		  /* vel */
	    gg->nbs[i][1]   = i+1;		  /* next pos */
	    gg->nbs[i][2]   = i+1+n;		  /* next vel */

	    gg->nnbs[i+n]   = 3;
	    gg->nbs[i+n][0] = i;		  /* pos */
	    gg->nbs[i+n][1] = i+1;		  /* next pos */
	    gg->nbs[i+n][2] = i+1+n;		  /* next vel */
	}
	else if (i == n-1)
	{
	    gg->nnbs[i]     = 3;
	    gg->nbs[i][0]   = i+n;		  /* vel */
	    gg->nbs[i][1]   = i-1;		  /* prev pos */
	    gg->nbs[i][2]   = i-1+n;		  /* prev vel */

	    gg->nnbs[i+n]   = 3;
	    gg->nbs[i+n][0] = i;		  /* pos */
	    gg->nbs[i+n][1] = i-1;		  /* prev pos */
	    gg->nbs[i+n][2] = i-1+n;		  /* prev vel */
	}
	else
	{
	    gg->nnbs[i]   = 5;
	    gg->nbs[i][0] = i+n;		  /* vel */
	    gg->nbs[i][1] = i+1;		  /* next pos */
	    gg->nbs[i][2] = i+1+n;		  /* next vel */
	    gg->nbs[i][3] = i-1;		  /* prev pos */
	    gg->nbs[i][4] = i-1+n;		  /* prev vel */

	    gg->nnbs[i+n]   = 5;
	    gg->nbs[i+n][0] = i;		  /* pos */
	    gg->nbs[i+n][1] = i+1;		  /* next pos */
	    gg->nbs[i+n][2] = i+1+n;		  /* next vel */
	    gg->nbs[i+n][3] = i-1;		  /* prev pos */
	    gg->nbs[i+n][4] = i-1+n;		  /* prev vel */
	}
    }
    GMRFLib_prepare_graph(gg);
    *graph = gg;

    return GMRFLib_SUCCESS;
}

/* NOT DOCUMENTED
   
  \brief Internal function in \c GMRFLib_nfold_graph().

  Returns a new graph created by expanding one graph \em g with another 
  graph \em  gg. The two graphs must have the same number of nodes.

  \param[out] ng The new expanded graph. The pointer <em>(*ng)</em> will be 
  allocated within the routine.
  \param[in] g, gg The graphs to be combined.

  \remarks Given two graphs \em g and \em gg, the new graph \em ng is 
  created by expanding the neighbourhood of \em g by the neighbourhood 
  of \em gg. That is, for each neighbour \f$ j_g \f$ of a node \f$ i_g \f$ 
  of \em g, a node \f$ k_{gg} \f$ in \em gg is included in the 
  neighbourhood of node \f$ i_g \f$ in \em g if \f$ k_{gg} \f$ and 
  \f$ j_g \f$ are neighbours in \em gg. The new graph \em ng is equal to the
  expanded version of of \em g.

  \sa GMRFLib_nfold_graph
 */
int GMRFLib_fold_graph(GMRFLib_graph_tp **ng, GMRFLib_graph_tp *g, GMRFLib_graph_tp *gg)
{
    /* 
       return ng = g * gg, meaning that ng is g expanded by gg
     */
    int i, ii, j, jj, k, kk, nneig, nnb, *hold, indx;
    GMRFLib_graph_tp *newg = NULL;

    /* 
       first the easy cases
     */
    if (!g && !gg)
    {
	*ng = NULL;
	return GMRFLib_SUCCESS;
    }
    if ((g && !gg) || (gg && !g))
    {
	GMRFLib_copy_graph(ng, (g ? g : gg));
	return GMRFLib_SUCCESS;
    }

    GMRFLib_ASSERT(g->n == gg->n, GMRFLib_EPARAMETER); /* of same size? */

    GMRFLib_make_empty_graph(&newg);
    newg->n    = g->n;
    newg->nnbs = Calloc(g->n, int);   MEMCHK(newg->nnbs);
    newg->nbs  = Calloc(g->n, int *); MEMCHK(newg->nbs);

    for(i=0;i<newg->n;i++)
    {
	/* 
	   count number of neighbors neighbors
	 */
	for(j=0, nneig=g->nnbs[i];j<g->nnbs[i];j++) nneig += gg->nnbs[ g->nbs[i][j] ];
	if (nneig)
	{
	    newg->nbs[i]  = Calloc(nneig, int); MEMCHK(newg->nbs[i]);
	    newg->nnbs[i] = 0;
	    for(j=0, k=0;j<g->nnbs[i];j++)
	    {
		jj = newg->nbs[i][k++] = g->nbs[i][j];
		for(ii=0;ii<gg->nnbs[jj];ii++) 
		{
		    kk = gg->nbs[jj][ii];
		    if (kk != i) newg->nbs[i][k++] = kk;
		}
	    }
	    newg->nnbs[i] = k;

	    /* 
	       make them unique
	     */
	    qsort((void *)newg->nbs[i], (size_t)newg->nnbs[i], sizeof(int), GMRFLib_intcmp);
	    for(j=k=0;j<newg->nnbs[i];j++)
		if (j == 0 || newg->nbs[i][j] != newg->nbs[i][k-1]) newg->nbs[i][k++] = newg->nbs[i][j];
	    newg->nnbs[i] = k;
	}
	else
	{
	    newg->nnbs[i] = 0;
	    newg->nbs[i]  = NULL;
	}
    }

    /* 	    
       rearrange into linear storage and free temporary storage
     */
    for(i=0, nnb=0;i<newg->n;i++) nnb += newg->nnbs[i];

    if (nnb)
    {
	hold = Calloc(nnb, int); MEMCHK(hold);
    }
    else
	hold = NULL;

    for(i=0, indx=0; i<newg->n;i++) 
    {
	if (newg->nnbs[i])
	{
	    memcpy(&hold[indx], newg->nbs[i], newg->nnbs[i]*sizeof(int));
	    FREE(newg->nbs[i]);
	    newg->nbs[i] = &hold[indx]; 
	}
	else
	{
	    FREE(newg->nbs[i]);
	}
	indx += newg->nnbs[i];
    }
    GMRFLib_prepare_graph(newg);		  
    *ng = newg;

    return GMRFLib_SUCCESS;
}

/*!
  \brief Returns a new graph created by multiple expanding of one graph by itself.

  <tt>nfold = 0: ng = I \n
  nfold = 1: ng = og \n
  nfold = 2: ng = og*og \n
  nfold = 3: ng = og*og*og \n
  etc </tt>

  \param[out] ng The new expanded graph. The pointer <em>(*ng)</em> will be 
  allocated within the routine.
  \param[in] og The graph to be expanded.
  \param[in] nfold The number of times to expand the graph.

  \remarks Given the old graph \em og, the new graph \em ng is created 
  by expanding the neighbours of \em og by it's own neighbourhood 
  <em>nfold-1</em> times. The function GMRFLib_fold_graph is called 
  successively <em>nfold-1</em> times, the first time using 
  <em>g = gg = og</em>, and in successive runs letting <em>g = gprev</em> 
  and <em>gg = og</em>, where \em gprev is the resulting graph of 
  the previous run.  For completeness, <em>nfold = 0</em> is allowed for, 
  meaning that a graph with \em n independent nodes is returned as 
  \em ng, where \em n is the number of nodes in the old graph \em og. If
  <em>nfold = 1</em>, the new graph is a copy of the old graph.

  \par Example
  See \ref ex_graph
 */
int GMRFLib_nfold_graph(GMRFLib_graph_tp **ng, GMRFLib_graph_tp *og, int nfold)
{
    /* 
       make new graph, 'ng', that is 'nfold' of 'og'

       nfold = 0: ng = I
       nfold = 1: ng = og
       nfold = 2: ng = og*og
       nfold = 3: ng = og*og*og
       etc

     */
    int i;
    GMRFLib_graph_tp *newg, *oldg;
    
    GMRFLib_ASSERT(nfold >= 0, GMRFLib_EPARAMETER);

    if (nfold == 0)
    {
	GMRFLib_make_empty_graph(&newg);
	newg->n    = og->n;
	newg->nnbs = Calloc(newg->n, int); MEMCHK(newg->nnbs);
	newg->nbs  = Calloc(newg->n, int *); MEMCHK(newg->nbs);
	GMRFLib_prepare_graph(newg);
    }
    else if (nfold == 1)
    {
	GMRFLib_copy_graph(&newg, og);
    }
    else
    {
	for(i=0, oldg = newg = NULL; i<nfold; i++)
	{
	    GMRFLib_fold_graph(&newg, oldg, og);
	    GMRFLib_free_graph(oldg);
	    if (i < nfold -1)
	    {
		oldg = newg;
		newg = NULL;
	    }
	}
    }

    *ng = newg;
    return GMRFLib_SUCCESS;
}

/*! 
  \brief Return a new graph which is the union of n_graphs graphs.

  <em>i~j</em> in union_graph, if <em>i~j</em> in
  <em>graph_array[0]...graph_array[n_graphs-1]</em> \n\n
*/
int GMRFLib_union_graph(GMRFLib_graph_tp **union_graph, GMRFLib_graph_tp **graph_array, int n_graphs)
{
    /* 
       return a new graph which is the union of n_graphs graphs: i~j in union_graph, if i~j in
       graph_array[0]...graph_array[n_graphs-1] 
    */
    
    int i, k, node, nnbs, idx, *hold, hold_idx;
    GMRFLib_graph_tp *tmp_graph;
    
    if (!graph_array || n_graphs <= 0) { *union_graph = NULL; return GMRFLib_SUCCESS; }
    for(k=1;k<n_graphs;k++) GMRFLib_ASSERT(graph_array[0]->n == graph_array[k]->n,  GMRFLib_EPARAMETER);
    
    GMRFLib_make_empty_graph(union_graph);
    (*union_graph)->n    = graph_array[0]->n;
    (*union_graph)->nnbs = Calloc((*union_graph)->n, int);  MEMCHK((*union_graph)->nnbs);
    (*union_graph)->nbs  = Calloc((*union_graph)->n, int*); MEMCHK((*union_graph)->nbs);

    for(node=0, nnbs=0;node<(*union_graph)->n;node++)
	for(k=0; k<n_graphs; k++) nnbs += graph_array[k]->nnbs[node];

    if (nnbs)
    {
	hold = Calloc(nnbs, int); MEMCHK(hold);
    }
    else
	hold = NULL;
    
    for(node=0, hold_idx=0;node<(*union_graph)->n;node++)
    {
	for(k=0, nnbs=0; k<n_graphs; k++) nnbs += graph_array[k]->nnbs[node];
	if (nnbs)
	{
	    (*union_graph)->nnbs[node] = nnbs;		  /* will include multiple counts */
	    (*union_graph)->nbs[node]  = &hold[hold_idx];

	    for(k=0, idx=0;k<n_graphs;k++)
		for(i=0;i<graph_array[k]->nnbs[node];i++)
		    (*union_graph)->nbs[node][idx++] = graph_array[k]->nbs[node][i];

	    qsort((*union_graph)->nbs[node], (size_t)(*union_graph)->nnbs[node], sizeof(int), GMRFLib_intcmp);
	    for(i=1, nnbs=0;i<(*union_graph)->nnbs[node];i++)
		if ((*union_graph)->nbs[node][i] != (*union_graph)->nbs[node][nnbs])
		    (*union_graph)->nbs[node][++nnbs] = (*union_graph)->nbs[node][i];
	    (*union_graph)->nnbs[node] = nnbs+1;
	    hold_idx += (*union_graph)->nnbs[node];
	}
	else
	{
	    (*union_graph)->nnbs[node] = 0;
	    (*union_graph)->nbs[node]  = NULL;
	}
    }
    GMRFLib_prepare_graph(*union_graph);	  /* this is required */
    
    /* 
       the union_graph is now (probably) to large as it acounts for multiple counts. the easiest way
       out of this, is to make (and use) a new copy, then free the current one.
     */
    GMRFLib_copy_graph(&tmp_graph, *union_graph);
    GMRFLib_free_graph(*union_graph);
    *union_graph = tmp_graph;
    
    return GMRFLib_SUCCESS;
}

/*!
  \brief Add missing neighbors to a uncomplete or malformed graph; if <em>i~j</em> but <em>j!~i</em>, then add <em>j~i</em>.

  \param[out] n_graph A pointer to a new completed graph, buildt
    from \em graph.
  \param[in] graph The incomplete graph.

  \note If \em graph is ok in the sense that <em>i~j</em>
    then <em>j~i</em>, then this routine just returns a copy of the
    graph, allthough it will be a bit slow.
 */
int GMRFLib_complete_graph(GMRFLib_graph_tp **n_graph, GMRFLib_graph_tp *graph)
{
    /* 
       return a new graph that is complete: if i~j but j!~i, then add j~i
    */
    
    int i, ii, j, k, *neigh_size, n_neig_tot, *storage=NULL, storage_idx;

    if (!graph) 
    {
	*n_graph = NULL;
	return GMRFLib_SUCCESS;
    }

    /* 
       setup new graph
     */
    GMRFLib_make_empty_graph(n_graph);
    neigh_size       = Calloc(graph->n, int);   MEMCHK(neigh_size);
    (*n_graph)->n    = graph->n;
    (*n_graph)->nbs  = Calloc(graph->n, int *); MEMCHK((*n_graph)->nbs);
    (*n_graph)->nnbs = Calloc(graph->n, int);   MEMCHK((*n_graph)->nnbs);

    for(i=0;i<(*n_graph)->n;i++)
    {
	(*n_graph)->nnbs[i] = graph->nnbs[i];
	neigh_size[i]       = MAX(16, 2*(*n_graph)->nnbs[i]);
	(*n_graph)->nbs[i]  = Calloc(neigh_size[i], int);
	if ((*n_graph)->nnbs[i])
	    memcpy((*n_graph)->nbs[i], graph->nbs[i], graph->nnbs[i]*sizeof(int));
    }

    /* 
       if i ~ j, then add j ~ i
     */
    for(i=0;i<(*n_graph)->n;i++)
    {
	for(j=0;j<(*n_graph)->nnbs[i];j++)
	{
	    ii = (*n_graph)->nbs[i][j];
	    if ((*n_graph)->nnbs[ii] == neigh_size[ii])
	    {
		neigh_size[ii] *= 2;
		(*n_graph)->nbs[ii] = realloc((*n_graph)->nbs[ii], neigh_size[ii]*sizeof(int));
	    }
	    (*n_graph)->nbs[ii][ (*n_graph)->nnbs[ii]++ ] = i;
	}
    }

    /* 
       some may appear more than once, count number of unique neighbors
     */
    for(i=0, n_neig_tot = 0;i<(*n_graph)->n;i++)
	if ((*n_graph)->nnbs[i])
	{
	    qsort((*n_graph)->nbs[i], (size_t)(*n_graph)->nnbs[i], sizeof(int), GMRFLib_intcmp);
	    for(j=1, k=0; j<(*n_graph)->nnbs[i];j++)
		if ((*n_graph)->nbs[i][j] != (*n_graph)->nbs[i][k]) 
		    (*n_graph)->nbs[i][++k] = (*n_graph)->nbs[i][j];
	    (*n_graph)->nnbs[i] = k+1;
	    n_neig_tot += (*n_graph)->nnbs[i];
	}

    if (n_neig_tot) { storage = Calloc(n_neig_tot, int); MEMCHK(storage); }
    for(i=0, storage_idx=0;i<(*n_graph)->n;i++)
    {
	if ((*n_graph)->nnbs[i])
	{
	    memcpy(&storage[storage_idx], (*n_graph)->nbs[i], sizeof(int)*(*n_graph)->nnbs[i]);
	    FREE((*n_graph)->nbs[i]);
	    (*n_graph)->nbs[i] = &storage[storage_idx];
	    storage_idx += (*n_graph)->nnbs[i];
	}
	else
	{
	    FREE((*n_graph)->nbs[i]);
	}
    }
    GMRFLib_prepare_graph(*n_graph);

    FREE(neigh_size);

    return GMRFLib_SUCCESS;
}
int GMRFLib_offset_graph(GMRFLib_graph_tp **new_graph, int n_new,  int offset, GMRFLib_graph_tp *graph)
{
    /* 
       insert graph into a larger graph with n_new nodes such that node `i' in graph corresponds to
       node `offset +i' in the larger graph. n_new must be larger than graph->n, and `offset' such
       that `0<= offset', and `offset + graph->n < new_graph->n'
    */
    
    int i, ii, j, n_neig, *hold, hold_idx;
    GMRFLib_graph_tp *g;
    

    GMRFLib_ASSERT(graph,  GMRFLib_EPARAMETER);
    GMRFLib_ASSERT(n_new >= graph->n,  GMRFLib_EPARAMETER);
    GMRFLib_ASSERT(offset >= 0 && offset + graph->n -1 < n_new,  GMRFLib_EPARAMETER);

    for(i=0, n_neig = 0; i < graph->n; i++) n_neig += graph->nnbs[i];

    GMRFLib_make_empty_graph(&g);
    g->n    = n_new;
    g->nnbs = Calloc(n_new,  int);   MEMCHK(g->nnbs);
    g->nbs  = Calloc(n_new,  int *); MEMCHK(g->nbs);

    if (n_neig)
    {
	hold = Calloc(n_neig, int);   MEMCHK(hold);
    }
    else
	hold = NULL;
  
    for(i=0, hold_idx=0;i<graph->n;i++)
    {
	ii          = i + offset;
	g->nnbs[ii] = graph->nnbs[i];
	g->nbs[ii]  = &hold[hold_idx];
	hold_idx   += graph->nnbs[i];

	for(j=0;j<graph->nnbs[i];j++) g->nbs[ii][j] = graph->nbs[i][j] + offset;
    }
    GMRFLib_prepare_graph(g);

    *new_graph = g;

    return GMRFLib_SUCCESS;
}
double GMRFLib_offset_Qfunc(int node, int nnode,  char *arg)
{
    GMRFLib_offset_arg_tp *a;

    a = (GMRFLib_offset_arg_tp *)arg;

    if (MIN(node, nnode) <  a->offset || 
	MAX(node, nnode) >= a->offset + a->n) return 0.0;

    return (*a->Qfunc)(node - a->offset,  nnode - a->offset, a->Qfunc_arg);
}
int GMRFLib_offset(GMRFLib_offset_tp **off, int n_new,  int offset,
		   GMRFLib_graph_tp *graph, GMRFLib_Qfunc_tp *Qfunc, char *Qfunc_arg) 
{
    /* 
       make an object containg the new graph, Qfunc and Qfunc_arg when an graph if inserted into a
       graph.
    */
    
    GMRFLib_offset_tp *val;
    GMRFLib_offset_arg_tp *offset_arg;
    
    GMRFLib_ASSERT(graph,  GMRFLib_EPARAMETER);
    GMRFLib_ASSERT(n_new >= graph->n,  GMRFLib_EPARAMETER);
    GMRFLib_ASSERT(offset >= 0 && offset + graph->n -1 < n_new,  GMRFLib_EPARAMETER);

    offset_arg            = Calloc(1, GMRFLib_offset_arg_tp); MEMCHK(offset_arg);
    offset_arg->Qfunc     = Qfunc;
    offset_arg->Qfunc_arg = Qfunc_arg;
    offset_arg->offset    = offset;
    offset_arg->n         = graph->n;

    val = Calloc(1, GMRFLib_offset_tp); MEMCHK(val);
    GMRFLib_offset_graph(&(val->graph), n_new, offset, graph);
    
    val->Qfunc     = GMRFLib_offset_Qfunc;
    val->Qfunc_arg = (char *) offset_arg;

    *off = val;

    return GMRFLib_SUCCESS;
}

/*
  Example for manual
 */

/*! \page ex_graph Graph specification and handling. 
  This page includes the examples
  \ref ex_graph_sec1 and \ref ex_graph_sec2 

  \section ex_graph_sec1 Reading and printing graphs

  \par Description:

  This program performs the following tasks:
  
  - Reads a graph from the file \c graph1.dat, given below,
  of the form required by the function \c GMRFLib_read_graph(). 
  The specified graph has 10 nodes, with each of them having between 1 and 5 neighbors. 
  \verbinclude example_doxygen_graph_2.txt
  
  - Prints the graph specification to standard output.

  - Computes a subgraph, by specifying nodes to be removed from
  the graph. More specifically, the first half of the nodes,
  that is node 0, ..., 4, are specified to be removed.  To
  compute the subgraph, the function \c GMRFLib_compute_subgraph() is called.
  
  - Frees allocated memory.

  \par Program code:

  \verbinclude example_doxygen_graph_1.txt

  \par Output:

  \verbinclude example_doxygen_graph_3.txt

  \section ex_graph_sec2 Creating lattice graphs, linear graphs and folded graphs 

  \par Description:

  This program describes how to create a lattice graph and a linear
  graph, and how to expand the neighbourhood of a graph. The lattice
  graph is defined on a 6 \c x 7 lattice, such that \c nr = 6 and
  \c nc = 7. The neighbourhood is defined to be a 3 \c x
  3-neighbourhood, such that \c mr = \c mc = 1 in (ref 1).  The
  node numbers corresponding to the lattice indices are extracted, and
  the reverse operation is also illustrated.  The linear graph is
  defined by a cyclic AR(1)-model with 10 elements. A third graph is
  generated by expanding the neighbourhood of the linear graph twice.

  \par Program code:

  \verbinclude example_doxygen_graph_4.txt

  \par Output:

  \verbinclude example_doxygen_graph_5.txt

*/  

    
    
    
    
