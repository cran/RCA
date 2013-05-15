// This file is part of RCA-R.
// Written by Amir Goldberg with the help of Ryan Turner.
// This library includes a set of altered igraph functions based on code
// provided by Gábor Csárdi (see http://igraph.sourceforge.net/ for details)

// RCA-R is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free
// Software Foundation; either version 2 of the License, or (at your
// option) any later version.

// RCA-R is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.

// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA


/*#include <igraph.h>
#include <igraph_community.h>
#include <igraph_constructors.h>
#include <igraph_memory.h>
#include <igraph_random.h>
#include <igraph_arpack.h>
#include <igraph_adjlist.h>
#include <igraph_interface.h>
#include <igraph_components.h>
#include <igraph_dqueue.h>
#include <igraph_progress.h>
#include <igraph_stack.h>
#include <igraph_spmatrix.h>
#include <igraph_statusbar.h>
#include <igraph_conversion.h>
#include <igraph_centrality.h>*/

// The header files of igraph-C package are provided for those who can not
// install the igraph-C package from source code. 

#include "igraph.h"
#include "igraph_community.h"
#include "igraph_constructors.h"
#include "igraph_memory.h"
#include "igraph_random.h"
#include "igraph_arpack.h"
#include "igraph_adjlist.h"
#include "igraph_interface.h"
#include "igraph_components.h"
#include "igraph_dqueue.h"
#include "igraph_progress.h"
#include "igraph_stack.h"
#include "igraph_spmatrix.h"
#include "igraph_statusbar.h"
#include "igraph_conversion.h"
#include "igraph_centrality.h"


#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/param.h>
#include <string.h>
#include <math.h>
#include "RCA.h"

void release_array(double **array, int size) {
    int i;
    
    for(i=0; i<size; i++) {
        free(array[i]);
    }
    free(array);
}

double **get_min_max(double **matrix, int numRows, int numCol) {
    double **m_array;
	double temp_min, temp_max;
	int i, j;
    
    m_array = (double **)malloc(2 * sizeof(double *));
    m_array[0] = (double *)malloc(numCol * sizeof(double));
    m_array[1] = (double *)malloc(numCol * sizeof(double));
    
    for (i=0; i<numCol; i++) {
		temp_min = matrix[0][i];
        temp_max = matrix[0][i];
		for (j=1; j<numRows; j++) {
			if (matrix[j][i] < temp_min) {
				temp_min = matrix[j][i];
			}
            if (matrix[j][i] > temp_max) {
				temp_max = matrix[j][i];
			}
		}
		m_array[0][i] = temp_min;
        m_array[1][i] = temp_max;
    }
    return m_array;
}

// Standardizes a matrix on a 0 to 1 range
double **standardize(double **matrix, int numObs, int numVars, double **min_max) {
	int i, j;
    
    double **s_matrix = (double **)malloc(numObs*sizeof(double *));
	for(i=0; i<numObs; i++) {
        s_matrix[i] = (double *)malloc(numVars*sizeof(double));
		for(j=0; j<numVars; j++) {
			s_matrix[i][j] = (matrix[i][j] - min_max[0][j])/(min_max[1][j] - min_max[0][j]);
		}
	}
    
    return s_matrix;
}

double **relatioanlity_matrix(double **input_matrix, int numObs, int numVars) {
	int i, j, k, l, m;
	double sigma;
	double delta;
    double lambda;
	double **output_matrix;
    double **diff_matrix;
	output_matrix = (double **)malloc(numObs*sizeof(double *));
	/*if (output_matrix == NULL) {
		fprintf(stderr, "Malloc");
		exit(1);
	}*/
    
	for (m=0; m<numObs; m++) {
		output_matrix[m] = (double *)malloc(numObs * sizeof(double));
		/*if(output_matrix[m] == NULL) {
			fprintf(stderr, "Malloc");
			exit(1);
		}*/
	}
    
    diff_matrix = (double **)malloc(numObs*sizeof(double *));
    /*if (diff_matrix == NULL) {
		fprintf(stderr, "Malloc");
		exit(1);
	}*/
    
    int diff_num = (numVars*(numVars-1))/2;
    
    for (m=0; m<numObs; m++) {
		diff_matrix[m] = (double *)malloc(diff_num * sizeof(double));
		/*if(diff_matrix[m] == NULL) {
			fprintf(stderr, "Malloc");
			exit(1);
		}*/
	}
    
    // caculate deltas within each observation 
    for(i=0; i<numObs; i++) {
        l = 0;
        for (j=0; j<numVars-1; j++) {
            for (k=j+1; k<numVars; k++) {
                diff_matrix[i][l] = input_matrix[i][j] - input_matrix[i][k];
                l++;
            }
        }
    }
    
    for(i=0; i<numObs-1; i++) {
		for(j=i+1; j<numObs; j++) {
            sigma = 0;
            if (i!=j) {
                for (k=0; k<diff_num; k++) {
                    delta = 1 - fabs(fabs(diff_matrix[i][k]) - fabs(diff_matrix[j][k]));
                    lambda = diff_matrix[i][k]*diff_matrix[j][k]<0 ? -1 : 1;
                    sigma += delta*lambda;
                }
            }
            output_matrix[i][j] = sigma/diff_num;
            output_matrix[j][i] = sigma/diff_num;
        }
    }
    
    release_array(diff_matrix,numObs);
    
	return output_matrix;
}

void bs_sample(double *sample,double **R,int dim) {
    int i, j, k, mult;
    int *indx = (int *)malloc(dim * sizeof(int));
    for (i=0; i<dim; i++) {
        mult = rand() % dim;
        indx[i] = (int)((mult+(rand()%dim)) % dim);
    }
        
    k = 0;
    for (i=0; i<(dim-1); i++) {
        for (j=i+1; j<dim; j++) {
            if(indx[i]==indx[j]) {
                sample[k] = 1;
            }
            else {
                sample[k] = R[indx[i]][indx[j]];
            }
            k++;
        }
    }
    
    free(indx);
}

double stats_mean (double * sample, int n){
    double total = 0;
    for(int i=0;i<n;i++){
        total += sample[i];
    }
    return total/n;
}

double stats_sd (double * sample, int n){
    double mean = stats_mean(sample, n);
    double varSum =0;
    for(int i=0;i<n;i++){
       varSum += (sample[i] - mean)*(sample[i] - mean);
    }
    varSum /=(n-1);
    return sqrt(varSum);
}

void r_boorstrap(double **R, int dim, double *mean, double *stdev, int bootstrap) {
    int i, j, k;
    int n = (dim*(dim-1))/2;
    double *means = NULL;
    double *stds = NULL;
    double *sample = (double *)malloc(n * sizeof(double));
    if(bootstrap) {
        means = (double *)malloc((bootstrap)*sizeof(double));
        stds = (double *)malloc(bootstrap*sizeof(double));
    }
    
    char text[256];

    srand(time(NULL));
        
    for (i=0; i<bootstrap; i++) {
        bs_sample(sample,R,dim);
        means[i] = stats_mean(sample,n);
        stds[i] = stats_sd(sample,n);
    }

    if (bootstrap) {
        *mean = (double)stats_mean(means,bootstrap);
        *stdev = (double)stats_mean(stds,bootstrap);
    }
    else {
        k = 0;
        for (i=0; i<(dim-1); i++) {
            for (j=i+1; j<dim; j++) {
                sample[k] = R[i][j];
                k++;
            }
        }
        *mean = (double)stats_mean(sample,n);
        *stdev = (double)stats_sd(sample,n);
    }
    free(sample);
    if(bootstrap) {
        free(means);
        free(stds);
    }
    
}

void get_significance(double **R, int dim, double mean, double stdev, double z_score) {
    
	int i, j;
    //double z_score =fabs(gsl_cdf_ugaussian_Pinv(p_value/2));
    
    double rangeVal = stdev*z_score;
    
	for(i=0; i<dim-1; i++) {
		for(j=i+1; j<dim; j++) {
			if (i != j) {
                if ((R[i][j]-mean)>(-rangeVal) && (R[i][j]-mean)<rangeVal) {
                    R[i][j] = 0;
                }
                else {
                    R[i][j] = fabs(R[i][j]-mean);
                }
                R[j][i] = R[i][j];
            }
		}
	}
}

void get_weights(igraph_t *g, igraph_vector_t *weights) {
    long int i, n;
    
    n = igraph_ecount(g);
    igraph_vector_init(weights, n);
    
    for (i=0; i<n; i++) {
        VECTOR(*weights)[i]=EAN(g, "weight", i);
    }
    
}

void generate_graph(double **R,int dim,igraph_t *g, igraph_matrix_t *weighted_adjacency) {
    int i,j;
    igraph_matrix_init(weighted_adjacency, dim, dim);
    
    for (i=0; i<dim; i++) {
        for (j=0; j<dim; j++) {
            igraph_matrix_set(weighted_adjacency, i, j,  R[i][j]);
        }
    }
    igraph_i_set_attribute_table(&igraph_cattribute_table);
    igraph_weighted_adjacency(g, weighted_adjacency, IGRAPH_ADJ_UPPER, 0, /*loops=*/ 0);
    
}

typedef struct RCA_igraph_i_community_leading_eigenvector_data_t {
    igraph_vector_t *idx;
    igraph_vector_t *idx2;
    igraph_adjlist_t *adjlist;
    igraph_inclist_t *inclist;
    igraph_vector_t *tmp;
    long int no_of_edges;
    igraph_vector_t *mymembership;
    long int comm;
    const igraph_vector_t *weights;
    const igraph_t *graph;
    igraph_vector_t *strength;
    igraph_real_t sumweights;
} RCA_igraph_i_community_leading_eigenvector_data_t;

void RCA_igraph_i_error_handler_none(const char *reason, const char *file,
                                 int line, int igraph_errno) {
    IGRAPH_UNUSED(reason);
    IGRAPH_UNUSED(file);
    IGRAPH_UNUSED(line);
    IGRAPH_UNUSED(igraph_errno);  
    /* do nothing */
}


int RCA_igraph_i_community_leading_eigenvector(igraph_real_t *to,
                                           const igraph_real_t *from,
                                           int n, void *extra) {
    
    RCA_igraph_i_community_leading_eigenvector_data_t *data=(RCA_igraph_i_community_leading_eigenvector_data_t*)extra;
    long int j, k, nlen, size=n;
    igraph_vector_t *idx=data->idx;
    igraph_vector_t *idx2=data->idx2;
    igraph_vector_t *tmp=data->tmp;
    igraph_adjlist_t *adjlist=data->adjlist;
    igraph_real_t ktx, ktx2;
    long int no_of_edges=data->no_of_edges;
    igraph_vector_t *mymembership=data->mymembership;
    long int comm=data->comm;
    
    /* Ax */
    for (j=0; j<size; j++) {
        long int oldid=VECTOR(*idx)[j];
        igraph_vector_t *neis=igraph_adjlist_get(adjlist, oldid);
        nlen=igraph_vector_size(neis);
        to[j]=0.0;
        VECTOR(*tmp)[j]=0.0;
        for (k=0; k<nlen; k++) {
            long int nei=VECTOR(*neis)[k];
            long int neimemb=VECTOR(*mymembership)[nei];
            if (neimemb==comm) {
                to[j] += from[ (long int) VECTOR(*idx2)[nei] ];
                VECTOR(*tmp)[j] += 1;
            }
        }
    }
    
    /* Now calculate k^Tx/2m */
    ktx=0.0; ktx2=0.0;
    for (j=0; j<size; j++) {
        long int oldid=VECTOR(*idx)[j];
        igraph_vector_t *neis=igraph_adjlist_get(adjlist, oldid);
        long int degree=igraph_vector_size(neis);
        ktx += from[j] * degree;
        ktx2 += degree;
    }
    ktx = ktx / no_of_edges/2.0;
    ktx2 = ktx2 / no_of_edges/2.0;
    
    /* Now calculate Bx */
    for (j=0; j<size; j++) {
        long int oldid=VECTOR(*idx)[j];
        igraph_vector_t *neis=igraph_adjlist_get(adjlist, oldid);
        igraph_real_t degree=igraph_vector_size(neis);
        to[j] = to[j] - ktx*degree;
        VECTOR(*tmp)[j] = VECTOR(*tmp)[j] - ktx2*degree;
    }
    
    /* -d_ij summa l in G B_il */
    for (j=0; j<size; j++) {
        to[j] -= VECTOR(*tmp)[j] * from[j];
    }
    
    return 0;
}

int RCA_igraph_i_community_leading_eigenvector2(igraph_real_t *to,
                                            const igraph_real_t *from,
                                            int n, void *extra) {
    
    RCA_igraph_i_community_leading_eigenvector_data_t *data=(RCA_igraph_i_community_leading_eigenvector_data_t*)extra;
    long int j, k, nlen, size=n;
    igraph_vector_t *idx=data->idx;
    igraph_vector_t *idx2=data->idx2;
    igraph_vector_t *tmp=data->tmp;
    igraph_adjlist_t *adjlist=data->adjlist;
    igraph_real_t ktx, ktx2;
    long int no_of_edges=data->no_of_edges;
    igraph_vector_t *mymembership=data->mymembership;
    long int comm=data->comm;
    
    /* Ax */
    for (j=0; j<size; j++) {
        long int oldid=VECTOR(*idx)[j];
        igraph_vector_t *neis=igraph_adjlist_get(adjlist, oldid);
        nlen=igraph_vector_size(neis);
        to[j]=0.0;
        VECTOR(*tmp)[j]=0.0;
        for (k=0; k<nlen; k++) {
            long int nei=VECTOR(*neis)[k];
            long int neimemb=VECTOR(*mymembership)[nei];
            if (neimemb==comm) {
                long int fi=VECTOR(*idx2)[nei];
                if (fi < size) {
                    to[j] += from[fi];
                }
                VECTOR(*tmp)[j] += 1;
            }
        }
    }
    
    /* Now calculate k^Tx/2m */
    ktx=0.0; ktx2=0.0;
    for (j=0; j<size+1; j++) {
        long int oldid=VECTOR(*idx)[j];
        igraph_vector_t *neis=igraph_adjlist_get(adjlist, oldid);
        long int degree=igraph_vector_size(neis);
        if (j<size) {
            ktx += from[j] * degree;
        }
        ktx2 += degree;
    }
    ktx = ktx / no_of_edges/2.0;
    ktx2 = ktx2 / no_of_edges/2.0;
    
    /* Now calculate Bx */
    for (j=0; j<size; j++) {
        long int oldid=VECTOR(*idx)[j];
        igraph_vector_t *neis=igraph_adjlist_get(adjlist, oldid);
        igraph_real_t degree=igraph_vector_size(neis);
        to[j] = to[j] - ktx*degree;
        VECTOR(*tmp)[j] = VECTOR(*tmp)[j] - ktx2*degree;
    }
    
    /* -d_ij summa l in G B_il */
    for (j=0; j<size; j++) {
        to[j] -= VECTOR(*tmp)[j] * from[j];
    }
    
    return 0;
}


int RCA_igraph_i_community_leading_eigenvector_weighted(igraph_real_t *to,
                                                    const igraph_real_t *from,
                                                    int n, void *extra) {
    
    RCA_igraph_i_community_leading_eigenvector_data_t *data=(RCA_igraph_i_community_leading_eigenvector_data_t*)extra;
    long int j, k, nlen, size=n;
    igraph_vector_t *idx=data->idx;
    igraph_vector_t *idx2=data->idx2;
    igraph_vector_t *tmp=data->tmp;
    igraph_inclist_t *inclist=data->inclist;
    igraph_real_t ktx, ktx2;
    igraph_vector_t *mymembership=data->mymembership;
    long int comm=data->comm;
    const igraph_vector_t *weights=data->weights;
    const igraph_t *graph=data->graph;
    igraph_vector_t *strength=data->strength;
    igraph_real_t sw=data->sumweights;
    
    /* Ax */
    for (j=0; j<size; j++) {
        long int oldid=VECTOR(*idx)[j];
        igraph_vector_t *inc=igraph_inclist_get(inclist, oldid);
        nlen=igraph_vector_size(inc);
        to[j]=0.0;
        VECTOR(*tmp)[j]=0.0;
        for (k=0; k<nlen; k++) {
            long int edge=VECTOR(*inc)[k];
            igraph_real_t w=VECTOR(*weights)[edge];
            long int nei=IGRAPH_OTHER(graph, edge, oldid);
            long int neimemb=VECTOR(*mymembership)[nei];
            if (neimemb==comm) {
                to[j] += from[ (long int) VECTOR(*idx2)[nei] ] * w;
                VECTOR(*tmp)[j] += w;
            }
        }
    }
    
    /* k^Tx/2m */
    ktx=0.0; ktx2=0.0;
    for (j=0; j<size; j++) {
        long int oldid=VECTOR(*idx)[j];
        igraph_real_t str=VECTOR(*strength)[oldid];
        ktx += from[j] * str;
        ktx2 += str;
    }
    ktx = ktx / sw / 2.0;
    ktx2 = ktx2 / sw / 2.0;
    
    /* Bx */
    for (j=0; j<size; j++) {
        long int oldid=VECTOR(*idx)[j];
        igraph_real_t str=VECTOR(*strength)[oldid];
        to[j] = to[j] - ktx * str;
        VECTOR(*tmp)[j] = VECTOR(*tmp)[j] - ktx2*str;
    }
    
    /* -d_ij summa l in G B_il */
    for (j=0; j<size; j++) {
        to[j] -= VECTOR(*tmp)[j] * from[j];
    }
    
    return 0;
}

int RCA_igraph_i_community_leading_eigenvector2_weighted(igraph_real_t *to,
                                                     const igraph_real_t *from,
                                                     int n, void *extra) {
    
    RCA_igraph_i_community_leading_eigenvector_data_t *data=(RCA_igraph_i_community_leading_eigenvector_data_t*)extra;
    long int j, k, nlen, size=n;
    igraph_vector_t *idx=data->idx;
    igraph_vector_t *idx2=data->idx2;
    igraph_vector_t *tmp=data->tmp;
    igraph_inclist_t *inclist=data->inclist;
    igraph_real_t ktx, ktx2;
    igraph_vector_t *mymembership=data->mymembership;
    long int comm=data->comm;
    const igraph_vector_t *weights=data->weights;
    const igraph_t *graph=data->graph;
    igraph_vector_t *strength=data->strength;
    igraph_real_t sw=data->sumweights;
    
    /* Ax */
    for (j=0; j<size; j++) {
        long int oldid=VECTOR(*idx)[j];
        igraph_vector_t *inc=igraph_inclist_get(inclist, oldid);
        nlen=igraph_vector_size(inc);
        to[j]=0.0;
        VECTOR(*tmp)[j]=0.0;
        for (k=0; k<nlen; k++) {
            long int edge=VECTOR(*inc)[k];
            igraph_real_t w=VECTOR(*weights)[edge];
            long int nei=IGRAPH_OTHER(graph, edge, oldid);
            long int neimemb=VECTOR(*mymembership)[nei];
            if (neimemb==comm) {
                long int fi=VECTOR(*idx2)[nei];
                if (fi < size) {
                    to[j] += from[fi] * w;
                }
                VECTOR(*tmp)[j] += w;
            }
        }
    }
    
    /* k^Tx/2m */
    ktx=0.0; ktx2=0.0;
    for (j=0; j<size+1; j++) {
        long int oldid=VECTOR(*idx)[j];
        igraph_real_t str=VECTOR(*strength)[oldid];
        if (j<size) { 
            ktx += from[j] * str;
        } 
        ktx2 += str;
    }
    ktx = ktx / sw / 2.0;
    ktx2 = ktx2 / sw / 2.0;
    
    /* Bx */
    for (j=0; j<size; j++) {
        long int oldid=VECTOR(*idx)[j];
        igraph_real_t str=VECTOR(*strength)[oldid];
        to[j] = to[j] - ktx * str;
        VECTOR(*tmp)[j] = VECTOR(*tmp)[j] - ktx2 * str;
    }
    
    /* -d_ij summa l in G B_il */
    for (j=0; j<size; j++) {
        to[j] -= VECTOR(*tmp)[j] * from[j];
    }
    
    return 0;
}

void rca_igraph_i_levc_free(igraph_vector_ptr_t *ptr) { 
    long int i, n=igraph_vector_ptr_size(ptr);
    for (i=0; i<n; i++) {
        igraph_vector_t *v=(igraph_vector_t*)VECTOR(*ptr)[i];
        if (v) {
            igraph_vector_destroy(v);
            igraph_free(v);
        }
    }
}


int RCA_igraph_community_leading_eigenvector(const igraph_t *graph,
                                             const igraph_vector_t *weights,
                                             igraph_matrix_t *merges,
                                             igraph_vector_t *membership,
                                             igraph_integer_t steps,
                                             igraph_arpack_options_t *options, 
                                             igraph_real_t *modularity,
                                             igraph_bool_t start,
                                             igraph_vector_t *eigenvalues,
                                             igraph_vector_ptr_t *eigenvectors,
                                             igraph_vector_t *history,
                                             igraph_community_leading_eigenvector_callback_t *callback,
                                             void *callback_extra) {
    
    long int no_of_nodes=igraph_vcount(graph);
    long int no_of_edges=igraph_ecount(graph);
    igraph_dqueue_t tosplit;
    igraph_vector_t idx, idx2, mymerges;
    igraph_vector_t strength, tmp;
    long int staken=0;
    igraph_adjlist_t adjlist;
    igraph_inclist_t inclist;
    long int i, j, k, l;
    long int communities;
    igraph_vector_t vmembership, *mymembership=membership;
    RCA_igraph_i_community_leading_eigenvector_data_t extra;
    igraph_arpack_storage_t storage;
    igraph_real_t mod=0;
    igraph_arpack_function_t *arpcb1 = 
    weights ? RCA_igraph_i_community_leading_eigenvector_weighted : 
    RCA_igraph_i_community_leading_eigenvector;
    igraph_arpack_function_t *arpcb2 = 
    weights ? RCA_igraph_i_community_leading_eigenvector2_weighted : 
    RCA_igraph_i_community_leading_eigenvector2;
    igraph_real_t sumweights=0.0;
    
    if (weights && no_of_edges != igraph_vector_size(weights)) {
        IGRAPH_ERROR("Invalid weight vector length", IGRAPH_EINVAL);
    }
    
    if (start && !membership) { 
        IGRAPH_ERROR("Cannot start from given configuration if memberships "
                     "missing", IGRAPH_EINVAL);
    }
    
    if (start && membership && 
        igraph_vector_size(membership) != no_of_nodes) {
        IGRAPH_ERROR("Wrong length for vector of predefined memberships", 
                     IGRAPH_EINVAL);
    }
    
    if (start && membership && igraph_vector_max(membership) >= no_of_nodes) {
        IGRAPH_WARNING("Too many communities in membership start vector");
    }
    
    if (igraph_is_directed(graph)) {
        IGRAPH_WARNING("This method was developed for undirected graphs");
    }
    
    if (steps < 0 || steps > no_of_nodes-1) {
        steps=no_of_nodes-1;
    }
    
    if (!membership) {
        mymembership=&vmembership;
        IGRAPH_VECTOR_INIT_FINALLY(mymembership, 0);
    }
    
    IGRAPH_VECTOR_INIT_FINALLY(&mymerges, 0);
    IGRAPH_CHECK(igraph_vector_reserve(&mymerges, steps*2));
    IGRAPH_VECTOR_INIT_FINALLY(&idx, 0);
    if (eigenvalues)  { igraph_vector_clear(eigenvalues);      }
    if (eigenvectors) { 
        igraph_vector_ptr_clear(eigenvectors); 
        IGRAPH_FINALLY(rca_igraph_i_levc_free, eigenvectors);
    }
    
    IGRAPH_STATUS("Starting leading eigenvector method.\n", 0);
    
    if (!start) {
        /* Calculate the weakly connected components in the graph and use them as
         * an initial split */
        IGRAPH_CHECK(igraph_clusters(graph, mymembership, &idx, 0, IGRAPH_WEAK));
        communities = igraph_vector_size(&idx);
        IGRAPH_STATUSF(("Starting from %li component(s).\n", 0, communities));
        if (history) { 
            IGRAPH_CHECK(igraph_vector_push_back(history, 
                                                 IGRAPH_LEVC_HIST_START_FULL));
        }
    } else {
        /* Just create the idx vector for the given membership vector */
        communities=igraph_vector_max(mymembership)+1;
        IGRAPH_STATUSF(("Starting from given membership vector with %li "
                        "communities.\n", 0, communities));
        if (history) { 
            IGRAPH_CHECK(igraph_vector_push_back(history, 
                                                 IGRAPH_LEVC_HIST_START_GIVEN));
            IGRAPH_CHECK(igraph_vector_push_back(history, communities));
        }
        IGRAPH_CHECK(igraph_vector_resize(&idx, communities));
        igraph_vector_null(&idx);
        for (i=0; i<no_of_nodes; i++) {
            int t=VECTOR(*mymembership)[i];
            VECTOR(idx)[t] += 1;
        }
    }
    
    IGRAPH_DQUEUE_INIT_FINALLY(&tosplit, 100);
    for (i = 0; i < communities; i++) {
        if (VECTOR(idx)[i] > 2) {
            igraph_dqueue_push(&tosplit, i);
        }
    }
    for (i=1; i<communities; i++) {
        /* Record merge */
        IGRAPH_CHECK(igraph_vector_push_back(&mymerges, i-1));
        IGRAPH_CHECK(igraph_vector_push_back(&mymerges, i));
        if (eigenvalues) { 
            IGRAPH_CHECK(igraph_vector_push_back(eigenvalues, IGRAPH_NAN));
        }
        if (eigenvectors) { 
            igraph_vector_t *v=igraph_Calloc(1, igraph_vector_t);
            if (!v) { 
                IGRAPH_ERROR("Cannot do leading eigenvector community detection", 
                             IGRAPH_ENOMEM); 
            }
            IGRAPH_FINALLY(igraph_free, v);
            IGRAPH_VECTOR_INIT_FINALLY(v, 0);
            IGRAPH_CHECK(igraph_vector_ptr_push_back(eigenvectors, v));
            IGRAPH_FINALLY_CLEAN(2);
        }
        if (history) {
            IGRAPH_CHECK(igraph_vector_push_back(history, IGRAPH_LEVC_HIST_SPLIT));
            IGRAPH_CHECK(igraph_vector_push_back(history, i-1));
        }
    }
    staken = communities - 1;
    
    IGRAPH_VECTOR_INIT_FINALLY(&tmp, no_of_nodes);
    IGRAPH_CHECK(igraph_vector_resize(&idx, no_of_nodes));
    igraph_vector_null(&idx);
    IGRAPH_VECTOR_INIT_FINALLY(&idx2, no_of_nodes);
    if (!weights) { 
        IGRAPH_CHECK(igraph_adjlist_init(graph, &adjlist, IGRAPH_ALL));
        IGRAPH_FINALLY(igraph_adjlist_destroy, &adjlist);
    } else {
        IGRAPH_CHECK(igraph_inclist_init(graph, &inclist, IGRAPH_ALL));
        IGRAPH_FINALLY(igraph_inclist_destroy, &inclist);
        IGRAPH_VECTOR_INIT_FINALLY(&strength, no_of_nodes);
        IGRAPH_CHECK(igraph_strength(graph, &strength, igraph_vss_all(), 
                                     IGRAPH_ALL, IGRAPH_LOOPS, weights));
        sumweights=igraph_vector_sum(weights);
    }
    
    options->ncv = 0;   /* 0 means "automatic" in igraph_arpack_rssolve */
    options->start = 0;
    options->which[0]='L'; options->which[1]='A';
    
    /* Memory for ARPACK */
    /* We are allocating memory for 20 eigenvectors since options->ncv won't be
     * larger than 20 when using automatic mode in igraph_arpack_rssolve */
    IGRAPH_CHECK(igraph_arpack_storage_init(&storage, no_of_nodes, 20, 
                                            no_of_nodes, 1));
    IGRAPH_FINALLY(igraph_arpack_storage_destroy, &storage);
    extra.idx=&idx;
    extra.idx2=&idx2;
    extra.tmp=&tmp;
    extra.adjlist=&adjlist;
    extra.inclist=&inclist;
    extra.weights=weights;
    extra.sumweights=sumweights;
    extra.graph=graph;
    extra.strength=&strength;
    extra.no_of_edges=no_of_edges;
    extra.mymembership=mymembership;
    
    while (!igraph_dqueue_empty(&tosplit) && staken < steps) {
        long int comm=igraph_dqueue_pop_back(&tosplit); /* depth first search */
        long int size=0;
        igraph_real_t tmpev;
        
        IGRAPH_STATUSF(("Trying to split community %li... ", 0, comm));
        //IGRAPH_ALLOW_INTERRUPTION();
        
        for (i=0; i<no_of_nodes; i++) {
            if (VECTOR(*mymembership)[i]==comm) {
                VECTOR(idx)[size]=i;
                VECTOR(idx2)[i]=size++;
            }
        }
        
        staken++;
        if (size<=2) {
            continue;
        }
        
        /* We solve two eigenproblems, one for the original modularity
         matrix, and one for the modularity matrix after deleting the
         last row and last column from it. This is a trick to find
         multiple leading eigenvalues, because ARPACK is sometimes
         unstable when the first two eigenvalues are requested, but it
         does much better for the single principal eigenvalue. */
        
        /* We start with the smaller eigenproblem. */
        
        options->n=size-1;
        options->info=0;
        options->nev=1;
        options->ldv=0;
        options->ncv = 0;   /* 0 means "automatic" in igraph_arpack_rssolve */
        options->nconv = 0;
        options->lworkl = 0;		/* we surely have enough space */
        extra.comm=comm;
        
        /* We try calling the solver twice, once from a random starting
         point, onece from a fixed one. This is because for some hard
         cases it tends to fail. We need to suppress error handling for
         the first call. */
        {
            int i;
            igraph_error_handler_t *errh=
            igraph_set_error_handler(RCA_igraph_i_error_handler_none);
            igraph_arpack_rssolve(arpcb2, &extra, options, &storage,
                                  /*values=*/ 0, /*vectors=*/ 0);
            igraph_set_error_handler(errh);
            if (options->nconv < 1) {
                /* Call again, from a fixed starting point */
                options->start=1;
                options->info=0;
                options->ncv=0;
                options->lworkl = 0;	/* we surely have enough space */
                for (i=0; i < options->n ; i++) {
                    storage.resid[i] = 1;
                }
                IGRAPH_CHECK(igraph_arpack_rssolve(arpcb2, &extra, options, &storage,
                                                   /*values=*/ 0, /*vectors=*/ 0));
                options->start=0;	
            }
        }
        
        if (options->nconv < 1) {
            IGRAPH_ERROR("ARPACK did not converge", IGRAPH_ARPACK_FAILED);
        }
        
        tmpev=storage.d[0];
        
        /* Now we do the original eigenproblem, again, twice if needed */
        
        options->n=size;
        options->info=0;
        options->nev=1;
        options->ldv=0;
        options->nconv=0;
        options->lworkl = 0;	/* we surely have enough space */
        options->ncv = 0;   /* 0 means "automatic" in igraph_arpack_rssolve */
        
        {
            int i;
            igraph_error_handler_t *errh=
            igraph_set_error_handler(RCA_igraph_i_error_handler_none);
            igraph_arpack_rssolve(arpcb1, &extra, options, &storage,
                                  /*values=*/ 0, /*vectors=*/ 0);
            igraph_set_error_handler(errh);
            if (options->nconv < 1) {
                /* Call again from a fixed starting point */
                options->start=1;
                options->info=0;
                options->ncv=0;
                options->lworkl = 0;	/* we surely have enough space */
                for (i=0; i < options->n; i++) { storage.resid[i] = 1; }
                IGRAPH_CHECK(igraph_arpack_rssolve(arpcb1, &extra, options, &storage, 
                                                   /*values=*/ 0, /*vectors=*/ 0));
                options->start=0;
            }
        }
        
        if (options->nconv < 1) {
            IGRAPH_ERROR("ARPACK did not converge", IGRAPH_ARPACK_FAILED);
        }
        
        /* Ok, we have the leading eigenvector of the modularity matrix*/
        
        /* ---------------------------------------------------------------*/
        /* To avoid numeric errors */
        if (fabs(storage.d[0]) < 1e-8) { 
            storage.d[0] = 0;
        }
        
        /* We replace very small (in absolute value) elements of the 
         leading eigenvector with zero, to get the same result, 
         consistently.*/
        for (i=0; i<size; i++) {
            if (fabs(storage.v[i]) < 1e-8) {
                storage.v[i]=0;
            }
        }
        
        /* Just to have the always the same result, we multiply by -1
         if the first (nonzero) element is not positive. */
        for (i=0; i<size; i++) {
            if (storage.v[i] != 0) { break; }
        }
        if (i<size && storage.v[i]<0) {
            for (i=0; i<size; i++) {
                storage.v[i] = - storage.v[i];
            }
        }
        /* ---------------------------------------------------------------*/
        
        if (callback) {
            igraph_vector_t vv;
            int ret;
            igraph_vector_view(&vv, storage.v, size);
            ret=callback(mymembership, comm, storage.d[0], &vv, 
                         arpcb1, &extra, callback_extra);
            if (ret) {
                break;
            }
        }
        
        if (eigenvalues) {
            IGRAPH_CHECK(igraph_vector_push_back(eigenvalues, storage.d[0]));
        }
        
        if (eigenvectors) { 
            igraph_vector_t *v=igraph_Calloc(1, igraph_vector_t);
            if (!v) { 
                IGRAPH_ERROR("Cannot do leading eigenvector community detection", 
                             IGRAPH_ENOMEM);
            }
            IGRAPH_FINALLY(igraph_free, v);
            IGRAPH_VECTOR_INIT_FINALLY(v, size);
            for (i=0; i<size; i++) {
                VECTOR(*v)[i]=storage.v[i];
            }
            IGRAPH_CHECK(igraph_vector_ptr_push_back(eigenvectors, v));
            IGRAPH_FINALLY_CLEAN(2);
        }
        
        if (storage.d[0] <= 0) {
            IGRAPH_STATUS("no split.\n", 0);
            if (history) { 
                IGRAPH_CHECK(igraph_vector_push_back(history, 
                                                     IGRAPH_LEVC_HIST_FAILED));
                IGRAPH_CHECK(igraph_vector_push_back(history, comm));
            }					     
            continue; 
        }
        
        /* Check for multiple leading eigenvalues */
        
        if (fabs(storage.d[0]-tmpev) < 1e-8) {
            IGRAPH_STATUS("multiple principal eigenvalue, no split.\n", 0);
            if (history) { 
                IGRAPH_CHECK(igraph_vector_push_back(history, 
                                                     IGRAPH_LEVC_HIST_FAILED));
                IGRAPH_CHECK(igraph_vector_push_back(history, comm));
            }					     
            continue;
        }
        
        /* Count the number of vertices in each community after the split */
        l=0;
        for (j=0; j<size; j++) {
            if (storage.v[j] < 0) {
                storage.v[j] = -1;
                l++;
            } else {
                storage.v[j] = 1;
            }
        }
        if (l==0 || l==size) {
            IGRAPH_STATUS("no split.\n", 0);
            if (history) { 
                IGRAPH_CHECK(igraph_vector_push_back(history, 
                                                     IGRAPH_LEVC_HIST_FAILED));
                IGRAPH_CHECK(igraph_vector_push_back(history, comm));
            }					     
            continue;
        }
        
        /* Check that Q increases with our choice of split */
        arpcb1(storage.v+size, storage.v, size, &extra);
        mod=0;
        for (i=0; i<size; i++) {
            mod += storage.v[size+i] * storage.v[i];
        }
        if (mod <= 1e-8) {
            IGRAPH_STATUS("no modularity increase, no split.\n", 0);
            if (history) { 
                IGRAPH_CHECK(igraph_vector_push_back(history, 
                                                     IGRAPH_LEVC_HIST_FAILED));
                IGRAPH_CHECK(igraph_vector_push_back(history, comm));
            }					     
            continue;      
        }
        
        communities++;
        IGRAPH_STATUS("split.\n", 0);
        
        /* Rewrite the mymembership vector */
        for (j=0; j<size; j++) {
            if (storage.v[j] < 0) {
                long int oldid=VECTOR(idx)[j];
                VECTOR(*mymembership)[oldid]=communities-1;
            }
        }
        
        /* Record merge */
        IGRAPH_CHECK(igraph_vector_push_back(&mymerges, comm));
        IGRAPH_CHECK(igraph_vector_push_back(&mymerges, communities-1));
        if (history) {
            IGRAPH_CHECK(igraph_vector_push_back(history, IGRAPH_LEVC_HIST_SPLIT));
            IGRAPH_CHECK(igraph_vector_push_back(history, comm));
        }
        
        /* Store the resulting communities in the queue if needed */
        if (l > 1) {
            IGRAPH_CHECK(igraph_dqueue_push(&tosplit, communities-1));
        }
        if (size-l > 1) {
            IGRAPH_CHECK(igraph_dqueue_push(&tosplit, comm));
        }
        
    }
    
    igraph_arpack_storage_destroy(&storage);
    IGRAPH_FINALLY_CLEAN(1);
    if (!weights) { 
        igraph_adjlist_destroy(&adjlist);
        IGRAPH_FINALLY_CLEAN(1);
    } else {
        igraph_inclist_destroy(&inclist);
        igraph_vector_destroy(&strength);
        IGRAPH_FINALLY_CLEAN(2);
    }
    igraph_dqueue_destroy(&tosplit);
    igraph_vector_destroy(&tmp);
    igraph_vector_destroy(&idx2);
    IGRAPH_FINALLY_CLEAN(3);
    
    IGRAPH_STATUS("Done.\n", 0);
    
    /* reform the mymerges vector */
    if (merges) {
        igraph_vector_null(&idx);
        l=igraph_vector_size(&mymerges);
        k=communities;
        j=0;
        IGRAPH_CHECK(igraph_matrix_resize(merges, l/2, 2));
        for (i=l; i>0; i-=2) {
            long int from=VECTOR(mymerges)[i-1];
            long int to=VECTOR(mymerges)[i-2];
            MATRIX(*merges, j, 0)=VECTOR(mymerges)[i-2];
            MATRIX(*merges, j, 1)=VECTOR(mymerges)[i-1];    
            if (VECTOR(idx)[from]!=0) {
                MATRIX(*merges, j, 1)=VECTOR(idx)[from]-1;
            }
            if (VECTOR(idx)[to]!=0) {
                MATRIX(*merges, j, 0)=VECTOR(idx)[to]-1;
            }
            VECTOR(idx)[to]=++k;
            j++;
        }      
    }
    
    if (eigenvectors) { IGRAPH_FINALLY_CLEAN(1); }
    igraph_vector_destroy(&idx);
    igraph_vector_destroy(&mymerges);
    IGRAPH_FINALLY_CLEAN(2);
    
    if (modularity) {
        IGRAPH_CHECK(igraph_modularity(graph, mymembership, modularity, 
                                       weights));
    }
    
    if (!membership) {
        igraph_vector_destroy(mymembership);
        IGRAPH_FINALLY_CLEAN(1);
    }
    
    return 0;
}


void partition(double **R, int numObs, int **membership, double *modularity, int ***merges, int merge_dim[]) {
    int i, j;
    igraph_t g;
    igraph_matrix_t weighted_adjacency;
    igraph_matrix_t i_merges;
    igraph_vector_t i_membership;
    igraph_vector_t eigenvalues;
    igraph_arpack_options_t options;
    igraph_vector_t weights;
    igraph_vector_ptr_t eigenvectors;
    
    igraph_matrix_init(&i_merges, 0, 0);
    igraph_vector_init(&i_membership, 0);
    igraph_vector_init(&eigenvalues, 0);
    igraph_vector_ptr_init(&eigenvectors, 0);
    igraph_arpack_options_init(&options);
    
    generate_graph(R,numObs,&g,&weighted_adjacency);
    get_weights(&g, &weights);
    RCA_igraph_community_leading_eigenvector(&g, &weights, &i_merges, &i_membership,
                                         numObs,
                                         &options, modularity, 
                                         /*start=*/ 0, &eigenvalues,
                                         &eigenvectors, /*history=*/ 0,
                                         /*callback=*/ 0, 
                                         /*callback_extra=*/ 0);
    
    if(membership) {
        *membership = (int *)malloc(igraph_vector_size(&i_membership)*sizeof(int));
        for (i=0; i<numObs; i++) {
            (*membership)[i]=(int)VECTOR(i_membership)[i];
        }
    }
    
    if(merges) {
        *merges = (int **)malloc(igraph_matrix_nrow(&i_merges)*sizeof(int *));
        for (i=0; i<igraph_matrix_nrow(&i_merges); i++) {
            (*merges)[i] = (int *)malloc(igraph_matrix_ncol(&i_merges)*sizeof(int));
            for (j=0; j<igraph_matrix_ncol(&i_merges); j++) {
                (*merges)[i][j]=(int)MATRIX(i_merges, i, j);
            }
        }
    }
    
    if(merge_dim) {
        merge_dim[0] = igraph_matrix_nrow(&i_merges);
        merge_dim[1] = igraph_matrix_ncol(&i_merges);
    }
    
    igraph_vector_ptr_destroy(&eigenvectors);
    igraph_vector_destroy(&i_membership);
    igraph_vector_destroy(&eigenvalues);
    igraph_matrix_destroy(&i_merges);
    igraph_matrix_destroy(&weighted_adjacency);
    igraph_destroy(&g);
}


void RCA(double **dataset, int numRows, int numCol, int **membership, double *modularity, int ***merges, int merge_dim[], double stats[], int bootstrap, double z_score) {
    double mean, stdev;
    double **min_max = get_min_max(dataset, numRows, numCol);
    double **s_data = standardize(dataset, numRows, numCol, min_max);
    double **R = relatioanlity_matrix(s_data, numRows, numCol);
    r_boorstrap(R,numRows,&mean,&stdev,bootstrap);
    get_significance(R,numRows,mean,stdev, z_score);

    partition(R,numRows,membership,modularity,merges,merge_dim);

    release_array(min_max,2);

    release_array(s_data,numRows);

    release_array(R,numRows);
    
    if(stats) {
        stats[0] = mean;
        stats[1] = stdev;
    }
    
}

