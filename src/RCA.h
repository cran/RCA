#ifndef RCA_H
#define RCA_H
void RCA(double **dataset, int numRows, int numCol, int **membership, double *modularity, int ***merges, int merge_dim[], double stats[], int bootstrap, double z_score);

#endif

#ifndef IGRAPH_UNUSED
#define IGRAPH_UNUSED(x) (void)(x)
#endif

#ifndef MAX_COL
#define MAX_COL 1000
#endif

