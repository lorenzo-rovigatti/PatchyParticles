#ifndef CLUSTER_H
#define CLUSTER_H


struct vertice {
	struct vertice *parent;  // struttura del grafo con i vari nodi
	unsigned int size;       // se il nodo è una radice questa variabile indica la dimensione del cluster
};

typedef struct _clusters {
	struct vertice* nodi;
	int nbonds;
} clusters;

typedef struct _cluster_distribution {
	int *size;
	int num;
	int *distribution;
	int num_particles;
} cluster_distribution;


clusters* getClusters(int ncolloids);

void freeClusters(clusters *c);

void resetClusters(clusters *c,int ncolloids);

int addBond(int particle1,int particle2,clusters *c);

int addNode(int particle,clusters *c);

int findRoot(clusters *c,int node);

int isRoot(clusters *c,int i);

int getClusterSize(int particle1,clusters *c);

int sameCluster(int particle1,int particle2,clusters *c);

cluster_distribution* getClusterDistribution(clusters *c,int num_particles);

void freeClusterDistribution(cluster_distribution *cd);

void calculateClusterDistribution(cluster_distribution *cd);

void printClusterDistribution(cluster_distribution *cd,FILE *pfile);


#endif
