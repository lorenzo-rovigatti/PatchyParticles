#include <stdio.h>
#include <stdlib.h>


#include "cluster.h"


clusters* getClusters(int ncolloids)
{
	int i;
	
	clusters *c=malloc(sizeof(clusters));
	
	c->nodi=calloc(ncolloids,sizeof(struct vertice));
	
	// initialize c->graph
	// tutti i nodi si autopuntano
	for (i=0;i<ncolloids;i++)
	{
		//(c->nodi)[i].parent=c->nodi+i;
		(c->nodi)[i].parent=NULL;
		//(c->nodi)[i].size=1;
		(c->nodi)[i].size=0;
	}
	
	c->nbonds=0;
	
	return c;
}

void freeClusters(clusters *c)
{
	free(c->nodi);
	free(c);
}

void resetClusters(clusters *c,int ncolloids)
{
	int i;
	
	for (i=0;i<ncolloids;i++)
	{
		//(c->nodi)[i].parent=c->nodi+i;
		(c->nodi)[i].parent=NULL;
		//(c->nodi)[i].size=1;
		(c->nodi)[i].size=0;
	}
	c->nbonds=0;
}

int addNode(int particle,clusters *c)
{
	if (c->nodi[particle].size==0)
	{
		c->nodi[particle].size=1;
		c->nodi[particle].parent=c->nodi+particle;
		return 1;
	}
	else
		return 0;
}

int addBond(int particle1,int particle2,clusters *c)
{
	c->nbonds++;
	
	int root1=findRoot(c,particle1);
	int root2=findRoot(c,particle2);
	
	if (root1!=root2)
	{
		if ( (c->nodi)[root1].size > (c->nodi)[root2].size )
		{
			(c->nodi)[root2].parent=c->nodi+root1;
			(c->nodi)[root1].size+=(c->nodi)[root2].size;
			
			return (c->nodi)[root1].size;
		}
		else
		{
			(c->nodi)[root1].parent=c->nodi+root2;
			(c->nodi)[root2].size+=(c->nodi)[root1].size;
			
			return (c->nodi)[root2].size;
		}
	}
	
	return (c->nodi)[root1].size;
}

int findRoot(clusters *c,int node)
{
	
	struct vertice *punt=c->nodi+node;
	
	while ( punt!=punt->parent )
	{
		punt=punt->parent;
	}
	
	// the position is given by punt-c->nodi;
	return punt-c->nodi;
}

int isRoot(clusters *c,int i)
{
	if ( (c->nodi)[i].parent==c->nodi+i )
		return 1;
	else
		return 0;
}

int getClusterSize(int particle1,clusters *c)
{
	
	int root1=findRoot(c,particle1);
	
	return (c->nodi)[root1].size;
}

int sameCluster(int particle1,int particle2,clusters *c)
{
	int root1=findRoot(c,particle1);
	int root2=findRoot(c,particle2);
	
	if (root1!=root2)
		return 0;
	else
		return 1;
}


cluster_distribution* getClusterDistribution(clusters *c,int num_particles)
{
	cluster_distribution* cd=malloc(sizeof(cluster_distribution));
	
	cd->size=calloc(num_particles,sizeof(int));
	
	cd->num=0;
	
	int i;
	
	for (i=0;i<num_particles;i++)
	{
		if ( (c->nodi)[i].parent==c->nodi+i )
		{
			(cd->size)[cd->num]=(c->nodi)[i].size;
			cd->num++;
		}
	}
	
	// va da 1 a num_particles
	cd->distribution=calloc(num_particles+1,sizeof(int));
	
	cd->num_particles=num_particles;
	
	return cd;
}

void freeClusterDistribution(cluster_distribution *cd)
{
	free(cd->size);
	free(cd->distribution);
	free(cd);
}

void calculateClusterDistribution(cluster_distribution *cd)
{
	
	int i;
	
	for (i=0;i<cd->num;i++)
	{
		(cd->distribution)[(cd->size)[i]]++;
	}
}


void printClusterDistribution(cluster_distribution *cd,FILE *pfile)
{
	int i;
	
	for (i=1;i<=cd->num_particles;i++)
	{
		if ( (cd->distribution)[i]>0 )
		{
			fprintf(pfile,"%d %d\n",i,(cd->distribution)[i]);
		}
	}
	
}

