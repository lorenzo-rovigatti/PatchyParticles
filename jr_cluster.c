#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "jr_vector.h"
#include "jr_bilista.h"
#include "jr_interaction_map.h"
#include "jr_cluster.h"

static bilista *List_solid;
static clusters *Cluster_solid;
static double *Csd_buffer;
static int Csd_size;



void clustersConstructor(int num)
{
	Cluster_solid=getClusters(num);
	List_solid=bilistaGet(num);

	Csd_size=num;
	Csd_buffer=calloc(Csd_size,sizeof(double));
}

void clustersFree()
{
	free(Csd_buffer);
	freeClusters(Cluster_solid);
	bilistaFree(List_solid);
}


int clustersGetLargetCluster(int ncolloids,int *ImSolid,interactionmap *ime,int *num_solid)
{
	resetClusters(Cluster_solid,ncolloids);
	bilistaReset(List_solid,ncolloids);

	*num_solid=0;
	int max_size=0;
	int particle1,particle2;

	for (particle1=0;particle1<ncolloids;particle1++)
	{

		if (ImSolid[particle1]==1)
		{
			bilistaInsert(List_solid,particle1);
			(*num_solid)++;

			addNode(particle1,Cluster_solid);


				// USUAL
				int j;
				int hm=(ime->howmany[particle1]<5?ime->howmany[particle1]:5);
				for (j=0;j<hm;j++)
				{
					particle2=ime->with[particle1][j];

					if (ImSolid[particle2]==1)
					{

						addNode(particle2,Cluster_solid);

						int size=addBond(particle1,particle2,Cluster_solid);

						if (size>max_size)
							max_size=size;

					}
				}
			}
		}

		if ( (max_size==0) && (*num_solid>0) )
			max_size=1;


		return max_size;
}



void saveClusterDistribution()
{
	int i;
	for (i=0;i<Csd_size;i++)
	{
		Csd_buffer[i]=0.;
	}

	// scorriamo sui nodi solidi
	i=0;
	int solid_particle=List_solid[-1].next;
	while (solid_particle!=-1)
	{
		i++;

#ifdef DEBUG
		// per definizione i nodi solidi appartengono al cluster
		assert((Cluster_solid->nodi)[solid_particle].size!=0);
		assert((Cluster_solid->nodi)[solid_particle].parent!=NULL);
#endif

		// controlliamo se e' una root ed in caso estraiamo la size
		if ( (Cluster_solid->nodi)[solid_particle].parent==Cluster_solid->nodi+solid_particle )
		{
			int size=(Cluster_solid->nodi)[solid_particle].size;

			// aggiorniamo hist
			Csd_buffer[size]+=1.0;
		}

		solid_particle=List_solid[solid_particle].next;

	}

	// in Csd_buffer[0] ci va il numero di particelle liquide
	Csd_buffer[0]=Csd_size-i;

#ifdef DEBUG
	assert(i==Num_solid);
#endif
}

double* clustersGetCsd(int *size)
{
	*size=Csd_size;

	return Csd_buffer;
}

////////////////////////////////////////
////////////////////////////////////////
/// GENERAL PURPOSE FUNCTIONS //////////
////////////////////////////////////////
////////////////////////////////////////


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
