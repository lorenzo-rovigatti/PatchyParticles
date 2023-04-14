#include <stdlib.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "vector.h"
#include "smart_allocator.h"
#include "interaction_map.h"



interactionmap* createInteractionMap(int max_elements,int max_neighbours)
{
	interactionmap *i=malloc(sizeof(interactionmap));;
	
	i->num=0;
	i->num_bonds=0;
	i->who=calloc(max_elements,sizeof(int));
	//Matrix2D(i->with,max_elements,max_elements,int);
	i->howmany=calloc(max_elements,sizeof(int));
	i->size=max_elements;
	
	// la creazione degli array bidimensionali con un massimo numero di vicini
	// e' delicata. Dobbiamo lasciare dello spazio in coda a questi array nel
	// caso ci sia uno sforamento nel numero di vicini dell'ultima particella.
	// lasciamo in fondo max_elements spazi vuoti
	Matrix2DSafe(i->with,max_elements,max_neighbours,max_elements,int);
	Matrix2DSafe(i->rij2,max_elements,max_neighbours,max_elements,double);
	Matrix2DSafe(i->rij,max_elements,max_neighbours,max_elements,vector);
	
	return i;
}

void resetInteractionMap(interactionmap *i)
{
	i->num=0;
	i->num_bonds=0;
	int j;
	for (j=0;j<i->size;j++)
	{
		i->howmany[j]=0;
	}
}

void freeInteractionMap(interactionmap *i)
{
	free(i->who);
	Free2D(i->with);
	Free2D(i->rij);
	Free2D(i->rij2);
	free(i->howmany);
	free(i);
}

void buildImeFromIm(interactionmap *im,interactionmap *ime)
{
	ime->num=im->num;
	ime->size=im->size;
	ime->num_bonds=im->num_bonds;
	
	int i,index,j;
	for (i=0;i<im->num;i++)
	{
		ime->who[i]=im->who[i];
		
		for (index=0;index<im->howmany[i];index++)
		{
			j=im->with[i][index];
			double rij2=im->rij2[i][index];
			vector rij=im->rij[i][index];
			
			ime->rij[i][ime->howmany[i]]=rij;
			(ime->rij[j][ime->howmany[j]]).x=-rij.x;
			(ime->rij[j][ime->howmany[j]]).y=-rij.y;
			(ime->rij[j][ime->howmany[j]]).z=-rij.z;
			
			ime->rij2[i][ime->howmany[i]]=rij2;
			ime->rij2[j][ime->howmany[j]]=rij2;
			
			ime->with[i][ime->howmany[i]++]=j;
			ime->with[j][ime->howmany[j]++]=i;
		}
	}
}