#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "jvector.h"
#include "interaction_map.h"
#include "global_definitions.h"
#include "secure_search.h"
#include "smart_allocator.h"
#include "log.h"
#include "cell_list.h"

#define SQR(x) ((x)*(x))

void celllistGetInteractionMap(args *argv,int argc);

void celllistConstructor(FILE *config_file,jvector *pos,int *ncolloids,jvector *box,double *cutoff,interactionmap *interactionList,double *cellsize,
			 listcell **list,void (**getInteractionMap)(args *argv,int argc),args **interactionMapArguments,int *numinteractionargs)
{
	SearchTable *s=searchNew();
	searchDouble("Cell_size",cellsize,s);
	searchFile(config_file,s);
	searchFree(s);

	*list=getList(*box,*cellsize,*ncolloids);

	*numinteractionargs=7;
	*interactionMapArguments=calloc(*numinteractionargs,sizeof(args));
	(*interactionMapArguments)[0].argv=(void*)pos;
	(*interactionMapArguments)[1].argv=(void*)ncolloids;
	(*interactionMapArguments)[2].argv=(void*)box;
	(*interactionMapArguments)[3].argv=(void*)interactionList;
	(*interactionMapArguments)[4].argv=(void*)cellsize;
	(*interactionMapArguments)[5].argv=(void*)(*list);
	(*interactionMapArguments)[6].argv=(void*)cutoff;

	*getInteractionMap=&celllistGetInteractionMap;

	logPrint("Cell size: %lf\n",*cellsize);
}


void celllistGetInteractionMap(args *argv,int argc)
{
	// lettura degli argomenti
	jvector *pos=(jvector*)argv[0].argv;
	int ncolloids=*((int*)argv[1].argv);
	jvector *box=(jvector*)argv[2].argv;
	interactionmap *im=(interactionmap *)argv[3].argv;
	double cellsize=*((double*)argv[4].argv);
	listcell *list=(listcell*)argv[5].argv;
	double cutoff=*((double*)argv[6].argv);

	fullUpdateList(list,pos,ncolloids,*box,cellsize);

	#ifdef NN
	memset(im->howmany,0,ncolloids*sizeof(int));
	calculateSystemInteractionMapOrdered(list,im,pos,box,cutoff);
	#else
	calculateSystemInteractionMap(list,im,pos,box,cutoff);
	#endif
}


// GENERIC FUNCTIONS

static int module(int n,int mo)
{
	n=n%mo;

	while (n<0)
	{
		n+=mo;
	}

	return n;
}

listcell* getList(jvector box_sides,double cutoff,int num_particles)
{
	listcell *l=malloc(sizeof(listcell));

	l->NumberCells_x=(int)(box_sides.x/cutoff);
	l->NumberCells_y=(int)(box_sides.y/cutoff);
	l->NumberCells_z=(int)(box_sides.z/cutoff);

	l->CellSize_x=box_sides.x/(double)l->NumberCells_x;
	l->CellSize_y=box_sides.y/(double)l->NumberCells_y;
	l->CellSize_z=box_sides.z/(double)l->NumberCells_z;

	l->HoC=(int*)calloc(l->NumberCells_x*l->NumberCells_y*l->NumberCells_z,sizeof(int));
	l->LinkedList=(int*)calloc(num_particles,sizeof(int));

	int i;
	int nnn=l->NumberCells_x*l->NumberCells_y*l->NumberCells_z;

	for (i=0;i<nnn;i++)
		(l->HoC)[i]=-1;

	return l;
}

void freeList(listcell *l)
{
	free(l->HoC);
	free(l->LinkedList);
	free(l);
}

void updateList(listcell *l,const jvector *pos,int num)
{
	int i;
	int ncell;                       // cell number
	int posx,posy,posz;              // cell coordinates
	int nnn;                         // total number of cells


	nnn=l->NumberCells_x*l->NumberCells_y*l->NumberCells_z;

	// HoC initialization
	for (i=0;i<nnn;i++)
		(l->HoC)[i]=-1;

	// colloids loop
	for (i=0;i<num;i++)
	{

		posx=(int)floor(pos[i].x/l->CellSize_x);
		posy=(int)floor(pos[i].y/l->CellSize_y);
		posz=(int)floor(pos[i].z/l->CellSize_z);

		posx=module(posx,l->NumberCells_x);
		posy=module(posy,l->NumberCells_y);
		posz=module(posz,l->NumberCells_z);

		ncell=posx+posy*l->NumberCells_x+posz*(l->NumberCells_x*l->NumberCells_y);

		(l->LinkedList)[i]=(l->HoC)[ncell];
		(l->HoC)[ncell]=i;
	}

}

void fullUpdateList(listcell *l,const jvector *pos,int num,jvector box_sides,double cutoff)
{
	int i;
	int ncell;                       // cell number
	int posx,posy,posz;              // cell coordinates

	int old=l->NumberCells_x*l->NumberCells_y*l->NumberCells_z;

	l->NumberCells_x=(int)(box_sides.x/cutoff);
	l->NumberCells_y=(int)(box_sides.y/cutoff);
	l->NumberCells_z=(int)(box_sides.z/cutoff);

	int nnn=l->NumberCells_x*l->NumberCells_y*l->NumberCells_z;

	if ( (nnn)>(old) )
		l->HoC=realloc(l->HoC,nnn*sizeof(int));

	l->CellSize_x=box_sides.x/(double)l->NumberCells_x;
	l->CellSize_y=box_sides.y/(double)l->NumberCells_y;
	l->CellSize_z=box_sides.z/(double)l->NumberCells_z;

	// HoC initialization
	for (i=0;i<nnn;i++)
		(l->HoC)[i]=-1;

	// colloids loop
	for (i=0;i<num;i++)
	{

		posx=(int)floor(pos[i].x/l->CellSize_x);
		posy=(int)floor(pos[i].y/l->CellSize_y);
		posz=(int)floor(pos[i].z/l->CellSize_z);

		posx=module(posx,l->NumberCells_x);
		posy=module(posy,l->NumberCells_y);
		posz=module(posz,l->NumberCells_z);

		ncell=posx+posy*l->NumberCells_x+posz*(l->NumberCells_x*l->NumberCells_y);

		(l->LinkedList)[i]=(l->HoC)[ncell];
		(l->HoC)[ncell]=i;
	}

}


void calculateInteractionMap(listcell *l,interactionmap *im)
{
	int i;
	int nn,nnn;                      // Number of cells on surface and in box
	int particle1,particle2;         // interacting particles
	int neighbour;                   // Neighbour cell
	int ix,iy,iz,alpha;              // present cell coordinates
	int jx,jy,jz;                    // next cell coordinates
	int kx,ky,kz;                    // previous cell coordinates


	nn=l->NumberCells_x*l->NumberCells_y;
	nnn=l->NumberCells_z*nn;


	int number_particles=0;

	// scan cells
	for (i=0;i<nnn;i++)
	{
		// reconstruct cartesian coordinates
		iz=i/nn;
		alpha=i%nn;
		ix=alpha%l->NumberCells_x;
		iy=alpha/l->NumberCells_x;

		jx=(ix+1)%l->NumberCells_x;
		jy=(iy+1)%l->NumberCells_y;
		jz=(iz+1)%l->NumberCells_z;

		kx=module(ix-1,l->NumberCells_x);
		ky=module(iy-1,l->NumberCells_y);
		kz=module(iz-1,l->NumberCells_z);


		particle1=(l->HoC)[i];

		while (particle1!=-1)
		{
			//interactions with particles in the same cell
			// and in neighbouring cells
			im->who[particle1]=particle1;
			number_particles++;
			int howmany_particle1=0;

			// 0 same cell
			particle2=(l->LinkedList)[particle1];
			while (particle2!=-1)
			{
				im->with[particle1][howmany_particle1++]=particle2;
				particle2=(l->LinkedList)[particle2];
			}

			// 1 cell y+1
			neighbour=(ix)+(jy)*l->NumberCells_x+(iz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				im->with[particle1][howmany_particle1++]=particle2;
				particle2=(l->LinkedList)[particle2];
			}

			// 2 cell x+1 y+1
			neighbour=(jx)+(jy)*l->NumberCells_x+(iz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				im->with[particle1][howmany_particle1++]=particle2;
				particle2=(l->LinkedList)[particle2];
			}

			// 3 cell x+1
			neighbour=(jx)+(iy)*l->NumberCells_x+(iz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				im->with[particle1][howmany_particle1++]=particle2;
				particle2=(l->LinkedList)[particle2];
			}

			// 4 cell x+1 y-1
			neighbour=(jx)+(ky)*l->NumberCells_x+(iz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				im->with[particle1][howmany_particle1++]=particle2;
				particle2=(l->LinkedList)[particle2];
			}

			// 5 cell z+1 y+1
			neighbour=(ix)+(jy)*l->NumberCells_x+(jz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				im->with[particle1][howmany_particle1++]=particle2;
				particle2=(l->LinkedList)[particle2];
			}

			// 6 cell z+1
			neighbour=(ix)+(iy)*l->NumberCells_x+(jz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				im->with[particle1][howmany_particle1++]=particle2;
				particle2=(l->LinkedList)[particle2];
			}

			// 7 cell z+1 x+1 y+1
			neighbour=(jx)+(jy)*l->NumberCells_x+(jz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				im->with[particle1][howmany_particle1++]=particle2;
				particle2=(l->LinkedList)[particle2];
			}

			// 8 cell z+1 x+1
			neighbour=(jx)+(iy)*l->NumberCells_x+(jz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				im->with[particle1][howmany_particle1++]=particle2;
				particle2=(l->LinkedList)[particle2];
			}

			// 9 cell z+1 x+1 y-1
			neighbour=(jx)+(ky)*l->NumberCells_x+(jz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				im->with[particle1][howmany_particle1++]=particle2;
				particle2=(l->LinkedList)[particle2];
			}

			// 10 cell z-1 y+1
			neighbour=(ix)+(jy)*l->NumberCells_x+(kz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				im->with[particle1][howmany_particle1++]=particle2;
				particle2=(l->LinkedList)[particle2];
			}

			// 11 cell x+1 y+1 z-1
			neighbour=(jx)+(jy)*l->NumberCells_x+(kz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				im->with[particle1][howmany_particle1++]=particle2;
				particle2=(l->LinkedList)[particle2];
			}

			// 12 cell x+1 z-1
			neighbour=(jx)+(iy)*l->NumberCells_x+(kz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				im->with[particle1][howmany_particle1++]=particle2;
				particle2=(l->LinkedList)[particle2];
			}

			// 13 cell z-1 x+1 y-1
			neighbour=(jx)+(ky)*l->NumberCells_x+(kz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				im->with[particle1][howmany_particle1++]=particle2;
				particle2=(l->LinkedList)[particle2];
			}

			im->howmany[particle1]=howmany_particle1;

			particle1=(l->LinkedList)[particle1];
		}
	}

	im->num=number_particles;

	assert(number_particles==im->size);
}

static void insertionSort_dist2dist(int *number,int num_el,double *dist2,double dist2_el,jvector *dist,jvector dist_el,int *length)
{
	int l=*length;
	number[l]=num_el;
	dist2[l]=dist2_el;
	dist[l]=dist_el;

	while ((l>0) && (dist2[l]<dist2[l-1]))
	{
		double buffer_dist2;
		buffer_dist2=dist2[l];
		dist2[l]=dist2[l-1];
		dist2[l-1]=buffer_dist2;

		int buffer_num;
		buffer_num=number[l];
		number[l]=number[l-1];
		number[l-1]=buffer_num;

		jvector buffer_dist;
		buffer_dist=dist[l];
		dist[l]=dist[l-1];
		dist[l-1]=buffer_dist;


		l--;
	}
	(*length)++;
}

int calculateSystemInteractionMap(listcell *l,interactionmap *im,jvector *pos,jvector *box,double cutoff)
{
	int i;
	int nn,nnn;                      // Number of cells on surface and in box
	int particle1,particle2;         // interacting particles
	int neighbour;                   // Neighbour cell
	int ix,iy,iz,alpha;              // present cell coordinates
	int jx,jy,jz;                    // next cell coordinates
	int kx,ky,kz;                    // previous cell coordinates


	nn=l->NumberCells_x*l->NumberCells_y;
	nnn=l->NumberCells_z*nn;


	int number_particles=0;

	double cutoff2=SQR(cutoff);

	im->num_bonds=0;

	// scan cells
	for (i=0;i<nnn;i++)
	{
		// reconstruct cartesian coordinates
		iz=i/nn;
		alpha=i%nn;
		ix=alpha%l->NumberCells_x;
		iy=alpha/l->NumberCells_x;

		jx=(ix+1)%l->NumberCells_x;
		jy=(iy+1)%l->NumberCells_y;
		jz=(iz+1)%l->NumberCells_z;

		kx=module(ix-1,l->NumberCells_x);
		ky=module(iy-1,l->NumberCells_y);
		kz=module(iz-1,l->NumberCells_z);


		particle1=(l->HoC)[i];

		while (particle1!=-1)
		{
			//interactions with particles in the same cell
			// and in neighbouring cells
			im->who[particle1]=particle1;
			number_particles++;
			int howmany_particle1=0;
			jvector dist;
			double dist2;

			// 0 same cell
			particle2=(l->LinkedList)[particle1];
			while (particle2!=-1)
			{
				dist.x=pos[particle1].x-pos[particle2].x;
				dist.y=pos[particle1].y-pos[particle2].y;
				dist.z=pos[particle1].z-pos[particle2].z;

				dist.x-=box->x*rint(dist.x/box->x);
				dist.y-=box->y*rint(dist.y/box->y);
				dist.z-=box->z*rint(dist.z/box->z);

				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);

				if (dist2<cutoff2)
				{
					im->rij2[particle1][howmany_particle1]=dist2;
					im->rij[particle1][howmany_particle1]=dist;
					im->with[particle1][howmany_particle1++]=particle2;
				}
				particle2=(l->LinkedList)[particle2];
			}

			// 1 cell y+1
			neighbour=(ix)+(jy)*l->NumberCells_x+(iz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				dist.x=pos[particle1].x-pos[particle2].x;
				dist.y=pos[particle1].y-pos[particle2].y;
				dist.z=pos[particle1].z-pos[particle2].z;

				dist.x-=box->x*rint(dist.x/box->x);
				dist.y-=box->y*rint(dist.y/box->y);
				dist.z-=box->z*rint(dist.z/box->z);

				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);

				if (dist2<cutoff2)
				{
					im->rij2[particle1][howmany_particle1]=dist2;
					im->rij[particle1][howmany_particle1]=dist;
					im->with[particle1][howmany_particle1++]=particle2;
				}
				particle2=(l->LinkedList)[particle2];
			}

			// 2 cell x+1 y+1
			neighbour=(jx)+(jy)*l->NumberCells_x+(iz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				dist.x=pos[particle1].x-pos[particle2].x;
				dist.y=pos[particle1].y-pos[particle2].y;
				dist.z=pos[particle1].z-pos[particle2].z;

				dist.x-=box->x*rint(dist.x/box->x);
				dist.y-=box->y*rint(dist.y/box->y);
				dist.z-=box->z*rint(dist.z/box->z);

				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);

				if (dist2<cutoff2)
				{
					im->rij2[particle1][howmany_particle1]=dist2;
					im->rij[particle1][howmany_particle1]=dist;
					im->with[particle1][howmany_particle1++]=particle2;
				}
				particle2=(l->LinkedList)[particle2];
			}

			// 3 cell x+1
			neighbour=(jx)+(iy)*l->NumberCells_x+(iz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				dist.x=pos[particle1].x-pos[particle2].x;
				dist.y=pos[particle1].y-pos[particle2].y;
				dist.z=pos[particle1].z-pos[particle2].z;

				dist.x-=box->x*rint(dist.x/box->x);
				dist.y-=box->y*rint(dist.y/box->y);
				dist.z-=box->z*rint(dist.z/box->z);

				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);

				if (dist2<cutoff2)
				{
					im->rij2[particle1][howmany_particle1]=dist2;
					im->rij[particle1][howmany_particle1]=dist;
					im->with[particle1][howmany_particle1++]=particle2;
				}
				particle2=(l->LinkedList)[particle2];
			}

			// 4 cell x+1 y-1
			neighbour=(jx)+(ky)*l->NumberCells_x+(iz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				dist.x=pos[particle1].x-pos[particle2].x;
				dist.y=pos[particle1].y-pos[particle2].y;
				dist.z=pos[particle1].z-pos[particle2].z;

				dist.x-=box->x*rint(dist.x/box->x);
				dist.y-=box->y*rint(dist.y/box->y);
				dist.z-=box->z*rint(dist.z/box->z);

				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);

				if (dist2<cutoff2)
				{
					im->rij2[particle1][howmany_particle1]=dist2;
					im->rij[particle1][howmany_particle1]=dist;
					im->with[particle1][howmany_particle1++]=particle2;
				}
				particle2=(l->LinkedList)[particle2];
			}

			// 5 cell z+1 y+1
			neighbour=(ix)+(jy)*l->NumberCells_x+(jz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				dist.x=pos[particle1].x-pos[particle2].x;
				dist.y=pos[particle1].y-pos[particle2].y;
				dist.z=pos[particle1].z-pos[particle2].z;

				dist.x-=box->x*rint(dist.x/box->x);
				dist.y-=box->y*rint(dist.y/box->y);
				dist.z-=box->z*rint(dist.z/box->z);

				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);

				if (dist2<cutoff2)
				{
					im->rij2[particle1][howmany_particle1]=dist2;
					im->rij[particle1][howmany_particle1]=dist;
					im->with[particle1][howmany_particle1++]=particle2;
				}
				particle2=(l->LinkedList)[particle2];
			}

			// 6 cell z+1
			neighbour=(ix)+(iy)*l->NumberCells_x+(jz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				dist.x=pos[particle1].x-pos[particle2].x;
				dist.y=pos[particle1].y-pos[particle2].y;
				dist.z=pos[particle1].z-pos[particle2].z;

				dist.x-=box->x*rint(dist.x/box->x);
				dist.y-=box->y*rint(dist.y/box->y);
				dist.z-=box->z*rint(dist.z/box->z);

				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);

				if (dist2<cutoff2)
				{
					im->rij2[particle1][howmany_particle1]=dist2;
					im->rij[particle1][howmany_particle1]=dist;
					im->with[particle1][howmany_particle1++]=particle2;
				}
				particle2=(l->LinkedList)[particle2];
			}

			// 7 cell z+1 x+1 y+1
			neighbour=(jx)+(jy)*l->NumberCells_x+(jz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				dist.x=pos[particle1].x-pos[particle2].x;
				dist.y=pos[particle1].y-pos[particle2].y;
				dist.z=pos[particle1].z-pos[particle2].z;

				dist.x-=box->x*rint(dist.x/box->x);
				dist.y-=box->y*rint(dist.y/box->y);
				dist.z-=box->z*rint(dist.z/box->z);

				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);

				if (dist2<cutoff2)
				{
					im->rij2[particle1][howmany_particle1]=dist2;
					im->rij[particle1][howmany_particle1]=dist;
					im->with[particle1][howmany_particle1++]=particle2;
				}
				particle2=(l->LinkedList)[particle2];
			}

			// 8 cell z+1 x+1
			neighbour=(jx)+(iy)*l->NumberCells_x+(jz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				dist.x=pos[particle1].x-pos[particle2].x;
				dist.y=pos[particle1].y-pos[particle2].y;
				dist.z=pos[particle1].z-pos[particle2].z;

				dist.x-=box->x*rint(dist.x/box->x);
				dist.y-=box->y*rint(dist.y/box->y);
				dist.z-=box->z*rint(dist.z/box->z);

				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);

				if (dist2<cutoff2)
				{
					im->rij2[particle1][howmany_particle1]=dist2;
					im->rij[particle1][howmany_particle1]=dist;
					im->with[particle1][howmany_particle1++]=particle2;
				}
				particle2=(l->LinkedList)[particle2];
			}

			// 9 cell z+1 x+1 y-1
			neighbour=(jx)+(ky)*l->NumberCells_x+(jz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				dist.x=pos[particle1].x-pos[particle2].x;
				dist.y=pos[particle1].y-pos[particle2].y;
				dist.z=pos[particle1].z-pos[particle2].z;

				dist.x-=box->x*rint(dist.x/box->x);
				dist.y-=box->y*rint(dist.y/box->y);
				dist.z-=box->z*rint(dist.z/box->z);

				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);

				if (dist2<cutoff2)
				{
					im->rij2[particle1][howmany_particle1]=dist2;
					im->rij[particle1][howmany_particle1]=dist;
					im->with[particle1][howmany_particle1++]=particle2;
				}
				particle2=(l->LinkedList)[particle2];
			}

			// 10 cell z-1 y+1
			neighbour=(ix)+(jy)*l->NumberCells_x+(kz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				dist.x=pos[particle1].x-pos[particle2].x;
				dist.y=pos[particle1].y-pos[particle2].y;
				dist.z=pos[particle1].z-pos[particle2].z;

				dist.x-=box->x*rint(dist.x/box->x);
				dist.y-=box->y*rint(dist.y/box->y);
				dist.z-=box->z*rint(dist.z/box->z);

				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);

				if (dist2<cutoff2)
				{
					im->rij2[particle1][howmany_particle1]=dist2;
					im->rij[particle1][howmany_particle1]=dist;
					im->with[particle1][howmany_particle1++]=particle2;
				}
				particle2=(l->LinkedList)[particle2];
			}

			// 11 cell x+1 y+1 z-1
			neighbour=(jx)+(jy)*l->NumberCells_x+(kz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				dist.x=pos[particle1].x-pos[particle2].x;
				dist.y=pos[particle1].y-pos[particle2].y;
				dist.z=pos[particle1].z-pos[particle2].z;

				dist.x-=box->x*rint(dist.x/box->x);
				dist.y-=box->y*rint(dist.y/box->y);
				dist.z-=box->z*rint(dist.z/box->z);

				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);

				if (dist2<cutoff2)
				{
					im->rij2[particle1][howmany_particle1]=dist2;
					im->rij[particle1][howmany_particle1]=dist;
					im->with[particle1][howmany_particle1++]=particle2;
				}
				particle2=(l->LinkedList)[particle2];
			}

			// 12 cell x+1 z-1
			neighbour=(jx)+(iy)*l->NumberCells_x+(kz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				dist.x=pos[particle1].x-pos[particle2].x;
				dist.y=pos[particle1].y-pos[particle2].y;
				dist.z=pos[particle1].z-pos[particle2].z;

				dist.x-=box->x*rint(dist.x/box->x);
				dist.y-=box->y*rint(dist.y/box->y);
				dist.z-=box->z*rint(dist.z/box->z);

				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);

				if (dist2<cutoff2)
				{
					im->rij2[particle1][howmany_particle1]=dist2;
					im->rij[particle1][howmany_particle1]=dist;
					im->with[particle1][howmany_particle1++]=particle2;
				}
				particle2=(l->LinkedList)[particle2];
			}

			// 13 cell z-1 x+1 y-1
			neighbour=(jx)+(ky)*l->NumberCells_x+(kz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				dist.x=pos[particle1].x-pos[particle2].x;
				dist.y=pos[particle1].y-pos[particle2].y;
				dist.z=pos[particle1].z-pos[particle2].z;

				dist.x-=box->x*rint(dist.x/box->x);
				dist.y-=box->y*rint(dist.y/box->y);
				dist.z-=box->z*rint(dist.z/box->z);

				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);

				if (dist2<cutoff2)
				{
					im->rij2[particle1][howmany_particle1]=dist2;
					im->rij[particle1][howmany_particle1]=dist;
					im->with[particle1][howmany_particle1++]=particle2;
				}
				particle2=(l->LinkedList)[particle2];
			}

			im->howmany[particle1]=howmany_particle1;
			particle1=(l->LinkedList)[particle1];
		}
	}

	im->num=number_particles;

	assert(number_particles==im->size);

	int max_neighbours=0;
	for (i=0;i<number_particles;i++)
	{
		if (im->howmany[i]>max_neighbours)
			max_neighbours=im->howmany[i];
	}

	return max_neighbours;
}

int calculateSystemInteractionMapOrdered(listcell *l,interactionmap *im,jvector *pos,jvector *box,double cutoff)
{
	int i;
	int nn,nnn;                      // Number of cells on surface and in box
	int particle1,particle2;         // interacting particles
	int neighbour;                   // Neighbour cell
	int ix,iy,iz,alpha;              // present cell coordinates
	int jx,jy,jz;                    // next cell coordinates
	int kx,ky,kz;                    // previous cell coordinates


	nn=l->NumberCells_x*l->NumberCells_y;
	nnn=l->NumberCells_z*nn;


	int number_particles=0;

	double cutoff2=SQR(cutoff);

	im->num_bonds=0;

	// scan cells
	for (i=0;i<nnn;i++)
	{
		// reconstruct cartesian coordinates
		iz=i/nn;
		alpha=i%nn;
		ix=alpha%l->NumberCells_x;
		iy=alpha/l->NumberCells_x;

		jx=(ix+1)%l->NumberCells_x;
		jy=(iy+1)%l->NumberCells_y;
		jz=(iz+1)%l->NumberCells_z;

		kx=module(ix-1,l->NumberCells_x);
		ky=module(iy-1,l->NumberCells_y);
		kz=module(iz-1,l->NumberCells_z);


		particle1=(l->HoC)[i];

		while (particle1!=-1)
		{
			//interactions with particles in the same cell
			// and in neighbouring cells
			im->who[particle1]=particle1;
			number_particles++;
			//int howmany_particle1=0;
			jvector dist;
			double dist2;

			// 0 same cell
			particle2=(l->LinkedList)[particle1];
			while (particle2!=-1)
			{
				dist.x=pos[particle1].x-pos[particle2].x;
				dist.y=pos[particle1].y-pos[particle2].y;
				dist.z=pos[particle1].z-pos[particle2].z;

				dist.x-=box->x*rint(dist.x/box->x);
				dist.y-=box->y*rint(dist.y/box->y);
				dist.z-=box->z*rint(dist.z/box->z);

				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);

				if (dist2<cutoff2)
				{
					if (particle1<particle2)
					{
						insertionSort_dist2dist(im->with[particle1],particle2,im->rij2[particle1],dist2,im->rij[particle1],dist,im->howmany+particle1);
					}
					else
					{
						dist.x*=-1.;
						dist.y*=-1.;
						dist.z*=-1.;
						insertionSort_dist2dist(im->with[particle2],particle1,im->rij2[particle2],dist2,im->rij[particle2],dist,im->howmany+particle2);
					}
				}
				particle2=(l->LinkedList)[particle2];
			}

			// 1 cell y+1
			neighbour=(ix)+(jy)*l->NumberCells_x+(iz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				dist.x=pos[particle1].x-pos[particle2].x;
				dist.y=pos[particle1].y-pos[particle2].y;
				dist.z=pos[particle1].z-pos[particle2].z;

				dist.x-=box->x*rint(dist.x/box->x);
				dist.y-=box->y*rint(dist.y/box->y);
				dist.z-=box->z*rint(dist.z/box->z);

				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);

				if (dist2<cutoff2)
				{
					if (particle1<particle2)
					{
						insertionSort_dist2dist(im->with[particle1],particle2,im->rij2[particle1],dist2,im->rij[particle1],dist,im->howmany+particle1);
					}
					else
					{
						dist.x*=-1.;
						dist.y*=-1.;
						dist.z*=-1.;
						insertionSort_dist2dist(im->with[particle2],particle1,im->rij2[particle2],dist2,im->rij[particle2],dist,im->howmany+particle2);
					}
				}
				particle2=(l->LinkedList)[particle2];
			}

			// 2 cell x+1 y+1
			neighbour=(jx)+(jy)*l->NumberCells_x+(iz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				dist.x=pos[particle1].x-pos[particle2].x;
				dist.y=pos[particle1].y-pos[particle2].y;
				dist.z=pos[particle1].z-pos[particle2].z;

				dist.x-=box->x*rint(dist.x/box->x);
				dist.y-=box->y*rint(dist.y/box->y);
				dist.z-=box->z*rint(dist.z/box->z);

				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);

				if (dist2<cutoff2)
				{
					if (particle1<particle2)
					{
						insertionSort_dist2dist(im->with[particle1],particle2,im->rij2[particle1],dist2,im->rij[particle1],dist,im->howmany+particle1);
					}
					else
					{
						dist.x*=-1.;
						dist.y*=-1.;
						dist.z*=-1.;
						insertionSort_dist2dist(im->with[particle2],particle1,im->rij2[particle2],dist2,im->rij[particle2],dist,im->howmany+particle2);
					}
				}
				particle2=(l->LinkedList)[particle2];
			}

			// 3 cell x+1
			neighbour=(jx)+(iy)*l->NumberCells_x+(iz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				dist.x=pos[particle1].x-pos[particle2].x;
				dist.y=pos[particle1].y-pos[particle2].y;
				dist.z=pos[particle1].z-pos[particle2].z;

				dist.x-=box->x*rint(dist.x/box->x);
				dist.y-=box->y*rint(dist.y/box->y);
				dist.z-=box->z*rint(dist.z/box->z);

				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);

				if (dist2<cutoff2)
				{
					if (particle1<particle2)
					{
						insertionSort_dist2dist(im->with[particle1],particle2,im->rij2[particle1],dist2,im->rij[particle1],dist,im->howmany+particle1);
					}
					else
					{
						dist.x*=-1.;
						dist.y*=-1.;
						dist.z*=-1.;
						insertionSort_dist2dist(im->with[particle2],particle1,im->rij2[particle2],dist2,im->rij[particle2],dist,im->howmany+particle2);
					}
				}
				particle2=(l->LinkedList)[particle2];
			}

			// 4 cell x+1 y-1
			neighbour=(jx)+(ky)*l->NumberCells_x+(iz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				dist.x=pos[particle1].x-pos[particle2].x;
				dist.y=pos[particle1].y-pos[particle2].y;
				dist.z=pos[particle1].z-pos[particle2].z;

				dist.x-=box->x*rint(dist.x/box->x);
				dist.y-=box->y*rint(dist.y/box->y);
				dist.z-=box->z*rint(dist.z/box->z);

				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);

				if (dist2<cutoff2)
				{
					if (particle1<particle2)
					{
						insertionSort_dist2dist(im->with[particle1],particle2,im->rij2[particle1],dist2,im->rij[particle1],dist,im->howmany+particle1);
					}
					else
					{
						dist.x*=-1.;
						dist.y*=-1.;
						dist.z*=-1.;
						insertionSort_dist2dist(im->with[particle2],particle1,im->rij2[particle2],dist2,im->rij[particle2],dist,im->howmany+particle2);
					}
				}
				particle2=(l->LinkedList)[particle2];
			}

			// 5 cell z+1 y+1
			neighbour=(ix)+(jy)*l->NumberCells_x+(jz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				dist.x=pos[particle1].x-pos[particle2].x;
				dist.y=pos[particle1].y-pos[particle2].y;
				dist.z=pos[particle1].z-pos[particle2].z;

				dist.x-=box->x*rint(dist.x/box->x);
				dist.y-=box->y*rint(dist.y/box->y);
				dist.z-=box->z*rint(dist.z/box->z);

				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);

				if (dist2<cutoff2)
				{
					if (particle1<particle2)
					{
						insertionSort_dist2dist(im->with[particle1],particle2,im->rij2[particle1],dist2,im->rij[particle1],dist,im->howmany+particle1);
					}
					else
					{
						dist.x*=-1.;
						dist.y*=-1.;
						dist.z*=-1.;
						insertionSort_dist2dist(im->with[particle2],particle1,im->rij2[particle2],dist2,im->rij[particle2],dist,im->howmany+particle2);
					}
				}
				particle2=(l->LinkedList)[particle2];
			}

			// 6 cell z+1
			neighbour=(ix)+(iy)*l->NumberCells_x+(jz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				dist.x=pos[particle1].x-pos[particle2].x;
				dist.y=pos[particle1].y-pos[particle2].y;
				dist.z=pos[particle1].z-pos[particle2].z;

				dist.x-=box->x*rint(dist.x/box->x);
				dist.y-=box->y*rint(dist.y/box->y);
				dist.z-=box->z*rint(dist.z/box->z);

				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);

				if (dist2<cutoff2)
				{
					if (particle1<particle2)
					{
						insertionSort_dist2dist(im->with[particle1],particle2,im->rij2[particle1],dist2,im->rij[particle1],dist,im->howmany+particle1);
					}
					else
					{
						dist.x*=-1.;
						dist.y*=-1.;
						dist.z*=-1.;
						insertionSort_dist2dist(im->with[particle2],particle1,im->rij2[particle2],dist2,im->rij[particle2],dist,im->howmany+particle2);
					}
				}
				particle2=(l->LinkedList)[particle2];
			}

			// 7 cell z+1 x+1 y+1
			neighbour=(jx)+(jy)*l->NumberCells_x+(jz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				dist.x=pos[particle1].x-pos[particle2].x;
				dist.y=pos[particle1].y-pos[particle2].y;
				dist.z=pos[particle1].z-pos[particle2].z;

				dist.x-=box->x*rint(dist.x/box->x);
				dist.y-=box->y*rint(dist.y/box->y);
				dist.z-=box->z*rint(dist.z/box->z);

				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);

				if (dist2<cutoff2)
				{
					if (particle1<particle2)
					{
						insertionSort_dist2dist(im->with[particle1],particle2,im->rij2[particle1],dist2,im->rij[particle1],dist,im->howmany+particle1);
					}
					else
					{
						dist.x*=-1.;
						dist.y*=-1.;
						dist.z*=-1.;
						insertionSort_dist2dist(im->with[particle2],particle1,im->rij2[particle2],dist2,im->rij[particle2],dist,im->howmany+particle2);
					}
				}
				particle2=(l->LinkedList)[particle2];
			}

			// 8 cell z+1 x+1
			neighbour=(jx)+(iy)*l->NumberCells_x+(jz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				dist.x=pos[particle1].x-pos[particle2].x;
				dist.y=pos[particle1].y-pos[particle2].y;
				dist.z=pos[particle1].z-pos[particle2].z;

				dist.x-=box->x*rint(dist.x/box->x);
				dist.y-=box->y*rint(dist.y/box->y);
				dist.z-=box->z*rint(dist.z/box->z);

				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);

				if (dist2<cutoff2)
				{
					if (particle1<particle2)
					{
						insertionSort_dist2dist(im->with[particle1],particle2,im->rij2[particle1],dist2,im->rij[particle1],dist,im->howmany+particle1);
					}
					else
					{
						dist.x*=-1.;
						dist.y*=-1.;
						dist.z*=-1.;
						insertionSort_dist2dist(im->with[particle2],particle1,im->rij2[particle2],dist2,im->rij[particle2],dist,im->howmany+particle2);
					}
				}
				particle2=(l->LinkedList)[particle2];
			}

			// 9 cell z+1 x+1 y-1
			neighbour=(jx)+(ky)*l->NumberCells_x+(jz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				dist.x=pos[particle1].x-pos[particle2].x;
				dist.y=pos[particle1].y-pos[particle2].y;
				dist.z=pos[particle1].z-pos[particle2].z;

				dist.x-=box->x*rint(dist.x/box->x);
				dist.y-=box->y*rint(dist.y/box->y);
				dist.z-=box->z*rint(dist.z/box->z);

				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);

				if (dist2<cutoff2)
				{
					if (particle1<particle2)
					{
						insertionSort_dist2dist(im->with[particle1],particle2,im->rij2[particle1],dist2,im->rij[particle1],dist,im->howmany+particle1);
					}
					else
					{
						dist.x*=-1.;
						dist.y*=-1.;
						dist.z*=-1.;
						insertionSort_dist2dist(im->with[particle2],particle1,im->rij2[particle2],dist2,im->rij[particle2],dist,im->howmany+particle2);
					}
				}
				particle2=(l->LinkedList)[particle2];
			}

			// 10 cell z-1 y+1
			neighbour=(ix)+(jy)*l->NumberCells_x+(kz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				dist.x=pos[particle1].x-pos[particle2].x;
				dist.y=pos[particle1].y-pos[particle2].y;
				dist.z=pos[particle1].z-pos[particle2].z;

				dist.x-=box->x*rint(dist.x/box->x);
				dist.y-=box->y*rint(dist.y/box->y);
				dist.z-=box->z*rint(dist.z/box->z);

				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);

				if (dist2<cutoff2)
				{
					if (particle1<particle2)
					{
						insertionSort_dist2dist(im->with[particle1],particle2,im->rij2[particle1],dist2,im->rij[particle1],dist,im->howmany+particle1);
					}
					else
					{
						dist.x*=-1.;
						dist.y*=-1.;
						dist.z*=-1.;
						insertionSort_dist2dist(im->with[particle2],particle1,im->rij2[particle2],dist2,im->rij[particle2],dist,im->howmany+particle2);
					}
				}
				particle2=(l->LinkedList)[particle2];
			}

			// 11 cell x+1 y+1 z-1
			neighbour=(jx)+(jy)*l->NumberCells_x+(kz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				dist.x=pos[particle1].x-pos[particle2].x;
				dist.y=pos[particle1].y-pos[particle2].y;
				dist.z=pos[particle1].z-pos[particle2].z;

				dist.x-=box->x*rint(dist.x/box->x);
				dist.y-=box->y*rint(dist.y/box->y);
				dist.z-=box->z*rint(dist.z/box->z);

				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);

				if (dist2<cutoff2)
				{
					if (particle1<particle2)
					{
						insertionSort_dist2dist(im->with[particle1],particle2,im->rij2[particle1],dist2,im->rij[particle1],dist,im->howmany+particle1);
					}
					else
					{
						dist.x*=-1.;
						dist.y*=-1.;
						dist.z*=-1.;
						insertionSort_dist2dist(im->with[particle2],particle1,im->rij2[particle2],dist2,im->rij[particle2],dist,im->howmany+particle2);
					}
				}
				particle2=(l->LinkedList)[particle2];
			}

			// 12 cell x+1 z-1
			neighbour=(jx)+(iy)*l->NumberCells_x+(kz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				dist.x=pos[particle1].x-pos[particle2].x;
				dist.y=pos[particle1].y-pos[particle2].y;
				dist.z=pos[particle1].z-pos[particle2].z;

				dist.x-=box->x*rint(dist.x/box->x);
				dist.y-=box->y*rint(dist.y/box->y);
				dist.z-=box->z*rint(dist.z/box->z);

				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);

				if (dist2<cutoff2)
				{
					if (particle1<particle2)
					{
						insertionSort_dist2dist(im->with[particle1],particle2,im->rij2[particle1],dist2,im->rij[particle1],dist,im->howmany+particle1);
					}
					else
					{
						dist.x*=-1.;
						dist.y*=-1.;
						dist.z*=-1.;
						insertionSort_dist2dist(im->with[particle2],particle1,im->rij2[particle2],dist2,im->rij[particle2],dist,im->howmany+particle2);
					}
				}
				particle2=(l->LinkedList)[particle2];
			}

			// 13 cell z-1 x+1 y-1
			neighbour=(jx)+(ky)*l->NumberCells_x+(kz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				dist.x=pos[particle1].x-pos[particle2].x;
				dist.y=pos[particle1].y-pos[particle2].y;
				dist.z=pos[particle1].z-pos[particle2].z;

				dist.x-=box->x*rint(dist.x/box->x);
				dist.y-=box->y*rint(dist.y/box->y);
				dist.z-=box->z*rint(dist.z/box->z);

				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);

				if (dist2<cutoff2)
				{
					if (particle1<particle2)
					{
						insertionSort_dist2dist(im->with[particle1],particle2,im->rij2[particle1],dist2,im->rij[particle1],dist,im->howmany+particle1);
					}
					else
					{
						dist.x*=-1.;
						dist.y*=-1.;
						dist.z*=-1.;
						insertionSort_dist2dist(im->with[particle2],particle1,im->rij2[particle2],dist2,im->rij[particle2],dist,im->howmany+particle2);
					}
				}
				particle2=(l->LinkedList)[particle2];
			}

			//im->howmany[particle1]=howmany_particle1;
			particle1=(l->LinkedList)[particle1];
		}
	}

	im->num=number_particles;

	assert(number_particles==im->size);

	int max_neighbours=0;
	for (i=0;i<number_particles;i++)
	{
		if (im->howmany[i]>max_neighbours)
			max_neighbours=im->howmany[i];
	}

	return max_neighbours;
}




void calculateInteractionMapWithCutoffDistance(listcell *l,interactionmap *im,jvector *pos,jvector *box,double cutoff)
{
	int i;
	int nn,nnn;                      // Number of cells on surface and in box
	int particle1,particle2;         // interacting particles
	int neighbour;                   // Neighbour cell
	int ix,iy,iz,alpha;              // present cell coordinates
	int jx,jy,jz;                    // next cell coordinates
	int kx,ky,kz;                    // previous cell coordinates


	nn=l->NumberCells_x*l->NumberCells_y;
	nnn=l->NumberCells_z*nn;


	int number_particles=0;

	double cutoff2=SQR(cutoff);

	// scan cells
	for (i=0;i<nnn;i++)
	{
		// reconstruct cartesian coordinates
		iz=i/nn;
		alpha=i%nn;
		ix=alpha%l->NumberCells_x;
		iy=alpha/l->NumberCells_x;

		jx=(ix+1)%l->NumberCells_x;
		jy=(iy+1)%l->NumberCells_y;
		jz=(iz+1)%l->NumberCells_z;

		kx=module(ix-1,l->NumberCells_x);
		ky=module(iy-1,l->NumberCells_y);
		kz=module(iz-1,l->NumberCells_z);

		particle1=(l->HoC)[i];

		while (particle1!=-1)
		{
			//interactions with particles in the same cell
			// and in neighbouring cells
			im->who[particle1]=particle1;
			number_particles++;
			int howmany_particle1=0;
			jvector dist;
			double dist2;

			// 0 same cell
			particle2=(l->LinkedList)[particle1];
			while (particle2!=-1)
			{
				dist.x=pos[particle1].x-pos[particle2].x;
				dist.y=pos[particle1].y-pos[particle2].y;
				dist.z=pos[particle1].z-pos[particle2].z;

				dist.x-=box->x*rint(dist.x/box->x);
				dist.y-=box->y*rint(dist.y/box->y);
				dist.z-=box->z*rint(dist.z/box->z);

				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);

				if (dist2<cutoff2)
				{
					im->with[particle1][howmany_particle1++]=particle2;
				}
				particle2=(l->LinkedList)[particle2];
			}

			// 1 cell y+1
			neighbour=(ix)+(jy)*l->NumberCells_x+(iz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				dist.x=pos[particle1].x-pos[particle2].x;
				dist.y=pos[particle1].y-pos[particle2].y;
				dist.z=pos[particle1].z-pos[particle2].z;

				dist.x-=box->x*rint(dist.x/box->x);
				dist.y-=box->y*rint(dist.y/box->y);
				dist.z-=box->z*rint(dist.z/box->z);

				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);

				if (dist2<cutoff2)
				{
					im->with[particle1][howmany_particle1++]=particle2;
				}

				particle2=(l->LinkedList)[particle2];
			}

			// 2 cell x+1 y+1
			neighbour=(jx)+(jy)*l->NumberCells_x+(iz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				dist.x=pos[particle1].x-pos[particle2].x;
				dist.y=pos[particle1].y-pos[particle2].y;
				dist.z=pos[particle1].z-pos[particle2].z;

				dist.x-=box->x*rint(dist.x/box->x);
				dist.y-=box->y*rint(dist.y/box->y);
				dist.z-=box->z*rint(dist.z/box->z);

				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);

				if (dist2<cutoff2)
				{
					im->with[particle1][howmany_particle1++]=particle2;
				}

				particle2=(l->LinkedList)[particle2];
			}

			// 3 cell x+1
			neighbour=(jx)+(iy)*l->NumberCells_x+(iz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				dist.x=pos[particle1].x-pos[particle2].x;
				dist.y=pos[particle1].y-pos[particle2].y;
				dist.z=pos[particle1].z-pos[particle2].z;

				dist.x-=box->x*rint(dist.x/box->x);
				dist.y-=box->y*rint(dist.y/box->y);
				dist.z-=box->z*rint(dist.z/box->z);

				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);

				if (dist2<cutoff2)
				{
					im->with[particle1][howmany_particle1++]=particle2;
				}

				particle2=(l->LinkedList)[particle2];
			}

			// 4 cell x+1 y-1
			neighbour=(jx)+(ky)*l->NumberCells_x+(iz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				dist.x=pos[particle1].x-pos[particle2].x;
				dist.y=pos[particle1].y-pos[particle2].y;
				dist.z=pos[particle1].z-pos[particle2].z;

				dist.x-=box->x*rint(dist.x/box->x);
				dist.y-=box->y*rint(dist.y/box->y);
				dist.z-=box->z*rint(dist.z/box->z);

				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);

				if (dist2<cutoff2)
				{
					im->with[particle1][howmany_particle1++]=particle2;
				}

				particle2=(l->LinkedList)[particle2];
			}

			// 5 cell z+1 y+1
			neighbour=(ix)+(jy)*l->NumberCells_x+(jz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				dist.x=pos[particle1].x-pos[particle2].x;
				dist.y=pos[particle1].y-pos[particle2].y;
				dist.z=pos[particle1].z-pos[particle2].z;

				dist.x-=box->x*rint(dist.x/box->x);
				dist.y-=box->y*rint(dist.y/box->y);
				dist.z-=box->z*rint(dist.z/box->z);

				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);

				if (dist2<cutoff2)
				{
					im->with[particle1][howmany_particle1++]=particle2;
				}

				particle2=(l->LinkedList)[particle2];
			}

			// 6 cell z+1
			neighbour=(ix)+(iy)*l->NumberCells_x+(jz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				dist.x=pos[particle1].x-pos[particle2].x;
				dist.y=pos[particle1].y-pos[particle2].y;
				dist.z=pos[particle1].z-pos[particle2].z;

				dist.x-=box->x*rint(dist.x/box->x);
				dist.y-=box->y*rint(dist.y/box->y);
				dist.z-=box->z*rint(dist.z/box->z);

				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);

				if (dist2<cutoff2)
				{
					im->with[particle1][howmany_particle1++]=particle2;
				}

				particle2=(l->LinkedList)[particle2];
			}

			// 7 cell z+1 x+1 y+1
			neighbour=(jx)+(jy)*l->NumberCells_x+(jz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				dist.x=pos[particle1].x-pos[particle2].x;
				dist.y=pos[particle1].y-pos[particle2].y;
				dist.z=pos[particle1].z-pos[particle2].z;

				dist.x-=box->x*rint(dist.x/box->x);
				dist.y-=box->y*rint(dist.y/box->y);
				dist.z-=box->z*rint(dist.z/box->z);

				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);

				if (dist2<cutoff2)
				{
					im->with[particle1][howmany_particle1++]=particle2;
				}

				particle2=(l->LinkedList)[particle2];
			}

			// 8 cell z+1 x+1
			neighbour=(jx)+(iy)*l->NumberCells_x+(jz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				dist.x=pos[particle1].x-pos[particle2].x;
				dist.y=pos[particle1].y-pos[particle2].y;
				dist.z=pos[particle1].z-pos[particle2].z;

				dist.x-=box->x*rint(dist.x/box->x);
				dist.y-=box->y*rint(dist.y/box->y);
				dist.z-=box->z*rint(dist.z/box->z);

				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);

				if (dist2<cutoff2)
				{
					im->with[particle1][howmany_particle1++]=particle2;
				}

				particle2=(l->LinkedList)[particle2];
			}

			// 9 cell z+1 x+1 y-1
			neighbour=(jx)+(ky)*l->NumberCells_x+(jz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				dist.x=pos[particle1].x-pos[particle2].x;
				dist.y=pos[particle1].y-pos[particle2].y;
				dist.z=pos[particle1].z-pos[particle2].z;

				dist.x-=box->x*rint(dist.x/box->x);
				dist.y-=box->y*rint(dist.y/box->y);
				dist.z-=box->z*rint(dist.z/box->z);

				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);

				if (dist2<cutoff2)
				{
					im->with[particle1][howmany_particle1++]=particle2;
				}

				particle2=(l->LinkedList)[particle2];
			}

			// 10 cell z-1 y+1
			neighbour=(ix)+(jy)*l->NumberCells_x+(kz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				dist.x=pos[particle1].x-pos[particle2].x;
				dist.y=pos[particle1].y-pos[particle2].y;
				dist.z=pos[particle1].z-pos[particle2].z;

				dist.x-=box->x*rint(dist.x/box->x);
				dist.y-=box->y*rint(dist.y/box->y);
				dist.z-=box->z*rint(dist.z/box->z);

				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);

				if (dist2<cutoff2)
				{
					im->with[particle1][howmany_particle1++]=particle2;
				}

				particle2=(l->LinkedList)[particle2];
			}

			// 11 cell x+1 y+1 z-1
			neighbour=(jx)+(jy)*l->NumberCells_x+(kz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				dist.x=pos[particle1].x-pos[particle2].x;
				dist.y=pos[particle1].y-pos[particle2].y;
				dist.z=pos[particle1].z-pos[particle2].z;

				dist.x-=box->x*rint(dist.x/box->x);
				dist.y-=box->y*rint(dist.y/box->y);
				dist.z-=box->z*rint(dist.z/box->z);

				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);

				if (dist2<cutoff2)
				{
					im->with[particle1][howmany_particle1++]=particle2;
				}

				particle2=(l->LinkedList)[particle2];
			}

			// 12 cell x+1 z-1
			neighbour=(jx)+(iy)*l->NumberCells_x+(kz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				dist.x=pos[particle1].x-pos[particle2].x;
				dist.y=pos[particle1].y-pos[particle2].y;
				dist.z=pos[particle1].z-pos[particle2].z;

				dist.x-=box->x*rint(dist.x/box->x);
				dist.y-=box->y*rint(dist.y/box->y);
				dist.z-=box->z*rint(dist.z/box->z);

				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);

				if (dist2<cutoff2)
				{
					im->with[particle1][howmany_particle1++]=particle2;
				}

				particle2=(l->LinkedList)[particle2];
			}

			// 13 cell z-1 x+1 y-1
			neighbour=(jx)+(ky)*l->NumberCells_x+(kz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				dist.x=pos[particle1].x-pos[particle2].x;
				dist.y=pos[particle1].y-pos[particle2].y;
				dist.z=pos[particle1].z-pos[particle2].z;

				dist.x-=box->x*rint(dist.x/box->x);
				dist.y-=box->y*rint(dist.y/box->y);
				dist.z-=box->z*rint(dist.z/box->z);

				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);

				if (dist2<cutoff2)
				{
					im->with[particle1][howmany_particle1++]=particle2;
				}

				particle2=(l->LinkedList)[particle2];
			}

			im->howmany[particle1]=howmany_particle1;

			particle1=(l->LinkedList)[particle1];
		}
	}

	im->num=number_particles;

	//assert(number_particles==im->size);
}


void calculateExtendedInteractionMapWithCutoffDistance(listcell *l,interactionmap *im,interactionmap *ime,jvector *pos,jvector *box,double cutoff)
{
	int i;
	int nn,nnn;                      // Number of cells on surface and in box
	int particle1,particle2;         // interacting particles
	int neighbour;                   // Neighbour cell
	int ix,iy,iz,alpha;              // present cell coordinates
	int jx,jy,jz;                    // next cell coordinates
	int kx,ky,kz;                    // previous cell coordinates


	nn=l->NumberCells_x*l->NumberCells_y;
	nnn=l->NumberCells_z*nn;


	int number_particles=0;

	double cutoff2=SQR(cutoff);

	im->num_bonds=0;

	// scan cells
	for (i=0;i<nnn;i++)
	{
		// reconstruct cartesian coordinates
		iz=i/nn;
		alpha=i%nn;
		ix=alpha%l->NumberCells_x;
		iy=alpha/l->NumberCells_x;

		jx=(ix+1)%l->NumberCells_x;
		jy=(iy+1)%l->NumberCells_y;
		jz=(iz+1)%l->NumberCells_z;

		kx=module(ix-1,l->NumberCells_x);
		ky=module(iy-1,l->NumberCells_y);
		kz=module(iz-1,l->NumberCells_z);


		particle1=(l->HoC)[i];

		while (particle1!=-1)
		{
			//interactions with particles in the same cell
			// and in neighbouring cells

			im->who[particle1]=particle1;
			ime->who[particle1]=particle1;
			int howmany_particle1=0;
			number_particles++;
			jvector dist;
			double dist2;

			// 0 same cell
			particle2=(l->LinkedList)[particle1];
			while (particle2!=-1)
			{
				dist.x=pos[particle1].x-pos[particle2].x;
				dist.y=pos[particle1].y-pos[particle2].y;
				dist.z=pos[particle1].z-pos[particle2].z;

				dist.x-=box->x*rint(dist.x/box->x);
				dist.y-=box->y*rint(dist.y/box->y);
				dist.z-=box->z*rint(dist.z/box->z);

				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);

				if (dist2<cutoff2)
				{

					im->with[particle1][howmany_particle1++]=particle2;

					ime->with[particle1][ime->howmany[particle1]++]=particle2;
					ime->with[particle2][ime->howmany[particle2]++]=particle1;

					im->num_bonds++;

				}
				particle2=(l->LinkedList)[particle2];
			}

			// 1 cell y+1
			neighbour=(ix)+(jy)*l->NumberCells_x+(iz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				dist.x=pos[particle1].x-pos[particle2].x;
				dist.y=pos[particle1].y-pos[particle2].y;
				dist.z=pos[particle1].z-pos[particle2].z;

				dist.x-=box->x*rint(dist.x/box->x);
				dist.y-=box->y*rint(dist.y/box->y);
				dist.z-=box->z*rint(dist.z/box->z);

				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);

				if (dist2<cutoff2)
				{
					im->with[particle1][howmany_particle1++]=particle2;

					ime->with[particle1][ime->howmany[particle1]++]=particle2;
					ime->with[particle2][ime->howmany[particle2]++]=particle1;

					im->num_bonds++;
				}
				particle2=(l->LinkedList)[particle2];
			}

			// 2 cell x+1 y+1
			neighbour=(jx)+(jy)*l->NumberCells_x+(iz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				dist.x=pos[particle1].x-pos[particle2].x;
				dist.y=pos[particle1].y-pos[particle2].y;
				dist.z=pos[particle1].z-pos[particle2].z;

				dist.x-=box->x*rint(dist.x/box->x);
				dist.y-=box->y*rint(dist.y/box->y);
				dist.z-=box->z*rint(dist.z/box->z);

				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);

				if (dist2<cutoff2)
				{
					im->with[particle1][howmany_particle1++]=particle2;

					ime->with[particle1][ime->howmany[particle1]++]=particle2;
					ime->with[particle2][ime->howmany[particle2]++]=particle1;

					im->num_bonds++;
				}
				particle2=(l->LinkedList)[particle2];
			}

			// 3 cell x+1
			neighbour=(jx)+(iy)*l->NumberCells_x+(iz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				dist.x=pos[particle1].x-pos[particle2].x;
				dist.y=pos[particle1].y-pos[particle2].y;
				dist.z=pos[particle1].z-pos[particle2].z;

				dist.x-=box->x*rint(dist.x/box->x);
				dist.y-=box->y*rint(dist.y/box->y);
				dist.z-=box->z*rint(dist.z/box->z);

				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);

				if (dist2<cutoff2)
				{
					im->with[particle1][howmany_particle1++]=particle2;

					ime->with[particle1][ime->howmany[particle1]++]=particle2;
					ime->with[particle2][ime->howmany[particle2]++]=particle1;

					im->num_bonds++;
				}
				particle2=(l->LinkedList)[particle2];
			}

			// 4 cell x+1 y-1
			neighbour=(jx)+(ky)*l->NumberCells_x+(iz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				dist.x=pos[particle1].x-pos[particle2].x;
				dist.y=pos[particle1].y-pos[particle2].y;
				dist.z=pos[particle1].z-pos[particle2].z;

				dist.x-=box->x*rint(dist.x/box->x);
				dist.y-=box->y*rint(dist.y/box->y);
				dist.z-=box->z*rint(dist.z/box->z);

				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);

				if (dist2<cutoff2)
				{
					im->with[particle1][howmany_particle1++]=particle2;

					ime->with[particle1][ime->howmany[particle1]++]=particle2;
					ime->with[particle2][ime->howmany[particle2]++]=particle1;

					im->num_bonds++;
				}
				particle2=(l->LinkedList)[particle2];
			}

			// 5 cell z+1 y+1
			neighbour=(ix)+(jy)*l->NumberCells_x+(jz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				dist.x=pos[particle1].x-pos[particle2].x;
				dist.y=pos[particle1].y-pos[particle2].y;
				dist.z=pos[particle1].z-pos[particle2].z;

				dist.x-=box->x*rint(dist.x/box->x);
				dist.y-=box->y*rint(dist.y/box->y);
				dist.z-=box->z*rint(dist.z/box->z);

				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);

				if (dist2<cutoff2)
				{
					im->with[particle1][howmany_particle1++]=particle2;

					ime->with[particle1][ime->howmany[particle1]++]=particle2;
					ime->with[particle2][ime->howmany[particle2]++]=particle1;

					im->num_bonds++;
				}
				particle2=(l->LinkedList)[particle2];
			}

			// 6 cell z+1
			neighbour=(ix)+(iy)*l->NumberCells_x+(jz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				dist.x=pos[particle1].x-pos[particle2].x;
				dist.y=pos[particle1].y-pos[particle2].y;
				dist.z=pos[particle1].z-pos[particle2].z;

				dist.x-=box->x*rint(dist.x/box->x);
				dist.y-=box->y*rint(dist.y/box->y);
				dist.z-=box->z*rint(dist.z/box->z);

				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);

				if (dist2<cutoff2)
				{
					im->with[particle1][howmany_particle1++]=particle2;

					ime->with[particle1][ime->howmany[particle1]++]=particle2;
					ime->with[particle2][ime->howmany[particle2]++]=particle1;

					im->num_bonds++;
				}
				particle2=(l->LinkedList)[particle2];
			}

			// 7 cell z+1 x+1 y+1
			neighbour=(jx)+(jy)*l->NumberCells_x+(jz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				dist.x=pos[particle1].x-pos[particle2].x;
				dist.y=pos[particle1].y-pos[particle2].y;
				dist.z=pos[particle1].z-pos[particle2].z;

				dist.x-=box->x*rint(dist.x/box->x);
				dist.y-=box->y*rint(dist.y/box->y);
				dist.z-=box->z*rint(dist.z/box->z);

				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);

				if (dist2<cutoff2)
				{
					im->with[particle1][howmany_particle1++]=particle2;

					ime->with[particle1][ime->howmany[particle1]++]=particle2;
					ime->with[particle2][ime->howmany[particle2]++]=particle1;

					im->num_bonds++;
				}
				particle2=(l->LinkedList)[particle2];
			}

			// 8 cell z+1 x+1
			neighbour=(jx)+(iy)*l->NumberCells_x+(jz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				dist.x=pos[particle1].x-pos[particle2].x;
				dist.y=pos[particle1].y-pos[particle2].y;
				dist.z=pos[particle1].z-pos[particle2].z;

				dist.x-=box->x*rint(dist.x/box->x);
				dist.y-=box->y*rint(dist.y/box->y);
				dist.z-=box->z*rint(dist.z/box->z);

				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);

				if (dist2<cutoff2)
				{
					im->with[particle1][howmany_particle1++]=particle2;

					ime->with[particle1][ime->howmany[particle1]++]=particle2;
					ime->with[particle2][ime->howmany[particle2]++]=particle1;

					im->num_bonds++;
				}
				particle2=(l->LinkedList)[particle2];
			}

			// 9 cell z+1 x+1 y-1
			neighbour=(jx)+(ky)*l->NumberCells_x+(jz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				dist.x=pos[particle1].x-pos[particle2].x;
				dist.y=pos[particle1].y-pos[particle2].y;
				dist.z=pos[particle1].z-pos[particle2].z;

				dist.x-=box->x*rint(dist.x/box->x);
				dist.y-=box->y*rint(dist.y/box->y);
				dist.z-=box->z*rint(dist.z/box->z);

				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);

				if (dist2<cutoff2)
				{
					im->with[particle1][howmany_particle1++]=particle2;

					ime->with[particle1][ime->howmany[particle1]++]=particle2;
					ime->with[particle2][ime->howmany[particle2]++]=particle1;

					im->num_bonds++;
				}
				particle2=(l->LinkedList)[particle2];
			}

			// 10 cell z-1 y+1
			neighbour=(ix)+(jy)*l->NumberCells_x+(kz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				dist.x=pos[particle1].x-pos[particle2].x;
				dist.y=pos[particle1].y-pos[particle2].y;
				dist.z=pos[particle1].z-pos[particle2].z;

				dist.x-=box->x*rint(dist.x/box->x);
				dist.y-=box->y*rint(dist.y/box->y);
				dist.z-=box->z*rint(dist.z/box->z);

				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);

				if (dist2<cutoff2)
				{
					im->with[particle1][howmany_particle1++]=particle2;

					ime->with[particle1][ime->howmany[particle1]++]=particle2;
					ime->with[particle2][ime->howmany[particle2]++]=particle1;

					im->num_bonds++;
				}
				particle2=(l->LinkedList)[particle2];
			}

			// 11 cell x+1 y+1 z-1
			neighbour=(jx)+(jy)*l->NumberCells_x+(kz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				dist.x=pos[particle1].x-pos[particle2].x;
				dist.y=pos[particle1].y-pos[particle2].y;
				dist.z=pos[particle1].z-pos[particle2].z;

				dist.x-=box->x*rint(dist.x/box->x);
				dist.y-=box->y*rint(dist.y/box->y);
				dist.z-=box->z*rint(dist.z/box->z);

				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);

				if (dist2<cutoff2)
				{
					im->with[particle1][howmany_particle1++]=particle2;

					ime->with[particle1][ime->howmany[particle1]++]=particle2;
					ime->with[particle2][ime->howmany[particle2]++]=particle1;

					im->num_bonds++;
				}
				particle2=(l->LinkedList)[particle2];
			}

			// 12 cell x+1 z-1
			neighbour=(jx)+(iy)*l->NumberCells_x+(kz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				dist.x=pos[particle1].x-pos[particle2].x;
				dist.y=pos[particle1].y-pos[particle2].y;
				dist.z=pos[particle1].z-pos[particle2].z;

				dist.x-=box->x*rint(dist.x/box->x);
				dist.y-=box->y*rint(dist.y/box->y);
				dist.z-=box->z*rint(dist.z/box->z);

				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);

				if (dist2<cutoff2)
				{
					im->with[particle1][howmany_particle1++]=particle2;

					ime->with[particle1][ime->howmany[particle1]++]=particle2;
					ime->with[particle2][ime->howmany[particle2]++]=particle1;

					im->num_bonds++;
				}
				particle2=(l->LinkedList)[particle2];
			}

			// 13 cell z-1 x+1 y-1
			neighbour=(jx)+(ky)*l->NumberCells_x+(kz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				dist.x=pos[particle1].x-pos[particle2].x;
				dist.y=pos[particle1].y-pos[particle2].y;
				dist.z=pos[particle1].z-pos[particle2].z;

				dist.x-=box->x*rint(dist.x/box->x);
				dist.y-=box->y*rint(dist.y/box->y);
				dist.z-=box->z*rint(dist.z/box->z);

				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);

				if (dist2<cutoff2)
				{
					im->with[particle1][howmany_particle1++]=particle2;

					ime->with[particle1][ime->howmany[particle1]++]=particle2;
					ime->with[particle2][ime->howmany[particle2]++]=particle1;

					im->num_bonds++;
				}
				particle2=(l->LinkedList)[particle2];
			}

			im->howmany[particle1]=howmany_particle1;

			particle1=(l->LinkedList)[particle1];
		}
	}

	im->num=number_particles;

	ime->num_bonds=im->num_bonds;

	assert(number_particles==im->size);
}




static void insertionSort_dist2(int *number,int num_el,double *dist2,double dist2_el,int *length)
{
	int l=*length;
	number[l]=num_el;
	dist2[l]=dist2_el;
	while ((l>0) && (dist2[l]<dist2[l-1]))
	{
		double buffer_dist2;
		buffer_dist2=dist2[l];
		dist2[l]=dist2[l-1];
		dist2[l-1]=buffer_dist2;

		int buffer_num;
		buffer_num=number[l];
		number[l]=number[l-1];
		number[l-1]=buffer_num;

		l--;
	}
	(*length)++;
}



void calculateInteractionMapWithCutoffDistanceOrdered(listcell *l,interactionmap *ime,jvector *pos,jvector *box,double cutoff)
{

	int i;
	int nn,nnn;                      // Number of cells on surface and in box
	int particle1,particle2;         // interacting particles
	int neighbour;                   // Neighbour cell
	int ix,iy,iz,alpha;              // present cell coordinates
	int jx,jy,jz;                    // next cell coordinates
	int kx,ky,kz;                    // previous cell coordinates


	nn=l->NumberCells_x*l->NumberCells_y;
	nnn=l->NumberCells_z*nn;


	int number_particles=0;

	double cutoff2=SQR(cutoff);

	ime->num_bonds=0;


	// abbiamo bisogno della mappa delle distanze
	int ncolloids=ime->size;
	double **dist2map;
	Matrix2D(dist2map,ncolloids,ncolloids,double);



	// scan cells
	for (i=0;i<nnn;i++)
	{
		// reconstruct cartesian coordinates
		iz=i/nn;
		alpha=i%nn;
		ix=alpha%l->NumberCells_x;
		iy=alpha/l->NumberCells_x;

		jx=(ix+1)%l->NumberCells_x;
		jy=(iy+1)%l->NumberCells_y;
		jz=(iz+1)%l->NumberCells_z;

		kx=module(ix-1,l->NumberCells_x);
		ky=module(iy-1,l->NumberCells_y);
		kz=module(iz-1,l->NumberCells_z);


		particle1=(l->HoC)[i];

		while (particle1!=-1)
		{
			//interactions with particles in the same cell
			// and in neighbouring cells
			ime->who[particle1]=particle1;
			number_particles++;
			jvector dist;
			double dist2;

			// 0 same cell
			particle2=(l->LinkedList)[particle1];
			while (particle2!=-1)
			{
				dist.x=pos[particle1].x-pos[particle2].x;
				dist.y=pos[particle1].y-pos[particle2].y;
				dist.z=pos[particle1].z-pos[particle2].z;

				dist.x-=box->x*rint(dist.x/box->x);
				dist.y-=box->y*rint(dist.y/box->y);
				dist.z-=box->z*rint(dist.z/box->z);

				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);

				if (dist2<cutoff2)
				{

					insertionSort_dist2(ime->with[particle1],particle2,dist2map[particle1],dist2,ime->howmany+particle1);
					insertionSort_dist2(ime->with[particle2],particle1,dist2map[particle2],dist2,ime->howmany+particle2);

					ime->num_bonds++;
				}
				particle2=(l->LinkedList)[particle2];
			}

			// 1 cell y+1
			neighbour=(ix)+(jy)*l->NumberCells_x+(iz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				dist.x=pos[particle1].x-pos[particle2].x;
				dist.y=pos[particle1].y-pos[particle2].y;
				dist.z=pos[particle1].z-pos[particle2].z;

				dist.x-=box->x*rint(dist.x/box->x);
				dist.y-=box->y*rint(dist.y/box->y);
				dist.z-=box->z*rint(dist.z/box->z);

				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);

				if (dist2<cutoff2)
				{
					insertionSort_dist2(ime->with[particle1],particle2,dist2map[particle1],dist2,ime->howmany+particle1);
					insertionSort_dist2(ime->with[particle2],particle1,dist2map[particle2],dist2,ime->howmany+particle2);

					ime->num_bonds++;
				}
				particle2=(l->LinkedList)[particle2];
			}

			// 2 cell x+1 y+1
			neighbour=(jx)+(jy)*l->NumberCells_x+(iz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				dist.x=pos[particle1].x-pos[particle2].x;
				dist.y=pos[particle1].y-pos[particle2].y;
				dist.z=pos[particle1].z-pos[particle2].z;

				dist.x-=box->x*rint(dist.x/box->x);
				dist.y-=box->y*rint(dist.y/box->y);
				dist.z-=box->z*rint(dist.z/box->z);

				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);

				if (dist2<cutoff2)
				{
					insertionSort_dist2(ime->with[particle1],particle2,dist2map[particle1],dist2,ime->howmany+particle1);
					insertionSort_dist2(ime->with[particle2],particle1,dist2map[particle2],dist2,ime->howmany+particle2);

					ime->num_bonds++;
				}
				particle2=(l->LinkedList)[particle2];
			}

			// 3 cell x+1
			neighbour=(jx)+(iy)*l->NumberCells_x+(iz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				dist.x=pos[particle1].x-pos[particle2].x;
				dist.y=pos[particle1].y-pos[particle2].y;
				dist.z=pos[particle1].z-pos[particle2].z;

				dist.x-=box->x*rint(dist.x/box->x);
				dist.y-=box->y*rint(dist.y/box->y);
				dist.z-=box->z*rint(dist.z/box->z);

				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);

				if (dist2<cutoff2)
				{
					insertionSort_dist2(ime->with[particle1],particle2,dist2map[particle1],dist2,ime->howmany+particle1);
					insertionSort_dist2(ime->with[particle2],particle1,dist2map[particle2],dist2,ime->howmany+particle2);

					ime->num_bonds++;
				}
				particle2=(l->LinkedList)[particle2];
			}

			// 4 cell x+1 y-1
			neighbour=(jx)+(ky)*l->NumberCells_x+(iz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				dist.x=pos[particle1].x-pos[particle2].x;
				dist.y=pos[particle1].y-pos[particle2].y;
				dist.z=pos[particle1].z-pos[particle2].z;

				dist.x-=box->x*rint(dist.x/box->x);
				dist.y-=box->y*rint(dist.y/box->y);
				dist.z-=box->z*rint(dist.z/box->z);

				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);

				if (dist2<cutoff2)
				{
					insertionSort_dist2(ime->with[particle1],particle2,dist2map[particle1],dist2,ime->howmany+particle1);
					insertionSort_dist2(ime->with[particle2],particle1,dist2map[particle2],dist2,ime->howmany+particle2);

					ime->num_bonds++;
				}
				particle2=(l->LinkedList)[particle2];
			}

			// 5 cell z+1 y+1
			neighbour=(ix)+(jy)*l->NumberCells_x+(jz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				dist.x=pos[particle1].x-pos[particle2].x;
				dist.y=pos[particle1].y-pos[particle2].y;
				dist.z=pos[particle1].z-pos[particle2].z;

				dist.x-=box->x*rint(dist.x/box->x);
				dist.y-=box->y*rint(dist.y/box->y);
				dist.z-=box->z*rint(dist.z/box->z);

				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);

				if (dist2<cutoff2)
				{
					insertionSort_dist2(ime->with[particle1],particle2,dist2map[particle1],dist2,ime->howmany+particle1);
					insertionSort_dist2(ime->with[particle2],particle1,dist2map[particle2],dist2,ime->howmany+particle2);

					ime->num_bonds++;
				}
				particle2=(l->LinkedList)[particle2];
			}

			// 6 cell z+1
			neighbour=(ix)+(iy)*l->NumberCells_x+(jz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				dist.x=pos[particle1].x-pos[particle2].x;
				dist.y=pos[particle1].y-pos[particle2].y;
				dist.z=pos[particle1].z-pos[particle2].z;

				dist.x-=box->x*rint(dist.x/box->x);
				dist.y-=box->y*rint(dist.y/box->y);
				dist.z-=box->z*rint(dist.z/box->z);

				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);

				if (dist2<cutoff2)
				{
					insertionSort_dist2(ime->with[particle1],particle2,dist2map[particle1],dist2,ime->howmany+particle1);
					insertionSort_dist2(ime->with[particle2],particle1,dist2map[particle2],dist2,ime->howmany+particle2);

					ime->num_bonds++;
				}
				particle2=(l->LinkedList)[particle2];
			}

			// 7 cell z+1 x+1 y+1
			neighbour=(jx)+(jy)*l->NumberCells_x+(jz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				dist.x=pos[particle1].x-pos[particle2].x;
				dist.y=pos[particle1].y-pos[particle2].y;
				dist.z=pos[particle1].z-pos[particle2].z;

				dist.x-=box->x*rint(dist.x/box->x);
				dist.y-=box->y*rint(dist.y/box->y);
				dist.z-=box->z*rint(dist.z/box->z);

				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);

				if (dist2<cutoff2)
				{
					insertionSort_dist2(ime->with[particle1],particle2,dist2map[particle1],dist2,ime->howmany+particle1);
					insertionSort_dist2(ime->with[particle2],particle1,dist2map[particle2],dist2,ime->howmany+particle2);

					ime->num_bonds++;
				}
				particle2=(l->LinkedList)[particle2];
			}

			// 8 cell z+1 x+1
			neighbour=(jx)+(iy)*l->NumberCells_x+(jz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				dist.x=pos[particle1].x-pos[particle2].x;
				dist.y=pos[particle1].y-pos[particle2].y;
				dist.z=pos[particle1].z-pos[particle2].z;

				dist.x-=box->x*rint(dist.x/box->x);
				dist.y-=box->y*rint(dist.y/box->y);
				dist.z-=box->z*rint(dist.z/box->z);

				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);

				if (dist2<cutoff2)
				{
					insertionSort_dist2(ime->with[particle1],particle2,dist2map[particle1],dist2,ime->howmany+particle1);
					insertionSort_dist2(ime->with[particle2],particle1,dist2map[particle2],dist2,ime->howmany+particle2);

					ime->num_bonds++;
				}
				particle2=(l->LinkedList)[particle2];
			}

			// 9 cell z+1 x+1 y-1
			neighbour=(jx)+(ky)*l->NumberCells_x+(jz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				dist.x=pos[particle1].x-pos[particle2].x;
				dist.y=pos[particle1].y-pos[particle2].y;
				dist.z=pos[particle1].z-pos[particle2].z;

				dist.x-=box->x*rint(dist.x/box->x);
				dist.y-=box->y*rint(dist.y/box->y);
				dist.z-=box->z*rint(dist.z/box->z);

				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);

				if (dist2<cutoff2)
				{
					insertionSort_dist2(ime->with[particle1],particle2,dist2map[particle1],dist2,ime->howmany+particle1);
					insertionSort_dist2(ime->with[particle2],particle1,dist2map[particle2],dist2,ime->howmany+particle2);

					ime->num_bonds++;
				}
				particle2=(l->LinkedList)[particle2];
			}

			// 10 cell z-1 y+1
			neighbour=(ix)+(jy)*l->NumberCells_x+(kz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				dist.x=pos[particle1].x-pos[particle2].x;
				dist.y=pos[particle1].y-pos[particle2].y;
				dist.z=pos[particle1].z-pos[particle2].z;

				dist.x-=box->x*rint(dist.x/box->x);
				dist.y-=box->y*rint(dist.y/box->y);
				dist.z-=box->z*rint(dist.z/box->z);

				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);

				if (dist2<cutoff2)
				{
					insertionSort_dist2(ime->with[particle1],particle2,dist2map[particle1],dist2,ime->howmany+particle1);
					insertionSort_dist2(ime->with[particle2],particle1,dist2map[particle2],dist2,ime->howmany+particle2);

					ime->num_bonds++;
				}
				particle2=(l->LinkedList)[particle2];
			}

			// 11 cell x+1 y+1 z-1
			neighbour=(jx)+(jy)*l->NumberCells_x+(kz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				dist.x=pos[particle1].x-pos[particle2].x;
				dist.y=pos[particle1].y-pos[particle2].y;
				dist.z=pos[particle1].z-pos[particle2].z;

				dist.x-=box->x*rint(dist.x/box->x);
				dist.y-=box->y*rint(dist.y/box->y);
				dist.z-=box->z*rint(dist.z/box->z);

				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);

				if (dist2<cutoff2)
				{
					insertionSort_dist2(ime->with[particle1],particle2,dist2map[particle1],dist2,ime->howmany+particle1);
					insertionSort_dist2(ime->with[particle2],particle1,dist2map[particle2],dist2,ime->howmany+particle2);

					ime->num_bonds++;
				}
				particle2=(l->LinkedList)[particle2];
			}

			// 12 cell x+1 z-1
			neighbour=(jx)+(iy)*l->NumberCells_x+(kz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				dist.x=pos[particle1].x-pos[particle2].x;
				dist.y=pos[particle1].y-pos[particle2].y;
				dist.z=pos[particle1].z-pos[particle2].z;

				dist.x-=box->x*rint(dist.x/box->x);
				dist.y-=box->y*rint(dist.y/box->y);
				dist.z-=box->z*rint(dist.z/box->z);

				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);

				if (dist2<cutoff2)
				{
					insertionSort_dist2(ime->with[particle1],particle2,dist2map[particle1],dist2,ime->howmany+particle1);
					insertionSort_dist2(ime->with[particle2],particle1,dist2map[particle2],dist2,ime->howmany+particle2);

					ime->num_bonds++;
				}
				particle2=(l->LinkedList)[particle2];
			}

			// 13 cell z-1 x+1 y-1
			neighbour=(jx)+(ky)*l->NumberCells_x+(kz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				dist.x=pos[particle1].x-pos[particle2].x;
				dist.y=pos[particle1].y-pos[particle2].y;
				dist.z=pos[particle1].z-pos[particle2].z;

				dist.x-=box->x*rint(dist.x/box->x);
				dist.y-=box->y*rint(dist.y/box->y);
				dist.z-=box->z*rint(dist.z/box->z);

				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);

				if (dist2<cutoff2)
				{
					insertionSort_dist2(ime->with[particle1],particle2,dist2map[particle1],dist2,ime->howmany+particle1);
					insertionSort_dist2(ime->with[particle2],particle1,dist2map[particle2],dist2,ime->howmany+particle2);

					ime->num_bonds++;
				}
				particle2=(l->LinkedList)[particle2];
			}

			particle1=(l->LinkedList)[particle1];
		}
	}

	ime->num=number_particles;

	assert(number_particles==ime->size);

	Free2D(dist2map);
}




void changeCell(listcell *l,const jvector *oldpos,const jvector *newpos,int num)
{
	int ncell_old,ncell_new;                       // cell number
	int posx,posy,posz;              // cell coordinates


	posx=(int)floor(oldpos->x/l->CellSize_x);
	posy=(int)floor(oldpos->y/l->CellSize_y);
	posz=(int)floor(oldpos->z/l->CellSize_z);

	posx=module(posx,l->NumberCells_x);
	posy=module(posy,l->NumberCells_y);
	posz=module(posz,l->NumberCells_z);

	ncell_old=posx+posy*l->NumberCells_x+posz*(l->NumberCells_x*l->NumberCells_y);

	posx=(int)floor(newpos->x/l->CellSize_x);
	posy=(int)floor(newpos->y/l->CellSize_y);
	posz=(int)floor(newpos->z/l->CellSize_z);

	posx=module(posx,l->NumberCells_x);
	posy=module(posy,l->NumberCells_y);
	posz=module(posz,l->NumberCells_z);

	ncell_new=posx+posy*l->NumberCells_x+posz*(l->NumberCells_x*l->NumberCells_y);


	if (ncell_old==ncell_new)
		return;


	int *old_item,*new_item;

	// delete the old position
	int j=(l->HoC)[ncell_old];

	old_item=l->HoC+ncell_old;

	while (j!=-1)
	{
		new_item=l->LinkedList+j;

		if (j==num)
		{
			*old_item=*new_item;
			break;
		}

		j=(l->LinkedList)[j];

		old_item=new_item;
	}

	// you should never be here because there must be a cancellation
	// in the list. This means that there are some problems in the list
	// assert(j!=-1);

	// add the new position
	(l->LinkedList)[num]=(l->HoC)[ncell_new];
	(l->HoC)[ncell_new]=num;


}

void addToList(listcell *l,const jvector *pos,int num)
{
	int ncell;                       // cell number
	int posx,posy,posz;              // cell coordinates


	posx=(int)floor(pos->x/l->CellSize_x);
	posy=(int)floor(pos->y/l->CellSize_y);
	posz=(int)floor(pos->z/l->CellSize_z);

	posx=module(posx,l->NumberCells_x);
	posy=module(posy,l->NumberCells_y);
	posz=module(posz,l->NumberCells_z);

	ncell=posx+posy*l->NumberCells_x+posz*(l->NumberCells_x*l->NumberCells_y);

	(l->LinkedList)[num]=(l->HoC)[ncell];
	(l->HoC)[ncell]=num;
}
